mod accumulator;
mod traverse;
mod database;

use anyhow::Result;
use bio::io::fasta::Reader;
use image::ImageBuffer;

use std::f64::consts::PI;
use std::path::PathBuf;
use std::time::Duration;
use std::time::Instant;
use std::{str, thread, time};

use clap::{Parser, Subcommand};
use std::sync::mpsc::channel;

use crate::accumulator::Accumulator;
use crate::traverse::make_points;

#[derive(Parser)]
#[clap(name = "hustar")]
#[clap(about = "builds subsequence frequency index of genomic data", long_about = None)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    #[clap(arg_required_else_help = true)]
    Build {
        fasta_file: String,
        index_file: String,
        sequence_length: usize,
    },
    #[clap(arg_required_else_help = true)]
    Visualize {
        index_file: String,
        sequence_length: usize,
        side_length: usize,
    },
}

fn filter_n(seq: &&[u8]) -> bool {
    for l in seq.iter() {
        let l = l.to_ascii_lowercase();
        match l {
            b'n' => return false,
            b'm' => return false,
            b'r' => return false,
            b'y' => return false,
            b'w' => return false,
            b'k' => return false,
            b'b' => return false,
            b's' => return false,
            _ => continue,
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use crate::*;

    #[test]
    fn test_min_max_thetas() {
        let width = 2;
        let height = 2;
        let x = 1;
        let y = 1;
        let (min, max) = min_max_thetas(width, height, x, y);

        assert_eq!(min, 0.0);
        assert_eq!(max, PI / 2.0);

        let width = 200;
        let height = 200;
        let x = 110;
        let y = 145;
        let (min, max) = min_max_thetas(width, height, x, y);

        assert_eq!(min, 1.33105321792444);
        assert_eq!(max, 1.356735643231075);
    }

    #[test]
    fn test_seqlen_from_position() {
        let width = 200;
        let height = 200;
        let x = 100;
        let y = 100;
        let seqlen = seqlen_from_position(width, height, x, y, 90);

        assert_eq!(seqlen, 1);

        let width = 200;
        let height = 200;
        let x = 200;
        let y = 100;
        let seqlen = seqlen_from_position(width, height, x, y, 90);

        assert_eq!(seqlen, 90);

        let width = 200;
        let height = 200;
        let x = 190;
        let y = 100;
        let seqlen = seqlen_from_position(width, height, x, y, 100);

        assert_eq!(seqlen, 95);
    }
}

fn min_max_thetas(width: usize, height: usize, x: u32, y: u32) -> (f64, f64) {
    let (min_theta, max_theta) = vec![(x, y), (x + 1, y), (x, y + 1), (x + 1, y + 1)]
        .iter()
        .fold((None, None), |(min, max), (x, y)| {
            let x: i32 = *x as i32 - (width / 2) as i32;
            let y: i32 = *y as i32 - (height / 2) as i32;
            let theta = {
                let mut theta = (y as f64 / x as f64).atan();
                if x == 0 {
                    theta = PI / 2.0;
                }
                if x < 0 {
                    theta += PI;
                }
                if theta < 0.0 {
                    theta = theta + PI + PI;
                }
                // rotate for effect; note, this complicates the assumptions on min/max for angles
                // however, a min-max on sequence could still be done.
                // theta = theta + PI / 7.23;
                // if theta > 2.0 * PI {
                //     theta = theta - 2.0 * PI;
                // }
                theta
            };
            let min = match min {
                Some(current) => {
                    if current < theta {
                        Some(current)
                    } else {
                        Some(theta)
                    }
                }
                None => Some(theta),
            };
            let max = match max {
                Some(current) => {
                    if current < theta {
                        Some(theta)
                    } else {
                        Some(current)
                    }
                }
                None => Some(theta),
            };
            (min, max)
        });
    let min_theta = min_theta.unwrap();
    let max_theta = max_theta.unwrap();
    (min_theta, max_theta)
}

fn seqlen_from_position(width: usize, height: usize, x: u32, y: u32, max_seqlen: usize) -> usize {
    let x: i32 = x as i32 - (width / 2) as i32;
    let y: i32 = y as i32 - (height / 2) as i32;
    let r = ((x * x + y * y) as f64).sqrt();
    let max = width as f64 / 2.0;
    let p = r / max;
    // println!("r={} max={} p={}", r, max, p);
    let p = p.min(1.0);
    let p = p.sqrt().sqrt(); // express smaller radii less
    let interp = p * max_seqlen as f64;
    // println!("p={} max_seqlen={} interp={}", p, max_seqlen, interp);
    interp.ceil().max(1.0) as usize
    // (p * max_seqlen as f64).ceil().min(1.0) as usize

    // let max = width as f64 / max_seqlen as f64 / 2.0;

    // ((r / max).floor() as usize + 1).min(max_seqlen)
}

fn main() {
    let args = Cli::parse();
    match &args.command {
        Commands::Build {
            fasta_file,
            index_file,
            sequence_length,
        } => {
            create(fasta_file, index_file, *sequence_length);
        }
        Commands::Visualize {
            index_file,
            sequence_length,
            side_length,
        } => {
            print(index_file, *sequence_length, *side_length).expect("while printing");
        }
        _ => panic!("oh no"),
    }
}

fn print(index_file: &str, seqlen: usize, side_length: usize) -> Result<()> {
    use crate::database::Database;

    let cpus = num_cpus::get();
    let thread_count = cpus;
    let thread_count = 1;

    println!("printing {} with thread_count={}", index_file, thread_count);
    // let size = buf_size_bytes(seqlen);
    let m = Database::open(index_file)?;

    let width = side_length;
    let height = side_length;

    let pixels = make_points(width as u32, height as u32, seqlen as u32);
    let chunk_size = pixels.len() / thread_count;
    let chunks = pixels.chunks(chunk_size);

    let (tx, rx) = channel();
    for worker_id in 0..thread_count {
        let tx = tx.clone();
        let m = m.clone();
        let pixels = pixels.clone();
        thread::spawn(move || {
            let mut acc = Accumulator::default();
            for p in pixels {
                let (gte, lt) = p.index_range();
                let c = acc.sum_to(&m, gte, lt);
                let val = (
                    p.w as usize,
                    p.h as usize,
                    (c, lt - gte),
                    p.seqlen,
                );
                let (gtes, lts) = p.seq_range();
                println!("sending point={:?} range=({}, {}) range=({}, {}) val={:?}", p, gte, lt, gtes, lts, val);
                tx.send(val).unwrap();
            }
            // let seqrange = make_coords_to_seq_range(width as usize, height as usize, seqlen);
            // for y in 0..height {
            //     for x in 0..width {
            //         if x % thread_count != worker_id {
            //             continue;
            //         }
            //         let testx: i32 = x as i32 - (width / 2) as i32;
            //         let testy: i32 = y as i32 - (height / 2) as i32;
            //         if (-3 < testx && testx < 3) || (-3 < testy && testy < 3) {
            //             // println!("{} {}", min, max);
            //             // TODO figure out why bright lines are forming on axis
            //             continue;
            //         }
            //         let (min, max) = seqrange(x as u32, y as u32);

            //         let len = min.len();
            //         let c = m.count(seq_to_idx(min.as_bytes()), seq_to_idx(max.as_bytes()));
            //         tx.send((x, y, c, len)).unwrap();
            //     }
            // }
            // tx.send((0, 0, (0, 0), 0)).unwrap();
            // println!("worker_id={} exiting", worker_id)
        });
    }
    println!("creating image");
    let mut buf: Vec<u64> = vec![0; width * height];
    let mut maxes: Vec<u64> = vec![0; seqlen + 1];
    let mut maxesbuf: Vec<u64> = vec![0; width * height];
    println!("creating count buffer");
    drop(tx);
    let mut last = Instant::now();
    let mut counter = 0;
    let mut count_sequences = 0;
    while let Ok((x, y, (count, count_unique_sequences), len)) = rx.recv() {
        buf[y * width + x] = count;
        if maxes[len] < count {
            maxes[len] = count;
        }
        counter = counter + 1;
        count_sequences += count_unique_sequences;
        let now = Instant::now();
        if now.duration_since(last) > Duration::from_secs_f64(2.3) {
            last = now;
            let total = width * height;
            let percentage = counter as f64 / total as f64 * 100.0;
            println!(
                "processed {} pixels of {} ({:.2}%)",
                counter, total, percentage
            );
        }
    }
    println!("visited {} sequences", count_sequences);
    println!("creating maxes buffer");
    for y in 0..height {
        for x in 0..width {
            let seqlen = seqlen_from_position(width, height, x as u32, y as u32, seqlen);
            maxesbuf[y * height + x] = maxes[seqlen];
        }
    }
    println!("creating image buffer");
    let mut img = ImageBuffer::from_fn(width as u32, height as u32, |x, y| {
        let x = x as usize;
        let y = y as usize;
        let c = buf[y * width + x];
        let t = maxesbuf[y * height + x];
        let p = c as f64 / t as f64;
        let p = p.sqrt().sqrt();
        let l = (p * 255.0) as u8;
        image::Luma([l])
    });
    img.save_with_format("out.png", image::ImageFormat::Png)
        .expect("while writing image");
    Ok(())
}

fn create<P: Into<PathBuf>>(fasta_file: &str, outpath: P, seqlen: usize) -> Result<()> {
    use crate::database;

    let mut b = database::DatabaseMut::create(outpath, seqlen)?;
    let r = Reader::from_file(fasta_file).unwrap();
    let mut records = r.records();
    while let Some(Ok(record)) = records.next() {
        if record.id() != "CM000667.2" {
            continue;
        }
        let t0 = Instant::now();
        println!("inserting sequences from {}", record.id());
        let log_modulus = 1_000_000;
        let i = record.seq().windows(seqlen).filter(filter_n);
        for (count, elem) in i.enumerate() {
            let str = String::from_utf8(elem.to_vec())?;
            b[str.as_str()] += 1;
            if count % log_modulus == 0 {
                println!(
                    "inserted {} elements elapsed={:.02}s",
                    count,
                    t0.elapsed().as_secs_f64()
                );
            }
        }
        println!(
            "inserted all records for {} elapsed={:.2}s",
            record.id(),
            t0.elapsed().as_secs_f64(),
        );
    }
    let t0 = Instant::now();
    // println!("writing to disk {}!", outpath);
    // serialize(&b, outpath);
    println!("wrote to disk in {:.2}s", t0.elapsed().as_secs_f64());
    Ok(())
}
