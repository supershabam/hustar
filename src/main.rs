mod accumulator;
mod database;
mod traverse;

use anyhow::Result;
use bio::io::fasta::Reader;
use image::ImageBuffer;

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
    use crate::traverse::filter_points;
    use crate::traverse::Point;
    use crossbeam::channel::bounded;

    let cpus = num_cpus::get();
    let thread_count = cpus;

    println!("printing {} with thread_count={}", index_file, thread_count);
    let m = Database::open(index_file)?;

    let width = side_length;
    let height = side_length;

    let pixels = make_points(width as u32, height as u32, seqlen as u32);
    let (work_tx, work_rx) = bounded(thread_count);
    thread::spawn(move || {
        let num_chunks = thread_count * 8;
        for chunk_id in 0..num_chunks {
            let pixels: Vec<Point> = pixels
                .iter()
                .filter(|&p| filter_points(num_chunks, chunk_id, p))
                .cloned()
                .collect();
            // let len = pixels.len();
            work_tx.send(pixels).expect("while sending chunk of pixels");
            // println!("sent chunk_id={chunk_id} of num_chunks={num_chunks} with len={}", len);
        }
    });
    let (tx, rx) = channel();
    for _ in 0..thread_count {
        let tx = tx.clone();
        let m = m.clone();
        let work_rx = work_rx.clone();
        thread::spawn(move || {
            for pixels in work_rx {
                let mut acc = Accumulator::default();
                for p in pixels {
                    let (gte, lt) = p.index_range();
                    let c = acc.sum_to(&m, gte, lt);
                    let val = (p.w as usize, p.h as usize, (c, lt - gte), p.seqlen);
                    tx.send(val).unwrap();
                }
            }
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
            maxesbuf[y * height + x] = maxes[seqlen];
        }
    }
    println!("maxes {:?}", maxes);
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
