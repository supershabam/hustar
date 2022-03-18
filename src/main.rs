mod encoding;
mod window;

mod mmap;
mod traverse;

use mmap::Mmap;

use serde::{Deserialize, Serialize};

use anyhow::Result;
use binwrite::BinWrite;
use bio::alphabets;
use bio::io::fasta::Reader;
use image::imageops::resize;
use image::ImageBuffer;
use show_image::{create_window, ImageInfo, ImageView};
use std::f64::consts::PI;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::path::PathBuf;
use std::rc::Rc;
use std::time::Duration;
use std::time::Instant;
use std::{str, thread, time};

use byteorder::BigEndian;
use zerocopy::byteorder::U64;
use zerocopy::LayoutVerified;

use clap::{Parser, Subcommand};
use std::sync::mpsc::channel;

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

fn seq_to_bits(seq: &[u8]) -> u64 {
    let mut b: u64 = 0;
    for l in seq.iter() {
        let l = l.to_ascii_lowercase();
        let p: u64 = match l {
            b'a' => 0b00,
            b'c' => 0b01,
            b'g' => 0b10,
            b't' => 0b11,
            _ => panic!("invalid letter {}", l as char),
        };
        b = b << 2;
        b = b | p;
    }
    b
}

fn bits_to_seq(bits: u64, bitsize: usize) -> String {
    let mut bits = bits;
    let mut seq = "".to_string();
    for _ in 0..bitsize {
        let tail = bits & 0b11;
        bits = bits >> 2;
        let s = match tail {
            0b00 => "a",
            0b01 => "c",
            0b10 => "g",
            0b11 => "t",
            _ => panic!("impossible tail"),
        };
        seq.push_str(s);
    }
    seq.chars().rev().collect::<String>()
}

/*
need to create a custom data structure for storing prefix hits
[
    b00: u64,
    b01: u64,
    b10: u64,
    b11: u64,
    b0000: u64,
    b0001: u64,
    b0010: u64,
    ...
    b1111: u64,
    ...
]

grows by 2 bits for each tier

len = 4^seqlen + 4^(seqlen-1) ...

can be modelled simply with a vector of u64 with custom indexing

get(seq: &[u8]) -> u64
  let base: u64 = 4^(len(seq)-1);
  let addr = seq_to_addr(seq);
  let idx = base + addr; // seq_to_idx
  buf[idx]

inc(seq: &[u8])
  for let l = (1..len(seq)) {
      let substr = seq.substr(l);
      let idx = seq_to_idx(substr);
      buf[idx] += 1;
  }
*/

fn make_buf<'a, P: Into<PathBuf>>(path: P, seqlen: usize) -> Result<Mmap> {
    let size = buf_size_bytes(seqlen);
    let m = Mmap::new(path, size)?;
    Ok(m)
}

fn buf_size_bytes(seqlen: usize) -> u64 {
    let base = 4_u64;
    let mut size = 0;
    for l in 1..=seqlen {
        size += base.pow(l as u32);
    }
    size * 8
}

fn inc(m: &mut Mmap, seq: &[u8]) {
    let unit: u64 = 4;
    let mut base: usize = 0;
    for l in 1..=seq.len() {
        let s = &seq[..l];
        let addr = seq_to_bits(s) as usize;
        let idx = base + addr;
        m[idx] += 1;
        base += unit.pow(l as u32) as usize;
    }
}

fn seq_to_idx(seq: &[u8]) -> usize {
    let unit: u64 = 4;
    let mut base: usize = 0;
    for l in 1..seq.len() {
        base += unit.pow(l as u32) as usize;
    }
    let addr = seq_to_bits(seq) as usize;
    let idx = base + addr;
    idx
}

fn get(m: &Mmap, seq: &[u8]) -> u64 {
    let unit: u64 = 4;
    let mut base: usize = 0;
    for l in 1..seq.len() {
        base += unit.pow(l as u32) as usize;
    }
    let addr = seq_to_bits(seq) as usize;
    let idx = base + addr;
    m[idx]
}

fn seqlen_to_bitmap_range(seqlen: usize) -> (usize, usize) {
    let mut start = 0;
    let mut range = 4;
    for _ in 1..seqlen {
        start += range;
        range = range * 4;
    }
    (start, start + range)
}
#[cfg(test)]
mod tests {
    use crate::*;

    #[test]
    fn test_seq_to_bits() {
        let seq = "actg".to_string();
        let bits: u64 = 0b00011110;
        assert_eq!(seq_to_bits(seq.as_bytes()), bits);
        assert_eq!(bits_to_seq(bits, 4), seq);
    }

    #[test]
    fn test_seq_to_angle() {
        assert_eq!(seq_to_angle(b"aaa"), 0.0);
        assert_eq!(seq_to_angle(b"aac"), 2.0 * PI / 4.0 / 4.0 / 4.0);
        assert_eq!(seq_to_angle(b"aag"), 2.0 * 2.0 * PI / 4.0 / 4.0 / 4.0);
        assert_eq!(seq_to_angle(b"aat"), 3.0 * 2.0 * PI / 4.0 / 4.0 / 4.0);
        assert_eq!(seq_to_angle(b"aca"), 4.0 * 2.0 * PI / 4.0 / 4.0 / 4.0);
    }

    #[test]
    fn test_coords() {
        let f = make_coords_to_seq(3.0);
        assert_eq!(f(seq_to_angle(b"aaa"), 6.3), "aaa".to_string());
        assert_eq!(f(seq_to_angle(b"aac"), 6.3), "aac".to_string());
        assert_eq!(f(seq_to_angle(b"cac"), 6.3), "cac".to_string());
        assert_eq!(f(seq_to_angle(b"tga"), 6.3), "tga".to_string());
    }

    #[test]
    fn test_seqlen_range() {
        assert_eq!(seqlen_to_bitmap_range(0), (0, 4));
        assert_eq!(seqlen_to_bitmap_range(1), (0, 4));
        assert_eq!(seqlen_to_bitmap_range(2), (4, 20));
        assert_eq!(seqlen_to_bitmap_range(3), (20, 84));
    }

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

    #[test]
    fn test_make_coords_to_seq_range() {
        let width = 200;
        let height = 200;
        let max_seqlen = 3;
        let x = 190;
        let y = 100;
        let f = make_coords_to_seq_range(width, height, max_seqlen);
        let (start, end) = f(x, y);

        assert_eq!(start, "aaa");
        assert_eq!(end, "aaa");

        let width = 200;
        let height = 200;
        let max_seqlen = 20;
        let x = 200;
        let y = 100;
        let f = make_coords_to_seq_range(width, height, max_seqlen);
        let (start, end) = f(x, y);

        assert_eq!(start, "aaaaaaaaaaaaaaaaaaaa");
        assert_eq!(end, "aaaacggacatatgaatggg");
    }
}

fn seq_to_angle(seq: &[u8]) -> f64 {
    let mut angle = 0.0;
    let mut step = 2.0 * PI / 4.0;
    for c in seq {
        let v = seq_to_bits(&[*c]) as f64;
        angle = angle + v * step;
        step = step / 4.0;
    }
    angle
}

fn make_coords_to_seq(r_step: f64) -> Box<dyn Fn(f64, f64) -> String> {
    Box::new(move |theta: f64, r: f64| -> String {
        let mut slice = 2.0 * PI / 4.0;
        let mut bits: u64 = 0;
        let mut bitsize: usize = 0;
        let n = (r / r_step).floor() as usize + 1;
        for _ in 0..n {
            let p = theta / slice;
            let tail = (p.floor() as u64) % 4;
            bits = bits << 2;
            bits = bits | tail;
            bitsize += 1;
            slice = slice / 4.0;
        }
        let seq = bits_to_seq(bits, bitsize);
        seq
    })
}

fn sequence_from_angle(theta: f64, seqlen: usize) -> String {
    let mut slice = 2.0 * PI / 4.0; // current width
    let mut bits: u64 = 0;
    let mut bitsize: usize = 0;
    for _ in 0..seqlen {
        let f = theta / slice;
        let tail = (f.floor() as u64) % 4;
        bits = bits << 2;
        bits = bits | tail;
        bitsize += 1;
        slice = slice / 4.0;
    }
    bits_to_seq(bits, bitsize)
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

fn make_coords_to_seq_range(
    width: usize,
    height: usize,
    max_seqlen: usize,
) -> Box<dyn Fn(u32, u32) -> (String, String)> {
    Box::new(move |x, y| -> (String, String) {
        // calculate all angles for the square and choose the min and max
        let (min_theta, max_theta) = min_max_thetas(width, height, x, y);
        let seqlen = seqlen_from_position(width, height, x, y, max_seqlen);
        (
            sequence_from_angle(min_theta, seqlen),
            sequence_from_angle(max_theta, seqlen),
        )
    })
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
    let cpus = num_cpus::get();
    let thread_count = cpus;

    println!("printing {} with thread_count={}", index_file, thread_count);
    // let size = buf_size_bytes(seqlen);
    let m = Mmap::open(index_file)?;

    let width = side_length;
    let height = side_length;

    let (tx, rx) = channel();
    for worker_id in 0..thread_count {
        let tx = tx.clone();
        let m = m.clone();
        thread::spawn(move || {
            let seqrange = make_coords_to_seq_range(width as usize, height as usize, seqlen);
            for y in 0..height {
                for x in 0..width {
                    if x % thread_count != worker_id {
                        continue;
                    }
                    let testx: i32 = x as i32 - (width / 2) as i32;
                    let testy: i32 = y as i32 - (height / 2) as i32;
                    if (-3 < testx && testx < 3) || (-3 < testy && testy < 3) {
                        // println!("{} {}", min, max);
                        // TODO figure out why bright lines are forming on axis
                        continue;
                    }
                    let (min, max) = seqrange(x as u32, y as u32);

                    let len = min.len();
                    let c = m.count(seq_to_idx(min.as_bytes()), seq_to_idx(max.as_bytes()));
                    tx.send((x, y, c, len)).unwrap();
                }
            }
            println!("worker_id={} exiting", worker_id)
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
            println!("processed {} pixels of {} ({:.2}%)", counter, total, percentage);
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
    let mut b = make_buf(outpath, seqlen)?;
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
            inc(&mut b, elem);
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
