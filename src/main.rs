mod encoding;
mod window;

use serde::{Deserialize, Serialize};

use binwrite::BinWrite;
use bio::alphabets;
use bio::io::fasta::Reader;
use image::imageops::resize;
use image::ImageBuffer;
use show_image::{create_window, ImageInfo, ImageView};
use std::f64::consts::PI;
use std::fs::File;
use std::io::prelude::*;
use std::path::PathBuf;
use std::fs::OpenOptions;
use std::rc::Rc;
use std::time::Instant;
use std::{str, time};
use anyhow::Result;
use memmap2::MmapMut;

use byteorder::BigEndian;
use zerocopy::byteorder::U64;
use zerocopy::LayoutVerified;

use clap::{Parser, Subcommand};

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

fn make_buf<'a>(path: &PathBuf, seqlen: usize) -> Result<MmapMut> {
    let base = 4_u64;
    let mut size = 0;
    for l in 1..=seqlen {
        size += base.pow(l as u32);
    }
    let mut file = OpenOptions::new()
                       .read(true)
                       .write(true)
                       .create(true)
                       .open(&path)?;
    file.set_len(size)?;
    let mut mmap = unsafe { MmapMut::map_mut(&file)? };
    // let (_, mut umap, _) = unsafe { mmap.align_to_mut::<u64>() };
    // Ok(Box::new(&mmap[..]))
    Ok(mmap)
}

fn inc(mmap: &mut MmapMut, seq: &[u8]) {
    let (_, mut b, _) = unsafe { mmap.align_to_mut::<u64>() };
    let unit: u64 = 4;
    let mut base: usize = 0;
    for l in 1..=seq.len() {
        let s = &seq[..l];
        let addr = seq_to_bits(s) as usize;
        let idx = base + addr;
        b[idx] += 1;
        base += unit.pow(l as u32) as usize;
    }
}

fn get(b: &[u64], seq: &[u8]) -> u64 {
    let unit: u64 = 4;
    let mut base: usize = 0;
    for l in 1..seq.len() {
        base += unit.pow(l as u32) as usize;
    }
    let addr = seq_to_bits(seq) as usize;
    let idx = base + addr;
    b[idx]
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

fn serialize(b: &[u64], path: &str) {
    let buf = bincode::serialize(b).expect("while serializing");
    let mut f = File::create(path).expect("while opening file");
    buf.write(&mut f).expect("while writing file");
}

fn deserialize(path: &str) -> Vec<u64> {
    let mut f = File::open(path).expect("while opening file");
    let mut encoded = Vec::new();
    f.read_to_end(&mut encoded).expect("while reading");
    let decoded: Vec<u64> = bincode::deserialize(&encoded[..]).unwrap();
    decoded
}

#[derive(Deserialize, Serialize, Debug, PartialEq, BinWrite)]
struct Bitmap {
    V: Vec<u64>,
    Seqlen: u32,
}

#[cfg(test)]
mod tests {
    use crate::*;
    #[test]
    fn test_serialization() {
        let seqlen = 3;
        let mut b = make_buf(&PathBuf::from(r"testing.bin"), seqlen).expect("while creating buf");
        inc(&mut b, b"acc");
        inc(&mut b, b"acg");
        // serialize(&b, "test.bin");
        // let b2 = deserialize("test.bin");
        // assert_eq!(b, b2);
    }

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

fn main() {
    let args = Cli::parse();
    match &args.command {
        Commands::Build { fasta_file, index_file, sequence_length} => {
            create(fasta_file, index_file, *sequence_length);
        },
        Commands::Visualize { index_file, sequence_length } => {
            print(index_file, *sequence_length);
        },
        _ => panic!("oh no"),
    }
}

fn read(index_file: &str, seqlen: usize, path: &str) {
    println!("opening {}", path);
    let bm = deserialize(path);
    let (start, end) = seqlen_to_bitmap_range(seqlen);
    let mut i = bm[start..end].iter().cloned();
    let e = encoding::Encoder::new(&mut i);
    let mut count = 0;
    for (idx, v) in e.enumerate() {
        println!("{:08b}", v);
        count = idx + 1;
    }
    println!("{} bytes", count);
}

fn print(index_file: &str, seqlen: usize) {
    println!("opening {}", index_file);
    let bm = deserialize(index_file);
    println!("creating image");
    let maxes: Vec<u64> = (1..=seqlen)
        .into_iter()
        .map(|l| {
            let (start, end) = seqlen_to_bitmap_range(l);
            bm[start..end]
                .iter()
                .fold(0, |acc, v| if acc > *v { acc } else { *v })
        })
        .collect();
    let width = 4000;
    let height = 4000;
    let circle_r = width as f64 / (1150.0);
    let coords_to_seq = make_coords_to_seq(circle_r);
    let mut img = ImageBuffer::from_fn(width, height, |px, py| {
        let x: f64 = (px as i32 - (width / 2) as i32) as f64;
        let y: f64 = (-(py as i32) + (height / 2) as i32) as f64;
        let theta = {
            let mut theta = (y / x).atan();
            if x == 0.0 {
                theta = 0.0;
            }
            if x < 0.0 {
                theta += PI;
            }
            if theta < 0.0 {
                theta = theta + PI + PI;
            }
            theta = theta + PI / 7.23;
            if theta > 2.0 * PI {
                theta = theta - 2.0 * PI;
            }
            theta
        };
        let r = (x * x + y * y).sqrt();
        
        let r = r.sqrt(); // make tiers closer to the origin smaller than extremities
        let r = r + theta / (2.0*PI); // add spiral between tiers
        let steps = (r / circle_r).ceil() as usize;
        let p = {
            if steps == 0 {
                0
            } else if steps > seqlen {
                0
            } else {
                let seq = coords_to_seq(theta, r);
                let v = get(&bm, seq.as_bytes());
                let normalized = v as f64 / maxes[seq.len() - 1] as f64;
                (normalized.sqrt().sqrt() * 255.0) as u8
            }
        };
        // println!(
        //     "pixel ({}, {}) xy ({}, {}) r={} theta={} seq={} p={}",
        //     px, py, x, y, r, theta, seq, p
        // );

        image::Luma([p])
    });
    img.save_with_format("out.png", image::ImageFormat::Png)
        .expect("while writing image");
}

fn create(fasta_file: &str, outpath: &str, seqlen: usize) -> Result<()> {
    let p = PathBuf::from(r"what.bin");
    let b = make_buf(&p, seqlen)?;
    let r = Reader::from_file(fasta_file).unwrap();
    let mut records = r.records();
    while let Some(Ok(record)) = records.next() {
        if record.id() != "CM000667.2" {
            continue
        }
        let t0 = Instant::now();
        println!("inserting sequences from {}", record.id());
        let log_modulus = 1_000_000;
        let i = record.seq().windows(seqlen).filter(filter_n);
        for (count, elem) in i.enumerate() {
            // inc(&mut b, elem);
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
    println!("writing to disk {}!", outpath);
    // serialize(&b, outpath);
    println!("wrote to disk in {:.2}s", t0.elapsed().as_secs_f64());
    Ok(())
}
