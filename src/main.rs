use serde::{Deserialize, Serialize};

use binwrite::BinWrite;
use bio::alphabets;
use bio::io::fasta::Reader;
use croaring::Treemap;
use image::imageops::resize;
use image::ImageBuffer;
use roaring::RoaringTreemap;
use show_image::{create_window, ImageInfo, ImageView};
use std::f64::consts::PI;
use std::fs::File;
use std::io::prelude::*;
use std::time::Instant;
use std::{str, time};

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

fn make_buf(seqlen: usize) -> Vec<u64> {
    let base = 4_u64;
    let mut size = 0;
    for l in 1..=seqlen {
        size += base.pow(l as u32);
    }
    vec![0_u64; size as usize]
}

fn inc(b: &mut Vec<u64>, seq: &[u8]) {
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
        let mut b = make_buf(seqlen);
        inc(&mut b, b"acc");
        inc(&mut b, b"acg");
        serialize(&b, "test.bin");
        let b2 = deserialize("test.bin");
        assert_eq!(b, b2);
    }

    #[test]
    fn test_seq_to_bits() {
        let seq = "actg".to_string();
        let bits: u64 = 0b00011110;
        assert_eq!(seq_to_bits(seq.as_bytes()), bits);
        assert_eq!(bits_to_seq(bits, 4), seq);
    }
}

fn main() {
    let path = "out.bin";
    let bm = deserialize(path);
    let max = bm.iter().fold(0, |acc, v| if acc > *v { acc } else { *v });
    let width = 25;
    let height = 25;
    let circle_r = width as f64 / 6.0;
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
            theta
        };
        let r = (x * x + y * y).sqrt();
        let slice = PI / 2.0; // slice size decreases over time
        let p = theta / slice;
        let bit = p.floor() as usize;

        let v = {
            if r > circle_r {
                0
            } else {
                bm[bit]
            }
        };
        let normalized = v as f64 / max as f64;
        let p = (normalized * 255.0) as u8;
        println!(
            "pixel ({}, {}) xy ({}, {}) r={} theta={} bit={} v={} p={}",
            px, py, x, y, r, theta, bit, v, p
        );

        image::Luma([(normalized * 255.0) as u8])
    });
    img.save_with_format("out.png", image::ImageFormat::Png)
        .expect("while writing image");
}

fn create() {
    let seqlen = 3;
    let mut b = make_buf(seqlen);
    let path = "files/ncbi-genomes-2022-02-23/GCA_000001405.29_GRCh38.p14_genomic.fna";
    let r = Reader::from_file(path).unwrap();
    let mut records = r.records();
    while let Some(Ok(record)) = records.next() {
        println!("inserting sequences from {}", record.id());
        if record.id() != "CM000663.2" {
            continue;
        }
        let i = record.seq().windows(seqlen).filter(filter_n);
        let modulus = 1_000_00;
        let mut last = Instant::now();
        for elem in i {
            // println!("inserting {}", str::from_utf8(elem).expect("while converting string"));
            inc(&mut b, elem);
        }
        println!("writing to disk!");
        let path = "out.bin";
        serialize(&b, path);
        println!("wrote {}", path);
    }
}
