use bio::alphabets;
use bio::io::fasta::Reader;
use croaring::Treemap;
use image::imageops::resize;
use image::ImageBuffer;
use roaring::RoaringTreemap;
use show_image::{create_window, ImageInfo, ImageView};
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
        println!("l={} s={:?} base={} addr={}",l, s, base, addr);
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

// #[show_image::main]
fn main() {
    let seqlen = 2;
    let mut b = make_buf(seqlen);
    inc(&mut b, b"gt");
    inc(&mut b, b"at");
    inc(&mut b, b"ga");
    println!("{}", get(&b, b"t"));
    println!("{:?}", b);
}

fn run() {
    let path = "files/ncbi-genomes-2022-02-23/GCA_000001405.29_GRCh38.p14_genomic.fna";
    let r = Reader::from_file(path).unwrap();
    let mut records = r.records();
    while let Some(Ok(record)) = records.next() {
        println!("inserting sequences from {}", record.id());
        if record.id() != "MU273354.1" {
            continue;
        }
        let i = record
            .seq()
            .windows(20)
            .filter(filter_n)
            .map(seq_to_bits);
        let mut count = 0_u64;
        let modulus = 1_000_000;
        let mut last = Instant::now();
        for (idx, elem) in i.enumerate() {
            count += 1;
            if idx % modulus == 0 {
                println!(
                    "processed total={} batch_size={} in {:.2}s",
                    count,
                    modulus,
                    last.elapsed().as_secs_f64()
                );
                last = Instant::now();
            }
        }
    }
}
