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
// #[show_image::main]
fn main() {
    let seq_width = 20;
    let path = "files/ncbi-genomes-2022-02-23/GCA_000001405.29_GRCh38.p14_genomic.fna";
    let r = Reader::from_file(path).unwrap();
    let mut records = r.records();
    while let Some(Ok(record)) = records.next() {
        println!("inserting sequences from {}", record.id());
        let i = record
            .seq()
            .windows(seq_width)
            .filter(filter_n)
            .map(seq_to_bits);
        let mut count = 0_u64;
        let modulus = 1_000_000;
        let mut last = Instant::now();
        for (idx, elem) in i.enumerate() {
            count += 1;
            if idx %  modulus== 0 {
                println!("processed total={} batch_size={} in {:.2}s", count, modulus, last.elapsed().as_secs_f64());
                last = Instant::now();
            }
        }
    }
}
