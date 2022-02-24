use bio::io::fasta::Reader;
use bio::alphabets;
use std::fs::File;
use std::io::prelude::*;
use std::str;
use roaring::RoaringTreemap;

fn filter_n(seq: &&[u8]) -> bool {
    for l in seq.iter() {
        let l = l.to_ascii_lowercase();
        match l {
            b'n' => return false,
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
            _ => panic!("invalid letter"),
        };
        b = b << 2;
        b = b | p;
    }
    b
}

fn main() {
    let seq_width = 20;
    let mut bm = RoaringTreemap::new();

    let path = "files/ncbi-genomes-2022-02-23/GCA_000001405.29_GRCh38.p14_genomic.fna";
    let r = Reader::from_file(path).unwrap();
    let mut records = r.records();

    while let Some(Ok(record)) = records.next() {
        println!("inserting sequences from {}", record.id());
        let i = record.seq().windows(seq_width).filter(filter_n).map(seq_to_bits);
        for (idx, elem) in i.enumerate() {
            bm.insert(elem);
            if idx%100000 == 0 {
                println!("idx={}; still inserting...", idx);
            }
        }
    }
}
