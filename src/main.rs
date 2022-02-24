use bio::alphabets;
use bio::io::fasta::Reader;
use croaring::Treemap;
use image::ImageBuffer;
use image::imageops::resize;
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
// #[show_image::main]
fn main() {
    let seq_width = 20;
    let mut bm = Treemap::create();

    let path = "files/ncbi-genomes-2022-02-23/GCA_000001405.29_GRCh38.p14_genomic.fna";
    let r = Reader::from_file(path).unwrap();
    let mut records = r.records();
    let mut img: [u8; 4096] = [0; 4096];
    // let window = create_window("image", Default::default()).unwrap();

    let buf = ImageBuffer::from_fn(4096, 4, |x, y| {
        image::Luma([0_u64])
    });
    let max = buf.iter().fold(0, |acc, p| {
        if *p > acc {
            *p
        } else {
            acc
        }
    });
    println!("{:?}", buf.dimensions());
    println!("max={}", max);
    let pixel = buf[(0, 0)];
    println!("{:?}", pixel);

    let buf2 = resize(&buf, 4096, 400, image::imageops::Lanczos3);
    println!("{:?}", buf2.dimensions());

    while let Some(Ok(record)) = records.next() {
        println!("inserting sequences from {}", record.id());
        let i = record
            .seq()
            .windows(seq_width)
            .filter(filter_n)
            .map(seq_to_bits);
        let mut last = Instant::now();
        for (idx, elem) in i.enumerate() {

            if idx % 100000 == 0 {
                let iidx: usize = (elem >> (40 - 12)) as usize; // 40 bits down to 12 (4096)
                img[iidx] = 1 + img[iidx];
                println!(
                    "idx={} elapsed={:.3}s; still inserting...",
                    idx,
                    last.elapsed().as_secs_f64()
                );
                // let image = ImageView::new(ImageInfo::rgb8(4096, 1), &img);
                // window.set_image("image-001", image).unwrap();
                last = Instant::now();
            }
        }
    }
}
