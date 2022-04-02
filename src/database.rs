use anyhow::Result;
use core::ops::{Index, IndexMut};
use memmap2::MmapMut;
use std::path::PathBuf;

pub struct DatabaseMut {
    mmap: MmapMut,
    seqlen: usize,
}

fn buf_size_bytes(seqlen: usize) -> u64 {
    let base = 4_u64;
    let mut size = 0;
    for l in 1..=seqlen {
        size += base.pow(l as u32);
    }
    size * 8
}

fn seq_to_index(seq: &str) -> usize {
    let unit: u64 = 4;
    let mut base: usize = 0;
    for l in 1..seq.len() {
        base += unit.pow(l as u32) as usize;
    }
    let addr = seq_to_addr(seq);
    let idx = base + addr;
    idx
}

fn seq_to_addr(seq: &str) -> usize {
    let mut b: u64 = 0;
    for l in seq.chars() {
        let l = l.to_ascii_lowercase();
        let p: u64 = match l {
            'a' => 0b00,
            'c' => 0b01,
            'g' => 0b10,
            't' => 0b11,
            _ => panic!("invalid letter {}", l as char),
        };
        b = b << 2;
        b = b | p;
    }
    b as usize
}

impl DatabaseMut {
    pub fn create<P: Into<PathBuf>>(path: P, seqlen: usize) -> Result<DatabaseMut> {
        use std::fs::OpenOptions;

        let size = buf_size_bytes(seqlen);
        let path = path.into();
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .open(&path)?;
        file.set_len(size)?;
        let mmap = unsafe { MmapMut::map_mut(&file)? };
        Ok(DatabaseMut { mmap: mmap, seqlen })
    }
}

impl Index<&str> for DatabaseMut {
    type Output = u64;
    fn index(&self, seq: &str) -> &Self::Output {
        let index = seq_to_index(seq);
        &self[index]
    }
}

impl IndexMut<&str> for DatabaseMut {
    fn index_mut(&mut self, seq: &str) -> &mut Self::Output {
        let index = seq_to_index(seq);
        &mut self[index]
    }
}

impl Index<usize> for DatabaseMut {
    type Output = u64;
    fn index(&self, index: usize) -> &Self::Output {
        let (_, buf, _) = unsafe { self.mmap.align_to::<u64>() };
        &buf[index]
    }
}

impl IndexMut<usize> for DatabaseMut {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        let (_, buf, _) = unsafe { self.mmap.align_to_mut::<u64>() };
        &mut buf[index]
    }
}
