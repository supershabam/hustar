use anyhow::Result;
use core::ops::{Index, IndexMut};
use memmap2::{Mmap, MmapMut};
use std::path::PathBuf;

pub struct DatabaseMut {
    mmap: MmapMut,
}

pub struct Database {
    mmap: Mmap,
    path: PathBuf,
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

pub fn index_to_seq(index: usize) -> String {
  let base = 4_usize;
  let mut offset = 0;
  let mut seqlen = 1;
  while index > offset + base.pow(seqlen) - 1{
    offset = offset + base.pow(seqlen);
    seqlen = seqlen+1;
  }
  let addr = index - offset;
  bits_to_seq(addr as u64, seqlen as usize)
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
        Ok(DatabaseMut { mmap })
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

impl Database {
    pub fn open<P: Into<PathBuf>>(path: P) -> Result<Database> {
        use std::fs::OpenOptions;

        let path = path.into();
        let file = OpenOptions::new()
            .read(true)
            .write(false)
            .create(false)
            .open(&path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(Database { mmap, path })
    }
}

impl Index<&str> for Database {
    type Output = u64;
    fn index(&self, seq: &str) -> &Self::Output {
        let index = seq_to_index(seq);
        &self[index]
    }
}

impl Index<usize> for Database {
    type Output = u64;
    fn index(&self, index: usize) -> &Self::Output {
        let (_, buf, _) = unsafe { self.mmap.align_to::<u64>() };
        &buf[index]
    }
}

impl Clone for Database {
    fn clone(&self) -> Self {
      let next = Database::open(self.path.clone()).expect("while opening file for database clone");
      next
    }
}

#[cfg(test)]
mod test {
  use super::*;

  #[test]
  fn test_index_seq_index() {
    for index in 0..255 {
      let seq = index_to_seq(index);
      println!("index={index} seq={seq}");
      let indexc = seq_to_index(&seq);
      assert_eq!(index, indexc);
    }
  }
}