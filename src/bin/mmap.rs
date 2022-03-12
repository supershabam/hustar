use std::fs::{File, OpenOptions};
use core::ops::{Index,IndexMut};
use std::io::Write;
use std::path::PathBuf;

use anyhow::Result;
use memmap2::MmapMut;

struct Buf {
    mmap: MmapMut,
}

impl Buf {
    fn new() -> Result<Buf> {
        let path: PathBuf = PathBuf::from(r"test.bin");
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .open(&path)?;
        file.set_len(4062)?;
        let mmap = unsafe { MmapMut::map_mut(&file)? };
        Ok(Buf { 
            mmap: mmap,
        })
    }
}

impl Index<usize> for Buf {
    type Output = u64;
    fn index(&self, index: usize) -> &Self::Output {
        let (_, buf, _) = unsafe { self.mmap.align_to::<u64>() };
        &buf[index]
    }
}

impl IndexMut<usize> for Buf {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        let (_, buf, _) = unsafe { self.mmap.align_to_mut::<u64>() };
        &mut buf[index]
    }
}

fn main() -> Result<()> {
    let mut buf = Buf::new()?;
    buf[5] = 0xdeadbeef;
    println!("{}", buf[5]);
    Ok(())
}
