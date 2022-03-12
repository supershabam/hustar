use std::fs::{File, OpenOptions};
use core::ops::{Index,IndexMut};
use std::path::{Path, PathBuf};

use anyhow::Result;
use memmap2::MmapMut;

pub struct Mmap {
    mmap: MmapMut,
}

impl Mmap {
    pub fn new<P: Into<PathBuf>>(path: P, size: u64) -> Result<Mmap> {
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .open(&path.into())?;
        file.set_len(size)?;
        let mmap = unsafe { MmapMut::map_mut(&file)? };
        Ok(Mmap { 
            mmap: mmap,
        })
    }
}

impl Index<usize> for Mmap {
    type Output = u64;
    fn index(&self, index: usize) -> &Self::Output {
        let (_, buf, _) = unsafe { self.mmap.align_to::<u64>() };
        &buf[index]
    }
}

impl IndexMut<usize> for Mmap {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        let (_, buf, _) = unsafe { self.mmap.align_to_mut::<u64>() };
        &mut buf[index]
    }
}
