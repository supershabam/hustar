use core::ops::{Index, IndexMut};
use std::fs::{File, OpenOptions};
use std::iter::Cloned;
use std::path::{Path, PathBuf};

use anyhow::Result;
use memmap2::MmapMut;

pub struct Mmap {
    mmap: MmapMut,
    path: String,
}

impl Mmap {
    pub fn new<P: Into<PathBuf>>(path: P, size: u64) -> Result<Mmap> {
        let path = path.into();
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .open(&path)?;
        file.set_len(size)?;
        let mmap = unsafe { MmapMut::map_mut(&file)? };
        Ok(Mmap { 
            mmap: mmap,
            path: path.into_os_string().into_string().expect("creating path"),
         })
    }

    pub fn open<P: Into<PathBuf>>(path: P) -> Result<Mmap> {
        let path = path.into();
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(false)
            .open(&path)?;
        let mmap = unsafe { MmapMut::map_mut(&file)? };
        Ok(Mmap {
            mmap: mmap,
            path: path.into_os_string().into_string().expect("creating path"),
         })
    }

    pub fn count(&self, min_inclusive: usize, max_inclusive: usize) -> u64 {
        let mut c = 0;
        for idx in min_inclusive..=max_inclusive {
            c = c + self[idx];
        }
        c
    }
}

impl Clone for Mmap {
    fn clone(&self) -> Self {
        let next = Self::open(self.path.clone()).expect("oops");
        next
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
