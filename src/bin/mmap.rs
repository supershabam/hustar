use std::fs::{OpenOptions, File};
use std::path::PathBuf;

use anyhow::Result;
use memmap2::MmapMut;

struct Buf {
    mmap: Box<MmapMut>,
}

impl Buf {
    fn new() -> Result<Self> {
        let path: PathBuf = PathBuf::from(r"test.bin");
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .open(&path)?;
        file.set_len(4062)?;
        let mut mmap = unsafe { Box::new(MmapMut::map_mut(&file)?) };
        Ok(Buf{
            mmap: mmap,
        })
    }
}

fn main() -> Result<()> {
    let mut buf = Buf::new()?;
    buf.mmap[3] = 8;
    Ok(())
}
