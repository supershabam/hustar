use std::fs::OpenOptions;
use std::path::PathBuf;

use anyhow::Result;
use memmap2::MmapMut;

struct Buf<'a> {
    mmap: &'a MmapMut,
    buf: &'a mut [u64],
}

impl Buf<'a> {
    fn new<'a>() -> Result<Self> {
        let path: PathBuf = PathBuf::from(r"test.bin");
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .open(&path)?;
        file.set_len(4062)?;
        let mut mmap = unsafe { MmapMut::map_mut(&file)? };
        let (_, mut buf, _) = unsafe { mmap.align_to_mut::<u64>() };
        Ok(Buf<'a>{
            mmap: &mmap,
            buf: buf,
        })
    }
}

fn main() -> Result<()> {
    let buf = Buf::new()?;
    Ok(())
}
