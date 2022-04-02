use std::cmp::Ordering;
use anyhow::Result;

use crate::mmap::Mmap;

#[derive(Default)]
pub struct Accumulator {
  gte: usize,
  lt: usize,
  sum: i64,
}

impl Accumulator {
  // sum_to mutates the accumulator internals to the provided range and returns the
  // sum under the range.
  //
  // TODO: it is possible to do less work. Say for the case where the acc is at [0,5)
  // and then is requested to move to [9,25).
  pub fn sum_to(&mut self, mmap: &Mmap, gte: usize, lt: usize) -> u64 {
    loop {
      match self.lt.cmp(&lt) {
        Ordering::Less => {
          let c = mmap[self.lt];
          self.sum = self.sum + c as i64;
          self.lt = self.lt + 1;
        },
        Ordering::Greater => {
          let c = mmap[self.lt];
          self.sum = self.sum - c as i64;
          self.lt = self.lt - 1;
        },
        Ordering::Equal => {
          break
        }
      }
    }
    loop {
      match self.gte.cmp(&gte) {
        Ordering::Less => {
          let c = mmap[self.gte];
          self.sum = self.sum - c as i64;
          self.gte = self.gte + 1;
        },
        Ordering::Greater => {
          let c = mmap[self.gte];
          self.sum = self.sum + c as i64;
          self.gte = self.gte - 1;
        },
        Ordering::Equal => {
          break
        }
      }
    }
    self.sum as u64
  }
}

#[cfg(test)]
mod test {
  use anyhow::Result;

  #[test]
  fn test_accumulator() -> Result<()> {
    use crate::accumulator::*;
    let mut a = Accumulator::default();
    let mmap = Mmap::open("seqlen=4.bin")?;
    let c = a.sum_to(&mmap, 2, 3);
    println!("count = {c}");
    let c = a.sum_to(&mmap, 1, 2);
    println!("count = {c}");
    Ok(())
  }
}