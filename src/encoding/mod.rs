use std::ops;

use bitvec::prelude::*;

struct Encoder<'a> {
  vals: Box<&'a dyn Iterator<Item=u64>>,
  bv: BitVec<u8, Msb0>,
}

impl<'a> Encoder<'a> {
  fn new<I>(vals: &'a I) -> Encoder<'a>
  where
    I: IntoIterator<Item = u64>,
  {
    Encoder{
      vals: Box::new(&vals.into_iter()),
      bv: BitVec::default(),
    }
  }
}

impl<'a> Iterator for Encoder<'a> {
  type Item = u8;

  fn next(&mut self) -> Option<Self::Item> {
    None
  }
}



#[cfg(test)]
mod tests {
    use crate::encoding::*;

    #[test]
    fn test() {
      let e = Encoder::new(vec![0_u64]);
      let v: Vec<u8> = e.into_iter().collect();
      assert_eq!(v, vec![0_u8, 0_u8, 0_u8, 0_u8]);
      let e = Encoder::new(0..90);

      let mut bv = bitvec![u8, Msb0;];
      bv.extend(&7u8.view_bits::<Msb0>()[..]);
      bv.push(false);
      bv.push(true);
      assert_eq!(bv.as_raw_slice(), &[0b0000_0111, 0b0100_0000]);
      assert_eq!(bv.len(), 10);
    }
}
