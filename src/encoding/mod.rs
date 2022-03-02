
// Example lazily consumes an iterator of values to act itself
// as an iterator.
// My intent is to use this pattern to flexibly digest the human genome
// which is 3.1e9 characters long.
struct Example<'a> {
    i: &'a mut dyn Iterator<Item = u64>,
    last: u64,
}

impl<'a> Example<'a> {
    fn new(i: &'a mut impl Iterator<Item = u64>) -> Example<'a> {
        Example { 
          i: i,
          last: 0,
        }
    }
}

impl<'a> Iterator for Example<'a> {
    type Item = u64;
    fn next(self: &mut Example<'a>) -> Option<Self::Item> {
      match self.i.next() {
        None => None,
        Some(i) => {
          self.last = i ^ self.last;
          Some(self.last)
        }
      }
    }
}

#[cfg(test)]
mod tests {
    use crate::encoding::*;

    #[test]
    fn test() {
      let v = vec![0, 1, 2, 3];
      let mut i = v.into_iter();
      let e = Example::new(&mut i);
      for v in e {
        println!("{:#016b}", v);
      }
    }
}

// use std::ops;

// use bitvec::prelude::*;

// struct Encoder<'a> {
//   vals: Box<&'a dyn Iterator<Item=u64>>,
//   bv: BitVec<u8, Msb0>,
// }

// impl<'a> Encoder<'a> {
//   fn new<I>(vals: &'a I) -> Encoder<'a>
//   where
//     I: IntoIterator<Item = u64>,
//   {
//     Encoder{
//       vals: Box::new(&vals.into_iter()),
//       bv: BitVec::default(),
//     }
//   }
// }

// impl<'a> Iterator for Encoder<'a> {
//   type Item = u8;

//   fn next(&mut self) -> Option<Self::Item> {
//     None
//   }
// }

// #[cfg(test)]
// mod tests {
//     use crate::encoding::*;

//     #[test]
//     fn test() {
//       let v = &vec![0_64];
//       let e = Encoder::new(v);
//       let v: Vec<u8> = e.into_iter().collect();
//       assert_eq!(v, vec![0_u8, 0_u8, 0_u8, 0_u8]);
//       let e = Encoder::new(&(0..90));

//       let mut bv = bitvec![u8, Msb0;];
//       bv.extend(&7u8.view_bits::<Msb0>()[..]);
//       bv.push(false);
//       bv.push(true);
//       assert_eq!(bv.as_raw_slice(), &[0b0000_0111, 0b0100_0000]);
//       assert_eq!(bv.len(), 10);
//     }
// }
