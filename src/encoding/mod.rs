use bitvec::prelude::*;

pub struct Encoder<'a> {
    i: &'a mut dyn Iterator<Item = u64>,
    buf: BitVec<u8, Msb0>,
    anchor: Option<i64>,
}

impl<'a> Encoder<'a> {
    pub fn new(i: &'a mut impl Iterator<Item = u64>) -> Self {
        Self {
            i: i,
            buf: BitVec::new(),
            anchor: None,
        }
    }
}

impl<'a> Iterator for Encoder<'a> {
    type Item = u8;
    fn next(self: &mut Encoder<'a>) -> Option<Self::Item> {
        loop {
            if self.buf.len() < 8 {
                let next = self.i.next();
                match next {
                    None => return None,
                    Some(c) => {
                        let delta = match self.anchor {
                            None => {
                                c as i64
                            }
                            Some(anchor) => c as i64 - anchor,
                        };
                        self.anchor = Some(delta);
                        let delta = vint64::signed::encode(delta);
                        self.buf.extend(delta.as_ref());
                        println!("extended {} bytes", delta.as_ref().len());
                    }
                }
            }
            let head = self.buf.drain(0..8).collect::<BitVec<u8, Msb0>>();
            let head = head.load_be::<u8>();
            return Some(head)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::encoding::*;

    #[test]
    fn test() {
        let v = vec![0xff_00_fe, 0x00];
        let mut i = v.into_iter();
        let e = Encoder::new(&mut i);
        for v in e {
            println!("{:064b}", 0xfe_u64);
            println!("{:08b}", v);
        }
    }

    #[test]
    fn test_vint() {
        let original = 0x7fff_ffff_ffff_ffff_i64;
        let delta = vint64::signed::encode(original);
        let mut bv: BitVec<u8, Msb0> = BitVec::new();
        bv.extend(delta.as_ref());
        let len = vint64::decoded_len(bv[0..8].load_be::<u8>());
        assert_eq!(delta.as_ref().len(), len);
        let buf: BitVec<u8, Msb0> = bv.drain(0..len * 8).collect();
        let mut buf = buf.as_raw_slice();
        let decoded = vint64::signed::decode(&mut buf).expect("while decoding");
        assert_eq!(decoded, original);
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
