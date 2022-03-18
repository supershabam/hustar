// traverse helps to iterate the mmap in a minimally querying way
// to populate a desired set of pixel points.

use core::f64::consts::PI;
use itertools::Itertools;

struct Point {
    x: i32,
    y: i32,
    seqlen: usize,
}

fn delta_theta(t1: f64, t2: f64) -> (f64, f64) {
  // assert t1 <= 2PI && t2 <= 2PI
  let mut counterclockwise = t2 - t1;
  if counterclockwise < 0.0 {
    counterclockwise += 2.0 * PI;
  }
  let mut clockwise = t1 - t2;
  if clockwise < 0.0 {
    clockwise += 2.0 * PI;
  }
  // return minimum angle
  if clockwise < counterclockwise {
    (t2, clockwise)
  } else {
    (t1, counterclockwise)
  }
}

impl Point {
    fn thetas(&self) -> (f64, f64) {
        let (x, y) = (self.x, self.y);
        let mut thetas: Vec<f64> = vec![(x, y), (x + 1, y), (x, y + 1), (x + 1, y + 1)]
            .iter()
            .filter_map(|(x, y)| {
                if *x == 0 && *y == 0 {
                    return None;
                }
                let theta = {
                    if *x == 0 && *y > 0 {
                        PI / 2.0
                    } else if *x == 0 && *y < 0 {
                        -PI / 2.0
                    } else if *x < 0 {
                        (*y as f64 / *x as f64).atan() + PI
                    } else {
                        (*y as f64 / *x as f64).atan()
                    }
                };
                let theta = (2.0 * PI + theta) % (2.0 * PI);
                Some(theta)
            })
            .collect();
        thetas.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let td = thetas.iter().permutations(2).fold((0.0, 0.0), |acc, t| {
          match t[..] {
            [&a, &b] => {
              let dt = delta_theta(a, b);
              // select maximum distance
              if dt.1 > acc.1 {
                dt
              } else {
                acc
              }
            },
            _ => panic!("permutations produced bad combination"),
          }
        });
        (td.0, td.0+td.1)
    }

    fn index_range(&self) -> (usize, usize) {
        let (t1, t2) = self.thetas();
        (self.index(t1), self.index(t2))
    }

    fn index(&self, theta: f64) -> usize {
        let base: usize = (1..self.seqlen).fold(0, |acc, i| acc + (1 << (2 * i)));
        let max: usize = 1 << self.seqlen * 2;
        let percentage = theta / (2.0 * PI);
        let addr = max as f64 * percentage;
        let addr = addr as usize;
        base + addr
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_indexes() {
        let p = Point {
            x: 0,
            y: 0,
            seqlen: 1,
        };
        assert_eq!((0, 1), p.index_range());

        let p = Point {
            x: 0,
            y: 0,
            seqlen: 2,
        };
        assert_eq!((4, 8), p.index_range());

        let p = Point {
            x: -1,
            y: 0,
            seqlen: 2,
        };
        assert_eq!((8, 12), p.index_range());

        let p = Point {
            x: 0,
            y: 0,
            seqlen: 3,
        };
        assert_eq!((20, 36), p.index_range());

        let p = Point {
            x: -1,
            y: 0,
            seqlen: 3,
        };
        assert_eq!((36, 52), p.index_range());
    }

    #[test]
    fn test_points() {
        let p = Point {
            x: 0,
            y: 0,
            seqlen: 0,
        };
        assert_eq!((0.0, PI / 2.0), p.thetas());

        let p = Point {
            x: -1,
            y: 0,
            seqlen: 0,
        };
        assert_eq!((PI / 2.0, PI), p.thetas());

        let p = Point {
            x: 0,
            y: -1,
            seqlen: 0,
        };
        assert_eq!((3.0 / 4.0 * 2.0 * PI, 2.0 * PI), p.thetas());
    }

    #[test]
    fn test_modulus() {
        assert_eq!(PI, (-PI + 2.0 * PI) % (2.0 * PI));
    }

    #[test]
    fn test_delta_theta() {
      assert_eq!((0.0, 0.0), delta_theta(0.0, 0.0))
    }
}
