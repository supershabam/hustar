// traverse helps to iterate the mmap in a minimally querying way
// to populate a desired set of pixel points.

use core::f64::consts::PI;
use itertools::Itertools;

#[derive(Debug, Clone, Copy, Default)]
pub struct Point {
    pub x: i32,
    pub y: i32,
    pub w: u32,
    pub h: u32,
    pub seqlen: usize,
}

fn theta_delta(t1: f64, t2: f64) -> (f64, f64) {
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

pub fn make_points(width: u32, height: u32, max_length: u32) -> Vec<Point> {
    let mut points: Vec<Point> = Vec::with_capacity(width as usize * height as usize);
    let edge = width.max(height);
    let edge = edge / 2;
    let edge_squared = edge * edge;
    for w in 0..width {
        for h in 0..height {
            let x = width as i32 / 2 - w as i32;
            let y = height as i32 / 2 - h as i32;
            let r_squared = x * x + y * y;
            let p = r_squared as f64 / edge_squared as f64;
            let p = p.min(1.0);
            let seqlen = (p * (max_length - 1) as f64) as usize + 1;
            points.push(Point {
                x: x,
                y: y,
                w,
                h,
                seqlen,
            });
        }
    }
    points.sort_by(|a, b| match a.seqlen.cmp(&b.seqlen) {
        std::cmp::Ordering::Equal => {
            let ar = a.index_range();
            let br = b.index_range();
            match ar.0.cmp(&br.0) {
                std::cmp::Ordering::Equal => ar.1.cmp(&br.1),
                ord => ord,
            }
        }
        ord => ord,
    });
    points
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
                    let td = theta_delta(a, b);
                    // select maximum distance
                    if td.1 > acc.1 {
                        td
                    } else {
                        acc
                    }
                }
                _ => panic!("permutations produced bad combination"),
            }
        });
        (td.0, td.0 + td.1)
    }

    // index_range returns the [min, max) range over which the pixel
    // overlaps in mmap index space.
    pub fn index_range(&self) -> (usize, usize) {
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
            ..Point::default()
        };
        assert_eq!((0, 1), p.index_range());

        let p = Point {
            x: 0,
            y: 0,
            seqlen: 2,
            ..Point::default()
        };
        assert_eq!((4, 8), p.index_range());

        let p = Point {
            x: -1,
            y: 0,
            seqlen: 2,
            ..Point::default()
        };
        assert_eq!((8, 12), p.index_range());

        let p = Point {
            x: 0,
            y: 0,
            seqlen: 3,
            ..Point::default()
        };
        assert_eq!((20, 36), p.index_range());

        let p = Point {
            x: -1,
            y: 0,
            seqlen: 3,
            ..Point::default()
        };
        assert_eq!((36, 52), p.index_range());
    }

    #[test]
    fn test_points() {
        let p = Point {
            x: 0,
            y: 0,
            seqlen: 0,
            ..Point::default()
        };
        assert_eq!((0.0, PI / 2.0), p.thetas());

        let p = Point {
            x: -1,
            y: 0,
            seqlen: 0,
            ..Point::default()
        };
        assert_eq!((PI / 2.0, PI), p.thetas());

        let p = Point {
            x: 0,
            y: -1,
            seqlen: 0,
            ..Point::default()
        };
        assert_eq!((3.0 / 4.0 * 2.0 * PI, 2.0 * PI), p.thetas());
    }

    #[test]
    fn test_modulus() {
        assert_eq!(PI, (-PI + 2.0 * PI) % (2.0 * PI));
    }

    #[test]
    fn test_theta_delta() {
        assert_eq!((0.0, 0.0), theta_delta(0.0, 0.0))
    }

    #[test]
    fn test_make_points() {
        let points = make_points(50, 50, 8);
        for p in &points {
            println!("point={:?} range={:?}", p, p.index_range());
        }
        // println!("points={:?}", points);
        assert!(false);
    }
}
