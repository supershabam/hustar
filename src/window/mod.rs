struct Window<'a> {
    valid: &'a [u8],
    length: usize,
}

impl<'a> Window<'a> {
    fn new(valid: &'a [u8], length: usize) -> Self {
        Self {
            valid: valid,
            length: length,
        }
    }

    fn reiterate(i: &mut impl Iterator<Item = u8>) -> impl Iterator<Item = Vec<u8>> {
        let v = vec![0_u8];
        vec![v].into_iter()
    }
}

#[cfg(test)]
mod tests {
    use crate::window::*;

    #[test]
    fn test() {
        let w = Window::new(b"actg", 5);
    }
}
