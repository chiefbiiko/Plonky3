// NOTE: This is the one from Blendy

// use ark_ff::Field;
use p3_field::{Field, TwoAdicField};

use super::{hypercube::BinaryHypercubePoint, MultilinearPoint};

pub struct LagrangePolynomialIterator<F: Field> {
    last_position: Option<usize>,
    point: Vec<F>,
    point_negated: Vec<F>,
    stack: Vec<F>,
    num_variables: usize,
}

impl<F: Field> LagrangePolynomialIterator<F> {
    pub fn new(point: &MultilinearPoint<F>) -> Self {
        let num_variables = point.0.len();

        // Initialize a stack with capacity for messages/ message_hats and the identity element
        let mut stack: Vec<F> = Vec::with_capacity(point.0.len() + 1);
        stack.push(F::one());

        let mut point = point.0.clone();
        let mut point_negated: Vec<_> = point.iter().map(|x| F::one() - *x).collect();
        // Iterate over the message_hats, update the running product, and push it onto the stack
        let mut running_product: F = F::one();
        for point_neg in &point_negated {
            running_product *= *point_neg;
            stack.push(running_product);
        }

        point.reverse();
        point_negated.reverse();

        // Return
        Self {
            num_variables,
            point,
            point_negated,
            stack,
            last_position: None,
        }
    }
}

impl<F: Field> Iterator for LagrangePolynomialIterator<F> {
    type Item = (BinaryHypercubePoint, F);
    // Iterator implementation for the struct
    fn next(&mut self) -> Option<Self::Item> {
        // a) Check if this is the first iteration
        if self.last_position == None {
            // Initialize last position
            self.last_position = Some(0);
            // Return the top of the stack
            return Some((BinaryHypercubePoint(0), *self.stack.last().unwrap()));
        }

        // b) Check if in the last iteration we finished iterating
        if self.last_position.unwrap() + 1 >= 1 << self.num_variables {
            return None;
        }

        // c) Everything else, first get bit diff
        let last_position = self.last_position.unwrap();
        let next_position = last_position + 1;
        let bit_diff = last_position ^ next_position;

        // Determine the shared prefix of the most significant bits
        let low_index_of_prefix = (bit_diff + 1).trailing_zeros() as usize;

        // Discard any stack values outside of this prefix
        self.stack.truncate(self.stack.len() - low_index_of_prefix);

        // Iterate up to this prefix computing lag poly correctly
        for bit_index in (0..low_index_of_prefix).rev() {
            let last_element = self.stack.last().unwrap();
            let next_bit: bool = (next_position & (1 << bit_index)) != 0;
            self.stack.push(match next_bit {
                true => *last_element * self.point[bit_index],
                false => *last_element * self.point_negated[bit_index],
            });
        }

        // Don't forget to update the last position
        self.last_position = Some(next_position);

        // Return the top of the stack
        Some((
            BinaryHypercubePoint(next_position),
            *self.stack.last().unwrap(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        // crypto::fields::Field64,
        poly_utils::{eq_poly, hypercube::BinaryHypercubePoint, MultilinearPoint},
    };
    use p3_mersenne_31::Mersenne31;
    use p3_field::AbstractField;

    use super::LagrangePolynomialIterator;

    type F = Mersenne31;

    #[test]
    fn test_blendy() {
        let one = F::from_canonical_u64(1);
        let (a, b) = (F::from_canonical_u64(2), F::from_canonical_u64(3));
        let point_1 = MultilinearPoint(vec![a, b]);

        let mut lag_iterator = LagrangePolynomialIterator::new(&point_1);

        assert_eq!(
            lag_iterator.next().unwrap(),
            (BinaryHypercubePoint(0), (one - a) * (one - b))
        );
        assert_eq!(
            lag_iterator.next().unwrap(),
            (BinaryHypercubePoint(1), (one - a) * b)
        );
        assert_eq!(
            lag_iterator.next().unwrap(),
            (BinaryHypercubePoint(2), a * (one - b))
        );
        assert_eq!(
            lag_iterator.next().unwrap(),
            (BinaryHypercubePoint(3), a * b)
        );
        assert_eq!(lag_iterator.next(), None);
    }

    #[test]
    fn test_blendy_2() {
        let point = MultilinearPoint(vec![F::from_canonical_u64(12), F::from_canonical_u64(13), F::from_canonical_u64(32)]);

        let mut last_b = None;
        for (b, lag) in LagrangePolynomialIterator::new(&point) {
            assert_eq!(eq_poly(&point, b), lag);
            assert!(b.0 < 1 << 3);
            last_b = Some(b);
        }
        assert_eq!(last_b, Some(BinaryHypercubePoint(7)));
    }

    #[test]
    fn test_blendy_3() {
        let point = MultilinearPoint(vec![
            F::from_canonical_u64(414151),
            F::from_canonical_u64(109849018),
            F::from_canonical_u64(033184190),
            F::from_canonical_u64(033184190),
            F::from_canonical_u64(033184190),
        ]);

        let mut last_b = None;
        for (b, lag) in LagrangePolynomialIterator::new(&point) {
            assert_eq!(eq_poly(&point, b), lag);
            last_b = Some(b);
        }
        assert_eq!(last_b, Some(BinaryHypercubePoint(31)));
    }
}