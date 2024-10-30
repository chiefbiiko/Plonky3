// use ark_ff::Field;
use p3_field::Field;
use rand::{
    distributions::{Distribution, Standard},
    Rng, RngCore,
};

use crate::utils::to_binary;

use hypercube::BinaryHypercubePoint;

//TODO
pub mod coeffs;
pub mod evals;
pub mod fold;
// pub mod gray_lag_poly;
pub mod hypercube;
pub mod sequential_lag_poly;
// pub mod streaming_evaluation_helper;

pub mod helpers;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultilinearPoint<F>(pub Vec<F>);

impl<F> MultilinearPoint<F>
where
    F: Field,
{
    pub fn n_variables(&self) -> usize {
        self.0.len()
    }

    pub fn from_binary_hypercube_point(point: BinaryHypercubePoint, num_variables: usize) -> Self {
        Self(
            to_binary(point.0, num_variables)
                .into_iter()
                .map(|x| if x { F::one() } else { F::zero() })
                .collect(),
        )
    }

    pub fn to_hypercube(&self) -> Option<BinaryHypercubePoint> {
        let mut counter = 0;
        for &coord in &self.0 {
            if coord == F::zero() {
                counter <<= 1;
            } else if coord == F::one() {
                counter = (counter << 1) + 1;
            } else {
                return None;
            }
        }

        Some(BinaryHypercubePoint(counter))
    }

    pub fn expand_from_univariate(point: F, num_variables: usize) -> Self {
        let mut res = Vec::with_capacity(num_variables);
        let mut cur = point;
        for _ in 0..num_variables {
            res.push(cur);
            cur = cur * cur;
        }

        // Reverse so higher power is first
        res.reverse();

        MultilinearPoint(res)
    }
}

impl<F> MultilinearPoint<F>
where
    Standard: Distribution<F>,
{
    pub fn rand(rng: &mut impl RngCore, num_variables: usize) -> Self {
        MultilinearPoint((0..num_variables).map(|_| rng.gen()).collect())
    }
}

impl<F> From<F> for MultilinearPoint<F> {
    fn from(value: F) -> Self {
        MultilinearPoint(vec![value])
    }
}

pub fn eq_poly<F>(coords: &MultilinearPoint<F>, point: BinaryHypercubePoint) -> F
where
    F: Field,
{
    let mut point = point.0;
    let n_variables = coords.n_variables();
    assert!(point < (1 << n_variables));

    let mut acc = F::one();

    for val in coords.0.iter().rev() {
        let b = point % 2;
        acc *= if b == 1 { *val } else { F::one() - *val };
        point >>= 1;
    }

    acc
}

pub fn eq_poly_outside<F>(coords: &MultilinearPoint<F>, point: &MultilinearPoint<F>) -> F
where
    F: Field,
{
    assert_eq!(coords.n_variables(), point.n_variables());

    let mut acc = F::one();

    for (&l, &r) in coords.0.iter().zip(&point.0) {
        acc *= l * r + (F::one() - l) * (F::one() - r);
    }

    acc
}

pub fn eq_poly3<F>(coords: &MultilinearPoint<F>, mut point: usize) -> F
where
    F: Field,
{
    let two = F::one() + F::one();
    let two_inv = two.inverse();//.unwrap();

    let n_variables = coords.n_variables();
    assert!(point < 3usize.pow(n_variables as u32));

    let mut acc = F::one();

    for &val in coords.0.iter().rev() {
        let b = point % 3;
        acc *= match b {
            0 => (val - F::one()) * (val - two) * two_inv,
            1 => val * (val - two) * (-F::one()),
            2 => val * (val - F::one()) * two_inv,
            _ => unreachable!(),
        };
        point /= 3;
    }

    acc
}

#[cfg(test)]
mod tests {
    use crate::poly_utils::eq_poly3;
    use crate::poly_utils::hypercube::BinaryHypercube;
    use crate::poly_utils::eq_poly;
    use p3_mersenne_31::Mersenne31;
    use p3_field::AbstractField;

    use super::coeffs::CoefficientList;
    use super::BinaryHypercubePoint;
    use super::MultilinearPoint;

    type F = Mersenne31;

    #[test]
    fn test_equality() {
        let point = MultilinearPoint(vec![F::from_canonical_u64(0), F::from_canonical_u64(0)]);
        assert_eq!(eq_poly(&point, BinaryHypercubePoint(0)), F::from_canonical_u64(1));
        assert_eq!(eq_poly(&point, BinaryHypercubePoint(1)), F::from_canonical_u64(0));
        assert_eq!(eq_poly(&point, BinaryHypercubePoint(2)), F::from_canonical_u64(0));
        assert_eq!(eq_poly(&point, BinaryHypercubePoint(3)), F::from_canonical_u64(0));

        let point = MultilinearPoint(vec![F::from_canonical_u64(1), F::from_canonical_u64(0)]);
        assert_eq!(eq_poly(&point, BinaryHypercubePoint(0)), F::from_canonical_u64(0));
        assert_eq!(eq_poly(&point, BinaryHypercubePoint(1)), F::from_canonical_u64(0));
        assert_eq!(eq_poly(&point, BinaryHypercubePoint(2)), F::from_canonical_u64(1));
        assert_eq!(eq_poly(&point, BinaryHypercubePoint(3)), F::from_canonical_u64(0));
    }

    #[test]
    fn test_equality_again() {
        let poly = CoefficientList::new(vec![F::from_canonical_u64(35), F::from_canonical_u64(97), F::from_canonical_u64(10), F::from_canonical_u64(32)]);
        let point = MultilinearPoint(vec![F::from_canonical_u64(42), F::from_canonical_u64(36)]);
        let eval = poly.evaluate(&point);

        assert_eq!(
            eval,
            BinaryHypercube::new(2)
                .map(
                    |i| poly.evaluate(&MultilinearPoint::from_binary_hypercube_point(i, 2))
                        * eq_poly(&point, i)
                )
                .sum()
        );
    }

    #[test]
    fn test_equality3() {
        let point = MultilinearPoint(vec![F::from_canonical_u64(0), F::from_canonical_u64(0)]);

        assert_eq!(eq_poly3(&point, 0), F::from_canonical_u64(1));
        assert_eq!(eq_poly3(&point, 1), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 2), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 3), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 4), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 5), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 6), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 7), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 8), F::from_canonical_u64(0));

        let point = MultilinearPoint(vec![F::from_canonical_u64(1), F::from_canonical_u64(0)]);

        assert_eq!(eq_poly3(&point, 0), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 1), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 2), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 3), F::from_canonical_u64(1));
        assert_eq!(eq_poly3(&point, 4), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 5), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 6), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 7), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 8), F::from_canonical_u64(0));

        let point = MultilinearPoint(vec![F::from_canonical_u64(0), F::from_canonical_u64(2)]);

        assert_eq!(eq_poly3(&point, 0), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 1), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 2), F::from_canonical_u64(1));
        assert_eq!(eq_poly3(&point, 3), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 4), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 5), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 6), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 7), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 8), F::from_canonical_u64(0));

        let point = MultilinearPoint(vec![F::from_canonical_u64(2), F::from_canonical_u64(2)]);

        assert_eq!(eq_poly3(&point, 0), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 1), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 2), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 3), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 4), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 5), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 6), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 7), F::from_canonical_u64(0));
        assert_eq!(eq_poly3(&point, 8), F::from_canonical_u64(1));
    }

    #[test]
    #[should_panic]
    fn test_equality_2() {
        let point = MultilinearPoint(vec![F::from_canonical_u64(0), F::from_canonical_u64(0)]);

        let _x = eq_poly(&point, BinaryHypercubePoint(4));
    }

    #[test]
    fn expand_from_univariate() {
        let num_variables = 4;

        let point0 = MultilinearPoint::expand_from_univariate(F::from_canonical_u64(0), num_variables);
        let point1 = MultilinearPoint::expand_from_univariate(F::from_canonical_u64(1), num_variables);
        let point2 = MultilinearPoint::expand_from_univariate(F::from_canonical_u64(2), num_variables);

        assert_eq!(point0.n_variables(), num_variables);
        assert_eq!(point1.n_variables(), num_variables);
        assert_eq!(point2.n_variables(), num_variables);

        assert_eq!(
            MultilinearPoint::from_binary_hypercube_point(BinaryHypercubePoint(0), num_variables),
            point0
        );

        assert_eq!(
            MultilinearPoint::from_binary_hypercube_point(
                BinaryHypercubePoint((1 << num_variables) - 1),
                num_variables
            ),
            point1
        );

        assert_eq!(
            MultilinearPoint(vec![F::from_canonical_u64(256), F::from_canonical_u64(16), F::from_canonical_u64(4), F::from_canonical_u64(2)]),
            point2
        );
    }

    #[test]
    fn from_hypercube_and_back() {
        let hypercube_point = BinaryHypercubePoint(24);
        assert_eq!(
            Some(hypercube_point),
            MultilinearPoint::<F>::from_binary_hypercube_point(hypercube_point, 5).to_hypercube()
        );
    }
}