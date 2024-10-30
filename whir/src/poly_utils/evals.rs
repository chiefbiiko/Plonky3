use std::ops::Index;

use p3_field::Field;

use super::sequential_lag_poly::LagrangePolynomialIterator;
use super::MultilinearPoint;

#[derive(Debug)]
pub struct EvaluationsList<F> {
    evals: Vec<F>,
    num_variables: usize,
}

impl<F> EvaluationsList<F>
where
    F: Field,
{
    pub fn new(evals: Vec<F>) -> Self {
        let len = evals.len();
        assert!(len.is_power_of_two());
        let num_variables = len.ilog2();

        EvaluationsList {
            evals,
            num_variables: num_variables as usize,
        }
    }

    pub fn evaluate(&self, point: &MultilinearPoint<F>) -> F {
        if let Some(point) = point.to_hypercube() {
            return self.evals[point.0];
        }

        let mut sum = F::zero();
        for (b, lag) in LagrangePolynomialIterator::new(point) {
            sum += lag * self.evals[b.0]
        }

        sum
    }

    pub fn evals(&self) -> &[F] {
        &self.evals
    }

    pub fn evals_mut(&mut self) -> &mut [F] {
        &mut self.evals
    }

    pub fn num_evals(&self) -> usize {
        self.evals.len()
    }

    pub fn num_variables(&self) -> usize {
        self.num_variables
    }
}

impl<F> Index<usize> for EvaluationsList<F> {
    type Output = F;
    fn index(&self, index: usize) -> &Self::Output {
        &self.evals[index]
    }
}

#[cfg(test)]
mod tests {
    use p3_field::AbstractField;
    use p3_mersenne_31::Mersenne31;

    use super::*;
    use crate::poly_utils::hypercube::BinaryHypercube;

    type F = Mersenne31;

    #[test]
    fn test_evaluation() {
        let evaluations_vec = vec![F::zero(), F::one(), F::zero(), F::one()];
        let evals = EvaluationsList::new(evaluations_vec.clone());

        for i in BinaryHypercube::new(2) {
            assert_eq!(
                evaluations_vec[i.0],
                evals.evaluate(&MultilinearPoint::from_binary_hypercube_point(i, 2))
            );
        }
    }
}
