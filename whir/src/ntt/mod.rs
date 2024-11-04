//! NTT and related algorithms.

mod matrix;
mod ntt;
mod transpose;
mod utils;
mod wavelet;

use self::matrix::MatrixMut;
// use ark_ff::FftField;
use p3_field::TwoAdicField;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub use self::{
    ntt::{intt, intt_batch, ntt, ntt_batch},
    transpose::transpose,
    wavelet::wavelet_transform,
};
use crate::utils::pow;

/// RS encode at a rate 1/`expansion`.
pub fn expand_from_coeff<F: TwoAdicField>(coeffs: &[F], expansion: usize) -> Vec<F> {
    let engine = ntt::NttEngine::<F>::new_from_cache();
    let expanded_size = coeffs.len() * expansion;
    let mut result = Vec::with_capacity(expanded_size);
    // Note: We can also zero-extend the coefficients and do a larger NTT.
    // But this is more efficient.

    // Do coset NTT.
    let root = engine.root(expanded_size);
    result.extend_from_slice(coeffs);
    #[cfg(not(feature = "parallel"))]
    for i in 1..expansion {
        // let root = root.exp_u64(i as u64);//.pow([i as u64]);
        let root = pow(root, i as u64);//root.exp_u64(i as u64);//.pow([i as u64]);
        let mut offset = F::one();
        result.extend(coeffs.iter().map(|x| {
            let val = *x * offset;
            offset *= root;
            val
        }));
    }
    #[cfg(feature = "parallel")]
    result.par_extend((1..expansion).into_par_iter().flat_map(|i| {
        let root_i = root.pow([i as u64]);
        coeffs
            .par_iter()
            .enumerate()
            .map_with(F::ZERO, move |root_j, (j, coeff)| {
                if root_j.is_zero() {
                    *root_j = root_i.pow([j as u64]);
                } else {
                    *root_j *= root_i;
                }
                *coeff * *root_j
            })
    }));

    ntt_batch(&mut result, coeffs.len());
    transpose(&mut result, expansion, coeffs.len());
    result
}