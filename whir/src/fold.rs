// use crate::ntt::intt_batch;
// use crate::parameters::FoldType;
// use ark_ff::{FftField, Field};
use p3_field::Field;

// Given the evaluation of f on the coset specified by coset_offset * <coset_gen>
// Compute the fold on that point
pub fn compute_fold<F: Field>(
    answers: &[F],
    folding_randomness: &[F],
    mut coset_offset_inv: F,
    mut coset_gen_inv: F,
    two_inv: F,
    folding_factor: usize,
) -> F {
    let mut answers = answers.to_vec();

    // We recursively compute the fold, rec is where it is
    for rec in 0..folding_factor {
        let offset = answers.len() / 2;
        let mut new_answers = vec![F::ZERO; offset];
        let mut coset_index_inv = F::ONE;
        for i in 0..offset {
            let f_value_0 = answers[i];
            let f_value_1 = answers[i + offset];
            let point_inv = coset_offset_inv * coset_index_inv;

            let left = f_value_0 + f_value_1;
            let right = point_inv * (f_value_0 - f_value_1);

            new_answers[i] =
                two_inv * (left + folding_randomness[folding_randomness.len() - 1 - rec] * right);
            coset_index_inv *= coset_gen_inv;
        }
        answers = new_answers;

        // Update for next one
        coset_offset_inv = coset_offset_inv * coset_offset_inv;
        coset_gen_inv = coset_gen_inv * coset_gen_inv;
    }

    answers[0]
}