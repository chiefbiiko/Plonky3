use crate::ntt::intt_batch;
use crate::parameters::FoldType;
// use ark_ff::{FftField, Field};
use p3_field::{Field, TwoAdicField};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

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
        let mut new_answers = vec![F::zero(); offset];
        let mut coset_index_inv = F::one();
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

pub fn restructure_evaluations<F: TwoAdicField>(
    mut stacked_evaluations: Vec<F>,
    fold_type: FoldType,
    _domain_gen: F,
    domain_gen_inv: F,
    folding_factor: usize,
) -> Vec<F> {
    let folding_size = 1_u64 << folding_factor;
    assert_eq!(stacked_evaluations.len() % (folding_size as usize), 0);
    match fold_type {
        FoldType::Naive => stacked_evaluations,
        FoldType::ProverHelps => {
            // TODO: This partially undoes the NTT transform from tne encoding.
            // Maybe there is a way to not do the full transform in the first place.

            // Batch inverse NTTs
            intt_batch(&mut stacked_evaluations, folding_size as usize);

            // Apply coset and size correction.
            // Stacked evaluation at i is f(B_l) where B_l = w^i * <w^n/k>
            let size_inv = F::from_canonical_u64(folding_size).inverse();//.unwrap();
            #[cfg(not(feature = "parallel"))]
            {
                let mut coset_offset_inv = F::one();
                for answers in stacked_evaluations.chunks_exact_mut(folding_size as usize) {
                    let mut scale = size_inv;
                    for v in answers.iter_mut() {
                        *v *= scale;
                        scale *= coset_offset_inv;
                    }
                    coset_offset_inv *= domain_gen_inv;
                }
            }
            #[cfg(feature = "parallel")]
            stacked_evaluations
                .par_chunks_exact_mut(folding_size as usize)
                .enumerate()
                .for_each_with(F::zero(), |offset, (i, answers)| {
                    if *offset == F::zero() {
                        // *offset = domain_gen_inv.pow([i as u64]);
                        *offset = domain_gen_inv.exp_u64(i as u64);
                    } else {
                        *offset *= domain_gen_inv;
                    }
                    let mut scale = size_inv;
                    for v in answers.iter_mut() {
                        *v *= scale;
                        scale *= &*offset;
                    }
                });

            stacked_evaluations
        }
    }
}

#[cfg(test)]
mod tests {
    // use p3_field::extension::BinomialExtensionField;
    // use ark_ff::{FftField, Field};
    use p3_field::AbstractField;
    use p3_field::{Field, TwoAdicField};
    use p3_mersenne_31::Mersenne31;
    use p3_bn254_fr::Bn254Fr;
    use p3_baby_bear::BabyBear;
    use p3_goldilocks::Goldilocks;

    use crate::{
        poly_utils::{coeffs::CoefficientList, MultilinearPoint, helpers::get_root_of_unity},
        utils::stack_evaluations,
    };

    use super::{compute_fold, restructure_evaluations};

    use p3_field::{extension::{HasTwoAdicBionmialExtension, BinomialExtensionField}, PrimeField64};

    // type F = Field64;
    // type F = Mersenne31;
    // type F = BinomialExtensionField<Goldilocks, 2>;
    type F = Goldilocks;

    #[test]
    fn test_folding_wip() {
        let num_variables = 5;
        let num_coeffs = 1 << num_variables;

        let domain_size = 32;
        let folding_factor = 3; // We fold in 8
        let folding_factor_exp = 1 << folding_factor;

        let poly = CoefficientList::new((0..num_coeffs).map(F::from_canonical_u64).collect());

        // let root_of_unity = F::get_root_of_unity(domain_size).unwrap();
        //
        //WRONG need get_root_of_unity(bits)
        // let root_of_unity = F::two_adic_generator(domain_size);
        let root_of_unity: F = get_root_of_unity(domain_size).expect("nth root of unity");
        println!("root_of_unity\t{:?}", root_of_unity);
        let index = 15;
        let folding_randomness: Vec<_> = (0..folding_factor).map(|i| F::from_canonical_u64(i as u64)).collect();

        let coset_offset = root_of_unity.exp_u64(index);// % F::from_canonical_u64(F::ORDER_U64);
        println!("coset_offset\t{:?}", coset_offset);
        let coset_gen = root_of_unity.exp_u64((domain_size / folding_factor_exp) as u64);// % F::ORDER_U64;
        println!("coset_gen\t{:?}", coset_gen);

        // Evaluate the polynomial on the coset
        let poly_eval: Vec<_> = (0..folding_factor_exp)
            .map(|i| {
                poly.evaluate(&MultilinearPoint::expand_from_univariate(
                    coset_offset * coset_gen.exp_u64(i as u64),
                    num_variables,
                ))
            })
            .collect();

        let fold_value = compute_fold(
            &poly_eval,
            &folding_randomness,
            coset_offset.inverse(),//.unwrap(),
            coset_gen.inverse(),//.unwrap(),
            F::from_canonical_u64(2).inverse(),//.unwrap(),
            folding_factor,
        );

        let truth_value = poly.fold(&MultilinearPoint(folding_randomness)).evaluate(
            &MultilinearPoint::expand_from_univariate(
                root_of_unity.exp_u64((folding_factor_exp as u64) * index),
                2,
            ),
        );

        assert_eq!(fold_value, truth_value);
    }

    #[test]
    fn test_folding_optimised() {
        let num_variables = 5;
        let num_coeffs = 1 << num_variables;

        let domain_size = 32;
        let folding_factor = 3; // We fold in 8
        let folding_factor_exp = 1 << folding_factor;

        let poly = CoefficientList::new((0..num_coeffs).map(F::from_canonical_u64).collect());

        // let root_of_unity = F::get_root_of_unity(domain_size).unwrap();
        // let root_of_unity = <F as TwoAdicField>::two_adic_generator(domain_size);
        // let root_of_unity = <F as TwoAdicField>::two_adic_generator(domain_size);
        let root_of_unity: F = get_root_of_unity(domain_size).expect("nth root of unity");
        println!("root_of_unity\t{:?}", root_of_unity);
        let root_of_unity_inv = root_of_unity.inverse();//.unwrap();

        let folding_randomness: Vec<_> = (0..folding_factor).map(|i| F::from_canonical_u64(i as u64)).collect();

        // Evaluate the polynomial on the domain
        let domain_evaluations: Vec<_> = (0..domain_size)
            // .map(|w| root_of_unity.pow([w as u64]))
            .map(|w| root_of_unity.exp_u64(w as u64))
            .map(|point| {
                poly.evaluate(&MultilinearPoint::expand_from_univariate(
                    point,
                    num_variables,
                ))
            })
            .collect();

        let unprocessed = stack_evaluations(domain_evaluations, folding_factor);

        let processed = restructure_evaluations(
            unprocessed.clone(),
            crate::parameters::FoldType::ProverHelps,
            root_of_unity,
            root_of_unity_inv,
            folding_factor,
        );

        let num = domain_size / folding_factor_exp;
        // let coset_gen_inv = root_of_unity_inv.pow(&[num]);
        let coset_gen_inv = root_of_unity_inv.exp_u64(num as u64);

        for index in 0..num {
            // let offset_inv = root_of_unity_inv.pow(&[index]);
            let offset_inv = root_of_unity_inv.exp_u64(index as u64);
            let span =
                (index * folding_factor_exp) as usize..((index + 1) * folding_factor_exp) as usize;

            let answer_unprocessed = compute_fold(
                &unprocessed[span.clone()],
                &folding_randomness,
                offset_inv,
                coset_gen_inv,
                F::from_canonical_u64(2).inverse(),//.unwrap(),
                folding_factor,
            );

            let answer_processed = CoefficientList::new(processed[span].to_vec())
                .evaluate(&MultilinearPoint(folding_randomness.clone()));

            assert_eq!(answer_processed, answer_unprocessed);
        }
    }
}