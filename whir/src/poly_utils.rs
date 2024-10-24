use p3_field::Field;


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