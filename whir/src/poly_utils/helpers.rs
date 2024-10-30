use p3_field::AbstractField;

pub fn eval_poly<AF: AbstractField>(poly: &[AF], x: AF) -> AF {
    let mut acc = AF::zero();
    for coeff in poly.iter().rev() {
        acc *= x.clone();
        acc += coeff.clone();
    }
    acc
}