use ark_ff::FftField;
use ark_poly::{domain::DomainCoeff, EvaluationDomain, GeneralEvaluationDomain};
use std::{fmt::Debug, marker::PhantomData};

use crate::is_pow_2;

pub struct Circulant<F: FftField, D: DomainCoeff<F> + Debug> {
    _f: PhantomData<F>,
    _d: PhantomData<D>,
}

impl<F: FftField, D: DomainCoeff<F> + Debug> Circulant<F, D> {
    pub fn mul_by_vec<T: DomainCoeff<F> + std::ops::MulAssign<D>>(repr: &[D], x: &[T]) -> Vec<T> {
        assert!(is_pow_2(repr.len()));
        let domain = GeneralEvaluationDomain::new(repr.len()).unwrap();
        let v = domain.fft(repr);

        let mut res = domain.fft(x);
        for (i, _) in x.iter().enumerate() {
            res[i] *= v[i]
        }

        domain.ifft(&res)
    }
}
