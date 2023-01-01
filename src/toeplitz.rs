use ark_ff::{FftField, Zero};
use ark_poly::{domain::DomainCoeff, univariate::DensePolynomial, Polynomial, UVPolynomial};
use ark_std::log2;

use crate::{
    circulant::{is_pow_2, Circulant},
};

pub fn next_pow2(n: usize) -> usize {
    let two: u32 = 2;
    let a: u32 = log2(n);

    if two.pow(a - 1) == n as u32 {
        return n;
    }

    two.pow(a).try_into().unwrap()
}


/*
    fm f(m-1) ... f1
    0   fm    ... f2
    .   ...   ... f3
    .   ... fm f(m-1)
    0   ...   ... fm 
*/
/// Succinct representation of Toeplitz matrix that is instantiated from polynomial 
/// on which mul by vector can be run efficiently
pub struct UpperToeplitz<F: FftField> {
    pub(crate) repr: Vec<F>,
}

impl<F: FftField> UpperToeplitz<F> {
    pub fn from_poly(poly: &DensePolynomial<F>) -> Self {
        let mut repr = poly.coeffs()[1..].to_vec();
        let next_pow2_degree = next_pow2(poly.degree());
        let to_extend = vec![F::zero(); next_pow2_degree - poly.degree()];
        repr.extend_from_slice(&to_extend);
        assert!(is_pow_2(repr.len()));
        Self { repr }
    }

    pub fn mul_by_vec<T: DomainCoeff<F> + std::ops::MulAssign<F> + Zero>(&self, x: &[T]) -> Vec<T> {
        let circulant_repr = self.to_circulant_repr();
        let zeroes = vec![T::zero(); x.len()];
        Circulant::mul_by_vec(&circulant_repr, &[x, zeroes.as_slice()].concat())
    }

    fn to_circulant_repr(&self) -> Vec<F> {
        let fm = *self.repr.last().unwrap();
        let mut circulant_repr = vec![F::zero(); self.repr.len() + 1];

        circulant_repr[0] = fm;
        circulant_repr[self.repr.len()] = fm;

        circulant_repr.extend_from_slice(&self.repr[..self.repr.len() - 1]);
        circulant_repr
    }
}