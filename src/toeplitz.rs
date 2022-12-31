use std::iter;

use ark_ff::{FftField, Zero};
use ark_poly::{domain::DomainCoeff, univariate::DensePolynomial, Polynomial, UVPolynomial};

use crate::{circulant::{Circulant, is_pow_2}, next_pow2};

pub fn is_toeplitz<T: PartialEq>(elems: &Vec<Vec<T>>) {
    let n = elems.len();
    let m = elems[0].len();

    for j in 0..m {
        check_diagonal(elems, 0, j);
    }

    for i in 1..n {
        check_diagonal(elems, i, 0);
    }
}

fn check_diagonal<T: PartialEq>(elems: &Vec<Vec<T>>, i_start: usize, j_start: usize) {
    let n = elems.len();
    let m = elems[0].len();

    for (i, j) in (i_start + 1..n).zip(j_start + 1..m) {
        if elems[i][j] != elems[i_start][j_start] {
            panic!("Not Toeplitz");
        }
    }
}

pub struct UpperToeplitz<F: FftField> {
    pub(crate) repr: Vec<F>
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
        let c = self.to_circulant();
        let zeroes = vec![T::zero(); x.len()];
        c.mul_by_vec(&[x, zeroes.as_slice()].concat())
    }

    pub fn to_circulant(&self) -> Circulant<F, F> {
        let fm = self.repr.last().unwrap().clone();
        let mut circulant_repr = vec![F::zero(); self.repr.len() + 1];

        circulant_repr[0] = fm; 
        circulant_repr[self.repr.len()] = fm; 

        circulant_repr.extend_from_slice(&self.repr[..self.repr.len() - 1]);
        assert_eq!(circulant_repr.len(), self.repr.len() * 2);
        
        Circulant::new(circulant_repr)
    }
}

pub struct Toeplitz<F: FftField> {
    pub(crate) elems: Vec<Vec<F>>,
}

impl<F: FftField> Toeplitz<F> {
    pub fn new(elems: Vec<Vec<F>>) -> Self {
        Self { elems }
    }

    pub fn from_poly(poly: &DensePolynomial<F>) -> Self {
        let mut coeffs_rev = poly.coeffs[..].to_vec().clone();
        let next_pow2_degree = next_pow2(poly.degree());
        let to_extend = vec![F::zero(); next_pow2_degree - poly.degree()];
        coeffs_rev.extend_from_slice(&to_extend);

        coeffs_rev.reverse();

        let elems: Vec<Vec<_>> = (0..coeffs_rev.len() - 1)
            .map(|i| {
                let mut row = vec![F::zero(); i];
                coeffs_rev.pop();
                row.extend_from_slice(&coeffs_rev);
                row
            })
            .collect();

        Self { elems }
    }

    pub fn mul_by_vec<T: DomainCoeff<F> + std::ops::MulAssign<F> + Zero>(&self, x: &[T]) -> Vec<T> {
        let c = self.to_circulant();
        let zeroes = vec![T::zero(); x.len()];
        c.mul_by_vec(&[x, zeroes.as_slice()].concat())
    }

    pub fn to_circulant(&self) -> Circulant<F, F> {
        // C_2n = [a0, ... an, a0, a(-n), ... a(-1)]
        let n = self.elems.len();

        let mut first_half = self.elems[n - 1][..].to_vec();
        first_half.reverse();

        let second_half_reversed = self.elems[0][..].to_vec();
        let second_half: Vec<_> = iter::once(second_half_reversed[0])
            .chain(second_half_reversed.into_iter().skip(1).rev())
            .collect();

        Circulant::new([first_half, second_half].concat())
    }
}

#[cfg(test)]
mod toeplitz_test {
    use ark_bn254::Fr;
    use ark_poly::{univariate::DensePolynomial, UVPolynomial};
    use ark_std::test_rng;
    use ndarray::arr2;

    use super::{is_toeplitz, Toeplitz, UpperToeplitz};

    #[test]
    fn test_upper_tp() {
        let mut rng = test_rng();
        let d = 17;

        let poly = DensePolynomial::<Fr>::rand(d, &mut rng);

        let t = UpperToeplitz::from_poly(&poly);
        let c = t.to_circulant();

        let t = Toeplitz::from_poly(&poly); 
        let c2 = t.to_circulant();
        assert_eq!(c.repr, c2.repr);
    }

    #[test]
    fn test_check_correctness() {
        let elems = vec![
            vec![6, 7, 8, 9],
            vec![4, 6, 7, 8],
            vec![1, 4, 6, 7],
            vec![0, 1, 4, 6],
            vec![2, 0, 1, 4],
        ];

        is_toeplitz(&elems);
    }

    #[test]
    fn test_vector_mul() {
        let poly = DensePolynomial::from_coefficients_slice(&[Fr::from(1u64), Fr::from(2u64), Fr::from(3u64), Fr::from(4u64), Fr::from(5u64)]);

        let t = Toeplitz::from_poly(&poly);
        is_toeplitz(&t.elems);

        let elems = [
            [Fr::from(5), Fr::from(4), Fr::from(3), Fr::from(2)],
            [Fr::from(0), Fr::from(5), Fr::from(4), Fr::from(3)],
            [Fr::from(0), Fr::from(0), Fr::from(5), Fr::from(4)],
            [Fr::from(0), Fr::from(0), Fr::from(0), Fr::from(5)],
        ];

        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(elems[i][j], t.elems[i][j]);
            }
        }

        let x = vec![Fr::from(8), Fr::from(4), Fr::from(2), Fr::from(1)];
        let m = arr2(&elems); 
        let v = arr2(&[[Fr::from(8)], [Fr::from(4)], [Fr::from(2)], [Fr::from(1)]]);

        /*
            5 4 3 2 | 8
            0 5 4 3 | 4   
            0 0 5 4 | 2
            0 0 0 5 | 1
        */

        let res_slow = m.dot(&v);
        let res_fast = t.mul_by_vec(&x);

        let t_dot_x_slow: Vec<_> = res_slow.into_iter().collect();
        let t_dot_x_fast = res_fast[..4].to_vec();
        assert_eq!(t_dot_x_slow, t_dot_x_fast);
    }
}
