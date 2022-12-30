use std::{fmt::{Display, Formatter, Result, Debug}};

use ark_ff::FftField;
use ark_poly::{domain::DomainCoeff, GeneralEvaluationDomain, EvaluationDomain};

pub struct Circulant<F: FftField, D: DomainCoeff<F> + Debug> {
    pub(crate) domain: GeneralEvaluationDomain<F>,
    pub(crate) repr: Vec<D>,
    pub(crate) v: Vec<D>,
    pub(crate) elems: Vec<Vec<D>>,
}

pub fn is_pow_2(x: usize) -> bool {
    (x & (x - 1)) == 0
}

impl<F: FftField, D: DomainCoeff<F> + Debug> Circulant<F, D> {
    pub fn new(repr: Vec<D>) -> Self {
        assert!(is_pow_2(repr.len()));
        let n: i32 = repr.len() as i32;
        let mut indices: Vec<i32> = (0..n).map(|i| -i).collect();

        let resolve_index = |i: i32| {
            if i >= 0 {
                repr[i as usize]
            } else {
                repr[repr.len() - (-i as usize)]
            }
        };

        let mut elems = vec![Vec::<D>::with_capacity(repr.len()); repr.len()];
        for i in 0..repr.len() {
            for j in 0..repr.len() {
                elems[i].push(resolve_index(indices[j]));
                indices[j] += 1;
            }
        }

        let domain = GeneralEvaluationDomain::new(repr.len()).unwrap();
        let v = domain.fft(&repr);

        Self {
            domain,
            v,
            repr, 
            elems
        }
    }

    pub fn mul_by_vec<T: DomainCoeff<F> + std::ops::MulAssign<D>>(&self, x: &Vec<T>) -> Vec<T> {
        assert_eq!(x.len(), self.repr.len());
        let mut res = self.domain.fft(x);
        for i in 0..x.len() {
            res[i] *= self.v[i]
        }
        self.domain.ifft(&res)
    }
}

impl<F: FftField, D: DomainCoeff<F> + Debug> Display for Circulant<F, D> {
    fn fmt(&self, f: &mut Formatter) -> Result {
        let n = self.repr.len(); 

        for i in 0..n {
            for j in 0..n {
                write!(f, "{:?}", self.elems[i][j])?;
                write!(f, "{}", " ")?;
            }

            write!(f, "{}", "\n")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod circulant_tests {
    use ark_bn254::{Fr, G1Projective, G1Affine};
    use ark_ec::{AffineCurve, ProjectiveCurve};
    use ark_ff::{UniformRand, FromBytes};
    use ark_poly::{GeneralEvaluationDomain, EvaluationDomain};
    use ark_std::test_rng;

    use super::Circulant;

    #[test]
    fn test_mul_by_vector() {
        let n = 32; 
        let mut rng = test_rng();

        let repr: Vec<_> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let c = Circulant::<Fr, Fr>::new(repr);

        let x: Vec<_> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        // mul by scalars   
        let _ = c.mul_by_vec(&x);

        let g_1 = G1Affine::prime_subgroup_generator();
        let x: Vec<_> = (0..n).map(|_| g_1.mul(Fr::rand(&mut rng))).collect(); 
        
        // mul by ec points
        let _ = c.mul_by_vec(&x);
    }
}
