use std::fmt::{Debug, Display, Formatter, Result};

use ark_ff::FftField;
use ark_poly::{domain::DomainCoeff, EvaluationDomain, GeneralEvaluationDomain};

pub struct Circulant<F: FftField, D: DomainCoeff<F> + Debug> {
    pub(crate) domain: GeneralEvaluationDomain<F>,
    pub(crate) repr: Vec<D>,
    pub(crate) v: Vec<D>,
}

pub fn is_pow_2(x: usize) -> bool {
    (x & (x - 1)) == 0
}

pub fn resolve_index<T: Copy>(i: i32, x: &Vec<T>) -> T {
    if i >= 0 {
        x[i as usize]
    } else {
        x[x.len() - (-i as usize)]
    }
}

impl<F: FftField, D: DomainCoeff<F> + Debug> Circulant<F, D> {
    pub fn new(repr: Vec<D>) -> Self {
        assert!(is_pow_2(repr.len()));

        let domain = GeneralEvaluationDomain::new(repr.len()).unwrap();
        let v = domain.fft(&repr);

        Self { domain, v, repr }
    }

    pub fn mul_by_vec<T: DomainCoeff<F> + std::ops::MulAssign<D>>(&self, x: &Vec<T>) -> Vec<T> {
        assert_eq!(x.len(), self.repr.len());
        let mut res = self.domain.fft(x);
        for (i, _) in x.iter().enumerate() {
            res[i] *= self.v[i]
        }
        self.domain.ifft(&res)
    }
}

impl<F: FftField, D: DomainCoeff<F> + Debug> Display for Circulant<F, D> {
    fn fmt(&self, f: &mut Formatter) -> Result {
        let mut indices: Vec<i32> = (0..self.repr.len() as i32).map(|i| -i).collect();

        for _ in 0..self.repr.len() {
            for index in indices.iter_mut().take(self.repr.len()) {
                write!(f, "a{}", index)?;
                *index += 1;
                write!(f, " ")?;
            }

            writeln!(f)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod circulant_tests {
    use ark_bn254::{Fr, G1Affine};
    use ark_ec::AffineCurve;
    use ark_ff::UniformRand;
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
