pub mod circulant;
pub mod toeplitz;

use ark_std::log2;
pub use toeplitz::Toeplitz;

pub fn next_pow2(n: usize) -> usize {
    let two: u32 = 2; 
    let a: u32 = log2(n);

    if two.pow(a - 1) == n as u32 {
        return n;
    }

    return two.pow(a).try_into().unwrap()
}

#[cfg(test)]
mod tests {
    use std::iter;

    use ark_bn254::{Fr, G1Affine, G1Projective, Bn254};
    use ark_ec::{msm::VariableBaseMSM, AffineCurve, PairingEngine};
    use ark_ff::{One, PrimeField, Field, Zero, UniformRand};
    use ark_poly::{
        univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
        UVPolynomial,
    };
    use ark_std::test_rng;

    use crate::{next_pow2, toeplitz::UpperToeplitz};

    pub fn commit<E: PairingEngine>(
        srs: &[E::G1Affine],
        poly: &DensePolynomial<E::Fr>,
    ) -> E::G1Projective {
        if srs.len() - 1 < poly.degree() {
            panic!(
                "SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
                poly.degree(),
                srs.len()
            );
        }
        let coeff_scalars: Vec<_> = poly.coeffs.iter().map(|c| c.into_repr()).collect();
        VariableBaseMSM::multi_scalar_mul(srs, &coeff_scalars)
    }

    pub fn open<E: PairingEngine>(
        srs: &[E::G1Affine],
        poly: &DensePolynomial<E::Fr>,
        challenge: E::Fr,
    ) -> (E::Fr, E::G1Affine) {
        let q = poly / &DensePolynomial::from_coefficients_slice(&[-challenge, E::Fr::one()]);
        if srs.len() - 1 < q.degree() {
            panic!(
                "Open g1: SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
                q.degree(),
                srs.len()
            );
        }
        let proof: E::G1Affine = commit::<E>(srs, &q).into();
        (poly.evaluate(&challenge), proof.into())
    }

    fn commit_in_each_omega_i<E: PairingEngine>(
        srs: &[E::G1Affine],
        domain: &GeneralEvaluationDomain<E::Fr>,
        poly: &DensePolynomial<E::Fr>,
    ) -> Vec<E::G1Affine> {
        domain.elements()
            .map(|omega_pow_i| open::<E>(srs, poly, omega_pow_i).1)
            .collect()
    }

    #[test]
    fn test_multipoint_commitment() {
        let n = 64;
        let mut rng = test_rng();

        let tau = Fr::rand(&mut rng);

        let powers_of_tau: Vec<Fr> = iter::successors(Some(Fr::one()), |p| Some(*p * tau))
            .take(n)
            .collect();

        let g1_gen = G1Affine::prime_subgroup_generator();

        let srs: Vec<G1Affine> = powers_of_tau
            .iter()
            .take(n)
            .map(|tp| g1_gen.mul(tp.into_repr()).into())
            .collect();

        let mut srs_proj: Vec<G1Projective> = srs.iter().map(|t| t.into_projective()).collect();
        srs_proj.reverse();

        let poly = DensePolynomial::<Fr>::rand(n, &mut rng);
        let t = UpperToeplitz::from_poly(&poly);

        let h_commitments = t.mul_by_vec(&srs_proj)[..n].to_vec();

        let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();

        let qs_fast = domain.fft(&h_commitments);
        let qs_slow = commit_in_each_omega_i::<Bn254>(&srs, &domain, &poly);
        assert_eq!(qs_fast, qs_slow);
    }

    #[test]
    fn test_smaller_degree() {
        let n = 32;
        let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
        let mut rng = test_rng();

        let tau = Fr::rand(&mut rng);

        let powers_of_tau: Vec<Fr> = iter::successors(Some(Fr::one()), |p| Some(*p * tau))
            .take(n)
            .collect();

        let g1_gen = G1Affine::prime_subgroup_generator();

        let srs: Vec<G1Affine> = powers_of_tau
            .iter()
            .take(n)
            .map(|tp| g1_gen.mul(tp.into_repr()).into())
            .collect();

        let d = 5;
        let next_pow_2_deg = next_pow2(d);
        let mut srs_proj: Vec<G1Projective> = srs.iter().take(next_pow_2_deg).map(|t| t.into_projective()).collect();
        srs_proj.reverse();

        let poly = DensePolynomial::<Fr>::rand(d, &mut rng);
        let poly_clone = poly.clone();

        let t = UpperToeplitz::from_poly(&poly);

        let mut h_commitments = t.mul_by_vec(&srs_proj)[..next_pow_2_deg].to_vec();
        let zero_cms = vec![G1Projective::zero(); n - next_pow_2_deg];
        h_commitments.extend_from_slice(&zero_cms);

        let qs_fast = domain.fft(&h_commitments);
        let qs_slow = commit_in_each_omega_i::<Bn254>(&srs, &domain, &poly_clone);
        assert_eq!(qs_fast, qs_slow);
    }

    #[test]
    fn test_poly_correctness() {
        let mut rng = test_rng(); 
        let n = 4; 

        let p = DensePolynomial::<Fr>::rand(3, &mut rng); 
        let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();

        let omega = domain.element(1);
        let omega_squared = domain.element(2);


        let q1 = &p / &DensePolynomial::from_coefficients_slice(&[-omega, Fr::one()]);

        let q1_by_hand = DensePolynomial::from_coefficients_slice(&[
            p.coeffs[1] + omega * p.coeffs[2] + omega_squared*p.coeffs[3],
            omega * p.coeffs[3] + p.coeffs[2],
            p.coeffs[3], 
        ]);

        assert_eq!(q1, q1_by_hand);
    }

    #[test]
    fn test_h_polys() {
        let mut rng = test_rng(); 
        let n = 64; 

        let p = DensePolynomial::<Fr>::rand(n, &mut rng); 
        let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();

        let omega = domain.element(1);
        let q1 = &p / &DensePolynomial::from_coefficients_slice(&[-omega, Fr::one()]);

        let tau = Fr::from(15); 
        let q1_eval = q1.evaluate(&tau);

        let mut h_coeffs = p.coeffs()[..].to_vec();

        let hs: Vec<DensePolynomial<_>> = (0..n).map(|_| {
            h_coeffs.remove(0);
            DensePolynomial::from_coefficients_slice(&h_coeffs)
        }).collect();

        let h_sum: Fr = hs.iter().enumerate().map(|(i, hi)| {
            hi.evaluate(&tau) * omega.pow([i as u64])
        }).sum();

        assert_eq!(q1_eval, h_sum);
    }
}
