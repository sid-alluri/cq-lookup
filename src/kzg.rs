use ark_std::{test_rng, UniformRand, One};
use ark_ec::pairing::Pairing;
use ark_ec::scalar_mul::variable_base::*;
use ark_ec::Group;
#[derive(Debug,Clone)]
pub struct Srs<P: Pairing>{
    pub srs_g1: Vec<P::G1Affine>,
    pub srs_g2: Vec<P::G2Affine>
}
impl<P: Pairing> Srs<P>{
    pub fn setup(big_n: usize)->Self{
        let mut r = test_rng();
        let tau = P::ScalarField::rand(&mut r); //-+-unsafe-+-\\

        let mut srs_g1 = vec![];
        let mut srs_g2 = vec![];
        let g1 = P::G1::generator();
        let g2 = P::G2::generator();
        let mut buf = P::ScalarField::one();
        for _i in 0..big_n{
            srs_g1.push((g1*buf).into());
            srs_g2.push((g2*buf).into());
            buf = buf*tau;
        }
        srs_g2.push((g2*buf).into()); // Observation: from point 1 in protocol, till N for srsg2 but only till N-1 for srsg1
        Srs{srs_g1, srs_g2}
    }

    pub fn commit_g1(srs_g1: &Vec<P::G1Affine>, poly_coeff: &[P::ScalarField]) -> P::G1
    {
        P::G1::msm(&srs_g1[..poly_coeff.len()], poly_coeff).unwrap()
    }
    pub fn commit_g2(srs_g2: &Vec<P::G2Affine>, poly_coeff: &[P::ScalarField]) -> P::G2
    {
        P::G2::msm(&srs_g2[..poly_coeff.len()], poly_coeff).unwrap()
    }

}

#[cfg(test)]
mod tests{
    use crate::kzg::*;
    use ark_bls12_381::{Fr, Bls12_381};

    #[test]
    fn com_add(){
        let srs = Srs::<Bls12_381>::setup(8);
        let a = Fr::from(4);
        let b = Fr::from(6);
        let c = a + b;
        let com_a = Srs::<Bls12_381>::commit_g1(&srs.srs_g1, &[a]); 
        let com_b = Srs::<Bls12_381>::commit_g1(&srs.srs_g1, &[b]); 
        let com_c = Srs::<Bls12_381>::commit_g1(&srs.srs_g1, &[c]); 
        assert_eq!(com_c, com_a + com_b);
    }
}