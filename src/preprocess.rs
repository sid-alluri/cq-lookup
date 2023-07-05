use crate::kzg::*;
use std::ops::Mul;
use ark_ec::{AffineRepr, Group, pairing::Pairing};
use ark_ff::Field;
use ark_poly::{univariate::DensePolynomial, GeneralEvaluationDomain, EvaluationDomain, DenseUVPolynomial,Polynomial};
use std::collections::BTreeMap; // Preserves the order, ordered hash map
use ark_std::{Zero, One};
use fk::UpperToeplitz;

#[derive(Debug,Clone)]

pub struct PreprocessedElements<P:Pairing>{
    pub srs: Srs<P>,
    pub big_n: usize,
    pub com2_vanishingpoly_big_n: P::G2, 
    pub com2_t_x: P::G2, 
    pub li_x : BTreeMap<u64, Vec<P::ScalarField>>,
    pub com1_li_x: BTreeMap<u64, P::G1>, 
    pub com1_li_0: BTreeMap<u64, P::G1>,
    pub com1_qi: BTreeMap<u64, P::G1>
}

impl<P: Pairing> PreprocessedElements<P>{
    
    pub fn gen(big_n: usize, t: Vec<P::ScalarField>) -> Self{    //-+-gen(N,t)-+-\\
       
        println!("Preprocessing");
        let srs = Srs::<P>::setup(big_n); 
        let fftdomain = GeneralEvaluationDomain::<P::ScalarField>::new(big_n).unwrap();
        
        // Z(x) = X^N-1, Com2_Z = g2*µ^N-g2
        let binding = DensePolynomial::from(fftdomain.vanishing_polynomial());
        let zv = binding.coeffs();
        // println!("Len {:?}", srs.srs_g1.len());
        let com2_vanishingpoly_big_n: P::G2 = srs.srs_g2[big_n].into_group() - P::G2::generator(); // return 
        assert_eq!(com2_vanishingpoly_big_n, Srs::<P>::commit_g2(&srs.srs_g2, &zv)); // in-line test
    
        // t: eval of T(x); convert to coeff form; then commit
        let t_coeff  = fftdomain.ifft(&t[..]);
        let com2_t_x = Srs::<P>::commit_g2(&srs.srs_g2, &t_coeff); // return
        
        // Lagrange Basis, L_i(W^j) = 1 when i=j else 0
        // i.e Given roots of unity, the evals of L_i is {1 at i = j else 0}
        let mut li_x = BTreeMap::new();
        let mut com1_li_x = BTreeMap::new();
        let mut com1_li_0 = BTreeMap::new();
        for i in 0..big_n{
            let mut li_eval = vec![P::ScalarField::zero() ; big_n];
            li_eval[i] = P::ScalarField::one();
            li_x.insert(i as u64, fftdomain.ifft(&li_eval[..])); // return
            let buf = li_x.get(&(i as u64)).unwrap();
            let li_x_poly = DensePolynomial::from_coefficients_slice(&buf[..]);
            assert_eq!(li_x_poly.evaluate(&fftdomain.element(i)), P::ScalarField::one());   // in-line test
            com1_li_x.insert(i as u64, Srs::<P>::commit_g1(&srs.srs_g1, &li_x.get(&(i as u64)).unwrap()[..])); // return
            com1_li_0.insert(i as u64, Srs::<P>::commit_g1(&srs.srs_g1, &li_x.get(&(i as u64)).unwrap()[1..])); // return
            // Explanation:
            // Suppose B(x) = (b_0 + b_1*x + ...); B(0) = b_0 + 0 + )......
            // Now, B(x)-B(0)/x = b_1 + b_2*x + .... | Hence considering coeff starting from index = 1 and ignoting b_0
        }
    
        // Based on Th-3.1 in cq and [FK23]
        let t_x = fftdomain.ifft(&t[..]);
        let tx_poly = DensePolynomial::from_coefficients_slice(&t_x[..]);
        let t_matrix = UpperToeplitz::from_poly(&tx_poly); // Ref: Implementation of fk by geometry research
        let mut srs1rev: Vec<P::G1> = srs.srs_g1.clone().iter().map(|a| a.into_group()).collect();
        srs1rev.reverse();
        let k = t_matrix.mul_by_vec(&srs1rev[..]); 
        let kdash = fftdomain.fft(&k[..big_n]); // only first N are material
    
        let domain_size = fftdomain.size();
        assert_eq!(domain_size, big_n);
        let n_inv = P::ScalarField::from(big_n as u64).inverse().unwrap();
        let scale = fftdomain.elements().map(|w| w*n_inv);
        let buf: Vec<P::G1> = kdash.iter().zip(scale).map(|(k,s)| k.mul(s)).collect();
        
        let mut com1_qi: BTreeMap<u64, P::G1> = BTreeMap::new(); // return
        for i in 0..big_n{
            com1_qi.insert(i as u64, buf[i]);
        }
    
        return PreprocessedElements 
                        {   srs,
                            big_n,
                            com2_vanishingpoly_big_n, 
                            com2_t_x, 
                            li_x, 
                            com1_li_x, 
                            com1_li_0, 
                            com1_qi };
    
    }

}


// #[cfg(test)]
// mod tests{
//     use crate::kzg::*;
//     use ark_bls12_381::{Fr, Bls12_381};

//     use super::PreprocessedElements;

//     #[test]
//     fn test_gen(){
//         let _srs = Srs::<Bls12_381>::setup(8);
//         let t = vec![1,2,3,4 as u64,5,6,7,8];
//         let tf: Vec<Fr> = t.iter().map(|ti| Fr::from(*ti)).collect();
//         let _pre = PreprocessedElements::<Bls12_381>::gen(8, tf);
//     }
// }
