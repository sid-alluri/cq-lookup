use crate::{transcript::Transcript, preprocess::PreprocessedElements, kzg::Srs};
use ark_ec::{pairing::Pairing};
use ark_poly::{univariate::DensePolynomial, GeneralEvaluationDomain, EvaluationDomain, DenseUVPolynomial, Polynomial};
use ark_std:: {Zero, One};
use ark_ff::Field;
use std::{collections::BTreeMap, ops::Mul, iter};

pub struct Phi<P:Pairing>{
    pub com1_m: P::G1,
    pub com1_f: P::G1,
    pub com1_a: P::G1,
    pub com1_qa: P::G1,
    pub com1_b_0: P::G1,
    pub com1_qb: P::G1,
    pub com_p: P::G1,
    pub a0: P::ScalarField,
    pub b0_gamma: P::ScalarField,
    pub q_b_gamma: P::ScalarField,
    pub f_gamma: P::ScalarField,
    pub com1_h: P::G1,
    pub com1_a0: P::G1,
 }

impl<P:Pairing> Phi<P>{
    pub fn prover(transcript: &Transcript<P>,pre: &PreprocessedElements<P>, 
                                t: Vec<P::ScalarField>, f: Vec<P::ScalarField>) -> Self
{    // first round
    println!("Proof Generation");

    let m = get_m::<P>(t.clone(),f.clone());
    let fftdomain = GeneralEvaluationDomain::<P::ScalarField>::new(f.clone().len()).unwrap();
    let mut com1_m = P::G1::zero(); // return
    for (i, mi) in m.iter(){
        com1_m += pre.com1_li_x.get(&(*i as u64)).unwrap().mul(P::ScalarField::from(*mi));
    }
    let f_coeff = fftdomain.fft(&f[..]);
    let com1_f = Srs::<P>::commit_g1(&pre.srs.srs_g1, &f_coeff[..]); // return
    
    // second round

    let beta = transcript.beta.clone();
    let mut a_eval = vec![]; 
    let mut com1_a: <P as Pairing>::G1 = P::G1::zero(); //return
    let mut com1_qa = P::G1::zero(); //return
    for (i, mi) in m.iter(){
        let mut a_i = P::ScalarField::from(*mi);
        a_i /= t[*i].clone()+beta;
        a_eval.push(a_i);
        com1_a += pre.com1_li_x.get(&(*i as u64)).unwrap().mul(a_i);
        com1_qa += pre.com1_qi.get(&(*i as u64)).unwrap().mul(a_i);
    }

    let a_coeff = fftdomain.ifft(&a_eval[..]);
    let mut b_eval = vec![];
    for f_i in f.clone(){
        b_eval.push((f_i+beta).inverse().unwrap());
    }

    let b_coeff = fftdomain.fft(&b_eval[..]);
    let com1_b_0  = Srs::<P>::commit_g1(&pre.srs.srs_g1, &b_coeff[1..]); // return
    let mut qb_eval: Vec<P::ScalarField> = vec![];
    let binding = DensePolynomial::from(fftdomain.vanishing_polynomial());
    let zv_coeff = binding.coeffs();
    let zv_eval = fftdomain.ifft(zv_coeff);
    for i in 0..b_eval.len(){
        let buf = ((b_eval[i]*(f[i]+beta))-P::ScalarField::one())*zv_eval[i].inverse().unwrap();
        qb_eval.push(buf);
    }

    let qb_coeff = fftdomain.fft(&qb_eval[..]);
    let com1_qb = Srs::<P>::commit_g1(&pre.srs.srs_g1, &qb_coeff[..]); // return

    let mut shifted_coeffs = vec![P::ScalarField::zero(); t.len() - 1 - (f.len() - 2)];
    shifted_coeffs.extend_from_slice(&b_coeff[1..]); //b0(x)
    let com_p = Srs::<P>::commit_g1(&pre.srs.srs_g1, &shifted_coeffs[..]); // return
  

    // third round

    let gamma = transcript.gamma.clone();
    let eta = transcript.eta.clone();
    let b0_x = DensePolynomial::from_coefficients_slice(&b_coeff[1..]);
    let f_x = DensePolynomial::from_coefficients_slice(&f_coeff[..]);
    let b0_gamma = b0_x.evaluate(&gamma);
    let f_gamma = f_x.evaluate(&gamma);

    let mut a0 = P::ScalarField::zero();
    for (i, m_i) in m.iter() {
        a0 += P::ScalarField::from(*m_i) / &(t[*i as usize] + beta) * pre.li_x.get(&(*i as u64)).unwrap()[0];
    }

    let q_b_x = DensePolynomial::from_coefficients_slice(&qb_coeff[..]); //Tweak
    let q_b_gamma = q_b_x.evaluate(&gamma);

    let com1_h = batch_open_g1::<P>( 
        &pre.srs.srs_g1,
        &[b0_x, f_x,q_b_x],
        gamma,
        eta,
    );

    
    let com1_a0 = Srs::<P>::commit_g1(&pre.srs.srs_g1, &a_coeff[1..]);  //return
 
    Phi{
        com1_m,
        com1_f,
        com1_a,
        com1_qa,
        com1_b_0,
        com1_qb,
        com_p,
        a0,
        b0_gamma,
        q_b_gamma,
        f_gamma,
        com1_h,
        com1_a0,
    }
}

}


pub fn get_m<P: Pairing>(table:Vec<P::ScalarField>, input: Vec<P::ScalarField>) -> BTreeMap<usize, u64>
{
    let mut table_map = BTreeMap::new();
    for i in 0..table.len(){
        let _ti = table_map.insert(table[i], i);
    }
    let mut m_map = BTreeMap::new();
    for j in &input{
        let i = table_map.get(j).ok_or("Error, can't do it").unwrap();
        let mi = m_map.entry(*i).or_insert(0 as u64);
        *mi += 1;
    }
    m_map
}


// The below function is from https://github.com/geometryresearch/cq/blob/c0e499cdf866631b5079a2ae6837e26df784d0eb/src/kzg.rs#L54

pub fn batch_open_g1<P: Pairing>(
    srs: &Vec<P::G1Affine>,
    polys: &[DensePolynomial<P::ScalarField>],
    opening_challenge: P::ScalarField,
    separation_challenge: P::ScalarField,
) -> P::G1 

{
    let powers_of_gamma = iter::successors(Some(separation_challenge), 
                        |p| {Some(*p * separation_challenge)});

    let mut batched = polys[0].clone();
    for (p_i, gamma_pow_i) in polys.iter().skip(1).zip(powers_of_gamma) 
    {
        batched += (gamma_pow_i, p_i);
    }

    let q = &batched
        / &DensePolynomial::from_coefficients_slice(&[-opening_challenge, P::ScalarField::one()]);

    if srs.len() - 1 < q.degree() {
        panic!(
            "Batch open g1: SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
            q.degree(),
            srs.len()
        );
    }

    Srs::<P>::commit_g1(srs, &q)
}



#[cfg(test)]
mod tests{
    use ark_bls12_381::{Bls12_381, Fr};

    use super::get_m;
    
    #[test]
    fn test_get_m(){
        let t = vec![1 as u64,2,3,4,5,6,7,8];
        let tf: Vec<Fr> = t.iter().map(|ti| Fr::from(*ti)).collect();
        let f = vec![2 as u64,2,2,2,4,4,3];
        let ff: Vec<Fr> = f.iter().map(|ti| Fr::from(*ti)).collect();
        let m  = get_m::<Bls12_381>(tf, ff);
        assert_eq!(*m.get(&1).unwrap(), 4);
        assert_eq!(*m.get(&3).unwrap(), 2);
        assert_eq!(*m.get(&2).unwrap(), 1);
    }
}