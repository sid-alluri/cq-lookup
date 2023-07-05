use crate::{transcript::Transcript, preprocess::PreprocessedElements, prover::Phi};
use ark_ec::{pairing::{Pairing}, Group};
use std::ops::Mul;

pub fn verifier<P: Pairing>(transcript: &Transcript<P>, phi: Phi<P>,
                            pre: &PreprocessedElements<P>, f: Vec<P::ScalarField>)
{   
    println!("Verification");

    let beta = transcript.beta;
    let gamma = transcript.gamma;
    let eta = transcript.eta;

    //check1 
    let lhs = P::pairing(phi.com1_a, pre.com2_t_x);
    let rhs = P::pairing(phi.com1_qa, pre.com2_vanishingpoly_big_n) + 
                P::pairing(phi.com1_m - phi.com1_a*beta, P::G2::generator());
    assert_eq!(lhs, rhs);

    //check2
    let small_n = f.len();
    let lhs = P::pairing(phi.com1_b_0, pre.srs.srs_g2[(pre.big_n-1-(small_n-2))]);
    let rhs = P::pairing(phi.com_p, P::G2::generator());
    assert_eq!(lhs, rhs);

    // round 3
    let q_b_gamma = phi.q_b_gamma; // Tweak
    let mu: P::ScalarField = phi.b0_gamma + eta*phi.f_gamma + eta*eta*q_b_gamma;
    let c: P::G1 = phi.com1_b_0 + phi.com1_f.mul(eta) + phi.com1_qb.mul(eta*eta);

    // check3
    let lhs = P::pairing(c - (P::G1::generator()*mu) + phi.com1_h*gamma, P::G2::generator()); 
    let rhs = P::pairing(phi.com1_h, pre.srs.srs_g2[1]); //SURE?
    assert_eq!(lhs, rhs);

    //check4
    let lhs = P::pairing(phi.com1_a - P::G1::generator().mul(phi.a0), P::G2::generator()); 
    let rhs = P::pairing(phi.com1_a0, pre.srs.srs_g2[1]); //SURE?
    assert_eq!(lhs, rhs);

}
