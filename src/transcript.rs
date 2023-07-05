use ark_ec::{pairing::Pairing};
use ark_std::{UniformRand, test_rng};

pub struct Transcript<P:Pairing>
{
    pub beta: P::ScalarField,
    pub gamma: P::ScalarField,
    pub eta: P::ScalarField,
}
impl<P:Pairing> Transcript<P>{
    pub fn initiate() -> Self{
    let mut r = test_rng();
    let beta = P::ScalarField::rand(&mut r);
    let gamma = P::ScalarField::rand(&mut r);
    let eta = P::ScalarField::rand(&mut r);
        Transcript{
            beta,
            gamma, 
            eta,
        }
    }
}