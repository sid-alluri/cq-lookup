pub mod kzg;
pub mod preprocess;
pub mod transcript;
pub mod prover;
pub mod verifier;
use ark_bls12_381::{Bls12_381, Fr};
use preprocess::PreprocessedElements;
use prover::Phi;
use transcript::Transcript;

pub fn cq(){
    let t = vec![1 as u64,2,3,4,5,6,7,8];
    let tf: Vec<Fr> = t.iter().map(|ti| Fr::from(*ti)).collect();
    let f = t.clone();
    let ff: Vec<Fr> = f.iter().map(|fi| Fr::from(*fi)).collect();
    let pre = PreprocessedElements::<Bls12_381>::gen(tf.len(), tf.clone());
    let transcript = Transcript::<Bls12_381>::initiate();
    let phi = Phi::<Bls12_381>::prover(&transcript, &pre, tf.clone(), ff.clone());
    verifier::verifier(&transcript, phi, &pre, ff); 
}

pub fn cqfalse()
{
    let t = vec![1 as u64,2,3,4,5,6,7,8];
    let tf: Vec<Fr> = t.iter().map(|ti| Fr::from(*ti)).collect();
    let f = t.clone();
    let f_false = vec![1 as u64,2,3,1];
    let ff: Vec<Fr> = f.iter().map(|fi| Fr::from(*fi)).collect();
    let ff_false: Vec<Fr> = f_false.iter().map(|fi| Fr::from(*fi)).collect(); // to check soundness
    let pre = PreprocessedElements::<Bls12_381>::gen(tf.len(), tf.clone());
    let transcript = Transcript::<Bls12_381>::initiate();
    let phi = Phi::<Bls12_381>::prover(&transcript, &pre, tf.clone(), ff.clone());
    verifier::verifier(&transcript, phi, &pre, ff_false); 
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn protocol() {
        cq();
        // cqfalse(); // It should fail
    }
}
