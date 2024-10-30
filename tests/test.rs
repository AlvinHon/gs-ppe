use ark_bls12_381::Bls12_381 as F;
use ark_ec::pairing::Pairing;
use ark_std::{test_rng, UniformRand};

use gs_ppe::{setup, CommitmentKeys, Matrix, ProofSystem, Variable};

type G1 = <F as Pairing>::G1;
type G2 = <F as Pairing>::G2;
type G1Affine = <F as Pairing>::G1Affine;
type G2Affine = <F as Pairing>::G2Affine;
type Fr = <F as Pairing>::ScalarField;

#[test]
fn test_proof() {
    let rng = &mut test_rng();
    let a = G1Affine::rand(rng);
    let b = G2Affine::rand(rng);

    let x_value = G1Affine::rand(rng);
    let x = Variable::<G1>::new(rng, x_value);
    let y_value = G2Affine::rand(rng);
    let y = Variable::<G2>::new(rng, y_value);
    let gamma = Matrix::<Fr>::rand(rng, 1, 1);

    let cks = CommitmentKeys::<F>::rand(rng);

    // Setup Proof System over Pairing Product Equation:
    // e(a, y) + e(x, b) + e(x, y)^gamma = T
    let proof_system = setup(rng, &cks, &[(a, y)], &[(x, b)], &gamma);

    let ProofSystem {
        equation,
        c,
        d,
        proof,
    } = proof_system;

    assert!(equation.verify(&cks.u, &cks.v, &c, &d, &proof));
}

#[test]
fn test_randomized_proof() {
    let rng = &mut test_rng();
    let a = G1Affine::rand(rng);
    let b = G2Affine::rand(rng);

    let x_value = G1Affine::rand(rng);
    let x = Variable::<G1>::new(rng, x_value);
    let y_value = G2Affine::rand(rng);
    let y = Variable::<G2>::new(rng, y_value);
    let gamma = Matrix::<Fr>::rand(rng, 1, 1);

    let cks = CommitmentKeys::<F>::rand(rng);

    // Setup Proof System over Pairing Product Equation:
    // e(a, y) + e(x, b) + e(x, y)^gamma = T
    let proof_system = setup(rng, &cks, &[(a, y)], &[(x, b)], &gamma);

    let ProofSystem {
        equation,
        mut c,
        mut d,
        mut proof,
    } = proof_system;

    let cr = c
        .iter_mut()
        .map(|c_i| c_i.randomize(rng, &cks.u))
        .collect::<Vec<_>>();
    let ds = d
        .iter_mut()
        .map(|d_j| d_j.randomize(rng, &cks.v))
        .collect::<Vec<_>>();

    proof.randomize(rng, &cks, &equation, &cr, &ds);

    assert!(equation.verify(&cks.u, &cks.v, &c, &d, &proof));
}
