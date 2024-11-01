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
fn test_proof_m_x_n() {
    let rng = &mut test_rng();
    let n = 3;
    let m = 2;
    let ay = (0..n)
        .map(|_| {
            let value = G2Affine::rand(rng);
            (G1Affine::rand(rng), Variable::<G2>::new(rng, value))
        })
        .collect::<Vec<_>>();
    let xb = (0..m)
        .map(|_| {
            let value = G1Affine::rand(rng);
            (Variable::<G1>::new(rng, value), G2Affine::rand(rng))
        })
        .collect::<Vec<_>>();

    let gamma = Matrix::<Fr>::rand(rng, m, n);

    let cks = CommitmentKeys::<F>::rand(rng);

    // Setup Proof System over Pairing Product Equation:
    // ∏e(a, y) + ∏e(x, b) + ∏e(x, y)^gamma = T
    let proof_system = setup(rng, &cks, &ay, &xb, &gamma);

    let ProofSystem {
        equation,
        c,
        d,
        proof,
    } = proof_system;

    assert!(equation.verify(&cks, &c, &d, &proof));
}

#[test]
fn test_randomized_proof_m_x_n() {
    let rng = &mut test_rng();
    let n = 3;
    let m = 2;
    let ay = (0..n)
        .map(|_| {
            let value = G2Affine::rand(rng);
            (G1Affine::rand(rng), Variable::<G2>::new(rng, value))
        })
        .collect::<Vec<_>>();
    let xb = (0..m)
        .map(|_| {
            let value = G1Affine::rand(rng);
            (Variable::<G1>::new(rng, value), G2Affine::rand(rng))
        })
        .collect::<Vec<_>>();

    let gamma = Matrix::<Fr>::rand(rng, m, n);

    let cks = CommitmentKeys::<F>::rand(rng);

    // Setup Proof System over Pairing Product Equation:
    // ∏e(a, y) + ∏e(x, b) + ∏e(x, y)^gamma = T
    let proof_system = setup(rng, &cks, &ay, &xb, &gamma);

    let ProofSystem {
        equation,
        c,
        d,
        proof,
    } = proof_system.randomize(rng, &cks);

    assert!(equation.verify(&cks, &c, &d, &proof));
}
