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
fn test_proof_m_zero() {
    let rng = &mut test_rng();
    let x_value = G1Affine::rand(rng);
    let x = Variable::new(rng, x_value);
    let b = G2Affine::rand(rng);
    let gamma = Matrix::<Fr>::new(&[[]]); // dim = (1, 0)

    let cks = CommitmentKeys::<F>::rand(rng);

    // Setup Proof System over Pairing Product Equation:
    // ∏e(x, b) = T
    let proof_system = setup(rng, &cks, &[], &[(x, b)], &gamma);

    let ProofSystem {
        equation,
        c,
        d,
        proof,
    } = proof_system;

    assert!(equation.verify(&cks, &c, &d, &proof));
}

#[test]
fn test_proof_n_zero() {
    let rng = &mut test_rng();
    let y_value = G2Affine::rand(rng);
    let y = Variable::new(rng, y_value);
    let a = G1Affine::rand(rng);
    let gamma = Matrix::<Fr>::zeros_column(1); // dim = (0, 1)

    let cks = CommitmentKeys::<F>::rand(rng);

    // Setup Proof System over Pairing Product Equation:
    // ∏e(a, y) = T
    let proof_system = setup(rng, &cks, &[(a, y)], &[], &gamma);

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

#[test]
fn test_homomorphic_proofs() {
    let rng = &mut test_rng();
    let cks = CommitmentKeys::<F>::rand(rng);

    // Setup a Proof System.
    let (a, b) = (G1Affine::rand(rng), G2Affine::rand(rng));
    let (x_value, y_value) = (G1Affine::rand(rng), G2Affine::rand(rng));
    let (x, y) = (
        Variable::<G1>::new(rng, x_value),
        Variable::<G2>::new(rng, y_value),
    );
    let gamma = Matrix::<Fr>::rand(rng, 1, 1);
    let proof_system_1 = setup(rng, &cks, &[(a, y)], &[(x, b)], &gamma);
    assert!(proof_system_1.equation.verify(
        &cks,
        &proof_system_1.c,
        &proof_system_1.d,
        &proof_system_1.proof
    ));

    // Setup another Proof System.
    let (a_p, b_p) = (G1Affine::rand(rng), G2Affine::rand(rng));
    let (x_p_value, y_p_value) = (G1Affine::rand(rng), G2Affine::rand(rng));
    let (x_p, y_p) = (
        Variable::<G1>::new(rng, x_p_value),
        Variable::<G2>::new(rng, y_p_value),
    );
    let gamma_p = Matrix::<Fr>::rand(rng, 1, 1);
    let proof_system_2 = setup(rng, &cks, &[(a_p, y_p)], &[(x_p, b_p)], &gamma_p);
    assert!(proof_system_2.equation.verify(
        &cks,
        &proof_system_2.c,
        &proof_system_2.d,
        &proof_system_2.proof
    ));

    // Homomorphic property of the commitment. It should be equavalent to the proof of the equation:
    // e(a, y) + a(a', y') + e(x, b) + e(x', b') + e(x, y)^gamma + e(x', y')^gamma' = T + T'
    let proof_system_sum = proof_system_1 + proof_system_2;
    assert!(proof_system_sum.equation.verify(
        &cks,
        &proof_system_sum.c,
        &proof_system_sum.d,
        &proof_system_sum.proof
    ));
}
