use ark_bls12_381::Bls12_381 as F;
use ark_ec::{
    pairing::{Pairing, PairingOutput},
    CurveGroup,
};
use ark_std::{rand::Rng, test_rng, UniformRand, Zero};
use criterion::{criterion_group, criterion_main, Criterion};
use gs_ppe::{CommitmentKeys, Equation, Matrix, Proof, Variable};
use std::ops::Mul;
use std::time::Duration;

type G1 = <F as Pairing>::G1;
type G2 = <F as Pairing>::G2;
type Fr = <F as Pairing>::ScalarField;

criterion_group! {
    name = gs_ppe;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(10));
    targets = bench_commit_g1, bench_commit_g2, bench_prove, bench_verify
}

criterion_main!(gs_ppe);

fn bench_commit_g1(c: &mut Criterion) {
    let rng = &mut test_rng();

    let mut group = c.benchmark_group("bench_commit_g1");

    for size in [5, 10, 20] {
        let cks = CommitmentKeys::<F>::rand(rng);
        let x_vec = x_variable_vec(rng, size);
        group.bench_with_input(format!("size: {}", size), &x_vec, |b, x_vec| {
            b.iter(|| {
                for x in x_vec {
                    cks.u.commit(&x);
                }
            })
        });
    }
}

fn bench_commit_g2(c: &mut Criterion) {
    let rng = &mut test_rng();

    let mut group = c.benchmark_group("bench_commit_g2");

    for size in [5, 10, 20] {
        let cks = CommitmentKeys::<F>::rand(rng);
        let y_vec = y_variable_vec(rng, size);
        group.bench_with_input(format!("size: {}", size), &y_vec, |b, y_vec| {
            b.iter(|| {
                for y in y_vec {
                    cks.v.commit(&y);
                }
            })
        });
    }
}

fn bench_prove(c: &mut Criterion) {
    let rng = &mut test_rng();

    let mut group = c.benchmark_group("bench_prove");

    for size in [5, 10, 20] {
        let cks = CommitmentKeys::<F>::rand(rng);
        let (equation, x, y) = prepare_prove(rng, size, size);
        group.bench_with_input(
            format!("size: {}", size),
            &(cks, equation, x, y),
            |b, (cks, equation, x, y)| {
                b.iter(|| {
                    Proof::new(rng, cks, equation, x, y);
                })
            },
        );
    }
}

fn bench_verify(c: &mut Criterion) {
    let rng = &mut test_rng();

    let mut group = c.benchmark_group("bench_verify");

    for size in [5, 10, 20] {
        let cks = CommitmentKeys::<F>::rand(rng);
        let (equation, x, y) = prepare_prove(rng, size, size);

        let c = x.iter().map(|x_i| cks.u.commit(x_i)).collect::<Vec<_>>();
        let d = y.iter().map(|y_i| cks.v.commit(y_i)).collect::<Vec<_>>();
        let proof = Proof::new(rng, &cks, &equation, &x, &y);

        group.bench_with_input(
            format!("size: {}", size),
            &(cks, equation, c, d, proof),
            |b, (cks, equation, c, d, proof)| {
                b.iter(|| {
                    equation.verify(cks, c, d, proof);
                })
            },
        );
    }
}

// ... utility functions ...

/// Returns a vector of `size` random `Variable<G1>`.
fn x_variable_vec(rng: &mut impl Rng, size: usize) -> Vec<Variable<G1>> {
    (0..size)
        .map(|_| {
            let value = G1::rand(rng);
            Variable::new(rng, value.into_affine())
        })
        .collect()
}

/// Returns a vector of `size` random `Variable<G2>`.
fn y_variable_vec(rng: &mut impl Rng, size: usize) -> Vec<Variable<G2>> {
    (0..size)
        .map(|_| {
            let value = G2::rand(rng);
            Variable::new(rng, value.into_affine())
        })
        .collect()
}

/// Prepare required input arguments to `prove` method, i.e. returns a random equation and
/// its corresponding `x` and `y` variables.
fn prepare_prove(
    rng: &mut impl Rng,
    m: usize,
    n: usize,
) -> (Equation<F>, Vec<Variable<G1>>, Vec<Variable<G2>>) {
    let gamma = Matrix::<Fr>::rand(rng, m, n);
    let x = x_variable_vec(rng, m);
    let y = y_variable_vec(rng, n);
    let a = (0..n)
        .map(|_| G1::rand(rng).into_affine())
        .collect::<Vec<_>>();
    let b = (0..m)
        .map(|_| G2::rand(rng).into_affine())
        .collect::<Vec<_>>();

    let ay_product = a
        .iter()
        .zip(y.iter())
        .fold(PairingOutput::zero(), |acc, (a, y)| {
            acc + F::pairing(a, y.value)
        });
    let xb_product = x
        .iter()
        .zip(b.iter())
        .fold(PairingOutput::zero(), |acc, (x, b)| {
            acc + F::pairing(x.value, b)
        });

    let mut xy_product = PairingOutput::zero();
    for (j, y_j) in y.iter().enumerate() {
        for (i, x_i) in x.iter().enumerate() {
            xy_product += F::pairing(x_i.value, y_j.value).mul(gamma[(i, j)]);
        }
    }
    let target = ay_product + xb_product + xy_product;
    let equation = Equation::<F>::new(a, b, gamma, target);

    (equation, x, y)
}
