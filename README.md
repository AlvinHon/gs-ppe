# SXDH Groth-Sahai Proofs for Pairing-Product Equations

Rust implementation of Groth-Sahai Proofs for Pairing-Product Equations defined in the paper [Commuting Signatures and Verifiable Encryption and an Application to Non-Interactively Delegatable Credentials](https://eprint.iacr.org/2010/233.pdf). The primary goal of this library is for being a building block in the `Commuting Signature` defined in the paper. It can be also used as a standalone library of a Zero Knowledge Proof.

The implementation of this work was inspired by the crate [groth-sahai-rs](https://github.com/jdwhite48/groth-sahai-rs) which provides proofs of a full set of equations. Although this work supports only pairing product equations, it allows **`randomization`** of commitment and proof.

This work uses cryptographic libraries provided by the rust ecosystem [Arkworks](https://github.com/arkworks-rs/). Similar to what they stated, this work has not been thoroughly audited for production use, please take your own risk to use it.

## Prove that I know `X` and `Y` in this equation:

```math
\prod_{j=1}^n e(A_j, Y_j)
\prod_{i=1}^m e(X_i, B_i)
\prod_{i=1}^m \prod_{j=1}^n e(X_i, Y_j)^{γ_{ij}}
= 
t_Τ
```

To do that, simply call the method like this:

```rust ignore
// Suppose `m` = 1 and `n` = 1. gamma is a matrix with dimension (1x1).
gs_ppe::setup(rng, &cks, &[(a, y)], &[(x, b)], &gamma);
```

It returns a result that contains components in the proof system:
- The specified pairing product `equation`.
- The commitments `c` and `d` which commit to the variables `x` and `y` respectively.
- The `proof` for proving `c` and `d` are committing to the variables `x` and `y` satisfying the `equation`.

```rust ignore
let cks = CommitmentKeys::<F>::rand(rng);
let ProofSystem {
    equation,
    c,
    d,
    proof,
} = setup(rng, &cks, &[(a, y)], &[(x, b)], &gamma);
assert!(equation.verify(&cks, &c, &d, &proof));
```

## Prove that you know `X` and `Y` in the equation

I only know the value of the constants `a`, `b`, `γ` and the output of the Pairing Product `t_Τ`. So I can construct the
`equation` like this:

```rust ignore
let equation = Equation::<F>::new(vec![a], vec![b], gamma, target);
```

We choose using the same commitment keys in the proof system. You then give me the `c` and `d`, the commitments to the hidden variables `X` and `Y`, together with the `proof`. Again, I can simply call the verification method like this:

```rust ignore
equation.verify(&cks, &c, &d, &proof);
```

## Randomization

To randomize the proof (as `RdProof()` stated in the paper), the commitments also need to be randomized (as `RdCom()` stated in the paper). The struct `ProofSystem` provides a single method `randomize` to do both for convenience. Once you get the result in type of `ProofSystem`, simply call the method `randomize`:

```rust ignore
let cks = CommitmentKeys::<F>::rand(rng);
let proof_system = setup(rng, &cks, &[(a, y)], &[(x, b)], &gamma);

let ProofSystem {
    equation,
    c, // a randomized `c`
    d, // a randomized  `d`
    proof, // a randomized `proof`
} = proof_system.randomize(rng, &cks);

assert!(equation.verify(&cks, &c, &d, &proof));
```

## Homomorphism 

Groth-Sahai Proofs are homomorphic. Given the equations `E` and `E'` under the same commitment key, the combined commitments `(c, d) + (c', d')` and the proofs `π + π'` give you a proof system with equation `E" = E + E'` like this:

```math
\prod_{j=1}^n e(A_j, Y_j)
\prod_{j=1}^{n'} e(A_j', Y_j')
\prod_{i=1}^m e(X_i, B_i)
\prod_{i=1}^{m'} e(X_i', B_i')
\prod_{i=1}^m \prod_{j=1}^n e(X_i, Y_j)^{γ_{ij}}
\prod_{i=1}^{m'} \prod_{j=1}^{n'} e(X_i', Y_j')^{γ'_{ij}}
= 
t_Τ + t_Τ'
```

To do that, simply add up those components. For convenience, the struct `ProofSystem` implements the trait `std::ops::Add` to add up those components at once.

```rust ignore
let cks = CommitmentKeys::<F>::rand(rng);
let proof_system_1 = setup(rng, &cks, &[(a, y)], &[(x, b)], &gamma);
let proof_system_2 = setup(rng, &cks, &[(a_p, y_p)], &[(x_p, b_p)], &gamma_p);

let  ProofSystem {
    equation, // E"
    c, // [c, c']
    d, // [d, d']
    proof, // π + π'
} = proof_system_1 + proof_system_2;
assert!(equation.verify(&cks, &c, &d, &proof));
```

## Reference:
- [Commuting Signatures and Verifiable Encryption and an Application to Non-Interactively Delegatable Credentials](https://eprint.iacr.org/2010/233.pdf), Georg Fuchsbauer
- [Commuting Signatures and Verifiable Encryption](https://www.iacr.org/archive/eurocrypt2011/66320227/66320227.pdf), Georg Fuchsbauer
- [Efficient Non-interactive Proof Systems for Bilinear Groups](https://eprint.iacr.org/archive/2007/155/1460357433.pdf), Jens Groth, Amit Sahai