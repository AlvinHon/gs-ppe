pub mod arith;
pub use arith::matrix::Matrix;

pub mod com;
pub use com::Com;

pub mod commit;
pub use commit::CommitmentKeys;

pub mod equation;
pub use equation::Equation;

pub mod extract;
pub use extract::ExtractionKey;

pub mod prove;
pub use prove::Proof;

pub mod variable;
pub use variable::Variable;

use ark_ec::pairing::Pairing;

pub fn setup<E: Pairing>(
    cks: &CommitmentKeys<E>,
    ay: &[(<E as Pairing>::G1Affine, Variable<<E as Pairing>::G1>)],
    xb: &[(Variable<<E as Pairing>::G2>, <E as Pairing>::G2Affine)],
    gamma: &Matrix<E::ScalarField>,
) -> ProofSystem<E> {
    todo!()
}

pub struct ProofSystem<E: Pairing> {
    pub equation: Equation<E>,
    pub c: Vec<Com<<E as Pairing>::G1>>,
    pub d: Vec<Com<<E as Pairing>::G2>>,
    pub proof: Proof<E>,
}
