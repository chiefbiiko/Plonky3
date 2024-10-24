// use std::{fmt::Display, marker::PhantomData, str::FromStr};
use core::{fmt::Display, marker::PhantomData, str::FromStr};

// use ark_crypto_primitives::merkle_tree::{Config, LeafParam, TwoToOneParam};
use serde::Serialize;

use p3_symmetric::{SerializingHasher32, CryptographicHasher, CompressionFunctionFromHasher, PseudoCompressionFunction};

//from https://github.com/arkworks-rs/crypto-primitives/blob/6a770ebf87c13d2e0aef0237d1f4669d94eb58a7/crypto-primitives/src/merkle_tree/mod.rs#L124
// pub type TwoToOneParam<P> = <<P as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters;
// pub type LeafParam<P> = <<P as Config>::LeafHash as CRHScheme>::Parameters;

//from https://github.com/WizardOfMenlo/whir/blob/901f92eee95d46f3a7fbe1a36eecc582812eeccf/src/crypto/merkle_tree/keccak.rs#L42C1-L43C36
// pub struct KeccakLeafHash<F>(PhantomData<F>);
// pub struct KeccakTwoToOneCRHScheme;

pub fn default_max_pow(num_variables: usize, log_inv_rate: usize) -> usize {
    num_variables + log_inv_rate - 3
}

#[derive(Debug, Clone, Copy, Serialize)]
pub enum SoundnessType {
    UniqueDecoding,
    ProvableList,
    ConjectureList,
}

impl Display for SoundnessType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match &self {
                SoundnessType::ProvableList => "ProvableList",
                SoundnessType::ConjectureList => "ConjectureList",
                SoundnessType::UniqueDecoding => "UniqueDecoding",
            }
        )
    }
}

impl FromStr for SoundnessType {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s == "ProvableList" {
            Ok(SoundnessType::ProvableList)
        } else if s == "ConjectureList" {
            Ok(SoundnessType::ConjectureList)
        } else if s == "UniqueDecoding" {
            Ok(SoundnessType::UniqueDecoding)
        } else {
            Err(format!("Invalid soundness specification: {}", s))
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct MultivariateParameters<F> {
    pub(crate) num_variables: usize,
    _field: PhantomData<F>,
}

impl<F> MultivariateParameters<F> {
    pub fn new(num_variables: usize) -> Self {
        Self {
            num_variables,
            _field: PhantomData,
        }
    }
}

impl<F> Display for MultivariateParameters<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Number of variables: {}", self.num_variables)
    }
}

#[derive(Debug, Clone, Copy)]
pub enum FoldType {
    Naive,
    ProverHelps,
}

impl FromStr for FoldType {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s == "Naive" {
            Ok(FoldType::Naive)
        } else if s == "ProverHelps" {
            Ok(FoldType::ProverHelps)
        } else {
            Err(format!("Invalid fold type specification: {}", s))
        }
    }
}

impl Display for FoldType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                FoldType::Naive => "Naive",
                FoldType::ProverHelps => "ProverHelps",
            }
        )
    }
}

#[derive(Clone)]
pub struct WhirParameters<H, C, PowStrategy>
where
    // MerkleConfig: Config,
    H: CryptographicHasher<P::Value, [PW::Value; DIGEST_ELEMS]>,
    H: CryptographicHasher<P, [PW; DIGEST_ELEMS]>,
    H: Sync,
    C: PseudoCompressionFunction<[PW::Value; DIGEST_ELEMS], 2>,
    C: PseudoCompressionFunction<[PW; DIGEST_ELEMS], 2>,
    C: Sync,
{
    pub starting_log_inv_rate: usize,
    pub folding_factor: usize,
    pub soundness_type: SoundnessType,
    pub security_level: usize,
    pub pow_bits: usize,

    pub fold_optimisation: FoldType,

    // PoW parameters
    pub _pow_parameters: PhantomData<PowStrategy>,

    // Merkle tree parameters
    // pub leaf_hash_params: LeafParam<MerkleConfig>,
    // pub two_to_one_params: TwoToOneParam<MerkleConfig>,
    // MMCS
    pub mmcs_h: SerializingHasher32<CryptographicHasher>,
    pub mmcs_c: CompressionFunctionFromHasher<CryptographicHashe, usize, usize>,
}

impl<H, C, PowStrategy> Display for WhirParameters<H, C, PowStrategy>
where
    // MerkleConfig: Config,
    H: CryptographicHasher<P::Value, [PW::Value; DIGEST_ELEMS]>,
    H: CryptographicHasher<P, [PW; DIGEST_ELEMS]>,
    H: Sync,
    C: PseudoCompressionFunction<[PW::Value; DIGEST_ELEMS], 2>,
    C: PseudoCompressionFunction<[PW; DIGEST_ELEMS], 2>,
    C: Sync,
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        writeln!(
            f,
            "Targeting {}-bits of security with {}-bits of PoW - soundness: {:?}",
            self.security_level, self.pow_bits, self.soundness_type
        )?;
        writeln!(
            f,
            "Starting rate: 2^-{}, folding_factor: {}, fold_opt_type: {}",
            self.starting_log_inv_rate, self.folding_factor, self.fold_optimisation,
        )
    }
}