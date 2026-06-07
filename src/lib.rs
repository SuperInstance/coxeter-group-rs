//! # coxeter-group-rs
//!
//! A library for computing with Coxeter groups: Coxeter matrices,
//! reflection representations, word reduction, and Bruhat order.
//!
//! ## Modules
//!
//! - [`matrix`] — Coxeter matrix construction and operations
//! - [`reflection`] — Reflection representation (geometric realization)
//! - [`word`] — Word reduction and normal forms
//! - [`bruhat`] — Bruhat order computation
//! - [`graph`] — Coxeter graph / Dynkin diagram utilities

pub mod bruhat;
pub mod graph;
pub mod matrix;
pub mod reflection;
pub mod word;

pub use matrix::CoxeterMatrix;
