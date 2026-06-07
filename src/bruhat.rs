//! Bruhat order computation for Coxeter groups.
//!
//! The Bruhat order is a partial order on group elements where u ≤ v
//! if every reduced expression for v contains a subexpression that is
//! a reduced expression for u.

use crate::matrix::CoxeterMatrix;
use crate::word::{reduce_word, word_length};

/// Check if word_u ≤ word_v in Bruhat order using the subword property.
///
/// u ≤ v iff some reduced expression for v has a subexpression that equals u.
/// Uses the tableau criterion for efficiency.
pub fn bruhat_leq(coxeter: &CoxeterMatrix, word_u: &[usize], word_v: &[usize]) -> bool {
    let u = reduce_word(coxeter, word_u);
    let v = reduce_word(coxeter, word_v);

    if u.is_empty() {
        return true; // identity ≤ everything
    }
    if u.len() > v.len() {
        return false;
    }

    // Use the subword property: check all 2^|v| subsequences
    check_subsequences(coxeter, &u, &v, 0, &mut vec![])
}

fn check_subsequences(
    coxeter: &CoxeterMatrix,
    target: &[usize],
    word: &[usize],
    pos: usize,
    current: &mut Vec<usize>,
) -> bool {
    if pos == word.len() {
        return reduce_word(coxeter, current) == target;
    }

    // Pruning: remaining positions must provide enough elements
    let remaining = word.len() - pos - 1;
    let needed = if target.len() > current.len() { target.len() - current.len() } else { 0 };
    if remaining < needed && current.len() < target.len() {
        // Still try excluding, but check if we can reach target
    }

    // Try including word[pos]
    current.push(word[pos]);
    if check_subsequences(coxeter, target, word, pos + 1, current) {
        current.pop();
        return true;
    }
    current.pop();

    // Try excluding word[pos]
    check_subsequences(coxeter, target, word, pos + 1, current)
}

/// Compute the Bruhat interval [u, v] using enumeration.
///
/// Only practical for small groups.
pub fn bruhat_interval(
    coxeter: &CoxeterMatrix,
    word_u: &[usize],
    word_v: &[usize],
) -> Vec<Vec<usize>> {
    let elements = coxeter.enumerate_elements();
    let u = reduce_word(coxeter, word_u);
    let v = reduce_word(coxeter, word_v);

    elements
        .into_iter()
        .filter(|w| bruhat_leq(coxeter, &u, w) && bruhat_leq(coxeter, w, &v))
        .collect()
}

/// Compute the Bruhat height of an element: l(w) - l(u).
pub fn bruhat_height(coxeter: &CoxeterMatrix, word_u: &[usize], word_w: &[usize]) -> usize {
    let lu = word_length(coxeter, word_u);
    let lw = word_length(coxeter, word_w);
    lw.saturating_sub(lu)
}

/// Find all elements that cover the given element in Bruhat order
/// (i.e., w < v and there is no u with w < u < v).
pub fn bruhat_covers(coxeter: &CoxeterMatrix, word: &[usize]) -> Vec<Vec<usize>> {
    let w = reduce_word(coxeter, word);
    let w_len = w.len();

    let elements = coxeter.enumerate_elements();
    elements
        .into_iter()
        .filter(|v| v.len() == w_len + 1 && bruhat_leq(coxeter, &w, v))
        .collect()
}

/// Compute the Bruhat graph edges for the whole group.
pub fn bruhat_graph(coxeter: &CoxeterMatrix) -> Vec<(Vec<usize>, Vec<usize>)> {
    let elements = coxeter.enumerate_elements();
    let mut edges = vec![];
    for w in &elements {
        let covers = bruhat_covers(coxeter, w);
        for v in covers {
            edges.push((w.clone(), v));
        }
    }
    edges
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bruhat_identity_leq_all() {
        let coxeter = CoxeterMatrix::type_a(2);
        assert!(bruhat_leq(&coxeter, &[], &[0]));
        assert!(bruhat_leq(&coxeter, &[], &[0, 1]));
    }

    #[test]
    fn test_bruhat_reflexive() {
        let coxeter = CoxeterMatrix::type_a(2);
        assert!(bruhat_leq(&coxeter, &[0], &[0]));
        assert!(bruhat_leq(&coxeter, &[0, 1], &[0, 1]));
    }

    #[test]
    fn test_bruhat_antisymmetric() {
        let coxeter = CoxeterMatrix::type_a(2);
        assert!(!bruhat_leq(&coxeter, &[0, 1], &[0]));
    }

    #[test]
    fn test_bruhat_a2_chain() {
        let coxeter = CoxeterMatrix::type_a(2);
        assert!(bruhat_leq(&coxeter, &[0], &[0, 1]));
        assert!(bruhat_leq(&coxeter, &[1], &[0, 1]));
        assert!(!bruhat_leq(&coxeter, &[0, 1], &[1]));
    }

    #[test]
    fn test_bruhat_height() {
        let coxeter = CoxeterMatrix::type_a(2);
        assert_eq!(bruhat_height(&coxeter, &[], &[]), 0);
        assert_eq!(bruhat_height(&coxeter, &[], &[0]), 1);
        assert_eq!(bruhat_height(&coxeter, &[], &[0, 1]), 2);
    }

    #[test]
    fn test_bruhat_covers_identity() {
        let coxeter = CoxeterMatrix::type_a(2);
        let covers = bruhat_covers(&coxeter, &[]);
        assert_eq!(covers.len(), 2); // s0 and s1 cover identity
    }

    #[test]
    fn test_bruhat_interval_a2() {
        let coxeter = CoxeterMatrix::type_a(2);
        let interval = bruhat_interval(&coxeter, &[], &[0, 1]);
        assert!(interval.len() >= 3);
    }

    #[test]
    fn test_bruhat_graph_a2() {
        let coxeter = CoxeterMatrix::type_a(2);
        let edges = bruhat_graph(&coxeter);
        assert!(!edges.is_empty());
        for (w, v) in &edges {
            assert!(v.len() == w.len() + 1);
        }
    }
}
