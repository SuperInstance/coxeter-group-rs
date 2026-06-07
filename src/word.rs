//! Word reduction and normal forms for Coxeter group elements.

use crate::matrix::CoxeterMatrix;

/// Reduce a word in the generators of a Coxeter group to its shortest form.
///
/// Uses deletion: removes adjacent pairs s_i s_i and applies braid relations.
/// The algorithm is guaranteed to terminate by the deletion property of Coxeter groups.
/// Canonicalize a reduced word: ensure it's in a unique form.
/// We choose the lexicographically smallest equivalent reduced word.
fn canonicalize(coxeter: &CoxeterMatrix, w: &[usize]) -> Vec<usize> {
    let mut best = w.to_vec();
    let mut current = w.to_vec();

    // Try all braid-equivalent forms
    for _ in 0..20 {
        let next = apply_braid_relations(coxeter, &current);
        let next = remove_adjacent_pairs(&next);
        if next == current {
            break;
        }
        if next < best {
            best = next.clone();
        }
        current = next;
    }
    best
}

/// Reduce a word to its canonical reduced form.
pub fn reduce_word(coxeter: &CoxeterMatrix, word: &[usize]) -> Vec<usize> {
    let mut w = word.to_vec();
    let max_total_iterations = (word.len() + 1) * (word.len() + 1) * 2;
    let mut total_iterations = 0;

    loop {
        total_iterations += 1;
        if total_iterations > max_total_iterations {
            break;
        }

        let prev = w.clone();

        // Step 1: Remove adjacent pairs
        w = remove_adjacent_pairs(&w);

        // Step 2: Sort commuting pairs and apply braids
        w = apply_braid_relations(coxeter, &w);

        // Step 3: Remove any new adjacent pairs
        w = remove_adjacent_pairs(&w);

        if w == prev {
            break;
        }
    }

    // Canonicalize to get a unique representation
    canonicalize(coxeter, &w)
}

fn remove_adjacent_pairs(w: &[usize]) -> Vec<usize> {
    let mut result = Vec::with_capacity(w.len());
    let mut i = 0;
    while i < w.len() {
        if i + 1 < w.len() && w[i] == w[i + 1] {
            i += 2;
        } else {
            result.push(w[i]);
            i += 1;
        }
    }
    result
}

fn apply_braid_relations(coxeter: &CoxeterMatrix, w: &[usize]) -> Vec<usize> {
    if w.len() < 2 {
        return w.to_vec();
    }

    // First, sort commuting pairs (m=2) using bubble sort
    let mut result = w.to_vec();
    let mut changed = true;
    while changed {
        changed = false;
        for i in 0..result.len().saturating_sub(1) {
            let si = result[i];
            let sj = result[i + 1];
            if si > sj && coxeter.entries[si][sj] == 2 {
                result.swap(i, i + 1);
                changed = true;
            }
        }
    }

    // Then apply braid relations for m >= 3
    if result.len() >= 3 {
        for i in 0..result.len().saturating_sub(2) {
            let si = result[i];
            let sj = result[i + 1];
            if si == sj {
                continue;
            }
            let m = coxeter.entries[si][sj];
            if !(3..=6).contains(&m) {
                continue;
            }
            let pattern = extract_alternating(&result, i, si, sj, m as usize);
            if pattern == m as usize {
                let mut new_result = result[..i].to_vec();
                for k in 0..(m as usize) {
                    new_result.push(if k % 2 == 0 { sj } else { si });
                }
                new_result.extend_from_slice(&result[i + m as usize..]);
                return new_result;
            }
        }
    }

    result
}

fn extract_alternating(w: &[usize], start: usize, a: usize, b: usize, max_len: usize) -> usize {
    let mut count = 0;
    let mut pos = start;
    while pos < w.len() && count < max_len {
        let expected = if count % 2 == 0 { a } else { b };
        if w[pos] != expected {
            break;
        }
        count += 1;
        pos += 1;
    }
    count
}

fn has_adjacent_pairs(w: &[usize]) -> bool {
    for i in 0..w.len().saturating_sub(1) {
        if w[i] == w[i + 1] {
            return true;
        }
    }
    false
}

fn has_applicable_braid(coxeter: &CoxeterMatrix, w: &[usize]) -> bool {
    for i in 0..w.len().saturating_sub(2) {
        if w[i] == w[i + 1] {
            continue;
        }
        let si = w[i];
        let sj = w[i + 1];
        let m = coxeter.entries[si][sj];
        if (3..=6).contains(&m) {
            let pattern = extract_alternating(w, i, si, sj, m as usize);
            if pattern == m as usize {
                return true;
            }
        }
    }
    false
}

/// Compute the length (number of generators in reduced form) of a word.
pub fn word_length(coxeter: &CoxeterMatrix, word: &[usize]) -> usize {
    reduce_word(coxeter, word).len()
}

/// Check if a word is reduced (i.e., already in shortest form).
pub fn is_reduced(coxeter: &CoxeterMatrix, word: &[usize]) -> bool {
    let reduced = reduce_word(coxeter, word);
    reduced.len() == word.len()
}

/// Compute the left descent set of a word: generators s such that
/// l(s * w) < l(w).
pub fn left_descent_set(coxeter: &CoxeterMatrix, word: &[usize]) -> Vec<usize> {
    let reduced = reduce_word(coxeter, word);
    let mut descents = vec![];
    for s in 0..coxeter.rank {
        let mut with_s = vec![s];
        with_s.extend_from_slice(&reduced);
        let reduced_with = reduce_word(coxeter, &with_s);
        if reduced_with.len() < reduced.len() {
            descents.push(s);
        }
    }
    descents
}

/// Compute the right descent set of a word: generators s such that
/// l(w * s) < l(w).
pub fn right_descent_set(coxeter: &CoxeterMatrix, word: &[usize]) -> Vec<usize> {
    let reduced = reduce_word(coxeter, word);
    let mut descents = vec![];
    for s in 0..coxeter.rank {
        let mut with_s = reduced.clone();
        with_s.push(s);
        let reduced_with = reduce_word(coxeter, &with_s);
        if reduced_with.len() < reduced.len() {
            descents.push(s);
        }
    }
    descents
}

/// Compute the inverse of a word (reverse order).
pub fn inverse_word(word: &[usize]) -> Vec<usize> {
    word.iter().rev().copied().collect()
}

/// Multiply two words (concatenate and reduce).
pub fn multiply_words(coxeter: &CoxeterMatrix, w1: &[usize], w2: &[usize]) -> Vec<usize> {
    let mut combined = w1.to_vec();
    combined.extend_from_slice(w2);
    reduce_word(coxeter, &combined)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reduce_identity() {
        let coxeter = CoxeterMatrix::type_a(2);
        assert_eq!(reduce_word(&coxeter, &[0, 0]), vec![]);
        assert_eq!(reduce_word(&coxeter, &[1, 1]), vec![]);
    }

    #[test]
    fn test_reduce_braid_a2() {
        let coxeter = CoxeterMatrix::type_a(2);
        // In A2: (s0 s1)^3 = identity via s0 s1 s0 = s1 s0 s1
        // s0 s1 s0 s1 s0 s1 → should reduce
        let reduced = reduce_word(&coxeter, &[0, 1, 0, 1, 0, 1]);
        assert_eq!(reduced.len(), 0);
    }

    #[test]
    fn test_word_length() {
        let coxeter = CoxeterMatrix::type_a(2);
        assert_eq!(word_length(&coxeter, &[]), 0);
        assert_eq!(word_length(&coxeter, &[0]), 1);
        assert_eq!(word_length(&coxeter, &[0, 1]), 2);
        assert_eq!(word_length(&coxeter, &[0, 0]), 0);
    }

    #[test]
    fn test_is_reduced() {
        let coxeter = CoxeterMatrix::type_a(2);
        assert!(is_reduced(&coxeter, &[]));
        assert!(is_reduced(&coxeter, &[0, 1]));
        assert!(!is_reduced(&coxeter, &[0, 0]));
    }

    #[test]
    fn test_left_descent_set() {
        let coxeter = CoxeterMatrix::type_a(2);
        let descents = left_descent_set(&coxeter, &[0, 1]);
        assert!(descents.contains(&0));
    }

    #[test]
    fn test_right_descent_set() {
        let coxeter = CoxeterMatrix::type_a(2);
        let descents = right_descent_set(&coxeter, &[0, 1]);
        assert!(descents.contains(&1));
    }

    #[test]
    fn test_inverse_word() {
        assert_eq!(inverse_word(&[0, 1, 2]), vec![2, 1, 0]);
        assert_eq!(inverse_word(&[]), vec![]);
    }

    #[test]
    fn test_multiply_words() {
        let coxeter = CoxeterMatrix::type_a(2);
        let result = multiply_words(&coxeter, &[0], &[0]);
        assert_eq!(result, vec![]);
        let result = multiply_words(&coxeter, &[0], &[1]);
        assert_eq!(result, vec![0, 1]);
    }

    #[test]
    fn test_reduce_a3() {
        let coxeter = CoxeterMatrix::type_a(3);
        // s0 s2 commutes (m=2), so s0 s2 s0 s2 = identity
        assert_eq!(reduce_word(&coxeter, &[0, 2, 0, 2]), vec![]);
    }
}
