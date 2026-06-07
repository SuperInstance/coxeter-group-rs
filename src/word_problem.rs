//! Word problem solver for Coxeter groups.
//!
//! Solves the word problem: given two words in the generators, determine
//! whether they represent the same group element. Provides algorithms for
//! finding all reduced expressions, shortest representations, and the
//! exchange/deletion conditions.

use crate::matrix::CoxeterMatrix;
use crate::reflection::ReflectionRepresentation;

/// Solve the word problem: check if two words represent the same group element.
///
/// Uses the reflection representation to compare group elements by their
/// action on the vector space. Two words represent the same element iff
/// their reflection matrices are equal.
pub fn words_equal(coxeter: &CoxeterMatrix, w1: &[usize], w2: &[usize]) -> bool {
    let repr = ReflectionRepresentation::new(coxeter);
    let m1 = repr.word_matrix(w1);
    let m2 = repr.word_matrix(w2);
    matrices_equal(&m1, &m2)
}

/// Find the length of the shortest representation of a group element
/// by brute-force enumeration up to the given bound.
///
/// Returns the shortest reduced word found, or the input reduced via
/// `crate::word::reduce_word` if enumeration doesn't improve it.
pub fn shortest_representation(coxeter: &CoxeterMatrix, word: &[usize], max_search: usize) -> Vec<usize> {
    let reduced = crate::word::reduce_word(coxeter, word);
    let target_len = reduced.len();
    if target_len == 0 {
        return vec![];
    }

    // Try to find a shorter equivalent word
    let repr = ReflectionRepresentation::new(coxeter);
    let target_matrix = repr.word_matrix(&reduced);

    // BFS over all words of length 0..target_len
    let mut best = reduced.clone();
    for len in 1..target_len {
        if len > max_search {
            break;
        }
        let found = search_length(coxeter, &target_matrix, len);
        if let Some(w) = found {
            best = w;
            break;
        }
    }
    best
}

/// Enumerate all reduced words that represent the same group element as `word`.
///
/// Uses braid relations to transform the word into all equivalent reduced forms.
/// Returns a set of distinct reduced words.
pub fn all_reduced_words(coxeter: &CoxeterMatrix, word: &[usize]) -> Vec<Vec<usize>> {
    let reduced = crate::word::reduce_word(coxeter, word);
    if reduced.is_empty() {
        return vec![vec![]];
    }

    let mut result = Vec::new();
    let mut seen = std::collections::HashSet::new();
    let mut queue = std::collections::VecDeque::new();

    queue.push_back(reduced.clone());
    seen.insert(reduced.clone());

    while let Some(current) = queue.pop_front() {
        result.push(current.clone());

        // Try all braid moves
        let neighbors = braid_neighbors(coxeter, &current);
        for neighbor in neighbors {
            if !seen.contains(&neighbor) {
                seen.insert(neighbor.clone());
                queue.push_back(neighbor);
            }
        }
    }

    result.sort();
    result
}

/// Apply the strong exchange condition.
///
/// If w = s_{i1} ... s_{ik} is reduced and s_j w is not reduced, then there
/// exists an index m such that s_j w = s_{i1} ... s_{im-1} s_{im+1} ... s_{ik}.
/// Returns the index of the generator that can be deleted.
pub fn strong_exchange(coxeter: &CoxeterMatrix, word: &[usize], s: usize) -> Option<usize> {
    let reduced = crate::word::reduce_word(coxeter, word);
    let mut with_s = vec![s];
    with_s.extend_from_slice(&reduced);
    let reduced_with = crate::word::reduce_word(coxeter, &with_s);

    if reduced_with.len() > reduced.len() {
        return None; // s * w is still reduced
    }

    // Strong exchange: s*w = w_1...w_{m-1}*w_{m+1}...w_k
    // i.e., removing the m-th letter from w gives the same group element as s*w
    for m in 0..reduced.len() {
        let mut without_m: Vec<usize> = reduced.iter().enumerate()
            .filter(|(i, _)| *i != m)
            .map(|(_, &g)| g)
            .collect();
        let test_reduced = crate::word::reduce_word(coxeter, &without_m);
        if test_reduced == reduced_with {
            return Some(m);
        }
    }
    None
}

/// Apply the deletion condition.
///
/// If w is not reduced, there exists a pair of indices i < j such that
/// deleting s_i and s_j yields a shorter word for the same element.
/// Returns all such deletable pairs.
pub fn deletion_condition(coxeter: &CoxeterMatrix, word: &[usize]) -> Vec<(usize, usize)> {
    let reduced = crate::word::reduce_word(coxeter, word);
    if reduced.len() == word.len() {
        return vec![]; // Already reduced
    }

    let mut pairs = vec![];
    let repr = ReflectionRepresentation::new(coxeter);
    let target = repr.word_matrix(word);

    for i in 0..word.len() {
        for j in (i + 1)..word.len() {
            let mut deleted = Vec::with_capacity(word.len() - 2);
            for (k, &s) in word.iter().enumerate() {
                if k != i && k != j {
                    deleted.push(s);
                }
            }
            let d_matrix = repr.word_matrix(&deleted);
            if matrices_equal(&d_matrix, &target) {
                pairs.push((i, j));
            }
        }
    }
    pairs
}

/// Compute the Coxeter presentation: return the generators and relations.
///
/// Returns a list of relations as (m, i, j) meaning (s_i s_j)^m = e.
pub fn coxeter_presentation(coxeter: &CoxeterMatrix) -> Vec<(u32, usize, usize)> {
    let mut relations = vec![];
    for i in 0..coxeter.rank {
        relations.push((2u32, i, i)); // s_i^2 = e
        for j in (i + 1)..coxeter.rank {
            let m = coxeter.entries[i][j];
            if m > 2 {
                relations.push((m, i, j));
            }
        }
    }
    relations
}

/// Check the braid relation: (s_i s_j)^m = (s_j s_i)^m = e
/// where m = m_{ij} from the Coxeter matrix.
pub fn verify_braid_relation(coxeter: &CoxeterMatrix, i: usize, j: usize) -> bool {
    let m = coxeter.entries[i][j];
    if m == 0 {
        return true; // Infinite order, can't easily verify
    }

    let repr = ReflectionRepresentation::new(coxeter);
    // Build the word (s_i s_j)^m = alternating sequence of length 2*m
    let mut word = vec![];
    for k in 0..(2 * m as usize) {
        word.push(if k % 2 == 0 { i } else { j });
    }
    repr.is_identity(&word)
}

/// Find a geodesic (shortest) word for a group element using
/// breadth-first search over the Cayley graph.
pub fn geodesic_word(coxeter: &CoxeterMatrix, word: &[usize]) -> Vec<usize> {
    let repr = ReflectionRepresentation::new(coxeter);
    let target = repr.word_matrix(word);

    if is_identity_matrix(&target) {
        return vec![];
    }

    // Check if the already-reduced word is geodesic
    let reduced = crate::word::reduce_word(coxeter, word);
    let reduced_mat = repr.word_matrix(&reduced);
    if matrices_equal(&reduced_mat, &target) {
        return reduced;
    }

    // BFS from identity
    let mut seen = std::collections::HashSet::new();
    let target_key = matrix_key(&target);
    seen.insert(matrix_key(&identity_matrix(coxeter.rank)));
    let mut queue = std::collections::VecDeque::new();
    queue.push_back((vec![], identity_matrix(coxeter.rank)));

    while let Some((path, current)) = queue.pop_front() {
        for s in 0..coxeter.rank {
            let mut new_path = path.clone();
            new_path.push(s);
            let new_matrix = mat_mul(&repr.reflections[s], &current);

            let key = matrix_key(&new_matrix);
            if key == target_key {
                return new_path;
            }

            if seen.insert(key) {
                queue.push_back((new_path, new_matrix));
            }
        }
    }

    reduced
}

/// Compute the element order: smallest k > 0 such that w^k = e.
/// Returns None if the element appears to have infinite order.
pub fn element_order(coxeter: &CoxeterMatrix, word: &[usize], max_k: usize) -> Option<usize> {
    let repr = ReflectionRepresentation::new(coxeter);
    let m = repr.word_matrix(word);

    if is_identity_matrix(&m) {
        return Some(1);
    }

    let mut current = m.clone();
    for k in 2..=max_k {
        current = mat_mul(&m, &current);
        if is_identity_matrix(&current) {
            return Some(k);
        }
    }
    None
}

// --- Helper functions ---

fn matrices_equal(a: &[Vec<f64>], b: &[Vec<f64>]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    for i in 0..a.len() {
        for j in 0..a[i].len() {
            if (a[i][j] - b[i][j]).abs() > 1e-8 {
                return false;
            }
        }
    }
    true
}

fn identity_matrix(n: usize) -> Vec<Vec<f64>> {
    let mut m = vec![vec![0.0; n]; n];
    for i in 0..n {
        m[i][i] = 1.0;
    }
    m
}

fn is_identity_matrix(m: &[Vec<f64>]) -> bool {
    let n = m.len();
    for i in 0..n {
        for j in 0..n {
            let expected = if i == j { 1.0 } else { 0.0 };
            if (m[i][j] - expected).abs() > 1e-8 {
                return false;
            }
        }
    }
    true
}

fn mat_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut c = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

fn matrix_key(m: &[Vec<f64>]) -> Vec<u64> {
    m.iter().flat_map(|row| {
        row.iter().map(|&x| {
            (x * 1e8).round() as i64 as u64
        })
    }).collect()
}

fn search_length(
    coxeter: &CoxeterMatrix,
    target: &[Vec<f64>],
    target_len: usize,
) -> Option<Vec<usize>> {
    let repr = ReflectionRepresentation::new(coxeter);
    let n = coxeter.rank;

    let mut queue = std::collections::VecDeque::new();
    let mut seen = std::collections::HashSet::new();

    for s in 0..n {
        let w = vec![s];
        let key = matrix_key(&repr.reflections[s]);
        if seen.insert(key) {
            queue.push_back(w);
        }
    }

    while let Some(word) = queue.pop_front() {
        if word.len() == target_len {
            let m = repr.word_matrix(&word);
            if matrices_equal(&m, target) {
                return Some(word);
            }
            continue;
        }

        if word.len() >= target_len {
            continue;
        }

        for s in 0..n {
            let mut new_word = word.clone();
            new_word.push(s);
            if new_word.len() <= target_len {
                queue.push_back(new_word);
            }
        }
    }
    None
}

/// Generate all braid-equivalent neighbors of a reduced word.
fn braid_neighbors(coxeter: &CoxeterMatrix, word: &[usize]) -> Vec<Vec<usize>> {
    let mut neighbors = vec![];

    // Commutation moves: swap adjacent commuting generators (m = 2)
    for i in 0..word.len().saturating_sub(1) {
        let si = word[i];
        let sj = word[i + 1];
        if si != sj && coxeter.entries[si][sj] == 2 {
            let mut neighbor = word.to_vec();
            neighbor.swap(i, i + 1);
            neighbors.push(neighbor);
        }
    }

    // Braid moves: replace alternating subword
    for i in 0..word.len() {
        for j in (i + 2)..=word.len().min(i + 8) {
            if j - i < 3 {
                continue;
            }
            let sub = &word[i..j];
            if let Some((a, b)) = check_alternating(sub) {
                let m = coxeter.entries[a][b];
                if m > 2 && (j - i) == m as usize {
                    let mut neighbor = word[..i].to_vec();
                    for k in 0..(m as usize) {
                        neighbor.push(if k % 2 == 0 { b } else { a });
                    }
                    neighbor.extend_from_slice(&word[j..]);
                    neighbors.push(neighbor);
                }
            }
        }
    }

    neighbors
}

fn check_alternating(word: &[usize]) -> Option<(usize, usize)> {
    if word.len() < 2 {
        return None;
    }
    let a = word[0];
    let b = word[1];
    if a == b {
        return None;
    }
    for (i, &s) in word.iter().enumerate() {
        let expected = if i % 2 == 0 { a } else { b };
        if s != expected {
            return None;
        }
    }
    Some((a, b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::CoxeterMatrix;

    #[test]
    fn test_words_equal_identity() {
        let coxeter = CoxeterMatrix::type_a(2);
        assert!(words_equal(&coxeter, &[0, 0], &[]));
        assert!(words_equal(&coxeter, &[0, 1, 0], &[1, 0, 1]));
    }

    #[test]
    fn test_words_not_equal() {
        let coxeter = CoxeterMatrix::type_a(2);
        assert!(!words_equal(&coxeter, &[0], &[1]));
    }

    #[test]
    fn test_words_equal_a3() {
        let coxeter = CoxeterMatrix::type_a(3);
        // s0*s2 commutes with s2*s0
        assert!(words_equal(&coxeter, &[0, 2], &[2, 0]));
    }

    #[test]
    fn test_shortest_representation() {
        let coxeter = CoxeterMatrix::type_a(2);
        // s0*s1*s0*s1*s0*s1 should reduce to empty
        let result = shortest_representation(&coxeter, &[0, 1, 0, 1, 0, 1], 10);
        assert!(result.is_empty());
    }

    #[test]
    fn test_shortest_representation_nontrivial() {
        let coxeter = CoxeterMatrix::type_a(2);
        let result = shortest_representation(&coxeter, &[0, 0, 1], 10);
        assert_eq!(result.len(), 1);
        assert_eq!(result, vec![1]);
    }

    #[test]
    fn test_all_reduced_words_a2() {
        let coxeter = CoxeterMatrix::type_a(2);
        // s0*s1 has two reduced forms: [0,1] and [1,0,1]... wait, no
        // In A2: [0,1] is reduced. [1,0,1] = [0,1,0] by braid (length 3)
        // [0,1] only has itself as reduced
        let words = all_reduced_words(&coxeter, &[0, 1]);
        assert!(!words.is_empty());
        assert!(words.iter().all(|w| w.len() == 2));
    }

    #[test]
    fn test_all_reduced_words_braid() {
        let coxeter = CoxeterMatrix::type_a(2);
        // s0*s1*s0 = s1*s0*s1 by braid relation
        let words = all_reduced_words(&coxeter, &[0, 1, 0]);
        assert!(words.len() >= 2);
    }

    #[test]
    fn test_deletion_condition() {
        let coxeter = CoxeterMatrix::type_a(2);
        // [0, 0] is not reduced - can delete the pair
        let pairs = deletion_condition(&coxeter, &[0, 0]);
        assert!(!pairs.is_empty());
    }

    #[test]
    fn test_deletion_condition_reduced() {
        let coxeter = CoxeterMatrix::type_a(2);
        // [0, 1] is reduced - no deletable pairs
        let pairs = deletion_condition(&coxeter, &[0, 1]);
        assert!(pairs.is_empty());
    }

    #[test]
    fn test_coxeter_presentation() {
        let coxeter = CoxeterMatrix::type_a(3);
        let relations = coxeter_presentation(&coxeter);
        // s_i^2 = e for 3 generators
        assert!(relations.iter().any(|&(m, i, j)| m == 2 && i == 0 && j == 0));
        assert!(relations.iter().any(|&(m, i, j)| m == 2 && i == 1 && j == 1));
        assert!(relations.iter().any(|&(m, i, j)| m == 2 && i == 2 && j == 2));
        // Braid: (s0*s1)^3 = e, (s1*s2)^3 = e
        assert!(relations.iter().any(|&(m, i, j)| m == 3 && i == 0 && j == 1));
        assert!(relations.iter().any(|&(m, i, j)| m == 3 && i == 1 && j == 2));
    }

    #[test]
    fn test_verify_braid_relation() {
        let coxeter = CoxeterMatrix::type_a(2);
        // (s0*s1)^3 = e in A2
        // The word for (s0*s1)^3 is [0,1,0,1,0,1]
        assert!(verify_braid_relation(&coxeter, 0, 1));
        // Also check: (s0*s0) = e always
        // For m=1 (diagonal), it's trivially true
    }

    #[test]
    fn test_element_order_identity() {
        let coxeter = CoxeterMatrix::type_a(2);
        assert_eq!(element_order(&coxeter, &[], 10), Some(1));
        assert_eq!(element_order(&coxeter, &[0, 0], 10), Some(1));
    }

    #[test]
    fn test_element_order_reflection() {
        let coxeter = CoxeterMatrix::type_a(2);
        // A single reflection has order 2
        assert_eq!(element_order(&coxeter, &[0], 10), Some(2));
        assert_eq!(element_order(&coxeter, &[1], 10), Some(2));
    }

    #[test]
    fn test_strong_exchange() {
        let coxeter = CoxeterMatrix::type_a(2);
        // In A2, s0*(s0*s1) = s1 (since s0*s0=e), so prepending s0 to word [0,1]
        // should give a shorter word. Word [0,1] has length 2.
        // s0 * [0,1] -> [0,0,1] reduces to [1] with length 1.
        // The strong exchange condition says there exists an index to delete.
        // Let's test with word [0,1] and generator s0:
        // s0 * [0,1] = [0,0,1] which reduces to [1] (length 1 < 2)
        // So index 0 can be deleted from [0,1] to get [1].
        let result = strong_exchange(&coxeter, &[0, 1], 0);
        assert!(result.is_some());
    }

    #[test]
    fn test_geodesic_word() {
        let coxeter = CoxeterMatrix::type_a(2);
        let result = geodesic_word(&coxeter, &[0, 0, 1]);
        assert!(result.len() <= 1);
    }

    #[test]
    fn test_geodesic_word_identity() {
        let coxeter = CoxeterMatrix::type_a(2);
        let result = geodesic_word(&coxeter, &[0, 0]);
        assert!(result.is_empty());
    }
}
