// Code from rust-bio, modified for case insensitivity and IUPAC alphabet support

// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Backward nondeterministic DAWG matching (BNDM).
//! Best-case complexity: O(n / m) with pattern of length m <= 64 and text of length n.
//! Worst case complexity: O(n * m).
//!

use bio::utils::{IntoTextIterator, TextSlice};
use iupac::expand_iupac;

// from shift_and.rs - only this function has been modified (added a test too)
pub fn masks<'a, I: IntoTextIterator<'a>>(pattern: I) -> ([u64; 256], u64) {
    let mut masks = [0; 256];
    let mut bit = 1;
    for &c in pattern {
        if let Some(nts) = expand_iupac(c) {
            for &nt in nts {
                masks[nt as usize] |= bit;
                masks[nt as usize - 32] |= bit; // Add uppercase version
            }
        } else {
            // not a valid IUPAC code
            if c > 64 && c < 91 {
                // add lowercase version if uppercase
                masks[c as usize + 32] |= bit;
            } else if c > 96 && c < 123 {
                // add uppercase version if lowercase
                masks[c as usize - 32] |= bit;
            }
            masks[c as usize] |= bit;
        }
        bit *= 2;
    }

    (masks, bit / 2)
}

/// BNDM algorithm.
pub struct BNDM {
    m: usize,
    masks: [u64; 256],
    accept: u64,
}

impl BNDM {
    /// Create a new instance for a given pattern.
    pub fn new<'a, P: IntoTextIterator<'a>>(pattern: P) -> Self
    where
        P::IntoIter: DoubleEndedIterator + ExactSizeIterator,
    {
        let pattern = pattern.into_iter();
        let m = pattern.len();
        assert!(m <= 64, "Expecting a pattern of at most 64 symbols.");
        // take the reverse pattern and build nondeterministic
        // suffix automaton
        let (masks, accept) = masks(pattern.rev());

        BNDM {
            m: m,
            masks: masks,
            accept: accept,
        }
    }

    /// Find all matches of pattern with a given text. Matches are returned as iterator over start positions.
    pub fn find_all<'a>(&'a self, text: TextSlice<'a>) -> Matches {
        Matches {
            bndm: self,
            window: self.m,
            text: text,
        }
    }
}

/// Iterator over start positions of matches.
pub struct Matches<'a> {
    bndm: &'a BNDM,
    window: usize,
    text: TextSlice<'a>,
}

impl<'a> Iterator for Matches<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        while self.window <= self.text.len() {
            let mut occ = None;
            // bit mask of ones, all states active
            let mut active = (1u64 << self.bndm.m) - 1;
            let (mut j, mut lastsuffix) = (1, 0);
            // while not in fail state
            while active != 0 {
                // process j-th symbol from right
                active &= self.bndm.masks[self.text[self.window - j] as usize];
                if active & self.bndm.accept != 0 {
                    // reached accepting state
                    if j == self.bndm.m {
                        occ = Some(self.window - self.bndm.m);
                        break;
                    } else {
                        // we reached the accepting state
                        // but not the end of the pattern
                        // hence, a suffix of the reverse pattern
                        // i.e. a prefix of the pattern of
                        // length j matches
                        // in case of a mismatch, we can shift
                        // to this prefix
                        lastsuffix = j;
                    }
                }
                j += 1;
                active <<= 1;
            }
            // shift the window
            self.window += self.bndm.m - lastsuffix;
            if occ.is_some() {
                return occ;
            }
        }

        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;

    #[test]
    fn test_find_all() {
        let text = b"dhjalkjwqtttattataflkjdklfj";
        let pattern = b"qtttattat";
        let bndm = BNDM::new(pattern);
        assert_eq!(bndm.find_all(text).collect_vec(), [8]);
    }
    #[test]
    fn test_find_all_mixed_case() {
        let text = b"dhjalkjwqtttattataflkjdklfj";
        let pattern = b"qTttAttAt";
        let bndm = BNDM::new(pattern);
        assert_eq!(bndm.find_all(text).collect_vec(), [8]);
    }
    #[test]
    fn test_find_all_overlapping() {
        let text = b"tagtagtagt";
        let pattern = b"tagt";
        let bndm = BNDM::new(pattern);
        assert_eq!(bndm.find_all(text).collect_vec(), [0, 3, 6]);
    }
    #[test]
    fn test_find_all_iupac() {
        let text = b"tagtagtagt";
        let pattern = b"nnnN";
        let bndm = BNDM::new(pattern);
        assert_eq!(bndm.find_all(text).collect_vec(), [0, 1, 2, 3, 4, 5, 6]);
    }
}
