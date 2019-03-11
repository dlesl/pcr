extern crate bio;
#[macro_use]
extern crate gb_io;
extern crate itertools;
#[macro_use]
extern crate log;
#[macro_use]
extern crate auto_impl;
use bio::alignment::sparse;
use bio::alphabets::dna::{self, revcomp};
use std::borrow::Cow;
use std::mem;
use std::str;

use gb_io::seq::{Feature, Location, QualifierKey, Seq};

pub mod bndm_iupac;
use bndm_iupac::BNDM;
pub mod iupac;
use iupac::nt_match;

#[cfg(feature = "parallel")]
extern crate rayon;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[cfg(not(feature = "parallel"))]
pub trait PrimerBounds: Clone {}

#[cfg(not(feature = "parallel"))]
impl<T> PrimerBounds for T where T: Clone {}

#[cfg(feature = "parallel")]
pub trait PrimerBounds: Clone + Sync + Send {}

#[cfg(feature = "parallel")]
impl<T> PrimerBounds for T where T: Clone + Sync + Send {}

/// A trait to allow client code to use its own type for primers,
/// e.g. to attach additional information. The trait is also implemented
/// for references and smart pointers to a type implementing `Primer`.
/// Must be `Clone`, since multiple matches can be found per primer.
#[auto_impl(&, Rc, Arc, Box)]
pub trait Primer: PrimerBounds {
    fn seq(&self) -> &[u8];
    fn seq_rc(&self) -> Cow<[u8]> {
        revcomp(self.seq()).into()
    }
    fn name(&self) -> &str;
    fn len(&self) -> i64 {
        self.seq().len() as i64
    }
}

impl Primer for &[u8] {
    fn seq(&self) -> &[u8] {
        self
    }
    fn seq_rc(&self) -> Cow<[u8]> {
        revcomp(*self).into()
    }
    fn name(&self) -> &str {
        str::from_utf8(self.seq()).unwrap_or("invalid_utf8")
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct Footprint<T: Primer> {
    pub start: i64,
    pub extent: i64, // - negative for other strand
    pub primer: T,
}

impl<T: Primer> Footprint<T> {
    pub fn len(&self) -> i64 {
        self.extent.abs()
    }
    pub fn is_forward(&self) -> bool {
        assert!(self.extent != 0);
        self.extent > 0
    }
    pub fn is_reverse(&self) -> bool {
        assert!(self.extent != 0);
        self.extent < 0
    }
}

#[derive(Clone, Copy, Debug)]
pub enum Method {
    /// Fastest with smaller numbers of primers, uses no extra memory. Supports
    /// IUPAC codes.
    Bndm,
    /// Index the template first - faster with large numbers of primers (hundreds),
    /// uses a LOT of memory. Doesn't support IUPAC codes.
    Index,
}

pub fn find_matches<'a, T: Primer + 'a>(
    template: &'a Seq,
    primers: impl IntoIterator<Item = T>,
    min_fp: i64,
    method: Method,
) -> Matches<T> {
    Annealer::new(template, min_fp, method).find_matches(primers)
}

struct Annealer<'a> {
    template: &'a Seq,
    min_fp: i64,
    method: Method,
}

#[derive(Clone, PartialEq, Debug)]
pub struct Matches<T: Primer> {
    pub fwd: Vec<Footprint<T>>,
    pub rev: Vec<Footprint<T>>,
}

impl<'a> Annealer<'a> {
    pub fn new(template: &'a Seq, min_fp: i64, method: Method) -> Annealer<'a> {
        Annealer {
            template,
            min_fp,
            method,
        }
    }
    pub fn find_matches<T: Primer + 'a>(self, primers: impl IntoIterator<Item = T>) -> Matches<T> {
        // Check that primers are appropriate:
        // - long enough
        // - no IUPAC codes in seed region if method == Index
        let primers: Vec<T> = primers
        .into_iter()
        .filter(|p| {
            if p.len() < self.min_fp {
                warn!("Primer '{}' too short, skipping", p.name());
                false
            } else {
                match self.method {
                    Method::Bndm => true,
                    Method::Index => {
                        if !dna::alphabet().is_word(p.seq()) {
                            warn!("Invalid bases in primer '{}', method 'Index' doesn't support IUPAC codes, skipping", p.name());
                            false
                        } else {
                            true
                        }
                    }
                }
            }
        }).collect();
        let seq = &self.template.seq;
        let mut res: (Vec<_>, Vec<_>) = match self.method {
            Method::Bndm => {
                #[cfg(feature = "parallel")]
                let res: (Vec<_>, Vec<_>) = rayon::join(
                    || {
                        primers
                            .par_iter()
                            .flat_map(|p| self.search_bndm_fwd(p, seq, 0).into_par_iter())
                            .collect()
                    },
                    || {
                        primers
                            .par_iter()
                            .flat_map(|p| self.search_bndm_rev(p, seq, 0).into_par_iter())
                            .collect()
                    },
                );
                #[cfg(not(feature = "parallel"))]
                let res: (Vec<_>, Vec<_>) = (
                    primers
                        .iter()
                        .flat_map(|p| self.search_bndm_fwd(p, seq, 0).into_iter())
                        .collect(),
                    primers
                        .iter()
                        .flat_map(|p| self.search_bndm_rev(p, seq, 0).into_iter())
                        .collect(),
                );
                res
            }
            Method::Index => {
                let seq_upper = seq.to_ascii_uppercase();
                #[cfg(not(feature = "parallel"))]
                let index = vec![(0, sparse::hash_kmers(&seq_upper, self.min_fp as usize))];
                #[cfg(feature = "parallel")]
                let index: Vec<_> = overlapping_chunks(
                    &seq_upper,
                    self.template.seq.len() / rayon::current_num_threads(),
                    self.min_fp as usize,
                )
                .into_par_iter()
                .map(|(offset, chunk)| (offset, sparse::hash_kmers(chunk, self.min_fp as usize)))
                .collect();
                let mut fwd = Vec::new();
                let mut rev = Vec::new();
                for (offset, index) in index {
                    fwd.extend(primers.iter().flat_map(|p| {
                        sparse::find_kmer_matches_seq1_hashed(
                            &index,
                            &p.seq()[p.seq().len() - self.min_fp as usize..].to_ascii_uppercase(),
                            self.min_fp as usize,
                        )
                        .into_iter()
                        .map(|(m, _)| self.extend_fwd(p.clone(), m as i64 + offset as i64))
                        .collect::<Vec<_>>()
                    }));
                    rev.extend(primers.iter().flat_map(|p| {
                        sparse::find_kmer_matches_seq1_hashed(
                            &index,
                            &p.seq_rc()[..self.min_fp as usize].to_ascii_uppercase(),
                            self.min_fp as usize,
                        )
                        .into_iter()
                        .map(|(m, _)| self.extend_rev(p.clone(), m as i64 + offset as i64))
                        .collect::<Vec<_>>()
                    }));
                }
                (fwd, rev)
            }
        };
        // collect matches in "origin" region
        if self.template.is_circular() {
            // one less so that we don't find the same match twice
            let start = seq.len() as i64 - self.min_fp + 1;
            let end = seq.len() as i64 + self.min_fp - 1;
            let origin = self.template.extract_range_seq(start, end);
            assert!(origin.len() as i64 == self.min_fp * 2 - 2);
            res.0.extend(
                primers
                    .iter()
                    .flat_map(|p| self.search_bndm_fwd(p, &origin, start)),
            );
            res.1.extend(
                primers
                    .iter()
                    .flat_map(|p| self.search_bndm_rev(p, &origin, start)),
            );
        }

        Matches {
            fwd: res.0,
            rev: res.1,
        }
    }

    // search for matches, for optimisation purposes, the slice to search and
    // the Location it starts at is provided so we can reuse it
    fn search_bndm_fwd<T: Primer>(
        &self,
        primer: &T,
        seq: &[u8],
        starts_at: i64,
    ) -> Vec<Footprint<T>> {
        let three_prime = &primer.seq()[primer.seq().len() - self.min_fp as usize..];
        BNDM::new(three_prime)
            .find_all(seq)
            .map(|m| self.extend_fwd(primer.clone(), m as i64 + starts_at))
            .collect()
    }

    fn search_bndm_rev<T: Primer>(
        &self,
        primer: &T,
        seq: &[u8],
        starts_at: i64,
    ) -> Vec<Footprint<T>> {
        let three_prime = &primer.seq_rc()[..self.min_fp as usize];
        BNDM::new(three_prime)
            .find_all(seq)
            .map(|m| self.extend_rev(primer.clone(), m as i64 + starts_at))
            .collect()
    }

    fn extend_fwd<T: Primer>(&self, primer: T, pos: i64) -> Footprint<T> {
        let primer_end = pos + self.min_fp;
        let mut first = primer_end - primer.seq().len() as i64;
        if !self.template.is_circular() {
            first = std::cmp::max(first, 0);
        }
        let template_ext = self.template.extract_range_seq(first, pos);
        let count = template_ext
            .iter()
            .rev()
            .zip(
                primer.seq()[..primer.seq().len() - self.min_fp as usize]
                    .iter()
                    .rev(),
            )
            .take_while(|&(&a, &b)| nt_match(a, b))
            .count() as i64;
        let start = pos - count;
        let start_wrapped = if start < 0 {
            assert!(self.template.is_circular());
            start + self.template.seq.len() as i64
        } else {
            start
        };
        Footprint {
            start: start_wrapped,
            extent: self.min_fp + count,
            primer,
        }
    }

    fn extend_rev<T: Primer>(&self, primer: T, pos: i64) -> Footprint<T> {
        let mut primer_end = pos + primer.seq().len() as i64;
        if !self.template.is_circular() {
            primer_end = std::cmp::min(primer_end, self.template.len());
        }
        let template_ext = self
            .template
            .extract_range_seq(pos + self.min_fp, primer_end);
        let count = template_ext
            .iter()
            .zip(primer.seq_rc()[self.min_fp as usize..].iter())
            .take_while(|&(&a, &b)| nt_match(a, b))
            .count() as i64;
        let start = pos - 1 + self.min_fp + count;
        let len = self.template.seq.len() as i64;
        let start_wrapped = if start >= len {
            assert!(self.template.is_circular());
            start - len
        } else {
            start
        };
        Footprint {
            start: start_wrapped,
            extent: -(self.min_fp + count),
            primer,
        }
    }
}

impl<T: Primer> Matches<T> {
    pub fn find_products<'a>(
        &'a self,
        template: &'a Seq,
        min_len: i64,
        max_len: i64,
    ) -> impl Iterator<Item = Product<T>> + 'a {
        self.fwd
            .iter()
            .flat_map(move |fwd| {
                self.rev
                    .iter()
                    .map(move |rev| Product(fwd.clone(), rev.clone()))
            })
            .filter(move |p| {
                p.len(template)
                    .map(|len| len >= min_len && len <= max_len)
                    .unwrap_or(false)
            })
    }
}

impl<T: Primer> Footprint<T> {
    pub fn annotate(&self, mut seq: Seq) -> Seq {
        let f = Feature {
            kind: gb_io::FeatureKind::from("primer_bind"),
            location: if self.is_forward() {
                seq.range_to_location(self.start, self.start + self.extent)
            } else {
                Location::Complement(Box::new(
                    seq.range_to_location(self.start + self.extent + 1, self.start + 1),
                ))
            },
            qualifiers: vec![(
                QualifierKey::from("PCR_primer"),
                Some(format!(
                    "name: {}, seq: {}",
                    self.primer.name(),
                    String::from_utf8_lossy(self.primer.seq()),
                )),
            )],
        };
        seq.features.push(f);
        seq
    }
}

impl<T: Primer> Product<T> {
    pub fn annotate(&self, mut seq: Seq, annotate_binding: bool) -> Seq {
        if annotate_binding {
            seq = self.1.annotate(self.0.annotate(seq));
        }
        let f = Feature {
            kind: feature_kind!("misc_feature"),
            location: seq.range_to_location(self.0.start, self.1.start + 1),
            qualifiers: vec![(
                QualifierKey::from("PCR_primers"),
                Some(format!(
                    "[fwd_name: {}, ]fwd_seq: {},\n\
                     [rev_name: {}, ]rev_seq: {}",
                    self.0.primer.name(),
                    String::from_utf8_lossy(self.0.primer.seq()),
                    self.1.primer.name(),
                    String::from_utf8_lossy(self.1.primer.seq())
                )),
            )],
        };
        seq.features.push(f);
        seq
    }
    pub fn extract(&self, template: &Seq) -> Seq {
        let Product(ref fwd, ref rev) = self;
        if fwd.start <= rev.start + rev.extent || template.is_circular() {
            let mut res = template.extract_range(fwd.start, rev.start + 1);
            let mut features = Vec::new();
            mem::swap(&mut features, &mut res.features);

            let features = features
                .into_iter()
                .map(|f| res.relocate_feature(f, fwd.primer.seq().len() as i64 - fwd.extent))
                .collect::<Result<Vec<_>, _>>()
                .unwrap();
            Seq {
                seq: fwd.primer.seq()[..(fwd.primer.seq().len() as i64 - fwd.extent) as usize]
                    .iter()
                    .chain(res.seq.iter())
                    .chain(rev.primer.seq_rc()[(-rev.extent) as usize..].iter())
                    .cloned()
                    .collect(),
                features,
                ..res
            }
        } else {
            panic!("invalid product");
        }
    }
    pub fn len(&self, template: &Seq) -> Option<i64> {
        let Product(ref fwd, ref rev) = self;
        let overhangs = (fwd.primer.len() - fwd.len()) + (rev.primer.len() - rev.len());
        if fwd.start <= rev.start + rev.extent {
            Some(rev.start - fwd.start + 1 + overhangs)
        } else if template.is_circular() {
            let rev_start_unwrapped = rev.start + template.len();
            let distance = rev_start_unwrapped - fwd.start + 1;
            Some(distance + overhangs)
        } else {
            None
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct Product<T: Primer>(pub Footprint<T>, pub Footprint<T>);

#[cfg(feature = "parallel")]
fn overlapping_chunks(v: &[u8], chunk_size: usize, overlap: usize) -> Vec<(usize, &[u8])> {
    use std::cmp;
    #[cfg(test)]
    const MIN_CHUNK_SIZE: usize = 1;
    #[cfg(not(test))]
    const MIN_CHUNK_SIZE: usize = 10_000;
    let chunk_size = cmp::max(cmp::max(chunk_size, overlap + 1), MIN_CHUNK_SIZE);
    if chunk_size >= v.len() {
        return vec![(0, v)];
    }
    let mut res = Vec::new();
    let mut last = overlap - 1;
    let mut next = cmp::min(v.len(), chunk_size);
    while last < v.len() {
        let offset = last - (overlap - 1);
        let slice = &v[offset..next];
        res.push((offset, slice));
        last = next;
        next = cmp::min(v.len(), last + chunk_size);
    }
    res
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_extend() {
        use gb_io::seq::Topology;
        let s = Seq {
            seq: b"0123456789".to_vec(),
            ..Seq::empty()
        };
        let primer = &b"245678"[..];
        let min_fp = 3;
        assert_eq!(
            Annealer::new(&s, min_fp, Method::Bndm).extend_fwd(primer, 6),
            Footprint {
                start: 4,
                extent: primer.len() as i64 - 1,
                primer
            }
        );
        //circular
        let s_circ = Seq {
            topology: Topology::Circular,
            ..s.clone()
        };
        let primer = &b"890123"[..];
        assert_eq!(
            Annealer::new(&s_circ, min_fp, Method::Bndm).extend_fwd(primer, 1),
            Footprint {
                start: 8,
                extent: primer.len() as i64,
                primer
            }
        );
        // revcomp leaves invalid characters alone so this is fine
        let primer = &b"665432"[..];
        assert_eq!(
            Annealer::new(&s, min_fp, Method::Bndm).extend_rev(primer, 2),
            Footprint {
                start: 6,
                extent: -(primer.len() as i64) + 1,
                primer,
            }
        );
        let primer = &b"65432109"[..];
        assert_eq!(
            Annealer::new(&s_circ, min_fp, Method::Bndm).extend_rev(primer, 9),
            Footprint {
                start: 6,
                extent: -(primer.len() as i64),
                primer,
            }
        );
    }

    #[test]
    fn test_search() {
        test_search_impl(Method::Bndm);
        test_search_impl(Method::Index);
    }
    fn test_search_impl(method: Method) {
        let seq = b"CGGAAATCCTCAAGCACCAGGTACGCTCATTGGTGCCAGCCGTGATGAAGACGAATTACCGGTCAAGGGC\
                 ATTTCCAATCTGAATAACATGGCAATGTTCAGCGTTTCTGGTCCGGGGATGAAAGGGATGGTCGGCATGG\
                 CGGCGCGCGTCTTTGCAGCGATGTCACGCGCCCGTATTTCCGTGGTGCTGATTACGCAATCATCTTCCGA";
        let s = Seq {
            seq: seq.to_vec(),
            ..Seq::empty()
        };

        let primers = &[&b"GGGAAATCCTCAAGCACCAG"[..], &b"TTGATTGCGTAATCAGCACCAC"[..]][..];
        let matches = Annealer::new(&s, 14, method).find_matches(primers);
        assert_eq!(
            (&matches.fwd, &matches.rev),
            (
                &vec![Footprint {
                    start: 1,
                    extent: 19,
                    primer: &primers[0],
                },],
                &vec![Footprint {
                    start: 201,
                    extent: -21,
                    primer: &primers[1],
                },],
            )
        );
    }

    #[test]
    fn test_search_circ() {
        test_search_circ_impl(Method::Bndm);
        test_search_circ_impl(Method::Index);
    }
    fn test_search_circ_impl(method: Method) {
        use gb_io::seq::Topology;
        let seq = b"GACTATTA";
        let s = Seq {
            seq: seq.to_vec(),
            topology: Topology::Circular,
            ..Seq::empty()
        };

        let primers = &[&b"AGACT"[..], &b"CTAAT"[..]][..];
        let matches = Annealer::new(&s, 4, method).find_matches(primers);
        assert_eq!(
            (&matches.fwd, &matches.rev),
            (
                &vec![Footprint {
                    start: 7,
                    extent: 5,
                    primer: &primers[0],
                },],
                &vec![Footprint {
                    start: 0,
                    extent: -5,
                    primer: &primers[1],
                },],
            )
        );
    }
    #[test]
    #[cfg(feature = "parallel")]
    fn test_overlapping_chunks() {
        let data: Vec<u8> = (0..50).collect();
        for i in 3..20 {
            for j in 1..i {
                println!("{:?}", overlapping_chunks(&data, i, j));
                assert_eq!(
                    data,
                    overlapping_chunks(&data, i, j)
                        .iter()
                        .map(|&(offset, c)| if offset == 0 { c } else { &c[(j - 1)..] })
                        .flatten()
                        .cloned()
                        .collect::<Vec<u8>>()
                );
            }
        }
        let text = b"abcdefghijkl";
        let chunks: Vec<_> = overlapping_chunks(&text[..], 5, 3)
            .into_iter()
            .map(|(_, c)| str::from_utf8(c).unwrap())
            .collect();
        for i in 0..text.len() - 3 {
            assert_eq!(
                chunks
                    .iter()
                    .filter(|c| c.contains(str::from_utf8(&text[i..i + 3]).unwrap()))
                    .count(),
                1
            );
        }
    }

    #[test]
    fn circular_footprint_rotation_test() {
        circular_footprint_rotation_test_impl(Method::Bndm);
        #[cfg(not(debug_assertions))] // too slow at this stage because we don't store the index
        circular_footprint_rotation_test_impl(Method::Index);
    }
    fn circular_footprint_rotation_test_impl(method: Method) {
        use gb_io::reader::SeqReader;
        use std::fs::File;
        let lac = SeqReader::new(File::open("tests/circ_test.gb").unwrap())
            .next()
            .unwrap()
            .unwrap();
        let mut lac_double = lac.seq.clone();
        lac_double.extend(lac.seq.iter());
        for n in 14..20 {
            for i in 0..(lac_double.len() - n) {
                let primer = &lac_double[i..i + n];
                let rprimer = revcomp(primer);
                let primers = [primer, &rprimer];
                let matches = Annealer::new(&lac, 14, method).find_matches(&primers);
                assert_eq!(matches.fwd.len(), 1);
                assert_eq!(matches.rev.len(), 1);
                let product = Product(matches.fwd[0].clone(), matches.rev[0].clone()).extract(&lac);
                assert_eq!(product.len() as usize, n);
                assert_eq!(&product.seq, &primer);
            }
        }
    }
}
