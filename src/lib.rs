extern crate bio;
#[macro_use]
extern crate gb_io;
extern crate itertools;
#[macro_use]
extern crate log;
use bio::alignment::sparse;
use bio::alphabets::dna::{self, revcomp};
use std::borrow::Borrow;
use std::fmt;
use std::mem;
use std::str;

use gb_io::seq::{Feature, QualifierKey, Seq};

mod bndm_iupac;
use bndm_iupac::BNDM;
pub mod iupac;
use iupac::nt_match;

#[cfg(feature = "parallel")]
extern crate rayon;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

// TODO: make this configurable
const MAX_PRODUCTS: usize = 10000;

#[derive(PartialEq, Eq, PartialOrd, Ord, Hash, Clone)]
pub struct Primer {
    seq: Vec<u8>,
    seq_rc: Vec<u8>,
    name: String,
}

impl Primer {
    pub fn new(name: String, seq: Vec<u8>) -> Primer {
        Primer {
            seq_rc: revcomp(&seq),
            seq: seq,
            name,
        }
    }
    pub fn seq(&self) -> &[u8] {
        &self.seq
    }
    pub fn seq_rc(&self) -> &[u8] {
        &self.seq_rc
    }
    pub fn name(&self) -> &str {
        self.name.as_str()
    }
}

impl fmt::Debug for Primer {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Primer {{ seq: \"{}\", seq_rc: \"{}\" }}",
            String::from_utf8_lossy(&self.seq),
            String::from_utf8_lossy(&self.seq_rc)
        )
    }
}

impl<'a> From<&'a [u8]> for Primer {
    fn from(i: &'a [u8]) -> Primer {
        Primer::new(String::new(), i.into())
    }
}

impl From<Vec<u8>> for Primer {
    fn from(i: Vec<u8>) -> Primer {
        Primer::new(String::new(), i)
    }
}

/// A trait to allow client code the option of using Primer, &Primer,
/// Rc<Primer>, ... when dealing with `Footprint`. Must be `Clone`, since
/// multiple `Footprint`s can be found for one primer
#[cfg(feature = "parallel")]
pub trait PrimerRef: Borrow<Primer> + Clone + Sync + Send {
    fn seq(&self) -> &[u8];
    fn seq_rc(&self) -> &[u8];
    fn name(&self) -> &str;
    fn len(&self) -> i64;
}

#[cfg(feature = "parallel")]
impl<T> PrimerRef for T
where
    T: Borrow<Primer> + Clone + Sync + Send,
{
    fn seq(&self) -> &[u8] {
        &self.borrow().seq()
    }
    fn seq_rc(&self) -> &[u8] {
        &self.borrow().seq_rc()
    }
    fn name(&self) -> &str {
        self.borrow().name()
    }
    fn len(&self) -> i64 {
        self.seq().len() as i64
    }
}

#[cfg(not(feature = "parallel"))]
pub trait PrimerRef: Borrow<Primer> + Clone {
    fn seq(&self) -> &[u8];
    fn seq_rc(&self) -> &[u8];
    fn name(&self) -> &str;
    fn len(&self) -> i64;
}

#[cfg(not(feature = "parallel"))]
impl<T> PrimerRef for T
where
    T: Borrow<Primer> + Clone,
{
    fn seq(&self) -> &[u8] {
        &self.borrow().seq()
    }
    fn seq_rc(&self) -> &[u8] {
        &self.borrow().seq_rc()
    }
    fn name(&self) -> &str {
        self.borrow().name()
    }
    fn len(&self) -> i64 {
        self.seq().len() as i64
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct Footprint<T: PrimerRef> {
    pub start: i64,
    pub extent: i64, // - negative for other strand
    pub primer: T,
}

impl<T: PrimerRef> Footprint<T> {
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

pub fn find_matches<T: PrimerRef>(
    template: &Seq,
    primers: &[T],
    min_fp: i64,
    method: Method,
) -> (Vec<Footprint<T>>, Vec<Footprint<T>>) {
    // Check that primers are appropriate:
    // - long enough
    // - no IUPAC codes in seed region if method == Index
    let primers: Vec<&T> = primers
        .iter()
        .filter(|p| {
            if p.seq().len() < min_fp as usize {
                warn!("Primer '{}' too short, skipping", p.name());
                false
            } else {
                match method {
                    Method::Bndm => true,
                    Method::Index => {
                        if !dna::alphabet().is_word(&p.seq()[p.seq().len() - min_fp as usize..]) {
                            warn!("Invalid bases in primer '{}', method 'Index' doesn't support IUPAC codes, skipping", p.name());
                            false
                        } else {
                            true
                        }
                    }
                }
            }
        })
        .collect();
    let seq = &template.seq;
    let mut res: (Vec<_>, Vec<_>) = match method {
        Method::Bndm => {
            #[cfg(feature = "parallel")]
            let res: (Vec<_>, Vec<_>) = rayon::join(
                || {
                    primers
                        .par_iter()
                        .flat_map(|&p| search_bndm_fwd(template, p, seq, 0, min_fp).into_par_iter())
                        .collect()
                },
                || {
                    primers
                        .par_iter()
                        .flat_map(|&p| search_bndm_rev(template, p, seq, 0, min_fp).into_par_iter())
                        .collect()
                },
            );
            #[cfg(not(feature = "parallel"))]
            let res: (Vec<_>, Vec<_>) = (
                primers
                    .iter()
                    .flat_map(|&p| search_bndm_fwd(template, p, seq, 0, min_fp).into_iter())
                    .collect(),
                primers
                    .iter()
                    .flat_map(|&p| search_bndm_rev(template, p, seq, 0, min_fp).into_iter())
                    .collect(),
            );
            res
        }
        Method::Index => {
            let seq_upper = template.seq.to_ascii_uppercase();
            #[cfg(not(feature = "parallel"))]
            let index = vec![(0, sparse::hash_kmers(&seq_upper, min_fp as usize))];
            #[cfg(feature = "parallel")]
            let index: Vec<_> = overlapping_chunks(
                &seq_upper,
                template.seq.len() / rayon::current_num_threads(),
                min_fp as usize,
            )
            .into_par_iter()
            .map(|(offset, chunk)| (offset, sparse::hash_kmers(chunk, min_fp as usize)))
            .collect();
            let mut fwd = Vec::new();
            let mut rev = Vec::new();
            for (offset, index) in index {
                fwd.extend(primers.iter().flat_map(|&p| {
                    sparse::find_kmer_matches_seq1_hashed(
                        &index,
                        &p.seq()[p.seq().len() - min_fp as usize..].to_ascii_uppercase(),
                        min_fp as usize,
                    )
                    .into_iter()
                    .map(|(m, _)| extend_fwd(template, p.clone(), m as i64 + offset as i64, min_fp))
                    .collect::<Vec<_>>()
                }));
                rev.extend(primers.iter().flat_map(|&p| {
                    sparse::find_kmer_matches_seq1_hashed(
                        &index,
                        &p.seq_rc()[..min_fp as usize].to_ascii_uppercase(),
                        min_fp as usize,
                    )
                    .into_iter()
                    .map(|(m, _)| extend_rev(template, p.clone(), m as i64 + offset as i64, min_fp))
                    .collect::<Vec<_>>()
                }));
            }
            (fwd, rev)
        }
    };
    // collect matches in "origin" region
    if template.is_circular() {
        // one less so that we don't find the same match twice
        let start = seq.len() as i64 - min_fp + 1;
        let end = seq.len() as i64 + min_fp - 1;
        let origin = template.extract_range_seq(start, end);
        assert!(origin.len() as i64 == min_fp * 2 - 2);
        res.0.extend(
            primers
                .iter()
                .flat_map(|&p| search_bndm_fwd(template, p, &origin, start, min_fp)),
        );
        res.1.extend(
            primers
                .iter()
                .flat_map(|&p| search_bndm_rev(template, p, &origin, start, min_fp)),
        );
    }
    res
}

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

// search for matches, for optimisation purposes, the slice to search and
// the position it starts at is provided so we can reuse it

fn search_bndm_fwd<T: PrimerRef>(
    template: &Seq,
    primer: &T,
    seq: &[u8],
    starts_at: i64,
    min_fp: i64,
) -> Vec<Footprint<T>> {
    let three_prime = &primer.seq()[primer.seq().len() - min_fp as usize..];
    BNDM::new(three_prime)
        .find_all(seq)
        .map(|m| extend_fwd(template, primer.clone(), m as i64 + starts_at, min_fp))
        .collect()
}

fn search_bndm_rev<T: PrimerRef>(
    template: &Seq,
    primer: &T,
    seq: &[u8],
    starts_at: i64,
    min_fp: i64,
) -> Vec<Footprint<T>> {
    let three_prime = &primer.seq_rc()[..min_fp as usize];
    BNDM::new(three_prime)
        .find_all(seq)
        .map(|m| extend_rev(template, primer.clone(), m as i64 + starts_at, min_fp))
        .collect()
}

fn extend_fwd<T: PrimerRef>(template: &Seq, primer: T, pos: i64, min_fp: i64) -> Footprint<T> {
    let primer_end = pos + min_fp;
    let mut first = primer_end - primer.seq().len() as i64;
    if !template.is_circular() {
        first = std::cmp::max(first, 0);
    }
    let template_ext = template.extract_range_seq(first, pos);
    let count = template_ext
        .iter()
        .rev()
        .zip(
            primer.seq()[..primer.seq().len() - min_fp as usize]
                .iter()
                .rev(),
        )
        .take_while(|&(&a, &b)| nt_match(a, b))
        .count() as i64;
    let start = pos - count;
    let start_wrapped = if start < 0 {
        assert!(template.is_circular());
        start + template.seq.len() as i64
    } else {
        start
    };
    Footprint {
        start: start_wrapped,
        extent: min_fp + count,
        primer,
    }
}

fn extend_rev<T: PrimerRef>(template: &Seq, primer: T, pos: i64, min_fp: i64) -> Footprint<T> {
    let mut primer_end = pos + primer.seq().len() as i64;
    if !template.is_circular() {
        primer_end = std::cmp::min(primer_end, template.len());
    }
    let template_ext = template.extract_range_seq(pos + min_fp, primer_end);
    let count = template_ext
        .iter()
        .zip(primer.seq_rc()[min_fp as usize..].iter())
        .take_while(|&(&a, &b)| nt_match(a, b))
        .count() as i64;
    let start = pos - 1 + min_fp + count;
    let len = template.seq.len() as i64;
    let start_wrapped = if start >= len {
        assert!(template.is_circular());
        start - len
    } else {
        start
    };
    Footprint {
        start: start_wrapped,
        extent: -(min_fp + count),
        primer,
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Product<T: PrimerRef>(pub Footprint<T>, pub Footprint<T>);

impl<T: PrimerRef> Product<T> {
    pub fn len(&self, template: &Seq) -> i64 {
        product_len(template, &self.0, &self.1).expect("Tried to get length of invalid Product")
    }
    pub fn extract(&self, template: &Seq) -> Seq {
        let &Product(ref fwd, ref rev) = self;
        if fwd.start + fwd.extent < rev.start + rev.extent || template.is_circular() {
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
            unimplemented!(); // product is really a primer dimer?
        }
    }
}

fn product_len<T: PrimerRef>(
    template: &Seq,
    fwd: &Footprint<T>,
    rev: &Footprint<T>,
) -> Option<i64> {
    let overhangs = (fwd.primer.len() - fwd.len()) + (rev.primer.len() - rev.len());
    if fwd.start <= rev.start {
        Some(rev.start - fwd.start + 1 + overhangs)
    } else if rev.start < fwd.start && template.is_circular() {
        let rev_start_unwrapped = rev.start + template.len();
        let distance = rev_start_unwrapped - fwd.start + 1;
        Some(distance + overhangs)
    } else {
        None
    }
}

pub fn find_products<T: PrimerRef>(
    template: &Seq,
    fwd_matches: &[Footprint<T>],
    rev_matches: &[Footprint<T>],
    min_len: i64,
    max_len: i64,
) -> Vec<Product<T>> {
    let mut products = Vec::new();
    for fwd in fwd_matches.iter() {
        for rev in rev_matches.iter() {
            if let Some(len) = product_len(template, fwd, rev) {
                if len >= min_len && len <= max_len {
                    products.push(Product(fwd.clone(), rev.clone())); //TODO: unwrap
                    assert!(
                        products.len() <= MAX_PRODUCTS,
                        "Maximum number of products ({}) exceeded!",
                        MAX_PRODUCTS
                    );
                }
            }
        }
    }
    products
}

pub fn annotate<T: PrimerRef>(mut res: Seq, products: Vec<Product<T>>) -> Seq {
    for Product(fwd, rev) in products {
        let f = Feature {
            kind: feature_kind!("misc_feature"),
            pos: res.range_to_position(fwd.start, rev.start + 1),
            qualifiers: vec![(
                QualifierKey::from("PCR_primers"),
                Some(format!(
                    "[fwd_name: {}, ]fwd_seq: {},\n\
                     [rev_name: {}, ]rev_seq: {}",
                    fwd.primer.name(),
                    String::from_utf8_lossy(fwd.primer.seq()),
                    rev.primer.name(),
                    String::from_utf8_lossy(rev.primer.seq())
                )),
            )],
        };
        res.features.push(f);
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
        let primer = Primer::from(&b"245678"[..]);
        let min_fp = 3;
        assert_eq!(
            extend_fwd(&s, &primer, 6, min_fp),
            Footprint {
                start: 4,
                extent: primer.seq.len() as i64 - 1,
                primer: &primer,
            }
        );
        //circular
        let s_circ = Seq {
            topology: Topology::Circular,
            ..s.clone()
        };
        let primer = Primer::from(&b"890123"[..]);
        assert_eq!(
            extend_fwd(&s_circ, &primer, 1, min_fp),
            Footprint {
                start: 8,
                extent: primer.seq.len() as i64,
                primer: &primer,
            }
        );
        // revcomp leaves invalid characters alone so this is fine
        let primer = Primer::from(&b"665432"[..]);
        assert_eq!(
            extend_rev(&s, &primer, 2, min_fp),
            Footprint {
                start: 6,
                extent: -(primer.seq.len() as i64) + 1,
                primer: &primer,
            }
        );
        let primer = Primer::from(&b"65432109"[..]);
        assert_eq!(
            extend_rev(&s_circ, &primer, 9, min_fp),
            Footprint {
                start: 6,
                extent: -(primer.seq.len() as i64),
                primer: &primer,
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

        let primers_owned = [
            Primer::from(&b"GGGAAATCCTCAAGCACCAG"[..]),
            Primer::from(&b"TTGATTGCGTAATCAGCACCAC"[..]),
        ];
        let primers: Vec<_> = primers_owned.iter().collect();
        assert_eq!(
            find_matches(&s, &primers, 14, method),
            (
                vec![Footprint {
                    start: 1,
                    extent: 19,
                    primer: primers[0],
                },],
                vec![Footprint {
                    start: 201,
                    extent: -21,
                    primer: primers[1],
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

        let primers_owned = [Primer::from(&b"AGACT"[..]), Primer::from(&b"CTAAT"[..])];
        let primers: Vec<_> = primers_owned.iter().collect();
        assert_eq!(
            find_matches(&s, &primers, 4, method),
            (
                vec![Footprint {
                    start: 7,
                    extent: 5,
                    primer: primers[0],
                },],
                vec![Footprint {
                    start: 0,
                    extent: -5,
                    primer: primers[1],
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
                let primer = Primer::from(&lac_double[i..i + n]);
                let rprimer = Primer::from(&primer.seq_rc[..]);
                let primers = [&primer, &rprimer];
                let (fwd, rev) = find_matches(&lac, &primers, 14, method);
                println!("{:?}, {:?}", fwd, rev);
                assert_eq!(fwd.len(), 1);
                assert_eq!(rev.len(), 1);
                let product = lac.extract_range_seq(fwd[0].start, rev[0].start + 1);
                assert_eq!(product.len(), n);
                assert_eq!(product, &primer.seq[..]);
            }
        }
    }
}
