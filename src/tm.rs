/* This is a partial port of `tm.py` from the `pydna` python package.
 * Its license:
 *
 * Copyright (c) 2013,2014,2015 Björn Johansson,
 * CBMA (Centro de Biologia Molecular e Ambiental)/Dept of Biology, University of
 * Minho, Braga, Portugal
 * All rights reserved.
 *
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the organizations:
 *       CBMA (Centro de Biologia Molecular e Ambiental) or
 *       Dept of Biology, University of Minho, Braga, Portugal
 *       nor the names of its contributors may be used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL Björn Johansson BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#[derive(PartialEq, Debug)]
pub struct InvalidPrimer;

/// Returns the melting temperature (Tm) of the primer using
/// the nearest neighbour algorithm. Formula and thermodynamic data
/// is taken from SantaLucia 1998 [#]_. This implementation gives the same
/// answer as the one provided by Biopython (See Examples).
/// Thermodynamic data used:
/// =====  ====  ====
/// pair   dH    dS
/// =====  ====  ====
/// AA/TT  7.9   22.2
/// AT/TA  7.2   20.4
/// TA/AT  7.2   21.3
/// CA/GT  8.5   22.7
/// GT/CA  8.4   22.4
/// CT/GA  7.8   21.0
/// GA/CT  8.2   22.2
/// CG/GC  10.6  27.2
/// GC/CG  9.8   24.4
/// GG/CC  8.0   19.9
/// =====  ====  ====
/// Parameters
/// ----------
/// primer : string
///     Primer sequence 5'-3'
/// Returns
/// -------
/// tm : float
///     tm of the primer
/// References
/// ----------
/// .. [#] SantaLucia J Jr. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. Proc Natl Acad Sci U S A 1998;95:1460–5.
///
//pub fn tmstaluc98(primer, dnac=50, saltc=50, **kwargs):
pub fn tmstaluc98(primer: &[u8], dnac: f32, saltc: f32) -> Result<f32, InvalidPrimer> {
    if primer.is_empty() {
        return Err(InvalidPrimer);
    }
    let nntermsl = |x: &[u8]| match x {
        b"AA" => Ok((7.9, 22.2)),
        b"TT" => Ok((7.9, 22.2)),
        b"AT" => Ok((7.2, 20.4)),
        b"TA" => Ok((7.2, 21.3)),
        b"CA" => Ok((8.5, 22.7)),
        b"TG" => Ok((8.5, 22.7)),
        b"GT" => Ok((8.4, 22.4)),
        b"AC" => Ok((8.4, 22.4)),
        b"CT" => Ok((7.8, 21.0)),
        b"AG" => Ok((7.8, 21.0)),
        b"GA" => Ok((8.2, 22.2)),
        b"TC" => Ok((8.2, 22.2)),
        b"CG" => Ok((10.6, 27.2)),
        b"GC" => Ok((9.8, 24.4)),
        b"GG" => Ok((8.0, 19.9)),
        b"CC" => Ok((8.0, 19.9)),
        b"A" => Ok((0.0, 0.0)),
        b"C" => Ok((0.0, 0.0)),
        b"G" => Ok((0.0, 0.0)),
        b"T" => Ok((0.0, 0.0)),
        _ => Err(InvalidPrimer),
    };

    let helixinit = |x| match x {
        b'G' => Ok((-0.1, 2.8)),
        b'C' => Ok((-0.1, 2.8)),
        b'A' => Ok((-2.3, -4.1)),
        b'T' => Ok((-2.3, -4.1)),
        _ => Err(InvalidPrimer),
    };
    let primer = primer.to_ascii_uppercase();
    let (mut dH, mut dS) = helixinit(primer[0])?;
    let (H, S) = helixinit(*primer.last().unwrap())?;
    dH += H;
    dS += S;
    for p in 0..primer.len() {
        let dn = &primer[p..std::cmp::min(p + 2, primer.len())];
        let (H, S) = nntermsl(dn)?;
        dH += H;
        dS += S;
    }
    let R = 1.987; // universal gas constant in Cal/degrees C*Mol
    let k = (dnac / 4.0) * 1e-9;
    dS -= 0.368 * ((primer.len() - 1) as f32) * (saltc / 1e3).ln();
    let tm = ((1000.0 * (-dH)) / (-dS + (R * k.ln()))) - 273.15;
    Ok(tm)
}


/// Returns the tm for a primer using a formula adapted to polymerases
/// with a DNA binding domain, such as the Phusion polymerase.
/// Parameters
/// ----------
/// primer : string
///     primer sequence 5'-3'
/// primerc : float
///    primer concentration in nM), set to 500.0 nm by default.
/// saltc : float, optional
///    Monovalent cation concentration in mM, set to 50.0 mM by default.
pub fn tmbresluc(primer: &[u8], primerc: f32, saltc: f32) -> Result<f32, InvalidPrimer> {
    let saltc = saltc / 1000.0;
    let pri = primerc / 10E7;
    let mut dS = -12.4;
    let mut dH = -3400.0;

    let s = primer.to_ascii_lowercase();

    for i in 0..(s.len() - 1) {
        let n1 = s[i];
        let n2 = s[i + 1];
        dH += dHBr(n1 - 97, n2 - 97)?;
        dS += dSBr(n1 - 97, n2 - 97)?;
    }
    let tm =
        (dH / (1.9872 * (pri / 1600.0).ln() + dS) + (16.6 * saltc.ln()) / (10f32).ln()) - 273.15;
    Ok(tm)
}

fn dSBr(a: u8, b: u8) -> Result<f32, InvalidPrimer> {
    match (a, b) {
        (0, 0) => Ok(-24.0),
        (0, 1) => Ok(-20.7),
        (0, 2) => Ok(-17.3),
        (0, 3) => Ok(-22.9),
        (0, 6) => Ok(-20.8),
        (0, 7) => Ok(-21.7),
        (0, 10) => Ok(-22.4),
        (0, 12) => Ok(-20.7),
        (0, 13) => Ok(-21.5),
        (0, 17) => Ok(-22.4),
        (0, 18) => Ok(-19.1),
        (0, 19) => Ok(-23.9),
        (0, 21) => Ok(-20.7),
        (0, 22) => Ok(-24.0),
        (0, 23) => Ok(-21.5),
        (0, 24) => Ok(-20.6),
        (1, 0) => Ok(-14.4),
        (1, 1) => Ok(-21.8),
        (1, 2) => Ok(-22.3),
        (1, 3) => Ok(-19.2),
        (1, 6) => Ok(-22.4),
        (1, 7) => Ok(-19.1),
        (1, 10) => Ok(-21.6),
        (1, 12) => Ok(-18.4),
        (1, 13) => Ok(-20.0),
        (1, 17) => Ok(-18.4),
        (1, 18) => Ok(-22.4),
        (1, 19) => Ok(-20.7),
        (1, 21) => Ok(-19.7),
        (1, 22) => Ok(-17.6),
        (1, 23) => Ok(-20.0),
        (1, 24) => Ok(-21.5),
        (2, 0) => Ok(-12.9),
        (2, 1) => Ok(-25.1),
        (2, 2) => Ok(-26.6),
        (2, 3) => Ok(-20.5),
        (2, 6) => Ok(-27.8),
        (2, 7) => Ok(-20.1),
        (2, 10) => Ok(-24.3),
        (2, 12) => Ok(-19.8),
        (2, 13) => Ok(-22.0),
        (2, 17) => Ok(-20.4),
        (2, 18) => Ok(-27.2),
        (2, 19) => Ok(-20.8),
        (2, 21) => Ok(-22.4),
        (2, 22) => Ok(-16.9),
        (2, 23) => Ok(-22.0),
        (2, 24) => Ok(-23.7),
        (3, 0) => Ok(-18.1),
        (3, 1) => Ok(-20.3),
        (3, 2) => Ok(-19.2),
        (3, 3) => Ok(-20.0),
        (3, 6) => Ok(-20.1),
        (3, 7) => Ok(-19.7),
        (3, 10) => Ok(-20.9),
        (3, 12) => Ok(-18.7),
        (3, 13) => Ok(-19.8),
        (3, 17) => Ok(-19.1),
        (3, 18) => Ok(-19.6),
        (3, 19) => Ok(-21.7),
        (3, 21) => Ok(-19.1),
        (3, 22) => Ok(-19.9),
        (3, 23) => Ok(-19.8),
        (3, 24) => Ok(-20.5),
        (6, 0) => Ok(-13.5),
        (6, 1) => Ok(-23.5),
        (6, 2) => Ok(-26.7),
        (6, 3) => Ok(-19.1),
        (6, 6) => Ok(-26.6),
        (6, 7) => Ok(-19.2),
        (6, 10) => Ok(-22.0),
        (6, 12) => Ok(-20.1),
        (6, 13) => Ok(-21.0),
        (6, 17) => Ok(-20.1),
        (6, 18) => Ok(-26.7),
        (6, 19) => Ok(-17.3),
        (6, 21) => Ok(-22.3),
        (6, 22) => Ok(-15.4),
        (6, 23) => Ok(-21.0),
        (6, 24) => Ok(-22.0),
        (7, 0) => Ok(-17.9),
        (7, 1) => Ok(-20.8),
        (7, 2) => Ok(-19.1),
        (7, 3) => Ok(-20.4),
        (7, 6) => Ok(-20.5),
        (7, 7) => Ok(-20.0),
        (7, 10) => Ok(-21.7),
        (7, 12) => Ok(-18.5),
        (7, 13) => Ok(-20.1),
        (7, 17) => Ok(-19.2),
        (7, 18) => Ok(-19.8),
        (7, 19) => Ok(-22.9),
        (7, 21) => Ok(-19.2),
        (7, 22) => Ok(-20.4),
        (7, 23) => Ok(-20.1),
        (7, 24) => Ok(-21.0),
        (10, 0) => Ok(-15.2),
        (10, 1) => Ok(-20.2),
        (10, 2) => Ok(-20.1),
        (10, 3) => Ok(-18.5),
        (10, 6) => Ok(-19.8),
        (10, 7) => Ok(-18.7),
        (10, 10) => Ok(-20.2),
        (10, 12) => Ok(-17.7),
        (10, 13) => Ok(-18.9),
        (10, 17) => Ok(-17.5),
        (10, 18) => Ok(-19.9),
        (10, 19) => Ok(-20.7),
        (10, 21) => Ok(-18.4),
        (10, 22) => Ok(-17.9),
        (10, 23) => Ok(-18.9),
        (10, 24) => Ok(-20.4),
        (12, 0) => Ok(-18.5),
        (12, 1) => Ok(-22.9),
        (12, 2) => Ok(-22.0),
        (12, 3) => Ok(-21.7),
        (12, 6) => Ok(-24.3),
        (12, 7) => Ok(-20.9),
        (12, 10) => Ok(-23.3),
        (12, 12) => Ok(-20.2),
        (12, 13) => Ok(-21.8),
        (12, 17) => Ok(-21.4),
        (12, 18) => Ok(-23.1),
        (12, 19) => Ok(-22.4),
        (12, 21) => Ok(-21.6),
        (12, 22) => Ok(-20.4),
        (12, 23) => Ok(-21.8),
        (12, 24) => Ok(-22.2),
        (13, 0) => Ok(-16.8),
        (13, 1) => Ok(-21.5),
        (13, 2) => Ok(-21.0),
        (13, 3) => Ok(-20.1),
        (13, 6) => Ok(-22.0),
        (13, 7) => Ok(-19.8),
        (13, 10) => Ok(-21.8),
        (13, 12) => Ok(-18.9),
        (13, 13) => Ok(-20.3),
        (13, 17) => Ok(-19.4),
        (13, 18) => Ok(-21.5),
        (13, 19) => Ok(-21.5),
        (13, 21) => Ok(-20.0),
        (13, 22) => Ok(-19.2),
        (13, 23) => Ok(-20.3),
        (13, 24) => Ok(-21.3),
        (17, 0) => Ok(-18.8),
        (17, 1) => Ok(-22.1),
        (17, 2) => Ok(-22.0),
        (17, 3) => Ok(-21.0),
        (17, 6) => Ok(-23.7),
        (17, 7) => Ok(-20.5),
        (17, 10) => Ok(-22.2),
        (17, 12) => Ok(-20.4),
        (17, 13) => Ok(-21.3),
        (17, 17) => Ok(-21.2),
        (17, 18) => Ok(-22.9),
        (17, 19) => Ok(-20.6),
        (17, 21) => Ok(-21.5),
        (17, 22) => Ok(-19.7),
        (17, 23) => Ok(-21.3),
        (17, 24) => Ok(-21.3),
        (18, 0) => Ok(-13.2),
        (18, 1) => Ok(-24.3),
        (18, 2) => Ok(-26.7),
        (18, 3) => Ok(-19.8),
        (18, 6) => Ok(-27.2),
        (18, 7) => Ok(-19.6),
        (18, 10) => Ok(-23.1),
        (18, 12) => Ok(-19.9),
        (18, 13) => Ok(-21.5),
        (18, 17) => Ok(-20.2),
        (18, 18) => Ok(-26.9),
        (18, 19) => Ok(-19.1),
        (18, 21) => Ok(-22.4),
        (18, 22) => Ok(-16.1),
        (18, 23) => Ok(-21.5),
        (18, 24) => Ok(-22.9),
        (19, 0) => Ok(-16.9),
        (19, 1) => Ok(-16.8),
        (19, 2) => Ok(-13.5),
        (19, 3) => Ok(-17.9),
        (19, 6) => Ok(-12.9),
        (19, 7) => Ok(-18.1),
        (19, 10) => Ok(-18.5),
        (19, 12) => Ok(-15.2),
        (19, 13) => Ok(-16.8),
        (19, 17) => Ok(-14.9),
        (19, 18) => Ok(-13.2),
        (19, 19) => Ok(-24.0),
        (19, 21) => Ok(-14.4),
        (19, 22) => Ok(-20.5),
        (19, 23) => Ok(-16.8),
        (19, 24) => Ok(-18.8),
        (21, 0) => Ok(-16.8),
        (21, 1) => Ok(-23.1),
        (21, 2) => Ok(-23.5),
        (21, 3) => Ok(-20.8),
        (21, 6) => Ok(-25.1),
        (21, 7) => Ok(-20.3),
        (21, 10) => Ok(-22.9),
        (21, 12) => Ok(-20.2),
        (21, 13) => Ok(-21.5),
        (21, 17) => Ok(-20.9),
        (21, 18) => Ok(-24.3),
        (21, 19) => Ok(-20.7),
        (21, 21) => Ok(-21.8),
        (21, 22) => Ok(-18.7),
        (21, 23) => Ok(-21.5),
        (21, 24) => Ok(-22.1),
        (22, 0) => Ok(-20.5),
        (22, 1) => Ok(-18.7),
        (22, 2) => Ok(-15.4),
        (22, 3) => Ok(-20.4),
        (22, 6) => Ok(-16.9),
        (22, 7) => Ok(-19.9),
        (22, 10) => Ok(-20.4),
        (22, 12) => Ok(-17.9),
        (22, 13) => Ok(-19.2),
        (22, 17) => Ok(-18.7),
        (22, 18) => Ok(-16.1),
        (22, 19) => Ok(-24.0),
        (22, 21) => Ok(-17.6),
        (22, 22) => Ok(-22.2),
        (22, 23) => Ok(-19.2),
        (22, 24) => Ok(-19.7),
        (23, 0) => Ok(-16.8),
        (23, 1) => Ok(-21.5),
        (23, 2) => Ok(-21.0),
        (23, 3) => Ok(-20.1),
        (23, 6) => Ok(-22.0),
        (23, 7) => Ok(-19.8),
        (23, 10) => Ok(-21.8),
        (23, 12) => Ok(-18.9),
        (23, 13) => Ok(-20.3),
        (23, 17) => Ok(-19.4),
        (23, 18) => Ok(-21.5),
        (23, 19) => Ok(-21.5),
        (23, 21) => Ok(-20.0),
        (23, 22) => Ok(-19.2),
        (23, 23) => Ok(-20.3),
        (23, 24) => Ok(-21.3),
        (24, 0) => Ok(-14.9),
        (24, 1) => Ok(-20.9),
        (24, 2) => Ok(-20.1),
        (24, 3) => Ok(-19.2),
        (24, 6) => Ok(-20.4),
        (24, 7) => Ok(-19.1),
        (24, 10) => Ok(-21.4),
        (24, 12) => Ok(-17.5),
        (24, 13) => Ok(-19.4),
        (24, 17) => Ok(-17.6),
        (24, 18) => Ok(-20.2),
        (24, 19) => Ok(-22.4),
        (24, 21) => Ok(-18.4),
        (24, 22) => Ok(-18.7),
        (24, 23) => Ok(-19.4),
        (24, 24) => Ok(-21.2),
        _ => Err(InvalidPrimer),
    }
}

fn dHBr(a: u8, b: u8) -> Result<f32, InvalidPrimer> {
    match (a, b) {
        (0, 0) => Ok(-9100.0),
        (0, 1) => Ok(-7633.3),
        (0, 2) => Ok(-6500.0),
        (0, 3) => Ok(-8500.0),
        (0, 6) => Ok(-7800.0),
        (0, 7) => Ok(-8066.7),
        (0, 10) => Ok(-8200.0),
        (0, 12) => Ok(-7800.0),
        (0, 13) => Ok(-8000.0),
        (0, 17) => Ok(-8450.0),
        (0, 18) => Ok(-7150.0),
        (0, 19) => Ok(-8600.0),
        (0, 21) => Ok(-7800.0),
        (0, 22) => Ok(-8850.0),
        (0, 23) => Ok(-8000.0),
        (0, 24) => Ok(-7550.0),
        (1, 0) => Ok(-5800.0),
        (1, 1) => Ok(-8866.7),
        (1, 2) => Ok(-9233.3),
        (1, 3) => Ok(-7722.2),
        (1, 6) => Ok(-9566.7),
        (1, 7) => Ok(-7611.1),
        (1, 10) => Ok(-8683.3),
        (1, 12) => Ok(-7516.7),
        (1, 13) => Ok(-8100.0),
        (1, 17) => Ok(-7683.3),
        (1, 18) => Ok(-9400.0),
        (1, 19) => Ok(-7800.0),
        (1, 21) => Ok(-8200.0),
        (1, 22) => Ok(-6800.0),
        (1, 23) => Ok(-8100.0),
        (1, 24) => Ok(-8516.7),
        (2, 0) => Ok(-5800.0),
        (2, 1) => Ok(-10233.3),
        (2, 2) => Ok(-11000.0),
        (2, 3) => Ok(-8500.0),
        (2, 6) => Ok(-11900.0),
        (2, 7) => Ok(-8200.0),
        (2, 10) => Ok(-9850.0),
        (2, 12) => Ok(-8400.0),
        (2, 13) => Ok(-9125.0),
        (2, 17) => Ok(-8850.0),
        (2, 18) => Ok(-11450.0),
        (2, 19) => Ok(-7800.0),
        (2, 21) => Ok(-9566.7),
        (2, 22) => Ok(-6800.0),
        (2, 23) => Ok(-9125.0),
        (2, 24) => Ok(-9400.0),
        (3, 0) => Ok(-6900.0),
        (3, 1) => Ok(-8000.0),
        (3, 2) => Ok(-7733.3),
        (3, 3) => Ok(-7722.2),
        (3, 6) => Ok(-8200.0),
        (3, 7) => Ok(-7566.7),
        (3, 10) => Ok(-8133.3),
        (3, 12) => Ok(-7316.7),
        (3, 13) => Ok(-7725.0),
        (3, 17) => Ok(-7550.0),
        (3, 18) => Ok(-7966.7),
        (3, 19) => Ok(-8066.7),
        (3, 21) => Ok(-7611.1),
        (3, 22) => Ok(-7483.3),
        (3, 23) => Ok(-7725.0),
        (3, 24) => Ok(-7900.0),
        (6, 0) => Ok(-5600.0),
        (6, 1) => Ok(-9533.3),
        (6, 2) => Ok(-11100.0),
        (6, 3) => Ok(-7700.0),
        (6, 6) => Ok(-11000.0),
        (6, 7) => Ok(-7733.3),
        (6, 10) => Ok(-8750.0),
        (6, 12) => Ok(-8350.0),
        (6, 13) => Ok(-8550.0),
        (6, 17) => Ok(-8300.0),
        (6, 18) => Ok(-11050.0),
        (6, 19) => Ok(-6500.0),
        (6, 21) => Ok(-9233.3),
        (6, 22) => Ok(-6050.0),
        (6, 23) => Ok(-8550.0),
        (6, 24) => Ok(-8800.0),
        (7, 0) => Ok(-6966.7),
        (7, 1) => Ok(-8233.3),
        (7, 2) => Ok(-7700.0),
        (7, 3) => Ok(-7988.9),
        (7, 6) => Ok(-8500.0),
        (7, 7) => Ok(-7722.2),
        (7, 10) => Ok(-8500.0),
        (7, 12) => Ok(-7333.3),
        (7, 13) => Ok(-7916.7),
        (7, 17) => Ok(-7733.3),
        (7, 18) => Ok(-8100.0),
        (7, 19) => Ok(-8500.0),
        (7, 21) => Ok(-7722.2),
        (7, 22) => Ok(-7733.3),
        (7, 23) => Ok(-7916.7),
        (7, 24) => Ok(-8100.0),
        (10, 0) => Ok(-5800.0),
        (10, 1) => Ok(-8183.3),
        (10, 2) => Ok(-8350.0),
        (10, 3) => Ok(-7333.3),
        (10, 6) => Ok(-8400.0),
        (10, 7) => Ok(-7316.7),
        (10, 10) => Ok(-8100.0),
        (10, 12) => Ok(-7075.0),
        (10, 13) => Ok(-7587.5),
        (10, 17) => Ok(-7100.0),
        (10, 18) => Ok(-8375.0),
        (10, 19) => Ok(-7800.0),
        (10, 21) => Ok(-7516.7),
        (10, 22) => Ok(-6800.0),
        (10, 23) => Ok(-7587.5),
        (10, 24) => Ok(-8075.0),
        (12, 0) => Ok(-7450.0),
        (12, 1) => Ok(-8933.3),
        (12, 2) => Ok(-8750.0),
        (12, 3) => Ok(-8500.0),
        (12, 6) => Ok(-9850.0),
        (12, 7) => Ok(-8133.3),
        (12, 10) => Ok(-9025.0),
        (12, 12) => Ok(-8100.0),
        (12, 13) => Ok(-8562.5),
        (12, 17) => Ok(-8650.0),
        (12, 18) => Ok(-9300.0),
        (12, 19) => Ok(-8200.0),
        (12, 21) => Ok(-8683.3),
        (12, 22) => Ok(-7825.0),
        (12, 23) => Ok(-8562.5),
        (12, 24) => Ok(-8475.0),
        (13, 0) => Ok(-6625.0),
        (13, 1) => Ok(-8558.3),
        (13, 2) => Ok(-8550.0),
        (13, 3) => Ok(-7916.7),
        (13, 6) => Ok(-9125.0),
        (13, 7) => Ok(-7725.0),
        (13, 10) => Ok(-8562.5),
        (13, 12) => Ok(-7587.5),
        (13, 13) => Ok(-8075.0),
        (13, 17) => Ok(-7875.0),
        (13, 18) => Ok(-8837.5),
        (13, 19) => Ok(-8000.0),
        (13, 21) => Ok(-8100.0),
        (13, 22) => Ok(-7312.5),
        (13, 23) => Ok(-8075.0),
        (13, 24) => Ok(-8275.0),
        (17, 0) => Ok(-7350.0),
        (17, 1) => Ok(-8583.3),
        (17, 2) => Ok(-8800.0),
        (17, 3) => Ok(-8100.0),
        (17, 6) => Ok(-9400.0),
        (17, 7) => Ok(-7900.0),
        (17, 10) => Ok(-8475.0),
        (17, 12) => Ok(-8075.0),
        (17, 13) => Ok(-8275.0),
        (17, 17) => Ok(-8375.0),
        (17, 18) => Ok(-9100.0),
        (17, 19) => Ok(-7550.0),
        (17, 21) => Ok(-8516.7),
        (17, 22) => Ok(-7450.0),
        (17, 23) => Ok(-8275.0),
        (17, 24) => Ok(-8175.0),
        (18, 0) => Ok(-5700.0),
        (18, 1) => Ok(-9883.3),
        (18, 2) => Ok(-11050.0),
        (18, 3) => Ok(-8100.0),
        (18, 6) => Ok(-11450.0),
        (18, 7) => Ok(-7966.7),
        (18, 10) => Ok(-9300.0),
        (18, 12) => Ok(-8375.0),
        (18, 13) => Ok(-8837.5),
        (18, 17) => Ok(-8575.0),
        (18, 18) => Ok(-11250.0),
        (18, 19) => Ok(-7150.0),
        (18, 21) => Ok(-9400.0),
        (18, 22) => Ok(-6425.0),
        (18, 23) => Ok(-8837.5),
        (18, 24) => Ok(-9100.0),
        (19, 0) => Ok(-6000.0),
        (19, 1) => Ok(-6833.3),
        (19, 2) => Ok(-5600.0),
        (19, 3) => Ok(-6966.7),
        (19, 6) => Ok(-5800.0),
        (19, 7) => Ok(-6900.0),
        (19, 10) => Ok(-7450.0),
        (19, 12) => Ok(-5800.0),
        (19, 13) => Ok(-6625.0),
        (19, 17) => Ok(-5900.0),
        (19, 18) => Ok(-5700.0),
        (19, 19) => Ok(-9100.0),
        (19, 21) => Ok(-5800.0),
        (19, 22) => Ok(-7550.0),
        (19, 23) => Ok(-6625.0),
        (19, 24) => Ok(-7350.0),
        (21, 0) => Ok(-6833.3),
        (21, 1) => Ok(-9133.3),
        (21, 2) => Ok(-9533.3),
        (21, 3) => Ok(-8233.3),
        (21, 6) => Ok(-10233.3),
        (21, 7) => Ok(-8000.0),
        (21, 10) => Ok(-8933.3),
        (21, 12) => Ok(-8183.3),
        (21, 13) => Ok(-8558.3),
        (21, 17) => Ok(-8533.3),
        (21, 18) => Ok(-9883.3),
        (21, 19) => Ok(-7633.3),
        (21, 21) => Ok(-8866.7),
        (21, 22) => Ok(-7233.3),
        (21, 23) => Ok(-8558.3),
        (21, 24) => Ok(-8583.3),
        (22, 0) => Ok(-7550.0),
        (22, 1) => Ok(-7233.3),
        (22, 2) => Ok(-6050.0),
        (22, 3) => Ok(-7733.3),
        (22, 6) => Ok(-6800.0),
        (22, 7) => Ok(-7483.3),
        (22, 10) => Ok(-7825.0),
        (22, 12) => Ok(-6800.0),
        (22, 13) => Ok(-7312.5),
        (22, 17) => Ok(-7175.0),
        (22, 18) => Ok(-6425.0),
        (22, 19) => Ok(-8850.0),
        (22, 21) => Ok(-6800.0),
        (22, 22) => Ok(-8200.0),
        (22, 23) => Ok(-7312.5),
        (22, 24) => Ok(-7450.0),
        (23, 0) => Ok(-6625.0),
        (23, 1) => Ok(-8558.3),
        (23, 2) => Ok(-8550.0),
        (23, 3) => Ok(-7916.7),
        (23, 6) => Ok(-9125.0),
        (23, 7) => Ok(-7725.0),
        (23, 10) => Ok(-8562.5),
        (23, 12) => Ok(-7587.5),
        (23, 13) => Ok(-8075.0),
        (23, 17) => Ok(-7875.0),
        (23, 18) => Ok(-8837.5),
        (23, 19) => Ok(-8000.0),
        (23, 21) => Ok(-8100.0),
        (23, 22) => Ok(-7312.5),
        (23, 23) => Ok(-8075.0),
        (23, 24) => Ok(-8275.0),
        (24, 0) => Ok(-5900.0),
        (24, 1) => Ok(-8533.3),
        (24, 2) => Ok(-8300.0),
        (24, 3) => Ok(-7733.3),
        (24, 6) => Ok(-8850.0),
        (24, 7) => Ok(-7550.0),
        (24, 10) => Ok(-8650.0),
        (24, 12) => Ok(-7100.0),
        (24, 13) => Ok(-7875.0),
        (24, 17) => Ok(-7375.0),
        (24, 18) => Ok(-8575.0),
        (24, 19) => Ok(-8450.0),
        (24, 21) => Ok(-7683.3),
        (24, 22) => Ok(-7175.0),
        (24, 23) => Ok(-7875.0),
        (24, 24) => Ok(-8375.0),
        _ => Err(InvalidPrimer),
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use assert_approx_eq::assert_approx_eq;
    #[test]
    fn tm() {
        let primer = b"agatcgactatctatcttatgcactatgtctat";
        assert_approx_eq!(tmbresluc(primer, 500.0, 50.0).unwrap(), 63.598602f32);
        assert_approx_eq!(tmstaluc98(primer, 50.0, 50.0).unwrap(), 55.08255f32);
    }
}