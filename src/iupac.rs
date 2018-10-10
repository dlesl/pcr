/// Returns a slice containing all of [agctu] which match a given IUPAC code.
/// Uppercase input is allowed, but the return value will be lowercase.
// Based on the table at https://www.bioinformatics.org/sms/iupac.html
pub fn expand_iupac(code: u8) -> Option<&'static [u8]> {
    match code.to_ascii_lowercase() {
        b'a' => Some(b"a"),
        b'c' => Some(b"c"),
        b'g' => Some(b"g"),
        b't' => Some(b"tu"),
        b'u' => Some(b"tu"),
        b'r' => Some(b"ag"),
        b'y' => Some(b"ctu"),
        b's' => Some(b"gc"),
        b'w' => Some(b"atu"),
        b'k' => Some(b"gtu"),
        b'm' => Some(b"ac"),
        b'b' => Some(b"cgtu"),
        b'd' => Some(b"agtu"),
        b'h' => Some(b"actu"),
        b'v' => Some(b"acg"),
        b'n' => Some(b"acgtu"),
        _ => {
            warn!(
                "Invalid IUPAC code '{}' encountered!",
                String::from_utf8_lossy(&[code][..])
            );
            None
        }
    }
}

/// Returns true when two NTs are equivalent. IUPAC codes in the primer will be
/// respected, however are not allowed in the template.
pub fn nt_match(template_nt: u8, primer_nt: u8) -> bool {
    if let Some(nts) = expand_iupac(primer_nt) {
        nts.iter().any(|&x| x == template_nt.to_ascii_lowercase())
    } else {
        primer_nt.to_ascii_lowercase() == template_nt.to_ascii_lowercase()
    }
}

#[test]
fn test_nt_match() {
    assert!(nt_match(b'a', b'a'));
    assert!(nt_match(b'a', b'A'));
    assert!(nt_match(b'A', b'a'));
    assert!(nt_match(b'u', b't'));
    assert!(nt_match(b't', b'u'));
    assert!(nt_match(b'a', b'n'));
    assert!(!nt_match(b'n', b'a'));
    assert!(nt_match(b'c', b'm'));
}
