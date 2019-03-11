/// Returns a slice containing all of [agctu] which match a given IUPAC code.
/// Input and output must be lowercase.
// Based on the table at https://www.bioinformatics.org/sms/iupac.html
pub fn expand_iupac(code: &u8) -> &[u8] {
    match *code {
        b'a' => b"a",
        b'c' => b"c",
        b'g' => b"g",
        b't' => b"tu",
        b'u' => b"tu",
        b'r' => b"ag",
        b'y' => b"ctu",
        b's' => b"gc",
        b'w' => b"atu",
        b'k' => b"gtu",
        b'm' => b"ac",
        b'b' => b"cgtu",
        b'd' => b"agtu",
        b'h' => b"actu",
        b'v' => b"acg",
        b'n' => b"acgtu",
        _ => {
            warn!(
                "Invalid IUPAC code '{}' encountered!",
                String::from_utf8_lossy(&[*code][..])
            );
            std::slice::from_ref(code)
        }
    }
}

/// Returns true when two NTs are equivalent. IUPAC codes in the primer will be
/// respected, however are not allowed in the template.
pub fn nt_match(template_nt: u8, primer_nt: u8) -> bool {
    let template_nt = template_nt.to_ascii_lowercase();
    expand_iupac(&primer_nt.to_ascii_lowercase())
        .iter()
        .any(|&x| x == template_nt)
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
