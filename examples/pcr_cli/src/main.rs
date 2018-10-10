#[macro_use]
extern crate clap;
extern crate env_logger;
#[macro_use]
extern crate log;
extern crate gb_io;
extern crate pcr;
#[macro_use]
extern crate failure;

use failure::Error;
use std::fs::File;
use std::io;
use std::io::Write;
use std::process::exit;

use gb_io::reader::SeqReader;
use pcr::{Method, Primer};

fn main() {
    match run() {
        Ok(_) => exit(0),
        Err(ref e) => {
            writeln!(&mut io::stderr(), "Error: {}", e).ok();
            exit(1)
        }
    }
}

fn run() -> Result<(), Error> {
    env_logger::init();
    let args = clap_app!(pcr =>
        (about: "Perform PCR on a genbank file using a list of primers.")
        (version: env!("CARGO_PKG_VERSION"))
        (@arg TEMPLATE: +required "Genbank file to use as template")
        (@arg PRIMERS: +required ... "Primer(s) to use")
        // (@arg outdir: -o --outdir +takes_value "Optional: output directory")
        (@arg matches_only: -m --("matches") "Print primer binding sites only")
        (@arg annotate: -a --("annotate") "Annotate the original sequence with the products")
        (@arg min_homology: -b --("min-homology") default_value("14") "Minimum homology required at 3' end of primer")
        (@arg min_len: -l --("min-len") default_value("100") "Minimum length of product")
        (@arg max_len: -L --("max-len") default_value("10000") "Maximum length of product")
        (@arg index: -i --("index") "Index template before searching (faster for large numbers of primers, uses (much) more memory)")
    ).get_matches();
    let template_file = value_t_or_exit!(args.value_of("TEMPLATE"), String);
    let min_homology = value_t_or_exit!(args.value_of("min_homology"), usize) as i64;
    if min_homology > 64 {
        return Err(format_err!("min-binding must be <= 64"));
    }
    let min_len = value_t_or_exit!(args.value_of("min_len"), usize) as i64;
    let max_len = value_t_or_exit!(args.value_of("max_len"), usize) as i64;
    let primers_in: Vec<_> = values_t_or_exit!(args.values_of("PRIMERS"), String);
    let mut primers = vec![];
    for p in primers_in {
        if p.len() < min_homology as usize {
            return Err(format_err!(
                "Primer `{}` shorter than minimum homology required ({} nt)",
                p,
                min_homology
            ));
        }
        primers.push(Primer::from(p.as_bytes()))
    }

    let method = if args.is_present("index") {
        Method::Index
    } else {
        Method::Bndm
    };

    debug!("loading...");
    let f = File::open(template_file)?;

    for record in SeqReader::new(f) {
        let mut record = record?;
        let primers: Vec<_> = primers.iter().collect();
        debug!(
            "searching {}...",
            record.name.as_ref().map_or("-", |x| &**x)
        );
        let (fwd, rev) = pcr::find_matches(&record, &primers, min_homology, method);
        if args.is_present("matches_only") {
            println!("Forward matches:");
            for m in fwd {
                println!("{}, {}", m.start, m.extent);
            }
            println!("Reverse matches:");
            for m in rev {
                println!("{}, {}", m.start, m.extent);
            }
        } else if args.is_present("annotate") {
            let products = pcr::find_products(&record, &fwd, &rev, min_len, max_len);
            let res = pcr::annotate(record, products);
            res.write(io::stdout())?;
        } else {
            let mut stdout = io::stdout();
            for p in pcr::find_products(&record, &fwd, &rev, min_len, max_len)
                .into_iter()
                .map(|p| p.extract(&record))
            {
                p.write(&mut stdout)?;
            }
        }
    }
    Ok(())
}
