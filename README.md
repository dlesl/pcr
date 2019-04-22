# About
This is a library and simple CLI for simulating PCR _in silico_,
written in Rust.

The main aim of this project was to learn Rust, use it at your own
risk! :)

At this stage only [gb_io](https://github.com/dlesl/gb-io/)'s `Seq`
format is supported for templates.

## Features
* Support for linear and circular templates
* Extraction of products preserving annotations

## Example
Amplify _lacZ_ from _E. coli_
```sh
$ pcr_cli mg1655.gb ATGACCATGATTACGGATTCAC TTATTTTTGACACCAGACCAAC
```

## Use it online
* [Simple PCR app](https://dlesl.github.io/clonifier/simple_pcr.html)
* [Clonifier](https://dlesl.github.io/clonifier/)
