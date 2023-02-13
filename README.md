# spASM

Tool for finding allele-specific methylation with the epiBED format.

## Install

spASM is written in Rust. If Rust is not installed, you can find instructions
[here](https://www.rust-lang.org/tools/install).

For now, spASM can be installed by cloning the repository and then built with `cargo` (installed when installing Rust).
```
$ git clone git@github.com:jamorrison/spASM.git
$ cd spASM
$ cargo build --release
```
By default, Rust builds in debug mode, so you will have to include `--release` to build in release mode with the
optimizations that are included therein. The default path the binary will be built at is
`(spasm top directory)/target/release/spasm`.

When development slows, static binaries will be made available on the GitHub releases page.

## Usage

```
spasm [options] <GENOME> <EPIBED>
```

### Required Arguments
|  Name  | Description                                                                   |
|:------:|:------------------------------------------------------------------------------|
| GENOME | Indexed FASTA (`samtools faidx`) file of your genome                          |
| EPIBED | Your epiBED file (if running with `-g`, then it must be bgzipped and tabixed) |

### Optional Arguments
| Short Option | Long Option      | Description | Default |
|:------------:|:----------------:|:------------|:--------|
|:    `-g`    :| `--region`       | region to extract (chr:start-end or chr)                                            | all                            |
|:    `-f`    :| `--fragment`     | collapse reads to fragment-level                                                    | left as individual fragments   |
|:    `-c`    :| `--fdr`          | type of false discovery rate correction to perform possibilities: BH (Benjamini-Hochberg), BY (Benjamini-Yekutieli), Bonferroni, Hochberg, Holm, No (do not apply false discovery correction) | BH |
|:    `-p`    :| `--pcutoff`      | p-value significance cutoff                                                         | 0.05                           |
|:    `-o`    :| `--output`       | output file name, compression level based on file name                              | stdout                         |
|:    `-O`    :| `--candidate`    | write only candidate locations, FDR-corrected p-value based on all locations probed | all locations printed          |
|:    `-N`    :| `--no-ambiguous` | write only locations with no ambiguous SNPs                                         | all SNPs are written           |
|:    `-b`    :| `--biscuit`      | write in BISCUIT ASM output format                                                  | output written in BEDPE format |
|:    `-v`    :| `--verbose`      | verbosity level (0: ERRORS ONLY | 1: WARNINGS + ERRORS | 2+: ALL)                   | 1                              |
|:    `-h`    :| `--help`         | Print help                                                                          |                                |
|:    `-V`    :| `--version`      | Print version                                                                       |                                |
