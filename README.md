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
spasm [options] <GENOME> <PATH>
```

### Required Arguments

|  Name  | Description                                                                   |
|:------:|:------------------------------------------------------------------------------|
| GENOME | Indexed FASTA (`samtools faidx`) file of your genome                          |
| PATH   | Path to epiBED file (if running with `-g`, then it must be bgzip'd + tabix'd) |

### Optional Arguments
| Short Option (`-`) | Long Option (`--`) | Description                                                                         | Default                        |
|:------------------:|:-------------------|:------------------------------------------------------------------------------------|:-------------------------------|
| `g`                | `region`           | region to extract (chr:start-end or chr)                                            | all                            |
| `n`                | `no-mate-merging`  | do not merge mate reads together into a single DNA fragment                         | mate reads are merged together |
| `c`                | `fdr`              | type of false discovery rate correction to perform possibilities: BH (Benjamini-Hochberg), BY (Benjamini-Yekutieli), Bonferroni, Hochberg, Holm, No (do not apply false discovery correction) | BH |
| `p`                | `pcutoff`          | p-value significance cutoff                                                         | 0.05                           |
| `o`                | `output`           | output file name, compression level based on file name                              | stdout                         |
| `O`                | `candidate`        | write only candidate locations, FDR-corrected p-value based on all locations probed | all locations printed          |
| `N`                | `no-ambiguous`     | write only locations with no ambiguous SNPs                                         | all SNPs are written           |
| `b`                | `biscuit`          | write in BISCUIT ASM output format                                                  | output written in BEDPE format |
| `v`                | `verbose`          | verbosity level (0: ERRORS ONLY, 1: WARNINGS + ERRORS, 2+: ALL)                     | 1                              |
| `h`                | `help`             | Print help                                                                          |                                |
| `V`                | `version`          | Print version                                                                       |                                |

## Merging Mates in Paired-End Data

In paired-end sequencing, the two read mates come from the same DNA fragment; therefore, they represent the same
"epi-haplotype." To recover this correlation, mate reads are merged into a single fragment.  On a qualitative level,
spASM will take the unique portions of reads 1 and 2, plus the read 1 portion of any locations that overlap between the
two reads. With respect to those overlapping regions, it should be noted:

  - spASM will use the data as it comes from `biscuit epiread`. By default, `biscuit epiread` will filter overlapping
  bases (including CpGs and SNPs) from read 2. This will then be passed to spASM, which will see these as filtered bases
  and won't include them in any calculations.
  - The only way to double count information from overlapping portions of reads 1 and 2 would be to include the `-d` in
  `biscuit epiread`.

## Output

spASM can produce output in one of two formats. The default method is a BEDPE-compliant format with the following
columns:

  1. SNP chromosome
  2. SNP start (0-based)
  3. SNP end (1-based, non-inclusive)
  4. CpG chromosome
  5. CpG start (0-based)
  6. CpG end (1-based, non-inclusive)
  7. "candidate" (if < p-value cutoff) or "non_candidate" (if >= p-value cutoff)
  8. p-value after applied false discovery correction
  9. SNP<sub>1</sub> (top row in contingency table)
  10. SNP<sub>2</sub> (bottom row in contigency table)
  11. CpG<sub>1</sub> (left column in contingency table)
  12. CpG<sub>2</sub> (right column in contigency table)
  13. SNP<sub>1</sub>-CpG<sub>1</sub> value in contingency table
  14. SNP<sub>1</sub>-CpG<sub>2</sub> value in contingency table
  15. SNP<sub>2</sub>-CpG<sub>1</sub> value in contingency table
  16. SNP<sub>2</sub>-CpG<sub>2</sub> value in contingency table

spASM can also produce output in the format generated by `biscuit asm`:

  1. Chromosome
  2. SNP position (0-based)
  3. CpG position (0-based)
  4. SNP<sub>1</sub> / SNP<sub>2</sub>
  5. CpG<sub>1</sub> / CpG<sub>2</sub>
  6. SNP<sub>1</sub>-CpG<sub>1</sub> value in contingency table
  7. SNP<sub>1</sub>-CpG<sub>2</sub> value in contingency table
  8. SNP<sub>2</sub>-CpG<sub>1</sub> value in contingency table
  9. SNP<sub>2</sub>-CpG<sub>2</sub> value in contingency table
  10. p-value after applied false discovery correction
  11. `.` (unused column included for consistency with biscuit format)
