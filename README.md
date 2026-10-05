# ldscore — LD score calculator

`ldscore` is a C++ program that computes LD scores from binary PLINK genotypes. For each retained SNP, the score is 1 plus the sum of unbiased *r*² estimates with earlier and later SNPs inside a genetic or physical window.

**Author:** Loic Yengo  
**License:** TBD

---

## Requirements

- C++23 compiler (`g++` recommended; `Makefile` uses `-std=c++23`)
- POSIX threads (`pthread`) and libm (`-lm`)

Tested on macOS with `make`. Genotypes are stored as one byte per call (dosage 0, 1, or 2, or missing). One autosome, or one `--mbfile` chunk, is loaded at a time.

---

## Build

From this directory:

```bash
make
```

This produces the `ldscore` executable here.

To rebuild from scratch:

```bash
make clean
make
```

Run `./ldscore --help` for the full flag list.

---

## Input formats

### Genotypes (`--bfile` or `--mbfile`)

PLINK SNP-major binary files: `[prefix].bed`, `[prefix].bim`, and `[prefix].fam`.

- `--bfile PREFIX` reads one fileset.
- `--mbfile FILE` reads one prefix per line. Blank lines and lines starting with `#` are ignored. Use either `--bfile` or `--mbfile`.

Every FAM must list the same individuals in the same order. A chromosome split across prefixes is scored separately in each prefix; the window does not cross files. A duplicate SNP id on the same chromosome across prefixes is an error.

Autosomes 1–22 are scored, in order. X, Y, MT, and any other chromosome are dropped. BIM column 3 (cM) is ignored. Base-pair positions on a chromosome must be non-decreasing.

Dosage is the count of allele A1.

### Sex-averaged map (`--map-average`)

Required unless `--mb-window` is set. Five columns, in order:

```text
chromosome  begin_bp  end_bp  cM_per_Mb  cumulative_cM
```

A header is optional and is not matched by name. Boundary points are (first Begin, 0 cM), then each (End, cumulative cM). Autosomes 1–22 only. Gaps, overlaps, a decreasing cM column, and a cumulative cM above 5000 are errors.

Without `--interpolate`, a SNP inside `[Begin, End)` takes the cM at Begin, and a SNP outside the span is dropped. With `--interpolate`, interior SNPs are placed by linear interpolation and SNPs outside the span are clamped to the end values.

### Individuals (`--keep`)

Two columns, FID and IID, separated by space or tab. Blank lines and lines starting with `#` are ignored. A header is optional. The first line is a header when it is `FID` and `IID`, or when that pair is not in the FAM and another individual follows.

This list is the eligible set for allele frequency and missingness. It is also the pool for `--nsample`. Without `--nsample`, LD scores use the same individuals. The default is every FAM row.

### LD-score sample (`--nsample`, `--write-sample`)

`--nsample N` draws N eligible individuals without replacement. The same people are used on every chromosome and every `--mbfile` prefix. If the pool has N people or fewer, everyone in the pool is used and `--seed` is not applied. Allele frequency and missingness still use the full eligible set. `--seed` is an integer from 0 through 4294967295; the default is 1, and it is used only when a subset is actually drawn.

`--write-sample` writes those LD-score individuals to `[prefix].sample` (`FID`, `IID`, with a header), in FAM order. That file can later be passed to `--keep`. A later `--keep` of the same people also applies them to allele frequency and missingness, which `--nsample` does not.

### SNP list (`--extract` or `--mextract`)

`--extract FILE` is one SNP identifier per line, matched to BIM column 2. The whole trimmed line is the identifier. Blank lines and lines starting with `#` are ignored. A repeated identifier is kept once.

`--mextract FILE` is one path per line. Each path is an `--extract` file, and the identifiers are pooled. A repeated path is an error. Use either `--extract` or `--mextract`.

The cut is applied after non-autosomal chromosomes are dropped, and before missingness, `--maf`, and the map or base-pair window. LD scores are computed only among SNPs that remain. A listed identifier absent from the BIM is reported in the log. A chromosome that loses every SNP is skipped. If no autosomal SNP remains, the run stops.

---

## Quick start

Genetic window of 1.5 cM, keeping SNPs with minor allele frequency at least 0.01:

```bash
./ldscore \
  --bfile myPLINK \
  --map-average myMap \
  --cm-window 1.5 \
  --maf 0.01 \
  --out PREFIX
```

Physical window of 1 Mb. A map is not read, and `--cm-window` is not used:

```bash
./ldscore \
  --bfile myPLINK \
  --mb-window 1 \
  --maf 0.01 \
  --out PREFIX
```

Several PLINK prefixes, a random LD sample, and a SNP list:

```bash
./ldscore \
  --mbfile prefixes.txt \
  --map-average myMap \
  --keep ids.txt \
  --nsample 1000 \
  --seed 1 \
  --write-sample \
  --extract snps.txt \
  --nthread 8 \
  --out PREFIX
```

---

## Options

| Flag | Default | Role |
|------|---------|------|
| `--bfile` | | One PLINK prefix |
| `--mbfile` | | File of PLINK prefixes. Not with `--bfile` |
| `--map-average` | required without `--mb-window` | Sex-averaged map. Ignored with `--mb-window` |
| `--cm-window` | 1 | Window in centiMorgans. Pairs at or beyond this distance are excluded. Ignored with `--mb-window` |
| `--mb-window` | off | Window in megabases. `--mb-window 1` keeps pairs less than 1,000,000 bp apart. The map, `--cm-window`, and `--interpolate` are not used |
| `--interpolate` | off | Linear map placement. Ignored with `--mb-window` |
| `--maf` | 0.01 | Minimum minor allele frequency, as a fraction in [0, 0.5]. SNPs below this value are dropped |
| `--max-missing` | 0.05 | Maximum missingness, as a fraction in [0, 1]. A SNP above this value is dropped. A SNP exactly at the threshold is kept |
| `--extract` / `--mextract` | off | SNP identifier list, or a file of such lists |
| `--keep` | all FAM rows | Eligible individuals for frequency, missingness, and the `--nsample` pool |
| `--nsample` | off | Number of eligible individuals drawn for the *r*² calculation |
| `--seed` | 1 | Seed for `--nsample`, from 0 through 4294967295 |
| `--write-sample` | off | Write the LD-score individuals to `[prefix].sample` |
| `--nthread` | 1 | POSIX threads for the SNP loop within each chromosome |
| `--out` | `ldscore` | Prefix for `[prefix].l2.ldscore` and `[prefix].log` |

`--maf` uses *p* = (sum of dosages / number of non-missing genotypes) / 2, and the minor allele frequency is min(*p*, 1−*p*). Both are computed on the eligible individuals, not on the `--nsample` draw.

---

## Output

| File | Contents |
|------|----------|
| `[prefix].l2.ldscore` | Header `CHR SNP BP L2`, then one row per retained SNP |
| `[prefix].log` | Command, filters, window, sample, and per-chromosome summaries |
| `[prefix].sample` | Written with `--write-sample`. Header `FID IID`, then the LD-score individuals in FAM order |

The chromosome label is the BIM chromosome string. Scores below 1e-8 are written as 1e-8.

---

## LD score

Each retained SNP starts at 1. A pair contributes the same unbiased *r*² to both SNPs when the distance is strictly below the window. A pair exactly at the window boundary is left out. Missing genotypes are skipped in that pair. With fewer than 3 observed individuals, or with no variance at either SNP, that pair contributes 0. Negative estimates are kept.

The unbiased estimate is

*r*² − (1 − *r*²) / (*n* − 2),

where *r*² is the squared Pearson correlation of the dosages and *n* is the number of individuals observed at both SNPs.

`--nthread 1` follows the single-thread sum, so those scores stay the same from run to run. The same thread count repeats the same scores. A different thread count can change the last bits, because the pairwise terms are added in a different order. Threads share the chromosome genotype matrix.

With fewer than 3 individuals in the LD sample, every pairwise term is 0 and every score stays at 1. The log says so.

---

## Citation

If you use `ldscore` in published work, please cite this software. (Formal citation details TBD.)

---

## Acknowledgments

Parts of this codebase were developed with assistance from [Cursor](https://cursor.com).

---

## Contact

Loic Yengo — for questions about the methods or software.
