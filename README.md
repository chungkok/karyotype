# cyto_auto_annotation

**Automated cytogenetic karyotype parsing and risk classification for AML and MDS**

This R script parses ISCN-formatted karyotype strings and applies guideline-based classification rules to produce structured, analysis-ready output. It is designed for use in research and diagnostic workflows involving Acute Myeloid Leukaemia (AML) and Myelodysplastic Syndromes (MDS).

---

## Features

- Parses free-text ISCN karyotype strings into structured abnormality tables
- Resolves multi-clone karyotypes with dependent clone expansion (`idem`, `sdl`, `sl`)
- Filters non-clonal cells (`[n=1]`) per ISCN convention
- Identifies the best clone per guideline rules (most abnormalities; preference for independent clones)
- Classifies each abnormality by type: deletion, trisomy, monosomy, derivative, dicentric, ring, marker, tetraploidy
- Detects complex and monosomal karyotypes
- Applies targeted locus flags: `del(5q)`, `del(7q)`, `-7`, `loss(17p)`, `MLL/KMT2A`
- Assigns disease-defining flags for AML and MDS cytogenetic patterns
- Risk stratification per:
  - **ELN 2022** (AML cytogenetic risk: Favorable / Intermediate / Adverse)
  - **IPSS-R** (MDS cytogenetic risk: Very good / Good / Intermediate / Poor / Very poor)
- Entity classification per **WHO 2022** (cytogenetics-inferrable entities)
- Outputs both a per-sample summary table and a per-abnormality detail table

---

## Dependencies

```r
install.packages(c("stringr", "dplyr", "purrr", "tidyr", "readr"))
```

| Package   | Role                              |
|-----------|-----------------------------------|
| `stringr` | Regex-based ISCN string parsing   |
| `dplyr`   | Table manipulation                |
| `purrr`   | Functional iteration over clones  |
| `tidyr`   | Data reshaping                    |
| `readr`   | CSV I/O                           |

---

## Input format

A CSV file with exactly two required columns:

| Column      | Type      | Description                              |
|-------------|-----------|------------------------------------------|
| `sample_id` | character | Unique sample identifier                 |
| `karyotype` | character | ISCN-formatted karyotype string          |

**Example (`test_1.csv`):**

```csv
sample_id,karyotype
p1,"46,XY[20]"
p2,"45,XY,-7[12]"
p3,"46~49,XX,-4,-5,add(7)(q22),+21,i(21)(?q10),+2~4mar"
```

- Multi-clone karyotypes separated by `/` are supported
- Clone cell counts in `[n]` or `[cpn]` brackets are parsed automatically
- Single-cell clones (`[1]`) are excluded (non-clonal per ISCN)
- Dependent clones using `idem` and `sdl`/`sl` keywords are fully expanded before scoring

---

## Usage

```r
source("cyto_auto_annotation_2026_v1.R")

df  <- read_csv("your_input.csv")   # must contain: sample_id, karyotype
out <- run_pipeline(df)

write_csv(out$summary, "summary.csv")
write_csv(out$detail,  "detail.csv")
```

---

## Output

The pipeline returns a named list with two tables: `$summary` and `$detail`.

### `summary.csv` — one row per sample

| Column                      | Type    | Description                                                          |
|-----------------------------|---------|----------------------------------------------------------------------|
| `sample_id`                 | chr     | Sample identifier (from input)                                       |
| `karyotype`                 | chr     | Original karyotype string (from input)                               |
| `total_abnormalities`       | int     | Aberration count from the best clone (using `n_count` weighting)     |
| `total_score`               | int     | Weighted abnormality score from the best clone                       |
| `clones`                    | int     | Number of abnormal clones detected                                   |
| `n_independent_clones`      | int     | Number of independent / stemline (`sl`) clones                       |
| `n_dependent_clones`        | int     | Number of dependent clones (`idem` / `sdl`)                          |
| `n_normal_clones`           | int     | Number of cytogenetically normal clones                              |
| `best_clone_type`           | chr     | Dependency class of the best clone (`independent`, `sl`, `idem`, `sdl`, `normal`) |
| `clone_dependency_summary`  | chr     | Human-readable clone structure summary                               |
| `complex_karyotype`         | lgl     | `TRUE` if ≥ 3 abnormalities in the best clone                        |
| `monosomal_karyotype`       | lgl     | `TRUE` if ≥ 2 monosomies in the best clone                           |
| `del5`                      | lgl     | Any `del(5...)` detected (any clone; matches `del(5q)`, `del(5p)`, `del(5)(q...)` etc.) |
| `del7q`                     | lgl     | `del(7q)` detected (any clone)                                       |
| `del7`                      | lgl     | `-7` (monosomy 7) detected (any clone)                               |
| `loss17p`                   | lgl     | `del(17p)` or `-17` detected (any clone)                             |
| `MLL`                       | lgl     | KMT2A/MLL rearrangement detected (any clone)                         |
| `MDS_defining_abnormality`  | lgl     | At least one MDS-defining cytogenetic pattern present                |
| `AML_defining_abnormality`  | lgl     | At least one AML-defining cytogenetic pattern present                |
| `MDS_defining_karyo`        | chr     | Matched MDS-defining abnormality strings (comma-separated)           |
| `AML_defining_karyo`        | chr     | Matched AML-defining abnormality strings (comma-separated)           |
| `classification`            | chr     | Karyotype classification (see categories below)                      |
| `eln_subtype`               | chr     | ELN cytogenetic subtype label                                        |
| `eln2022_cyto_risk`         | chr     | ELN 2022 AML cytogenetic risk tier                                   |
| `ipssr_cyto_risk`           | chr     | IPSS-R MDS cytogenetic risk tier                                     |
| `who2022_entity`            | chr     | WHO 2022 cytogenetics-inferrable disease entity                      |

**`classification` values:**

| Value                          | Meaning                                              |
|--------------------------------|------------------------------------------------------|
| `Normal`                       | No cytogenetic abnormalities detected                |
| `AML-defining abnormality`     | Matches a recognised AML cytogenetic pattern         |
| `MDS-defining abnormality`     | Matches a recognised MDS cytogenetic pattern         |
| `Monosomal complex karyotype`  | ≥ 3 abnormalities AND ≥ 2 monosomies                 |
| `Complex karyotype`            | ≥ 3 abnormalities                                    |
| `Monosomal karyotype`          | ≥ 2 monosomies                                       |
| `Single abnormality`           | Exactly 1 abnormality                                |
| `Two abnormalities`            | Exactly 2 abnormalities                              |
| `Other abnormal`               | Abnormal but not fitting above categories            |

---

### `detail.csv` — one row per abnormality per sample

| Column          | Type | Description                                                        |
|-----------------|------|--------------------------------------------------------------------|
| `sample_id`     | chr  | Sample identifier                                                  |
| `karyotype`     | chr  | Original karyotype string                                          |
| `clone`         | int  | Clone index (1-based, ordered as in input string)                  |
| `abnormality`   | chr  | Cleaned ISCN abnormality token                                     |
| `type`          | chr  | Abnormality type (see table below)                                 |
| `chromosomes`   | chr  | Chromosome(s) involved                                             |
| `breakpoint`    | chr  | Breakpoint notation extracted from parentheses                     |
| `weight`        | num  | Score weight (derivatives and markers with ranges are weighted >1) |
| `n_count`       | int  | Aberration unit count (derivatives count as 2 per R5)              |
| `is_best_clone` | lgl  | `TRUE` if this row belongs to the best/scoring clone               |
| `clone_type`    | chr  | Clone dependency class for this clone                              |

**`type` values:**

`deletion`, `trisomy`, `monosomy`, `derivative`, `dicentric`, `isodicentric`, `pseudodicentric`, `pseudoderivative`, `ring`, `marker`, `tetraploidy`, `other`

---

## Key implementation rules

| Rule | Description |
|------|-------------|
| R1   | Clones with cell count `[1]` are excluded as non-clonal |
| R3   | Best clone = clone with the highest `n_count`; ties broken in favour of independent/`sl` clones |
| R5   | Derivative chromosomes count as 2 aberration units (`n_count = 2`) |
| R7   | Chromosome count ≥ 69 triggers a `tetraploidy` flag |
| R8   | Total abnormality count and score are derived from the best clone only |
| R9   | Dependent clones (`idem`, `sdl`) are fully expanded by inheriting parent abnormalities before scoring |

---

## Citation

If you use this script in your research or clinical work, please cite:

> Kok, C.H. (2026). *cyto_auto_annotation: Automated cytogenetic karyotype parsing and risk classification for AML and MDS* [R script]. GitHub. https://github.com/[your-username]/[your-repo]

---

## Author

**Chung H Kok** (2026)

---

## Classification references

- Döhner H, Wei AH, Appelbaum FR, et al. *Diagnosis and management of AML in adults: 2022 recommendations from an international expert panel on behalf of the ELN.* Blood. 2022;140(12):1345–1377. doi:[10.1182/blood.2022016867](https://doi.org/10.1182/blood.2022016867)
- Greenberg PL, Tuechler H, Schanz J, et al. *Revised International Prognostic Scoring System (IPSS-R) for Myelodysplastic Syndromes.* Blood. 2012;120(12):2454–2465. doi:[10.1182/blood-2012-03-420489](https://doi.org/10.1182/blood-2012-03-420489)
- Khoury JD, Solary E, Abla O, et al. *The 5th edition of the World Health Organization Classification of Haematolymphoid Tumours: Myeloid and Histiocytic/Dendritic Neoplasms.* Leukemia. 2022;36(7):1703–1719. doi:[10.1038/s41375-022-01613-1](https://doi.org/10.1038/s41375-022-01613-1)
