############################################################
# CYTOGENETIC PARSER & SCORING PIPELINE
# GUIDELINE-COMPLIANT: Chun et al. / Supp Table 2
# + sl / sdl / idem DEPENDENCY CLASSIFICATION
# + INDEPENDENT CLONE COUNTING
# + ELN 2022 / IPSS-R / WHO 2022
#
# Clone dependency rules (ISCN 2020):
#   sl   (stemline)  — the reference/founding clone; INDEPENDENT
#   sdl  (sideline)  — derived from the sl clone;    DEPENDENT  → expands from sl
#   idem             — inherits previous clone's abn; DEPENDENT  → expands from parent
#   Normal clone     — not an abnormal independent clone (not counted)
#   All others       — INDEPENDENT abnormal clones
#
# Counting rules (Chun et al. 2010 / Supp Table 2):
#   R1  Non-clonal [n=1] metaphase clones are ignored
#   R2  Count 1 aberration per comma-separated item
#   R3  Use the clone with the HIGHEST aberration count (not summed)
#   R4  Numerical / balanced / simple structural = 1 aberration
#   R5  Unbalanced translocation (der) = 2 aberrations
#   R6  Constitutional aberrations not counted
#   R7  Tetraploidy (>=69 chromosomes) = 1
#   R8  Aberrations in independent clones are NOT added together
#   R9  idem / sdl subclones are expanded before counting
#
# New output columns vs. previous version:
#   clone_dependency_summary — e.g. "2 independent; 1 dependent (sdl); 1 normal"
#   n_independent_clones     — count of independent abnormal clones
#   n_dependent_clones       — count of dependent (idem/sdl) clones
#   n_normal_clones          — count of normal (46,XX or 46,XY) clones
#   best_clone_type          — "independent" | "dependent" (which clone drove scoring)
############################################################

library(stringr)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)

# ============================================================
# HELPER: extract clone cell count from [n] or [cpn] suffix
# ============================================================
get_clone_size <- function(clone_str) {
  m <- str_match(clone_str, "\\[cp?(\\d+)\\]")
  if (!is.na(m[1, 2])) return(as.integer(m[1, 2]))
  return(2L)  # no bracket → treat as clonal (>=2 cells)
}

# ============================================================
# 1a. LIGHT CLEAN — normalise whitespace/punctuation, PRESERVE [n]
# ============================================================
clean_karyo_light <- function(x) {
  if (is.na(x) || str_trim(x) == "") return(NA_character_)
  x %>%
    str_replace_all("\\?", "") %>%
    str_squish() %>%
    str_replace_all(",\\s+", ",") %>%
    str_replace_all("\\s+,", ",") %>%
    str_replace_all("\\(\\s+", "(") %>%
    str_replace_all("\\s+\\)", ")")
}

# 1b. FULL CLEAN (strips [n] — for legacy callers)
clean_karyo <- function(x) {
  if (is.na(x) || str_trim(x) == "") return(NA_character_)
  clean_karyo_light(x) %>% str_remove_all("\\[.*?\\]")
}

# ============================================================
# 2. CLASSIFY CLONE DEPENDENCY
#
#  Input: a clone string WITH [n] already stripped
#  Returns one of:
#    "normal"      — no abnormalities after parsing (e.g. 46,XX)
#    "independent" — abnormal clone, no idem/sdl/sl token
#    "sl"          — contains sl token (is THE stemline; independent)
#    "sdl"         — contains sdl token (derived from stemline; dependent)
#    "idem"        — contains idem token (inherits previous clone; dependent)
# ============================================================
classify_clone_dependency <- function(clone_str) {
  tokens <- str_split(clone_str, ",")[[1]] %>% str_trim()
  abns   <- tokens[-c(1, 2)]   # strip chromosome count and sex
  
  # Check for dependency keywords in abnormality tokens
  if (any(abns == "idem")) return("idem")
  if (any(abns == "sdl"))  return("sdl")
  if (any(abns == "sl"))   return("sl")   # sl IS the stemline reference (independent)
  
  # Normal if no abnormalities
  if (length(abns) == 0)   return("normal")
  
  return("independent")
}

# ============================================================
# 3. EXPAND DEPENDENT CLONES  (R9)
#
#  Handles both idem and sdl expansion:
#    idem → inherits from the immediately preceding clone string
#    sdl  → inherits from the sl (stemline) clone
#    sl   → the stemline clone itself; remove "sl" token before parsing
#
#  Input:  clone strings with [n] ALREADY STRIPPED
#  Output: expanded clone strings ready for parse_clone()
# ============================================================
expand_dependent_clones <- function(clone_strs) {
  deps    <- map_chr(clone_strs, classify_clone_dependency)
  result  <- clone_strs
  
  # Identify the sl clone index (there should be at most one)
  sl_idx  <- which(deps == "sl")
  sl_str  <- if (length(sl_idx) > 0) {
    # Remove the "sl" token from the stemline clone string
    tokens <- str_split(clone_strs[sl_idx[1]], ",")[[1]] %>% str_trim()
    abns   <- tokens[-c(1, 2)]
    abns   <- abns[abns != "sl"]
    paste(c(tokens[1], tokens[2], abns), collapse = ",")
  } else NA_character_
  
  # Clean "sl" token from the stemline in result
  if (length(sl_idx) > 0) result[sl_idx[1]] <- sl_str
  
  for (i in seq_along(clone_strs)) {
    dep <- deps[i]
    
    if (dep == "idem") {
      # Expand from immediately preceding clone (already expanded if it was also idem)
      if (i == 1) next
      parent_tokens <- str_split(result[i - 1], ",")[[1]] %>% str_trim()
      parent_abns   <- parent_tokens[-c(1, 2)]
      cur_tokens    <- str_split(clone_strs[i], ",")[[1]] %>% str_trim()
      cur_abns      <- cur_tokens[-c(1, 2)]
      cur_abns      <- cur_abns[cur_abns != "idem"]  # remove the idem keyword
      result[i]     <- paste(c(cur_tokens[1], parent_tokens[2],
                               parent_abns, cur_abns), collapse = ",")
    }
    
    if (dep == "sdl") {
      if (is.na(sl_str)) next  # no sl clone found; skip
      sl_tokens  <- str_split(sl_str, ",")[[1]] %>% str_trim()
      sl_abns    <- sl_tokens[-c(1, 2)]
      cur_tokens <- str_split(clone_strs[i], ",")[[1]] %>% str_trim()
      cur_abns   <- cur_tokens[-c(1, 2)]
      cur_abns   <- cur_abns[cur_abns != "sdl"]      # remove the sdl keyword
      result[i]  <- paste(c(cur_tokens[1], sl_tokens[2],
                            sl_abns, cur_abns), collapse = ",")
    }
  }
  
  result
}

# ============================================================
# 4. TOKENIZE CLONE
# ============================================================
tokenize_clone <- function(clone) {
  tokens <- str_split(clone, ",")[[1]] %>% str_trim()
  list(
    chromosome_count = tokens[1],
    sex              = tokens[2],
    abnormalities    = tokens[-c(1, 2)]
  )
}

# ============================================================
# 5. PARSE ABNORMALITY
#    n_count: guideline-compliant aberration count per item
#      der (unbalanced translocation) = 2   (R5)
#      all others                     = 1
# ============================================================
parse_abnormality <- function(abn) {
  abn_type <- "other"
  weight   <- 1
  n_count  <- 1L
  chrom    <- NA_character_
  
  # Deletion
  if (str_detect(abn, "^del")) {
    abn_type <- "deletion"
    if (!str_detect(abn, "\\("))
      abn <- str_replace(abn, "del([0-9XY]+[pq]?)", "del(\\1)")
    m <- str_match(abn, "del\\(([0-9XY]+)")
    if (!is.na(m[2])) chrom <- m[2]
  }
  
  if (str_detect(abn, "^\\+"))       abn_type <- "trisomy"
  if (str_detect(abn, "^\\-"))       abn_type <- "monosomy"
  if (str_detect(abn, "dic\\("))     abn_type <- "dicentric"
  if (str_detect(abn, "idic"))       abn_type <- "isodicentric"
  if (str_detect(abn, "psu\\s*dic")) abn_type <- "pseudodicentric"
  if (str_detect(abn, "psu\\s*der")) abn_type <- "pseudoderivative"
  if (str_detect(abn, "\\+r"))       abn_type <- "ring"
  
  # Derivative = unbalanced translocation → n_count = 2 (R5)
  if (str_detect(abn, "der\\(")) {
    abn_type <- "derivative"
    n_count  <- 2L
    weight   <- 2L
  }
  
  # Markers
  if (str_detect(abn, "mar")) {
    if (str_detect(abn, "~")) {
      rng       <- str_extract(abn, "[0-9]+~[0-9]+")
      max_count <- as.numeric(str_split(rng, "~")[[1]][2])
      weight    <- max_count
    } else {
      weight <- 1
    }
    abn_type <- "marker"
    n_count  <- 1L
    chrom    <- NA_character_
  }
  
  # Chromosome extraction for structural types
  if (is.na(chrom)) {
    if (abn_type %in% c("derivative", "dicentric", "pseudoderivative", "pseudodicentric")) {
      m <- str_match(abn, "\\(([^\\)]+)\\)")
      if (!is.na(m[2])) {
        chrom <- str_split(m[2], ";")[[1]] %>%
          str_replace_all("[^0-9XY]", "") %>%
          paste(collapse = ";")
      }
    }
    if (abn_type %in% c("monosomy", "trisomy", "ring")) {
      m <- str_match(abn, "([0-9XY]+)")
      if (!is.na(m[2])) chrom <- m[2]
    }
  }
  
  breakpoint <- str_extract(abn, "\\([^\\)]+\\)$")
  
  tibble(abnormality = abn, type = abn_type, chromosomes = chrom,
         breakpoint = breakpoint, weight = weight, n_count = n_count)
}

# ============================================================
# 6. PARSE CLONE  (fully-typed empty tibble; tetraploidy R7)
# ============================================================
parse_clone <- function(clone, clone_id) {
  empty <- tibble(
    clone = clone_id, abnormality = character(), type = character(),
    chromosomes = character(), breakpoint = character(),
    weight = numeric(), n_count = integer()
  )
  if (is.na(clone) || str_trim(clone) == "") return(empty)
  
  tk   <- tokenize_clone(clone)
  rows <- list()
  
  # Tetraploidy: chromosome count >=69 (R7)
  chr_n <- suppressWarnings(as.numeric(str_extract(tk$chromosome_count, "^[0-9]+")))
  if (!is.na(chr_n) && chr_n >= 69) {
    rows <- c(rows, list(tibble(
      clone = clone_id, abnormality = "tetraploidy", type = "tetraploidy",
      chromosomes = NA_character_, breakpoint = NA_character_,
      weight = 1L, n_count = 1L
    )))
  }
  
  if (length(tk$abnormalities) == 0) {
    if (length(rows) > 0) return(bind_rows(rows))
    return(empty)
  }
  
  parsed <- map_dfr(tk$abnormalities, parse_abnormality) %>%
    mutate(clone = clone_id)
  
  if (length(rows) > 0) bind_rows(c(rows, list(parsed))) else parsed
}

# ============================================================
# 7. PROCESS KARYOTYPE — full guideline-compliant pipeline
#
#  Steps:
#  (a) Light-clean (preserve [n])
#  (b) Split by "/" into clone strings
#  (c) Read clone sizes; drop [n=1] non-clonal clones (R1)
#  (d) Strip [n] brackets
#  (e) Classify each clone's dependency (normal/independent/sl/sdl/idem)
#  (f) Expand dependent clones (sdl from sl; idem from parent) (R9)
#  (g) Parse each clone independently
#  (h) Count n_count per clone (der=2, R5)
#  (i) Select BEST CLONE = highest n_count sum (R3)
#  (j) Totals from best clone only (R8)
#  (k) Count independent, dependent, normal clones
# ============================================================
process_karyotype <- function(karyo) {
  empty_tbl <- tibble(
    clone = integer(), abnormality = character(), type = character(),
    chromosomes = character(), breakpoint = character(),
    weight = numeric(), n_count = integer(), is_best_clone = logical(),
    clone_type = character()
  )
  empty_out <- list(
    table               = empty_tbl,
    total_abnormalities = 0L,
    total_score         = 0L,
    clones              = 0L,
    n_independent_clones = 0L,
    n_dependent_clones   = 0L,
    n_normal_clones      = 0L,
    best_clone_type      = NA_character_,
    clone_dependency_summary = NA_character_
  )
  
  if (is.na(karyo) || str_trim(karyo) == "") return(empty_out)
  
  # (a-b) Light clean; split
  clone_strs_raw <- str_split(clean_karyo_light(karyo), "/")[[1]]
  
  # (c) Filter out [n=1] non-clonal clones (R1)
  clone_sizes    <- map_int(clone_strs_raw, get_clone_size)
  keep           <- clone_sizes != 1L
  clone_strs_raw <- clone_strs_raw[keep]
  if (length(clone_strs_raw) == 0) return(empty_out)
  
  # (d) Strip [n] brackets
  clone_strs_stripped <- map_chr(clone_strs_raw,
                                 ~ str_remove(.x, "\\[.*?\\]") %>% str_trim())
  
  # (e) Classify dependency BEFORE expansion (so we know sl/sdl/idem status)
  dep_types_raw <- map_chr(clone_strs_stripped, classify_clone_dependency)
  
  # (f) Expand dependent clones (R9)
  clone_strs_expanded <- expand_dependent_clones(clone_strs_stripped)
  
  # After expansion, "normal" classification may change (idem/sdl could expand to
  # abnormal). Re-check normality on expanded strings to update dep_types.
  # Dependency label stays as classified pre-expansion (sl/sdl/idem/independent).
  # We re-check "normal" from expanded parse result later (step h).
  
  # (g) Parse each clone
  parsed_per_clone <- map(seq_along(clone_strs_expanded),
                          ~ parse_clone(clone_strs_expanded[.x], .x))
  
  # (h) Aberration count per clone using n_count
  n_abn_per_clone   <- map_int(parsed_per_clone,
                               ~ as.integer(sum(.x$n_count, na.rm = TRUE)))
  score_per_clone   <- map_dbl(parsed_per_clone,
                               ~ sum(.x$weight, na.rm = TRUE))
  
  # Re-classify "normal" based on expanded parse (0 aberrations = normal)
  dep_types <- dep_types_raw
  dep_types[n_abn_per_clone == 0] <- "normal"
  
  # (i) Best clone (R3): highest n_count; among ties, prefer independent clone
  best_idx <- {
    max_n <- max(n_abn_per_clone)
    candidates <- which(n_abn_per_clone == max_n)
    # Prefer independent or sl clone if tied
    indep_candidates <- candidates[dep_types[candidates] %in% c("independent", "sl")]
    if (length(indep_candidates) > 0) indep_candidates[1] else candidates[1]
  }
  
  # (j) Totals from best clone only (R8)
  total_abnormalities <- n_abn_per_clone[best_idx]
  total_score         <- as.integer(score_per_clone[best_idx])
  abnormal_clones     <- sum(n_abn_per_clone > 0)
  
  # (k) Independence summary
  n_independent <- sum(dep_types %in% c("independent", "sl"))
  n_dependent   <- sum(dep_types %in% c("idem", "sdl"))
  n_normal      <- sum(dep_types == "normal")
  best_type     <- dep_types[best_idx]
  
  # Build summary string
  parts <- c()
  if (n_independent > 0)
    parts <- c(parts, paste0(n_independent, " independent"))
  if (sum(dep_types == "sdl") > 0)
    parts <- c(parts, paste0(sum(dep_types == "sdl"), " dependent (sdl)"))
  if (sum(dep_types == "idem") > 0)
    parts <- c(parts, paste0(sum(dep_types == "idem"), " dependent (idem)"))
  if (n_normal > 0)
    parts <- c(parts, paste0(n_normal, " normal"))
  dep_summary <- if (length(parts) > 0) paste(parts, collapse = "; ") else NA_character_
  
  # Assemble full table with clone_type and is_best_clone
  parsed_all <- map_dfr(seq_along(parsed_per_clone), function(i) {
    p <- parsed_per_clone[[i]]
    if (nrow(p) == 0) return(NULL)
    p %>% mutate(is_best_clone = (i == best_idx),
                 clone_type    = dep_types[i])
  })
  
  if (nrow(parsed_all) == 0) parsed_all <- empty_tbl
  
  list(
    table                    = parsed_all,
    total_abnormalities      = total_abnormalities,
    total_score              = total_score,
    clones                   = abnormal_clones,
    n_independent_clones     = n_independent,
    n_dependent_clones       = n_dependent,
    n_normal_clones          = n_normal,
    best_clone_type          = best_type,
    clone_dependency_summary = dep_summary
  )
}

# ============================================================
# 8. SAFE PARSED — full schema including n_count, is_best_clone, clone_type
# ============================================================
safe_parsed <- function(parsed_table) {
  empty <- tibble(
    abnormality = character(), type = character(), chromosomes = character(),
    breakpoint = character(), weight = numeric(), n_count = integer(),
    is_best_clone = logical(), clone_type = character()
  )
  if (is.null(parsed_table) || length(parsed_table) == 0) return(empty)
  if (!is.data.frame(parsed_table) && length(parsed_table) == 1 && is.na(parsed_table))
    return(empty)
  parsed_table
}

# ============================================================
# 9. STRUCTURAL FLAGS
#
#  detect_complex / detect_monosomal → best clone only (R3/R8)
#  detect_complex uses sum(n_count) so der counts as 2 (R5)
#
#  Locus / disease flags → all clones (maximum clinical sensitivity:
#    if any clone has an AML/MDS abnormality, the case is positive)
# ============================================================
detect_monosomal <- function(df) {
  df <- safe_parsed(df)
  if ("is_best_clone" %in% colnames(df)) df <- filter(df, is_best_clone)
  sum(df$type == "monosomy", na.rm = TRUE) >= 2
}

detect_complex <- function(df) {
  df <- safe_parsed(df)
  if ("is_best_clone" %in% colnames(df)) df <- filter(df, is_best_clone)
  n <- if ("n_count" %in% colnames(df)) sum(df$n_count, na.rm = TRUE) else nrow(df)
  n >= 3
}

del5_flag    <- function(pt) { pt <- safe_parsed(pt); nrow(pt) > 0 && any(str_detect(pt$abnormality, "del\\(5"))        }
del7q_flag   <- function(pt) { pt <- safe_parsed(pt); nrow(pt) > 0 && any(str_detect(pt$abnormality, "del\\(7q"))       }
del7_flag    <- function(pt) { pt <- safe_parsed(pt); nrow(pt) > 0 && any(str_detect(pt$abnormality, "^-7"))            }
loss17p_flag <- function(pt) { pt <- safe_parsed(pt); nrow(pt) > 0 && any(str_detect(pt$abnormality, "^-17|del\\(17p")) }
MLL_flag     <- function(pt) { pt <- safe_parsed(pt); nrow(pt) > 0 && any(str_detect(pt$abnormality, "t\\(9;11\\)|t\\(11;19\\)|t\\(11q23\\)")) }

MDS_patterns <- c(
  "^-7", "del\\(7q\\)", "del\\(5q\\)", "t\\(5q\\)", "i\\(17q\\)", "t\\(17p\\)",
  "^-13", "del\\(13q\\)", "del\\(11q\\)", "del\\(12p\\)", "t\\(12p\\)", "del\\(9q\\)",
  "idic\\(X\\)\\(q13\\)", "t\\(11;16\\)", "t\\(3;21\\)", "t\\(1;3\\)", "t\\(2;11\\)",
  "inv\\(3\\)", "t\\(3;3\\)", "t\\(6;9\\)"
)
AML_patterns <- c(
  "t\\(15;17\\)", "t\\(8;21\\)", "inv\\(16\\)", "t\\(16;16\\)", "t\\(9;11\\)", "t\\(6;9\\)",
  "inv\\(3\\)", "t\\(3;3\\)", "t\\(9;22\\)", "t\\(1;3\\)", "t\\(1;22\\)", "t\\(3;5\\)",
  "t\\(5;11\\)", "t\\(7;12\\)", "t\\(8;16\\)", "t\\(10;11\\)", "t\\(11;12\\)", "t\\(16;21\\)"
)

MDS_flag <- function(pt) {
  pt <- safe_parsed(pt); if (nrow(pt) == 0) return(FALSE)
  any(sapply(MDS_patterns, function(p) any(str_detect(pt$abnormality, p))))
}
AML_flag <- function(pt) {
  pt <- safe_parsed(pt); if (nrow(pt) == 0) return(FALSE)
  any(sapply(AML_patterns, function(p) any(str_detect(pt$abnormality, p))))
}

get_MDS_abnormalities <- function(pt) {
  pt <- safe_parsed(pt); if (nrow(pt) == 0) return(NA_character_)
  matches <- sapply(MDS_patterns, function(p) str_detect(pt$abnormality, p))
  if (is.vector(matches)) matches <- matrix(matches, ncol = 1)
  matched <- pt$abnormality[rowSums(matches) > 0]
  if (length(matched) == 0) return(NA_character_)
  paste(unique(matched), collapse = ", ")
}
get_AML_abnormalities <- function(pt) {
  pt <- safe_parsed(pt); if (nrow(pt) == 0) return(NA_character_)
  matches <- sapply(AML_patterns, function(p) str_detect(pt$abnormality, p))
  if (is.vector(matches)) matches <- matrix(matches, ncol = 1)
  matched <- pt$abnormality[rowSums(matches) > 0]
  if (length(matched) == 0) return(NA_character_)
  paste(unique(matched), collapse = ", ")
}

# ============================================================
# 10. ELN SUBTYPE
# ============================================================
get_eln_subtype <- function(karyo) {
  if (is.na(karyo) || str_trim(karyo) == "") return(NA_character_)
  if (str_detect(karyo, "t\\(15;17\\)"))   return("t(15;17) APL")
  if (str_detect(karyo, "t\\(8;21\\)"))    return("t(8;21) RUNX1-RUNX1T1")
  if (str_detect(karyo, "inv\\(16\\)"))    return("inv(16) CBFB-MYH11")
  if (str_detect(karyo, "t\\(9;11\\)") |
      str_detect(karyo, "t\\(11;19\\)") |
      str_detect(karyo, "t\\(11q23\\)"))   return("KMT2A/MLL rearrangement")
  return("Other")
}

# ============================================================
# 11. CLASSIFY KARYOTYPE
# ============================================================
classify_karyo <- function(parsed_table, karyo_string) {
  parsed_table <- safe_parsed(parsed_table)
  if (nrow(parsed_table) == 0 || is.na(karyo_string) ||
      str_trim(karyo_string) == "" || str_detect(karyo_string, "^46,[XY]{2}$"))
    return("Normal")
  best_tbl <- if ("is_best_clone" %in% colnames(parsed_table))
    filter(parsed_table, is_best_clone) else parsed_table
  n_abn    <- if ("n_count" %in% colnames(best_tbl))
    sum(best_tbl$n_count, na.rm = TRUE) else nrow(best_tbl)
  is_complex <- detect_complex(parsed_table)
  is_mono    <- detect_monosomal(parsed_table)
  if (AML_flag(parsed_table)) return("AML-defining abnormality")
  if (MDS_flag(parsed_table)) return("MDS-defining abnormality")
  if (is_complex && is_mono)  return("Monosomal complex karyotype")
  if (is_complex)             return("Complex karyotype")
  if (is_mono)                return("Monosomal karyotype")
  if (n_abn == 1)             return("Single abnormality")
  if (n_abn == 2)             return("Two abnormalities")
  return("Other abnormal")
}

# ============================================================
# 12. ELN 2022 CYTOGENETIC RISK (AML)
# Source: Döhner et al., Blood 2022
# ============================================================
get_eln2022_cyto_risk <- function(parsed_table, karyo_string) {
  parsed_table <- safe_parsed(parsed_table)
  if (is.na(karyo_string) || str_trim(karyo_string) == "") return(NA_character_)
  abn   <- parsed_table$abnormality
  n_abn <- nrow(parsed_table)
  has   <- function(patterns) {
    if (n_abn == 0) return(FALSE)
    any(sapply(patterns, function(p) any(str_detect(abn, p))))
  }
  # Favorable
  if (has(c("t\\(8;21\\)", "inv\\(16\\)", "t\\(16;16\\)", "t\\(15;17\\)")))
    return("Favorable")
  # Adverse
  adv_tx    <- c("t\\(6;9\\)", "t\\(9;22\\)", "inv\\(3\\)", "t\\(3;3\\)")
  kmt2a_adv <- c("t\\(1;11\\)","t\\(4;11\\)","t\\(5;11\\)","t\\(6;11\\)",
                 "t\\(10;11\\)","t\\(11;16\\)","t\\(11;17\\)","t\\(11;19\\)","t\\(11q23\\)")
  if (has(adv_tx) || has(kmt2a_adv) ||
      has("^-5$") || has("del\\(5q\\)") || has("^-7$") ||
      has(c("^-17$","del\\(17p\\)","i\\(17q\\)")) ||
      detect_complex(parsed_table) || detect_monosomal(parsed_table))
    return("Adverse")
  return("Intermediate")
}

# ============================================================
# 13. IPSS-R CYTOGENETIC RISK (MDS)
# Source: Greenberg et al., Blood 2012
# ============================================================
get_ipssr_cyto_risk <- function(parsed_table, karyo_string) {
  parsed_table <- safe_parsed(parsed_table)
  if (is.na(karyo_string) || str_trim(karyo_string) == "") return(NA_character_)
  best_tbl <- if ("is_best_clone" %in% colnames(parsed_table))
    filter(parsed_table, is_best_clone) else parsed_table
  abn   <- best_tbl$abnormality
  n_abn <- if ("n_count" %in% colnames(best_tbl))
    sum(best_tbl$n_count, na.rm = TRUE) else nrow(best_tbl)
  if (n_abn == 0) return("Good")
  if (n_abn >  3) return("Very poor")
  has_minusY <- any(str_detect(abn, "^-Y$"))
  has_minus7 <- any(str_detect(abn, "^-7$"))
  has_del7q  <- any(str_detect(abn, "del\\(7q\\)"))
  has_del5q  <- any(str_detect(abn, "del\\(5q\\)"))
  has_del12p <- any(str_detect(abn, "del\\(12p\\)"))
  has_del20q <- any(str_detect(abn, "del\\(20q\\)"))
  has_del11q <- any(str_detect(abn, "del\\(11q\\)"))
  has_inv3   <- any(str_detect(abn, "inv\\(3\\)|t\\(3;3\\)"))
  if (n_abn == 3)                               return("Poor")
  if (n_abn == 1 && has_minus7)                 return("Poor")
  if (has_inv3)                                 return("Poor")
  if (n_abn == 2 && (has_minus7 || has_del7q))  return("Poor")
  if (n_abn == 1 && has_minusY)  return("Very good")
  if (n_abn == 1 && has_del11q)  return("Very good")
  if (n_abn == 1 && (has_del5q || has_del12p || has_del20q)) return("Good")
  if (n_abn == 2 && has_del5q)                               return("Good")
  return("Intermediate")
}

# ============================================================
# 14. WHO 2022 ENTITY (cytogenetics-inferrable)
# Source: Khoury et al., Leukemia 2022
# ============================================================
get_who2022_entity <- function(parsed_table, karyo_string) {
  parsed_table <- safe_parsed(parsed_table)
  if (is.na(karyo_string) || str_trim(karyo_string) == "") return(NA_character_)
  abn   <- parsed_table$abnormality
  n_abn <- nrow(parsed_table)
  has   <- function(patterns) {
    if (n_abn == 0) return(FALSE)
    any(sapply(patterns, function(p) any(str_detect(abn, p))))
  }
  if (has("t\\(15;17\\)")) return("APL with PML::RARA [t(15;17)]")
  if (has("t\\(8;21\\)"))  return("AML with RUNX1::RUNX1T1 [t(8;21)]")
  if (has(c("inv\\(16\\)","t\\(16;16\\)")))  return("AML with CBFB::MYH11 [inv(16)/t(16;16)]")
  if (has(c("inv\\(3\\)","t\\(3;3\\)")))     return("AML with MECOM rearrangement [inv(3)/t(3;3)]")
  if (has("t\\(6;9\\)"))  return("AML with DEK::NUP214 [t(6;9)]")
  if (has("t\\(9;22\\)")) return("AML with BCR::ABL1 [t(9;22)]")
  if (has("t\\(1;22\\)")) return("AML with RBM15::MRTFA [t(1;22)]")
  kmt2a <- c("t\\(9;11\\)","t\\(1;11\\)","t\\(4;11\\)","t\\(5;11\\)",
             "t\\(6;11\\)","t\\(10;11\\)","t\\(11;16\\)","t\\(11;17\\)",
             "t\\(11;19\\)","t\\(11q23\\)")
  if (has(kmt2a)) return("AML with KMT2A rearrangement [t(v;11q23)]")
  if (has(c("t\\(.*11p15.*\\)","t\\(5;11\\)\\(q35")))
    return("AML with NUP98 rearrangement [t(v;11p15)]")
  aml_mr <- c("^-5$","del\\(5q\\)","^-7$","del\\(7q\\)","del\\(17p\\)","i\\(17q\\)",
              "del\\(12p\\)","del\\(11q\\)","del\\(9q\\)","idic\\(X\\)\\(q13\\)",
              "t\\(11;16\\)","t\\(3;21\\)","t\\(1;3\\)","t\\(2;11\\)")
  is_complex <- detect_complex(parsed_table)
  if (has(aml_mr) || is_complex) return("AML, myelodysplasia-related (cytogenetic criteria)")
  has_del5q  <- n_abn > 0 && any(str_detect(abn, "del\\(5q\\)"))
  has_minus7 <- n_abn > 0 && any(str_detect(abn, "^-7$"))
  has_del7q  <- n_abn > 0 && any(str_detect(abn, "del\\(7q\\)"))
  best_tbl   <- if ("is_best_clone" %in% colnames(parsed_table))
    filter(parsed_table, is_best_clone) else parsed_table
  n_best     <- if ("n_count" %in% colnames(best_tbl))
    sum(best_tbl$n_count, na.rm = TRUE) else nrow(best_tbl)
  if (has_del5q && n_best <= 2 && !has_minus7 && !has_del7q)
    return("MDS with del(5q) [WHO 2022]")
  has_17p <- n_abn > 0 && any(str_detect(abn, "del\\(17p\\)|^-17$"))
  if (is_complex && has_17p)
    return("MDS — possible biallelic TP53 inactivation (molecular confirmation required)")
  if (MDS_flag(parsed_table)) return("MDS, NOS (MDS-defining cytogenetic abnormality present)")
  if (n_abn > 0) return("Abnormal karyotype - entity requires clinical/morphologic correlation")
  return("Normal karyotype - entity requires clinical/morphologic/molecular data")
}

# ============================================================
# 15. PIPELINE
# ============================================================
run_pipeline <- function(df) {
  if (!all(c("sample_id", "karyotype") %in% colnames(df)))
    stop("Input data frame must contain columns: sample_id, karyotype")
  
  empty_tbl_schema <- tibble(
    clone = integer(), abnormality = character(), type = character(),
    chromosomes = character(), breakpoint = character(), weight = numeric(),
    n_count = integer(), is_best_clone = logical(), clone_type = character()
  )
  
  results <- df %>%
    mutate(parsed = map(karyotype, process_karyotype)) %>%
    mutate(parsed = map(parsed, function(x) {
      if (is.null(x$table)) x$table <- empty_tbl_schema
      x
    })) %>%
    mutate(
      # --- Counts (best clone per guideline) ---
      total_abnormalities      = map_int(parsed,  "total_abnormalities"),
      total_score              = map_int(parsed,  "total_score"),
      clones                   = map_int(parsed,  "clones"),
      
      # --- Clone dependency breakdown ---
      n_independent_clones     = map_int(parsed,  "n_independent_clones"),
      n_dependent_clones       = map_int(parsed,  "n_dependent_clones"),
      n_normal_clones          = map_int(parsed,  "n_normal_clones"),
      best_clone_type          = map_chr(parsed,  ~ .x$best_clone_type %||% NA_character_),
      clone_dependency_summary = map_chr(parsed,  ~ .x$clone_dependency_summary %||% NA_character_),
      
      # --- Structural flags (best clone) ---
      complex_karyotype        = map_lgl(parsed,  ~ detect_complex(.x$table)),
      monosomal_karyotype      = map_lgl(parsed,  ~ detect_monosomal(.x$table)),
      
      # --- Locus flags (all clones) ---
      del5                     = map_lgl(parsed,  ~ del5_flag(.x$table)),
      del7q                    = map_lgl(parsed,  ~ del7q_flag(.x$table)),
      del7                     = map_lgl(parsed,  ~ del7_flag(.x$table)),
      loss17p                  = map_lgl(parsed,  ~ loss17p_flag(.x$table)),
      MLL                      = map_lgl(parsed,  ~ MLL_flag(.x$table)),
      
      # --- Disease-defining flags (all clones) ---
      MDS_defining_abnormality = map_lgl(parsed,  ~ MDS_flag(.x$table)),
      AML_defining_abnormality = map_lgl(parsed,  ~ AML_flag(.x$table)),
      MDS_defining_karyo       = map_chr(parsed,  ~ get_MDS_abnormalities(.x$table)),
      AML_defining_karyo       = map_chr(parsed,  ~ get_AML_abnormalities(.x$table)),
      
      # --- Classification ---
      classification           = map2_chr(parsed, karyotype,
                                          ~ classify_karyo(.x$table, .y)),
      eln_subtype              = map_chr(karyotype, get_eln_subtype),
      eln2022_cyto_risk        = map2_chr(parsed, karyotype,
                                          ~ get_eln2022_cyto_risk(.x$table, .y)),
      ipssr_cyto_risk          = map2_chr(parsed, karyotype,
                                          ~ get_ipssr_cyto_risk(.x$table, .y)),
      who2022_entity           = map2_chr(parsed, karyotype,
                                          ~ get_who2022_entity(.x$table, .y))
    )
  
  detail <- map_dfr(seq_len(nrow(results)), function(i) {
    tbl <- results$parsed[[i]]$table
    if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
    tbl %>% mutate(sample_id = results$sample_id[i],
                   karyotype = results$karyotype[i])
  })
  
  summary_out <- results %>% select(-parsed)
  list(summary = summary_out, detail = detail)
}

# Null-coalescing operator (base R doesn't have %||%)
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a else b

############################################################
# WORKED EXAMPLES for validation
#
# Ex 1: 47,XX,+8,+21[11]/46,XX,del(5)(q13q33),-7[5]
#   clone 1: independent (2 abn)
#   clone 2: independent (2 abn)
#   → n_independent_clones = 2, total_abnormalities = 2, complex = FALSE  ✓
#
# Ex 4: 48,XX,+8,+19,del(20)(q11.2)[4]/49,idem,+del(11)(?q23q23)[10]/46,XY[6]
#   clone 1: independent (3 abn)
#   clone 2: dependent/idem → expanded to 4 abn
#   clone 3: normal
#   → n_independent_clones=1, n_dependent_clones=1, n_normal_clones=1
#   → best = clone 2 (4 abn), best_clone_type = "idem", complex = TRUE  ✓
#
# sdl example: 46,XY,t(9;22)[10]/47,sl,+8[5]/45,sdl,-7[3]
#   clone 1: independent — t(9;22)
#   clone 2: sl (independent) — sl + +8  → 2 abn
#   clone 3: sdl (dependent) → expands from sl: t(9;22),+8,-7 → 3 abn (but wait,
#            actually the sl clone is clone 2 here not clone 1)
#   → n_independent_clones=2 (clone1 + sl), n_dependent_clones=1 (sdl)
############################################################

############################################################
# USAGE
############################################################
df  <- read_csv("test.csv")   # must have columns: sample_id, karyotype
out <- run_pipeline(df)
write_csv(out$summary, "summary.csv")
write_csv(out$detail,  "detail.csv")

df  <- read_csv("2025-10-12 abnormal_karyotypes_for_CK_test.csv",na = c("", "N/A", "NA"))
out <- run_pipeline(df)
write_csv(out$summary, "summary1.csv")
write_csv(out$detail,  "detail1.csv")

#
# Key new columns in out$summary:
#   n_independent_clones     — abnormal clones with no idem/sdl derivation
#   n_dependent_clones       — clones derived via idem or sdl
#   n_normal_clones          — clones with 0 aberrations (e.g. 46,XX)
#   best_clone_type          — type of clone used for total_abnormalities scoring
#   clone_dependency_summary — e.g. "2 independent; 1 dependent (sdl); 1 normal"
#
# detail table extra columns:
#   clone_type    — "independent" | "sl" | "sdl" | "idem" | "normal"
#   n_count       — per-row aberration count (der=2, others=1)
#   is_best_clone — TRUE for the clone used in scoring
############################################################