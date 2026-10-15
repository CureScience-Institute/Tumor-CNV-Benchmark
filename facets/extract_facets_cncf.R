#!/usr/bin/env Rscript
# Pulls the per-segment integer copy-number table out of run_facets.py's
# saved <patient>.rds (full facetsSuite fit object) and writes it as a plain
# TSV -- the only FACETS output today with integer copy number at all is this
# .rds; the sibling .seg file only has continuous seg.mean (log2 ratio), and
# the summary .txt has no per-segment breakdown, purity/ploidy only.
#
# Deliberately does NOT call facets::emcncf() itself and does not need the
# `facets` R package installed to run -- run-facets-wrapper.R (inside the
# facets-suite Docker image, see run_facets.py) already ran that during the
# original fit and saved its result as a plain list of data.frames in the
# .rds (confirmed live by reading a real saved object: names snps/segs/
# purity/ploidy/diplogr/alballogr/flags/em_flags/loglik -- no facets S4
# classes involved, base R readRDS() deserializes it fine on its own). This
# script only reads that list back out, so it can run anywhere plain R is
# available (e.g. bioinf-vm), not just inside the facets-suite container.
#
# segs' tcn.em/lcn.em (not the pre-EM tcn/lcn) are used as the final integer
# total/minor copy number -- tcn.em is FACETS' own EM-refined call, the same
# one facetsSuite reports as the definitive per-segment copy number. lcn.em
# is sometimes NA (e.g. a clonal-fraction-adjusted segment with no confident
# minor-allele call) even when tcn.em is a real integer -- that's expected
# and left as NA in the output; total CN is still usable on its own for the
# gene-CN report this feeds.
#
# Usage:
#   Rscript extract_facets_cncf.R <patient>.rds <out.tsv>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: extract_facets_cncf.R <patient>.rds <out.tsv>")
}
rds_path <- args[1]
out_tsv <- args[2]

fit <- readRDS(rds_path)
segs <- fit$segs

out <- data.frame(
  chrom = segs$chrom,
  start = segs$start,
  end = segs$end,
  num_mark = segs$num.mark,
  nhet = segs$nhet,
  tcn = segs$tcn.em,
  lcn = segs$lcn.em
)

write.table(out, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf(
  "%d segments -> %s (purity=%.3f, ploidy=%.3f, diplogr=%.3f)\n",
  nrow(out), out_tsv, fit$purity, fit$ploidy, fit$diplogr
))
