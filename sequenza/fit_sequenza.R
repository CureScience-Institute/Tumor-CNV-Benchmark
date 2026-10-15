#!/usr/bin/env Rscript
# Fits Sequenza on a binned seqz file -- the R half of pipeline/sequenza/
# (run_sequenza.py drives sequenza-utils bam2seqz/seqz_binning first; this
# script is everything after that). Unlike FACETS' facetsSuite, sequenza
# ships no ready-made CLI wrapper -- this is the whole wrapper, calling the
# three documented entry points in sequence: sequenza.extract() ->
# sequenza.fit() -> sequenza.results(), the last of which writes the
# standard output bundle (segments table, purity/ploidy confidence-interval
# table, alternative-solutions table, genome-view/model-fit/chromosome-view
# PDFs, a saved extract .RData) straight into --out-dir. Signatures verified
# directly against a real r-sequenza 3.0.0 install (args(sequenza.extract)
# etc.), not from memory or docs.
suppressPackageStartupMessages(library(sequenza))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = "Fit Sequenza on a binned seqz file")
parser$add_argument("--seqz-file", required = TRUE, help = "Binned seqz.gz from sequenza-utils seqz_binning")
parser$add_argument("--sample-id", required = TRUE)
parser$add_argument("--out-dir", required = TRUE)
# female is deliberately required, not defaulted -- sequenza.fit()/results()
# default to female=TRUE, which this wrapper does not want to silently
# inherit for a male patient's X/Y copy-number calls.
parser$add_argument("--female", required = TRUE, choices = c("true", "false"))
# sequenza.extract()'s own default (80) -- its dnacopy/PCF segmentation
# penalty, the closest thing here to FACETS' cval. Left as a plain float
# default here (not NULL) since sequenza.extract() itself requires a value.
parser$add_argument("--gamma", type = "double", default = 80)
args <- parser$parse_args()

female <- args$female == "true"
dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

# sequenza.extract() wraps copynumber's segmentation in try(), so an
# hg38-unaware copynumber (Bioconductor's build) only prints "assembly must be
# one of hg19..." per chromosome and yields one unsegmented segment each.
# Fail loudly instead -- run_sequenza.py installs the hg38-capable fork.
if (!exists("hg38", envir = asNamespace("copynumber"), inherits = FALSE)) {
  stop("installed copynumber package has no hg38 support -- see COPYNUMBER_HG38_URL in run_sequenza.py")
}
# chr1-22 + X only. Confirmed live (CS-26-0090, female): stray chrY
# positions from mismapped reads made sequenza.results() crash in
# seq.default(by = 1/CNn), since female=TRUE sets chrY's normal copy number
# to 0. Excluding Y is Sequenza's own documented setup; X stays, with its
# sex-aware normal copy number.
chromosomes <- paste0("chr", c(1:22, "X"))

# sequenza's find.breaks() calls copynumber::winsorize(logR) WITHOUT an
# assembly argument, so outlier winsorizing used hg19 arm boundaries even
# though segmentation itself gets assembly = "hg38" (r-sequenza 3.0.0
# segments.R; the fork's README says hg38 must be passed explicitly to each
# call). Make hg38 the default inside the copynumber namespace so every
# unqualified call uses the right centromeres.
hg38_winsorize <- copynumber::winsorize
formals(hg38_winsorize)$assembly <- "hg38"
utils::assignInNamespace("winsorize", hg38_winsorize, ns = "copynumber")
stopifnot(identical(formals(copynumber::winsorize)$assembly, "hg38"))

extract <- sequenza.extract(args$seqz_file, gamma = args$gamma, assembly = "hg38",
                            chromosome.list = chromosomes, verbose = TRUE)
# sequenza.extract() silently drops requested chromosomes that are absent
# from the seqz file (extract.R: chromosome.list[chromosome.list %in% chr.vect]),
# so a truncated or mis-ordered seqz would fit on a partial genome without
# any error. Fail instead.
missing_chr <- setdiff(chromosomes, extract$chromosomes)
if (length(missing_chr) > 0) {
  stop("seqz is missing chromosomes: ", paste(missing_chr, collapse = ", "))
}
# XY must name the chr-prefixed contigs: sequenza.fit()/results() default to
# XY = c(X = "X", Y = "Y") and pick sex chromosomes with `%in% XY`, so on
# GRCh38 names a male's chrX was silently modelled as a 2-copy autosome
# (confirmed in the installed r-sequenza 3.0.0 source, 2026-09-24; same fix
# as OICR's production wrapper, github.com/oicr-gsi/sequenza).
xy <- c(X = "chrX", Y = "chrY")
fit <- sequenza.fit(extract, female = female, XY = xy)
sequenza.results(extract, fit, sample.id = args$sample_id, out.dir = args$out_dir, female = female, XY = xy)
