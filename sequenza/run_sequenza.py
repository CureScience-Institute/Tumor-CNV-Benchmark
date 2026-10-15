#!/usr/bin/env python3
"""Standalone Sequenza (allele-specific tumor/normal copy-number +
purity/ploidy) run, invoked as a plain Python script from a Batch task --
no Nextflow, no nf-core/sarek integration (Sequenza has none). Same
relationship to sarek's output pipeline/facets/ already has: reads sarek's
own recalibrated tumor/normal CRAMs straight out of blob, as an
independent, separately-scheduled Batch job.

Engine: a throwaway per-task conda env (bioconda `sequenza-utils` +
`r-sequenza` + `r-argparse`), installed fresh into node-local scratch --
NOT the pool image's own baked Miniconda. That install lives at
$HOME/miniconda3 for the image-build admin user (see
pipeline/cloud-resources/provision-build-vm.sh), a fixed path tied to that
user's home directory; the Batch task itself runs as a distinct per-task
autoUser with its own $HOME, so the baked install's actual path and
permissions from a different user's context were never verified and
weren't worth depending on. A fresh conda env costs about the same as
facets' ~1GB Docker pull and matches this repo's established "install
fresh per task, no ACR/baked dependency" convention (cloud-neoantigen's
apt/pip installs, facets' Docker pull) instead of introducing a new one.

CONFIRMED LIVE (2026-09-17, not just read about): bioconda's current
`sequenza-utils` 3.0.0 crashes on every real subcommand invocation under
Python 3.12 (`argparse.ArgumentError: conflicting subparser`) -- traced to
sequenza's own commands.py double-registering each subcommand's parser.
Reproduced fixed by pinning `python=3.10` for the conda env; this is the
reason CONDA_PACKAGES below pins python explicitly rather than letting
conda solve to its newest available build. See pipeline/sequenza/README.md.

Pipeline:
  1. `samtools view -T <reference> -L <capture BED> -b` converts sarek's
     recalibrated CRAMs to BAM (same byte-identical-reference requirement
     as pipeline/facets/ -- reuses the same pipeline/resources/gatk-grch38-fasta/).
     Uses the conda env's own bundled samtools, not a separately apt-installed
     one -- confirmed sequenza-utils' bioconda recipe already pulls in samtools
     and tabix as dependencies, so nothing else needs installing.
  2. `sequenza-utils bam2seqz` merges both BAMs + the reference + a
     genome-wide GC-content wiggle track (pipeline/resources/sequenza-gc-hg38/,
     patient-agnostic, generated once) into a raw seqz file.
  3. `sequenza-utils seqz_binning -w 50` reduces it to a binned seqz file --
     written into --out-dir (not scratch): small enough to keep, and lets a
     later re-run of just fit_sequenza.R (different --gamma) skip re-doing
     steps 1/2.
  4. `fit_sequenza.R` (this repo's own script -- Sequenza ships no
     ready-made CLI wrapper the way facetsSuite does) calls
     sequenza.extract() -> sequenza.fit() -> sequenza.results(), writing the
     standard output bundle into --out-dir.

Deliberately does NOT feed its purity/ploidy estimate back into sarek's own
--sarek-ascat-ploidy -- see pipeline/facets/README.md's identical rationale,
which applies here unchanged.
"""
import argparse
import gzip
import hashlib
import os
import shutil
import subprocess
import tempfile
import urllib.request
from pathlib import Path

# Pinned (was "latest"): the build every validated run so far actually got
# -- "latest" resolved to exactly this file on 2026-09-24 (same sha256).
MINICONDA_URL = "https://repo.anaconda.com/miniconda/Miniconda3-py314_26.7.1-1-Linux-x86_64.sh"
MINICONDA_SHA256 = "e8b25b92b262499141c5bd57a98d3c008024185fa951494b9cd9b6d94e72338b"
CONDA_ENV_NAME = "sequenza"
# r-iotools/r-bit64 pinned: confirmed live that r-sequenza's compiled
# gc.sample.stats() calls iotools::chunk.apply(..., parallel=parallel) --
# a parameter iotools renamed to CH.PARALLEL in newer releases (conda's
# solved default, 0.3.5, crashes with "unused argument (parallel = 1)").
# r-iotools=0.3.2 is the newest release still exposing the old `parallel`
# name. That downgrade in turn needs r-bit64 pinned back too: bit64's
# latest build (4.8.2) calls a base-R S4 method (chooseOpsMethod) that
# doesn't exist before R 4.3, and conda's own r-sequenza recipe resolves
# R 4.2 -- 4.0.5 is the newest bit64 confirmed to load cleanly under that
# R version.
CONDA_PACKAGES = [
    # Versions pinned to what the validated runs resolved (2026-09-24):
    # sequenza-utils 3.0.0, r-sequenza 3.0.0, samtools/htslib 1.24.
    "python=3.10", "sequenza-utils=3.0.0", "r-sequenza=3.0.0", "samtools=1.24", "htslib=1.24", "r-argparse",
    "r-iotools=0.3.2", "r-bit64=4.0.5",
]

# hg38-capable drop-in for Bioconductor's `copynumber`, which sequenza.extract()
# uses for segmentation. Bioconductor's build only knows hg16-hg19: confirmed
# live (2026-09-23) that assembly="hg38" makes it fail with "assembly must be
# one of hg19, hg18, hg17 or hg16" on EVERY chromosome, after which Sequenza
# silently falls back to one unsegmented segment per chromosome. The fork
# (recommended by Sequenza's authors for hg38) only adds hg38 cytobands;
# pinned to a commit so the install is reproducible.
COPYNUMBER_HG38_URL = (
    "https://github.com/aroneklund/copynumber/archive/"
    "2f54c083a4a69a76fcf441a1e411085e0174b5d9.tar.gz"
)

PRIMARY_CHROMOSOMES = [f"chr{c}" for c in [*range(1, 23), "X"]]

REFERENCE_FASTA_REL = "pipeline/resources/gatk-grch38-fasta/Homo_sapiens_assembly38.fasta"
GC_WIGGLE_REL = "pipeline/resources/sequenza-gc-hg38/gc50Base.wig.gz"
FIT_SCRIPT_REL = "pipeline/sequenza/fit_sequenza.R"

OUTPUT_ROOT = Path("/data4/NewOutput")

SEQUENZA_ENV = Path(
    "/data4/local-bioinformatics/RNA-CNV/sequenza/env"
)

REFERENCE_FASTA = Path(
    "/data4/SarekReferences/gatk-grch38-fasta/"
    "Homo_sapiens_assembly38.fasta"
)

GC_WIGGLE = Path(
    "/data4/SarekReferences/sequenza-gc-hg38/"
    "gc50Base.wig.gz"
)

FIT_SCRIPT = Path(
    "/data4/local-bioinformatics/RNA-CNV/sequenza/"
    "fit_sequenza.R"
)

BG_CAPTURE_BED = Path(
    "/data4/cloud-bioinformatics_v1/BED/"
    "S33613271_Regions.cleaned.sorted.bed"
)

CEGAT_CAPTURE_BED = Path(
    "/data4/SarekReferences/CeGaT_bed/"
    "MANE_Select_v1.3_GRCh38_CDS.GATK_compatible.bed"
)


def run(cmd, **kwargs):
    print(f"+ {' '.join(str(c) for c in cmd)}", flush=True)
    subprocess.run(cmd, check=True, **kwargs)


def install_conda_env(scratch_dir):
    """Fresh Miniconda + a python=3.10-pinned conda env, entirely inside
    scratch_dir -- dies with the node like everything else here, no
    persistence needed since nothing amortizes it across runs yet."""
    installer = scratch_dir / "miniconda.sh"
    urllib.request.urlretrieve(MINICONDA_URL, installer)
    digest = hashlib.sha256(installer.read_bytes()).hexdigest()
    if digest != MINICONDA_SHA256:
        raise RuntimeError(f"Miniconda installer sha256 {digest} != pinned {MINICONDA_SHA256}")
    conda_root = scratch_dir / "miniconda3"
    run(["bash", str(installer), "-b", "-p", str(conda_root)])
    conda_bin = str(conda_root / "bin" / "conda")
    # Non-interactive ToS acceptance for the two default Anaconda channels --
    # conda create otherwise refuses to solve at all (CondaToSNonInteractiveError).
    for channel in ("https://repo.anaconda.com/pkgs/main", "https://repo.anaconda.com/pkgs/r"):
        run([conda_bin, "tos", "accept", "--override-channels", "--channel", channel])
    run([conda_bin, "create", "-y", "-n", CONDA_ENV_NAME, "-c", "bioconda", "-c", "conda-forge", *CONDA_PACKAGES])
    env_bin = conda_root / "envs" / CONDA_ENV_NAME / "bin"
    # Replace the env's Bioconductor copynumber with the hg38-capable fork
    # (see COPYNUMBER_HG38_URL). Pure R, and its dependencies are already in
    # the env as the original package's dependencies.
    tarball = scratch_dir / "copynumber-hg38.tar.gz"
    urllib.request.urlretrieve(COPYNUMBER_HG38_URL, tarball)
    # Explicit -l into the env's own library, and R_PROFILE_USER blanked so a
    # personal ~/.Rprofile can't redirect it: confirmed live on bioinf-vm that
    # one prepending ~/R/library to .libPaths() made R CMD INSTALL put the fork
    # there instead, where the --vanilla Rscript fit below never looks.
    env_r_library = env_bin.parent / "lib" / "R" / "library"
    run(
        [str(env_bin / "R"), "CMD", "INSTALL", "--no-test-load", "-l", str(env_r_library), str(tarball)],
        env={**os.environ, "R_PROFILE_USER": os.devnull, "R_LIBS_USER": os.devnull},
    )
    return env_bin


def _ordered_subset(names, wanted):
    """`wanted` names in the order they appear in `names`."""
    wanted_set = set(wanted)
    return [n for n in names if n in wanted_set]


def check_contig_order(gc_wiggle, fasta_fai, bam, chromosomes):
    """Fail unless `chromosomes` appear, all present, in the SAME relative
    order in the GC wiggle, the reference .fai and the BAM header.

    bam2seqz merges the pileup and the wiggle forward-only (sequenza-utils
    izip.py): if their contig orders differ it silently skips whole
    chromosomes -- the same failure class as the out-of-order FACETS SNP VCF
    (2026-09-24). Safe today only because all three share GATK order; this
    turns a future regenerated wiggle or reordered PRIMARY_CHROMOSOMES into
    a loud error instead of missing data.
    """
    with gzip.open(gc_wiggle, "rt") as fh:
        wig = [l.split("chrom=")[1].split()[0] for l in fh if l.startswith("variableStep")]
    with open(fasta_fai) as fh:
        fai = [l.split("\t")[0] for l in fh]
    header = subprocess.run(["samtools", "view", "-H", str(bam)], check=True,
                            capture_output=True, text=True).stdout
    sq = [f.split(":", 1)[1] for l in header.splitlines() if l.startswith("@SQ")
          for f in l.split("\t") if f.startswith("SN:")]
    for label, names in (("GC wiggle", wig), ("reference .fai", fai), ("BAM header", sq)):
        got = _ordered_subset(names, chromosomes)
        if got != list(chromosomes):
            missing = sorted(set(chromosomes) - set(got))
            raise RuntimeError(f"{label}: contigs missing {missing} or out of order vs {chromosomes}: {got}")


def cram_to_bam(cram_path, bam_path, reference_fasta, capture_bed, threads):
    """Restricted to the capture panel BED, same rationale as
    pipeline/facets/run_facets.py's identical helper: a WES CRAM's
    off-target reads are never looked at downstream."""
    run([
        "samtools", "view", "-@", str(threads), "-T", str(reference_fasta),
        "-L", str(capture_bed), "-b", "-o", str(bam_path), str(cram_path),
    ])
    run(["samtools", "index", str(bam_path)])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    ap.add_argument(
        "patient_id",
        nargs="?",
    )

    ap.add_argument(
        "--sex",
        choices=["male", "female"],
        default=None,
    )

    ap.add_argument(
        "--bed-type",
        choices=["1", "2"],
        default=None,
    )

    ap.add_argument("--threads", type=int, default=40)
    ap.add_argument("--gamma", type=float, default=None, help="Default per sequenza.extract(): 80")
    args = ap.parse_args()

    patient_id = args.patient_id
    if not patient_id:
        patient_id = input("Patient folder name: ").strip()

    sex = args.sex
    if not sex:
        sex_input = input("Sex (1 = Male, 2 = Female): ").strip()
        if sex_input == "1":
            sex = "male"
        elif sex_input == "2":
            sex = "female"
        else:
            raise SystemExit("Sex must be 1 (Male) or 2 (Female)")

    bed_type = args.bed_type
    if not bed_type:
        bed_type = input("BED type (1 = BostonGene, 2 = CeGaT): ").strip()

    if bed_type == "1":
        capture_bed = BG_CAPTURE_BED
    elif bed_type == "2":
        capture_bed = CEGAT_CAPTURE_BED
    else:
        raise SystemExit("BED type must be 1 (BostonGene) or 2 (CeGaT)")

    female = "true" if sex == "female" else "false"

    patient_root = OUTPUT_ROOT / patient_id
    out_dir = patient_root / "sequenza"
    out_dir.mkdir(parents=True, exist_ok=True)

    reference_fasta = REFERENCE_FASTA
    gc_wiggle = GC_WIGGLE
    fit_script = FIT_SCRIPT

    recal_dir = (
        patient_root
        / "sarek-results"
        / "preprocessing"
        / "recalibrated"
    )

    normal_cram = (
        recal_dir
        / f"{patient_id}_normal"
        / f"{patient_id}_normal.recal.cram"
    )

    tumor_cram = (
        recal_dir
        / f"{patient_id}_tumor"
        / f"{patient_id}_tumor.recal.cram"
    )

    print()
    print("============================================================")
    print("Sequenza local run")
    print("============================================================")
    print(f"Patient:       {patient_id}")
    print(f"Sex:           {sex}")
    print(f"BED type:      {'BostonGene' if bed_type == '1' else 'CeGaT'}")
    print(f"Capture BED:   {capture_bed}")
    print(f"Reference:     {reference_fasta}")
    print(f"GC wiggle:     {gc_wiggle}")
    print(f"Sequenza env:  {SEQUENZA_ENV}")
    print(f"Fit script:    {fit_script}")
    print(f"Normal CRAM:   {normal_cram}")
    print(f"Tumor CRAM:    {tumor_cram}")
    print(f"Output:        {out_dir}")
    print(f"Threads:       {args.threads}")
    print("============================================================")
    print()

    scratch_dir = Path(tempfile.mkdtemp(prefix="sequenza_"))
    try:
        env_bin = SEQUENZA_ENV / "bin"

        # Prepend the conda env so every bare command below (samtools,
        # tabix, sequenza-utils, Rscript, and sequenza-utils' OWN internal
        # shell-outs to samtools/tabix) resolves to this env's builds, not
        # whatever else might be on the node's PATH.
        os.environ["PATH"] = f"{env_bin}{os.pathsep}{os.environ['PATH']}"

        normal_bam = scratch_dir / "normal.bam"
        tumor_bam = scratch_dir / "tumor.bam"
        cram_to_bam(normal_cram, normal_bam, reference_fasta, capture_bed, args.threads)
        cram_to_bam(tumor_cram, tumor_bam, reference_fasta, capture_bed, args.threads)

        check_contig_order(gc_wiggle, Path(f"{reference_fasta}.fai"), normal_bam, PRIMARY_CHROMOSOMES)
        check_contig_order(gc_wiggle, Path(f"{reference_fasta}.fai"), tumor_bam, PRIMARY_CHROMOSOMES)

        raw_seqz = scratch_dir / f"{patient_id}.raw.seqz.gz"
        binned_seqz = out_dir / f"{patient_id}.binned.seqz.gz"
        run([
            "sequenza-utils", "bam2seqz",
            "-n", str(normal_bam), "-t", str(tumor_bam),
            "-F", str(reference_fasta), "-gc", str(gc_wiggle),
            # Primary chromosomes only: otherwise bam2seqz walks all 3,366
            # contigs of the GATK reference (alt/HLA/decoy) for no usable
            # signal. chrY is dropped here too -- fit_sequenza.R already
            # excludes it from the fit (female crash, see there).
            "-C", *PRIMARY_CHROMOSOMES,
            "-o", str(raw_seqz),
        ])
        run([
            "sequenza-utils", "seqz_binning",
            "-w", "50",
            "-s", str(raw_seqz),
            "-o", str(binned_seqz),
        ])

        fit_cmd = [
            # --vanilla (implies --no-init-file) skips any personal
            # .Rprofile -- confirmed live that a pre-existing
            # ~/.Rprofile on this dev VM (from unrelated prior R work)
            # unconditionally prepends ~/R/library to .libPaths(), which
            # shadowed this conda env's own correctly-pinned package
            # versions (bit64/iotools) with incompatible ones and broke
            # sequenza.extract() with cryptic "unused argument"/"object
            # not found" errors that had nothing to do with the actual
            # conda environment being wrong.
            "Rscript", "--vanilla", str(fit_script),
            "--seqz-file", str(binned_seqz),
            "--sample-id", patient_id,
            "--out-dir", str(out_dir),
            "--female", female,
        ]

        if args.gamma is not None:
            fit_cmd += ["--gamma", str(args.gamma)]

        run(fit_cmd)

    finally:
        # Node-local scratch (including the throwaway conda env itself)
        # dies with the VM anyway; cleaned up early so a retry doesn't trip
        # over a half-finished install.
        shutil.rmtree(scratch_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
