#!/usr/bin/env python3
"""FACETS (allele-specific tumor/normal copy-number + purity/ploidy) run,
invoked as a plain Python script -- no Nextflow, no nf-core/sarek
integration (FACETS has none). Chained directly into render_sarek()'s own
Batch task (see pipeline/catalog/generate_submit_scripts.py), running right
after sarek's own nextflow command exits, reading sarek's just-produced
recalibrated tumor/normal CRAMs straight off local disk via --recal-dir
(no blob round-trip needed, since they're already sitting in sarek's own
$OUTDIR). Can also be run as a fully separate job against blob-staged CRAMs
by omitting --recal-dir -- see pipeline/facets/README.md for why chaining
was chosen over that original standalone design.

Pipeline, all via MSKCC's public `philipjonsson/facets-suite:dev` Docker
image (bundles facets + facetsSuite + a prebuilt snp-pileup binary -- see
https://github.com/mskcc/facets-suite). Note the deployed image's actual
`run-facets-wrapper.R`/`snp-pileup-wrapper.R` (confirmed 2026-09-17 by
`docker run ... cat /usr/bin/run-facets-wrapper.R`) differ in real ways from
the current GitHub `master` source: `--cval` defaults to 100 (not 50),
`--purity-cval` defaults to NULL/single-pass (not 100), and there is no
`--facets-lib-path` flag at all -- verify against the live container, not
GitHub, before changing any flag here.

  1. `samtools view -T <reference> -b` converts sarek's recalibrated CRAMs to
     BAM. This step is the reason a byte-identical copy of the GATK.GRCh38
     reference sarek recalibrated against (pipeline/resources/gatk-grch38-fasta/)
     is required -- htslib validates a CRAM's embedded per-contig checksum
     against whatever reference it's handed, so any other GRCh38 FASTA (e.g.
     pipeline/resources/grch38-fasta/, a different GENCODE build used
     elsewhere in this repo) would mis-decode or hard-fail.
  2. `snp-pileup-wrapper.R` (inside the container) pileups both BAMs against
     an exome-restricted common-SNPs VCF (pipeline/resources/facets-wes-hg38/)
     -> a gzipped counts file. Normal is passed before tumor -- FACETS'
     column convention (File1 = normal, File2 = tumor) is order-dependent,
     not name-dependent.
  3. `run-facets-wrapper.R` (inside the container) fits FACETS on that counts
     file with `--genome hg38` (the wrapper's own default is hg19), writing
     the purity/ploidy summary (`<sample>.txt`), an IGV-style segment file,
     and a diagnostic PNG into --out-dir/<sample-id>/ (the wrapper nests its
     own output one directory deeper, named after --sample-id). Deliberately
     NOT passing --everything -- confirmed live that its NTAI calculation
     crashes unconditionally on this deployed image for a WES sample's
     sparse per-chromosome segmentation, and the core purity/ploidy output
     never depended on it anyway (only extra gene/arm-level/QC files did).

Deliberately does NOT feed its purity/ploidy estimate back into sarek's own
`--sarek-ascat-ploidy` -- see pipeline/facets/README.md. That stays a manual,
human-gated re-run, same as every other cross-pipeline hand-off in this repo.

The CRAM->BAM conversion writes to a node-local scratch dir, not --out-dir:
those BAMs are a disposable intermediate (the CRAM they came from already
made it to blob via sarek's own upload), so re-uploading them from every
facets run would just waste blob storage and upload time.
"""
import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

FACETS_SUITE_IMAGE = "philipjonsson/facets-suite:dev"

OUTPUT_ROOT = Path("/data4/NewOutput")

REFERENCE_FASTA = Path(
    "/data4/SarekReferences/gatk-grch38-fasta/"
    "Homo_sapiens_assembly38.fasta"
)

SNP_VCF = Path(
    "/data4/SarekReferences/facets-wes-hg38/"
    "common_snps_on_target_hg38.vcf.gz"
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


def cram_to_bam(cram_path, bam_path, reference_fasta, capture_bed, threads):
    """samtools view + index. Restricted to the capture panel BED (the same
    one sarek/ASCAT already use) since a WES CRAM has off-target reads that
    snp-pileup will never look at -- skipping them keeps the intermediate
    BAM and the pileup step both smaller and faster."""
    run([
        "samtools", "view", "-@", str(threads), "-T", str(reference_fasta),
        "-L", str(capture_bed), "-b", "-o", str(bam_path), str(cram_path),
    ])
    run(["samtools", "index", str(bam_path)])


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    ap.add_argument(
        "patient_id",
        nargs="?",
        help="Patient folder name under /data4/NewOutput",
    )

    ap.add_argument(
        "--bed-type",
        choices=["1", "2"],
        default=None,
        help="1 = BostonGene, 2 = CeGaT",
    )

    ap.add_argument("--threads", type=int, default=40)
    ap.add_argument("--cval", type=int, default=None, help="Deployed image's own default: 100")
    ap.add_argument("--purity-cval", type=int, default=None, help="Deployed image's own default: NULL (single-pass)")
    ap.add_argument("--seed", type=int, default=None, help="Deployed image's own default: 100")

    args = ap.parse_args()

    patient_id = args.patient_id
    if not patient_id:
        patient_id = input("Patient folder name: ").strip()

    bed_type = args.bed_type
    if not bed_type:
        bed_type = input("BED type (1 = BostonGene, 2 = CeGaT): ").strip()

    if bed_type == "1":
        capture_bed = BG_CAPTURE_BED
    elif bed_type == "2":
        capture_bed = CEGAT_CAPTURE_BED
    else:
        raise SystemExit("BED type must be 1 (BostonGene) or 2 (CeGaT)")

    patient_root = OUTPUT_ROOT / patient_id

    out_dir = patient_root / "facets"
    out_dir.mkdir(parents=True, exist_ok=True)

    reference_fasta = REFERENCE_FASTA
    snp_vcf = SNP_VCF

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
    print("FACETS local run")
    print("============================================================")
    print(f"Patient:       {patient_id}")
    print(f"BED type:      {'BostonGene' if bed_type == '1' else 'CeGaT'}")
    print(f"Capture BED:   {capture_bed}")
    print(f"Reference:     {reference_fasta}")
    print(f"SNP VCF:       {snp_vcf}")
    print(f"Normal CRAM:   {normal_cram}")
    print(f"Tumor CRAM:    {tumor_cram}")
    print(f"Output:        {out_dir}")
    print(f"Threads:       {args.threads}")
    print("============================================================")
    print()

    run(["docker", "pull", FACETS_SUITE_IMAGE])

    scratch_dir = Path(tempfile.mkdtemp(prefix="facets_"))
    try:
        normal_bam = scratch_dir / "normal.bam"
        tumor_bam = scratch_dir / "tumor.bam"

        cram_to_bam(
            normal_cram,
            normal_bam,
            reference_fasta,
            capture_bed,
            args.threads,
        )

        cram_to_bam(
            tumor_cram,
            tumor_bam,
            reference_fasta,
            capture_bed,
            args.threads,
        )

        docker_mounts = [
            "-u", f"{os.getuid()}:{os.getgid()}",
            "-v", f"{scratch_dir}:/bam:ro",
            "-v", f"{snp_vcf.parent}:/vcf:ro",
            "-v", f"{out_dir}:/out",
        ]

        run([
            "docker", "run", "--rm",
            *docker_mounts,
            FACETS_SUITE_IMAGE,
            "snp-pileup-wrapper.R",
            "--vcf-file", f"/vcf/{snp_vcf.name}",
            "--normal-bam", "/bam/normal.bam",
            "--tumor-bam", "/bam/tumor.bam",
            "--output-prefix", f"/out/{patient_id}",
        ])

    finally:
        # Node-local scratch, dies with the VM anyway -- cleaned up early so
        # a re-run of this script (retry) doesn't hit snp-pileup-wrapper.R's
        # own "output file already exists" guard from a half-finished attempt.
        shutil.rmtree(scratch_dir, ignore_errors=True)

    counts_file = out_dir / f"{patient_id}.snp_pileup.gz"

    facets_cmd = [
        "docker", "run", "--rm",
        "-u", f"{os.getuid()}:{os.getgid()}",
        "-v", f"{out_dir}:/out",
        FACETS_SUITE_IMAGE,
        "run-facets-wrapper.R",
        "--counts-file", f"/out/{counts_file.name}",
        "--sample-id", patient_id,
        "--directory", "/out",
        "--genome", "hg38",
        # NOT --everything: confirmed live (2026-09-17) that its NTAI/arm-
        # level/QC extras crash unconditionally on this deployed image for a
        # WES sample -- calculate_ntai() hits "missing value where TRUE/FALSE
        # needed" on this sample's sparse per-chromosome segmentation, after
        # the core purity/ploidy fit has already succeeded. The actual
        # purity/ploidy summary (print_run_details() -> <sample>.txt) never
        # depended on --everything in the first place -- only the additional
        # gene/arm-level/QC files did.
    ]

    if args.cval is not None:
        facets_cmd += ["--cval", str(args.cval)]

    if args.purity_cval is not None:
        facets_cmd += ["--purity-cval", str(args.purity_cval)]

    if args.seed is not None:
        facets_cmd += ["--seed", str(args.seed)]

    run(facets_cmd)


if __name__ == "__main__":
    main()
