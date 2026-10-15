#!/usr/bin/env python3
"""
======================================================================
SUPPORTED YAML FORMAT
======================================================================

Example 1: one tumor FASTQ pair + one normal FASTQ pair

samples:

  tumor_dna_r1: /path/to/tumor_R1.fastq.gz
  tumor_dna_r2: /path/to/tumor_R2.fastq.gz

  normal_dna_r1: /path/to/normal_R1.fastq.gz
  normal_dna_r2: /path/to/normal_R2.fastq.gz

  rna_r1: /path/to/rna_R1.fastq.gz
  rna_r2: /path/to/rna_R2.fastq.gz


Example 2: multiple tumor/normal FASTQ lanes

samples:

  tumor_dna_r1:
    - /path/to/tumor_L1_R1.fastq.gz
    - /path/to/tumor_L2_R1.fastq.gz

  tumor_dna_r2:
    - /path/to/tumor_L1_R2.fastq.gz
    - /path/to/tumor_L2_R2.fastq.gz

  normal_dna_r1:
    - /path/to/normal_L1_R1.fastq.gz
    - /path/to/normal_L2_R1.fastq.gz

  normal_dna_r2:
    - /path/to/normal_L1_R2.fastq.gz
    - /path/to/normal_L2_R2.fastq.gz


Notes:
    - GATK.GRCh38 is always used through nf-core/sarek iGenomes.
======================================================================
"""

import argparse
import csv
import re
import subprocess
from pathlib import Path

import yaml

# pipeline settings

NEXTFLOW = Path("/home/kevin/bin/nextflow-24.04.4")

YAML_DIR = Path("/data3/Code/HG19_code")

OUTPUT_ROOT = Path("/data4/NewOutput")

# shared Sarek work and cache
# Different patients use the same Nextflow work/cache location so reusable
# reference preparation and cached processes can be reused when possible.

SHARED_SAREK_DIR = Path("/data4/SarekShared")

SHARED_WORKDIR = SHARED_SAREK_DIR / "work"

LOCAL_RESOURCES_CONFIG = SHARED_SAREK_DIR / "local_resources.config"

MAX_CPUS = 25

BG_INTERVALS_BED = Path(
    "/data4/cloud-bioinformatics_v1/BED/"
    "S33613271_Regions.cleaned.sorted.bed"
)

CEGAT_INTERVALS_BED = Path(
    "/data4/SarekReferences/CeGaT_bed/"
    "MANE_Select_v1.3_GRCh38_CDS.GATK_compatible.bed"
)

# ASCAT's exome-restricted reference files (built against the same
# S33613271 kit as INTERVALS_BED above) -- actual paths are set via
# `params.ascat_*` in local_resources.config, not passed on the CLI. Listed
# here only so main() can fail fast if they're missing, same as INTERVALS_BED.
ASCAT_WES_DIR = Path(
    "/data4/SarekReferences/ASCAT_WES_hg38"
)

ASCAT_WES_FILES = [
    ASCAT_WES_DIR / "battenberg_alleles_on_target_hg38.zip",
    ASCAT_WES_DIR / "battenberg_loci_on_target_hg38.zip",
    ASCAT_WES_DIR / "GC_G1000_on_target_hg38.zip",
    ASCAT_WES_DIR / "RT_G1000_on_target_hg38.zip",
]

SAREK_REVISION = "3.5.1"
GENOME = "GATK.GRCh38"
IGENOMES_BASE = Path(
    "/data4/SarekReferences/igenomes"
)
ALIGNER = "bwa-mem"
TOOLS = "mutect2,snpeff,ascat"

# downstream mutation annotation settings

DNA_SNV_DIR = Path(
    "/data4/local-bioinformatics/DNA_SNV"
)

FILTER_MUTECT2_SCRIPT = (
    DNA_SNV_DIR / "filter_mutect2_calls.sh"
)

ANNOTATE_COSMIC_SCRIPT = (
    DNA_SNV_DIR / "annotate_cosmic.sh"
)

ANNOTATE_ONCOKB_SCRIPT = (
    DNA_SNV_DIR / "annotate_oncokb.sh"
)

ANNOTATE_CLINVAR_SCRIPT = (
    DNA_SNV_DIR / "annotate_clinvar.sh"
)

GENERATE_ABERRATIONS_SCRIPT = (
    DNA_SNV_DIR / "generate_aberrations_xml.py"
)

COSMIC_DIR = (
    DNA_SNV_DIR / "COSMIC"
)

ONCOKB_DIR = (
    DNA_SNV_DIR / "OncoKB"
)

CLINVAR_SOMATIC_VCF = (
    DNA_SNV_DIR
    / "ClinVar"
    / "clinvar.38.somatic.vcf.gz"
)

MAPPING_ENV = "mapping_env"


def sanitize_patient_id(patient_name: str):
    """
    Convert a patient name into a safe Sarek identifier.
    Example:
        Joseph Katz -> Joseph_Katz
        The output directory itself still keeps the original patient name or ID with space.
    """
    patient_id = patient_name.strip()

    patient_id = re.sub(r"\s+", "_", patient_id)
    patient_id = re.sub(r"[^A-Za-z0-9_.-]", "_", patient_id)
    patient_id = re.sub(r"_+", "_", patient_id)
    patient_id = patient_id.strip("_")

    if not patient_id:
        raise SystemExit("Invalid patient name")

    return patient_id


def normalize_fastq_entry(value, field_name):
    """
    Allow either:

        tumor_dna_r1: /path/file.fastq.gz

    or:

        tumor_dna_r1:
          - /path/lane1.fastq.gz
          - /path/lane2.fastq.gz

    Always return a list of Path objects.
    """

    if value is None:
        raise SystemExit(
            f"Required YAML field missing: samples.{field_name}"
        )

    if isinstance(value, str):
        values = [value]

    elif isinstance(value, list):
        values = value

    else:
        raise SystemExit(
            f"samples.{field_name} must be either a string "
            f"or a YAML list"
        )

    if not values:
        raise SystemExit(
            f"samples.{field_name} is empty"
        )

    paths = []

    for value in values:
        if not isinstance(value, str):
            raise SystemExit(
                f"Every entry in samples.{field_name} must be a string"
            )

        path = Path(value).expanduser()

        paths.append(path)

    return paths


def validate_fastqs(
    tumor_r1,
    tumor_r2,
    normal_r1,
    normal_r2,
):
    """
    Do Validation QC
    Confirm R1/R2 counts match and every FASTQ exists.
    """

    if len(tumor_r1) != len(tumor_r2):
        raise SystemExit(
            "Tumor FASTQ pair count mismatch: "
            f"{len(tumor_r1)} R1 file(s), "
            f"{len(tumor_r2)} R2 file(s)"
        )

    if len(normal_r1) != len(normal_r2):
        raise SystemExit(
            "Normal FASTQ pair count mismatch: "
            f"{len(normal_r1)} R1 file(s), "
            f"{len(normal_r2)} R2 file(s)"
        )

    all_files = (
        tumor_r1
        + tumor_r2
        + normal_r1
        + normal_r2
    )

    missing = [
        path
        for path in all_files
        if not path.exists()
    ]

    if missing:
        missing_text = "\n".join(
            str(path)
            for path in missing
        )

        raise SystemExit(
            "FASTQ file(s) not found:\n"
            f"{missing_text}"
        )


def load_fastqs_from_yaml(yaml_path: Path):
    """
    Read only tumor/normal DNA FASTQs.
    RNA and references sections are intentionally ignored now for mutect2.
    """

    if not yaml_path.exists():
        raise SystemExit(
            f"YAML file not found: {yaml_path}"
        )

    with yaml_path.open() as f:
        config = yaml.safe_load(f)

    if not isinstance(config, dict):
        raise SystemExit(
            f"Invalid YAML structure: {yaml_path}"
        )

    samples = config.get("samples")

    if not isinstance(samples, dict):
        raise SystemExit(
            f"'samples:' section missing or invalid in {yaml_path}"
        )

    tumor_r1 = normalize_fastq_entry(
        samples.get("tumor_dna_r1"),
        "tumor_dna_r1",
    )

    tumor_r2 = normalize_fastq_entry(
        samples.get("tumor_dna_r2"),
        "tumor_dna_r2",
    )

    normal_r1 = normalize_fastq_entry(
        samples.get("normal_dna_r1"),
        "normal_dna_r1",
    )

    normal_r2 = normalize_fastq_entry(
        samples.get("normal_dna_r2"),
        "normal_dna_r2",
    )

    validate_fastqs(
        tumor_r1,
        tumor_r2,
        normal_r1,
        normal_r2,
    )

    return (
        tumor_r1,
        tumor_r2,
        normal_r1,
        normal_r2,
    )


def build_samplesheet(
    samplesheet_path: Path,
    patient_id: str,
    sex: str,
    tumor_r1,
    tumor_r2,
    normal_r1,
    normal_r2,
):

    rows = []

    normal_sample = f"{patient_id}_normal"
    tumor_sample = f"{patient_id}_tumor"

    for lane_index, (r1, r2) in enumerate(
        zip(normal_r1, normal_r2),
        start=1,
    ):
        rows.append({
            "patient": patient_id,
            "sex": sex,
            "status": "0",
            "sample": normal_sample,
            "lane": f"L{lane_index}",
            "fastq_1": str(r1.resolve()),
            "fastq_2": str(r2.resolve()),
        })

    for lane_index, (r1, r2) in enumerate(
        zip(tumor_r1, tumor_r2),
        start=1,
    ):
        rows.append({
            "patient": patient_id,
            "sex": sex,
            "status": "1",
            "sample": tumor_sample,
            "lane": f"L{lane_index}",
            "fastq_1": str(r1.resolve()),
            "fastq_2": str(r2.resolve()),
        })

    with samplesheet_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "patient",
                "sex",
                "status",
                "sample",
                "lane",
                "fastq_1",
                "fastq_2",
            ],
        )

        writer.writeheader()
        writer.writerows(rows)

    print(
        f"Wrote {len(rows)} rows to {samplesheet_path}"
    )

    return rows


# shared Sarek resources

def build_local_resources_config(
    ascat_ploidy,
    ascat_purity,
):
    """
    Create a shared Nextflow resource config.
    The complete Sarek run can use up to 25 CPU cores in total.
    Individual Sarek processes keep their original nf-core CPU settings.
    """

    SHARED_SAREK_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    SHARED_WORKDIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    ascat_override = ""

    if (
        ascat_ploidy is not None
        and ascat_purity is not None
    ):
        ascat_override = f"""
    ascat_ploidy = {ascat_ploidy}
    ascat_purity = {ascat_purity}"""

    config_text = f"""executor {{
    cpus = {MAX_CPUS}
}}

params {{
    ascat_genome = 'hg38'
    ascat_alleles = '{ASCAT_WES_FILES[0]}'
    ascat_loci = '{ASCAT_WES_FILES[1]}'
    ascat_loci_gc = '{ASCAT_WES_FILES[2]}'
    ascat_loci_rt = '{ASCAT_WES_FILES[3]}'{ascat_override}
}}
"""

    LOCAL_RESOURCES_CONFIG.write_text(
        config_text
    )

    print(
        f"Shared Sarek work directory: {SHARED_WORKDIR}"
    )

    print(
        f"Maximum CPU cores: {MAX_CPUS}"
    )


# Sarek command

def build_command(
    samplesheet: Path,
    outdir: Path,
    workdir: Path,
    intervals_bed: Path,
):
    """
        FASTQ
          -> BWA-MEM
          -> MarkDuplicates
          -> BQSR
          -> Mutect2
          -> FilterMutectCalls
          -> snpEff

    Reference resources are resolved through:
        --genome GATK.GRCh38
    """

    return [
        str(NEXTFLOW),
        "run",
        "nf-core/sarek",

        "-r",
        SAREK_REVISION,

        "-profile",
        "docker",

        "-c",
        str(LOCAL_RESOURCES_CONFIG),

        "-work-dir",
        str(workdir),

        "--input",
        str(samplesheet),

        "--outdir",
        str(outdir),

        "--genome",
        GENOME,

        "--igenomes_base",
        str(IGENOMES_BASE),

        "--wes",

        # BWA-MEM
        # BWA-MEM the Sarek 3.5.1 default, but it is written
        "--aligner",
        ALIGNER,

        "--intervals",
        str(intervals_bed),
        "--tools",
        TOOLS,

        "-resume",
    ]


def run_mutation_downstream(
    outdir: Path,
):
    """
        Mutect2 + snpEff and then
          -> filter_mutect2_calls.sh
          -> annotate_cosmic.sh
          -> annotate_oncokb.sh
          -> annotate_clinvar.sh
          -> generate_aberrations_xml.py
    """

    required_paths = [
        FILTER_MUTECT2_SCRIPT,
        ANNOTATE_COSMIC_SCRIPT,
        ANNOTATE_ONCOKB_SCRIPT,
        ANNOTATE_CLINVAR_SCRIPT,
        GENERATE_ABERRATIONS_SCRIPT,
        COSMIC_DIR,
        ONCOKB_DIR,
        CLINVAR_SOMATIC_VCF,
    ]

    missing = [
        path
        for path in required_paths
        if not path.exists()
    ]

    if missing:
        missing_text = "\n".join(
            str(path)
            for path in missing
        )

        raise SystemExit(
            "Required downstream resource(s) not found:\n"
            f"{missing_text}"
        )

    mutation_dir = (
        outdir
        / "mutation-filtering"
    )

    subprocess.run(
        [
            "conda",
            "run",
            "-n",
            MAPPING_ENV,
            "bash",
            str(FILTER_MUTECT2_SCRIPT),
            str(outdir),
        ],
        check=True,
    )

    subprocess.run(
        [
            "conda",
            "run",
            "-n",
            MAPPING_ENV,
            "bash",
            str(ANNOTATE_COSMIC_SCRIPT),
            str(outdir),
            str(COSMIC_DIR),
        ],
        check=True,
    )

    cosmic_files = list(
        mutation_dir.glob(
            "*.annotated_calls.tsv"
        )
    )

    if len(cosmic_files) != 1:
        raise RuntimeError(
            "Expected exactly one COSMIC annotated_calls.tsv "
            f"under {mutation_dir}, found {len(cosmic_files)}"
        )

    cosmic_tsv = cosmic_files[0]

    subprocess.run(
        [
            "bash",
            str(ANNOTATE_ONCOKB_SCRIPT),
            str(cosmic_tsv),
            str(ONCOKB_DIR),
        ],
        check=True,
    )

    oncokb_tsv = Path(
        str(cosmic_tsv).replace(
            ".tsv",
            ".oncokb.tsv",
        )
    )

    if not oncokb_tsv.exists():
        raise RuntimeError(
            f"OncoKB output not found: {oncokb_tsv}"
        )

    subprocess.run(
        [
            "conda",
            "run",
            "-n",
            MAPPING_ENV,
            "bash",
            str(ANNOTATE_CLINVAR_SCRIPT),
            str(oncokb_tsv),
            str(CLINVAR_SOMATIC_VCF),
        ],
        check=True,
    )

    clinvar_tsv = Path(
        str(oncokb_tsv).replace(
            ".tsv",
            ".clinvar.tsv",
        )
    )

    if not clinvar_tsv.exists():
        raise RuntimeError(
            f"ClinVar output not found: {clinvar_tsv}"
        )

    aberrations_xml = (
        mutation_dir
        / "aberrations.xml"
    )

    subprocess.run(
        [
            "python3",
            str(GENERATE_ABERRATIONS_SCRIPT),
            str(clinvar_tsv),
            str(ONCOKB_DIR),
            str(aberrations_xml),
        ],
        check=True,
    )

    all_mutations_xml = (
        mutation_dir
        / "aberrations_all_mutations.xml"
    )


# main
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="only write samplesheet + print command",
    )

    args = parser.parse_args()

    patient_name = input(
        "Patient name: "
    ).strip()

    yaml_filename = input(
        "YAML filename: "
    ).strip()

    # patient sex
    # 0 = female (XX)
    # 1 = male (XY)

    bed_input = input(
        "BED type (1 = BostonGene, 2 = CeGaT): "
    ).strip()

    if bed_input == "1":
        intervals_bed = BG_INTERVALS_BED
        bed_name = "BostonGene"

    elif bed_input == "2":
        intervals_bed = CEGAT_INTERVALS_BED
        bed_name = "CeGaT"

    else:
        raise SystemExit(
            "Invalid BED selection. Please enter 1 for BostonGene or 2 for CeGaT."
        )

    print(
        f"Selected BED: {bed_name}"
    )

    print(
        f"Intervals: {intervals_bed}"
    )

    sex_input = input(
        "Sex (0 = female, 1 = male): "
    ).strip()

    if sex_input == "0":
        sex = "XX"

    elif sex_input == "1":
        sex = "XY"

    else:
        raise SystemExit(
            "Invalid sex input. Please enter 0 for female or 1 for male."
        )

    ascat_ploidy_input = input(
        "ASCAT ploidy (press Enter for automatic): "
    ).strip()

    ascat_purity_input = input(
        "ASCAT purity (0-1, press Enter for automatic): "
    ).strip()

    if not ascat_ploidy_input and not ascat_purity_input:
        ascat_ploidy = None
        ascat_purity = None

        print(
            "ASCAT purity/ploidy: automatic estimation"
        )

    elif ascat_ploidy_input and ascat_purity_input:
        try:
            ascat_ploidy = float(
                ascat_ploidy_input
            )
        except ValueError:
            raise SystemExit(
                "Invalid ASCAT ploidy. Please enter a numeric value "
                "or leave both fields blank for automatic estimation."
            )

        try:
            ascat_purity = float(
                ascat_purity_input
            )
        except ValueError:
            raise SystemExit(
                "Invalid ASCAT purity. Please enter a numeric value between "
                "0 and 1, or leave both fields blank for automatic estimation."
            )

        if ascat_ploidy <= 0:
            raise SystemExit(
                "ASCAT ploidy must be greater than 0."
            )

        if not 0 < ascat_purity <= 1:
            raise SystemExit(
                "ASCAT purity must be greater than 0 and less than or equal to 1."
            )

    else:
        raise SystemExit(
            "Please either provide BOTH ASCAT ploidy and purity, "
            "or leave BOTH blank for automatic estimation."
        )

    if not patient_name:
        raise SystemExit(
            "Patient name cannot be empty"
        )

    if not yaml_filename:
        raise SystemExit(
            "YAML filename cannot be empty"
        )

    if "/" in patient_name or "\\" in patient_name:
        raise SystemExit(
            "Patient name may contain spaces, but cannot contain / or \\"
        )

    patient_id = sanitize_patient_id(
        patient_name
    )

    yaml_path = (
        YAML_DIR
        / yaml_filename
    )

    patient_output_dir = (
        OUTPUT_ROOT
        / patient_name
    )

    launchdir = (
        patient_output_dir
        / "nextflow-launch"
    )

    outdir = (
        patient_output_dir
        / "sarek-results"
    )

    workdir = (
        SHARED_WORKDIR
    )

    samplesheet = (
        patient_output_dir
        / f"samplesheet_{patient_id}.csv"
    )

    if not NEXTFLOW.exists():
        raise SystemExit(
            f"Nextflow not found: {NEXTFLOW}"
        )

    if not intervals_bed.exists():
        raise SystemExit(
            f"INTERVALS_BED not found: {intervals_bed}"
        )

    missing = [
        f
        for f in ASCAT_WES_FILES
        if not f.exists()
    ]

    if missing:
        raise SystemExit(
            f"ASCAT WES reference file(s) not found: {missing}"
        )

    (
        tumor_r1,
        tumor_r2,
        normal_r1,
        normal_r2,
    ) = load_fastqs_from_yaml(
        yaml_path
    )

    patient_output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    build_samplesheet(
        samplesheet_path=samplesheet,
        patient_id=patient_id,
        sex=sex,
        tumor_r1=tumor_r1,
        tumor_r2=tumor_r2,
        normal_r1=normal_r1,
        normal_r2=normal_r2,
    )

    # build shared Sarek work/cache and resource settings

    build_local_resources_config(
        ascat_ploidy=ascat_ploidy,
        ascat_purity=ascat_purity,
    )

    cmd = build_command(
        samplesheet=samplesheet,
        outdir=outdir,
        workdir=workdir,
        intervals_bed=intervals_bed,
    )

    if args.dry_run:
        return

    # Run Sarek

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    workdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    launchdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    subprocess.run(
        cmd,
        cwd=launchdir,
        check=True,
    )

    run_mutation_downstream(
        outdir=outdir,
    )


if __name__ == "__main__":
    main()