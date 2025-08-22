from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

from flair_test_suite.stages.stage_utils import (
    collect_upstream_pairs,
    isoform_expected_outputs,
    build_flair_cmds,
)
from flair_test_suite.lib.paths import PathBuilder


def test_collect_upstream_pairs_standard(tmp_path):
    run_id = "run1"
    align_pb = PathBuilder(tmp_path, run_id, "align", "sigA")
    align_pb.stage_dir.mkdir(parents=True, exist_ok=True)
    bed = align_pb.stage_dir / f"{run_id}_flair.bed"
    bed.write_text("bed\n")

    pairs, sigs, mode = collect_upstream_pairs(
    "correct", {"align": align_pb}, run_id, "flair.bed", "{chrom}_{start}_{end}.bed"
    )
    assert pairs == [(bed, run_id)]
    assert sigs == [align_pb.signature]
    assert mode == "standard"


def test_collect_upstream_pairs_regionalized(tmp_path):
    run_id = "run1"
    reg_pb = PathBuilder(tmp_path, run_id, "regionalize", "sigR")
    corr_pb = PathBuilder(tmp_path, run_id, "correct", "sigC")
    reg_pb.stage_dir.mkdir(parents=True, exist_ok=True)
    corr_pb.stage_dir.mkdir(parents=True, exist_ok=True)
    # region details
    (reg_pb.stage_dir / "region_details.tsv").write_text("chrom\tstart\tend\nchr1\t0\t10\n")
    bed = corr_pb.stage_dir / "chr1_0_10_all_corrected.bed"
    bed.write_text("bed\n")

    pairs, sigs, mode = collect_upstream_pairs(
        "collapse",
        {"regionalize": reg_pb, "correct": corr_pb},
        run_id,
        "all_corrected.bed",
        "{chrom}_{start}_{end}_all_corrected.bed",
        file_check=lambda p: p.exists(),
    )
    assert pairs == [(bed, "chr1_0_10")]
    assert set(sigs) == {reg_pb.signature, corr_pb.signature}
    assert mode == "regionalized"


def test_isoform_expected_outputs():
    std = isoform_expected_outputs("run1", None, False)
    assert std["isoforms_bed"] == Path("run1.isoforms.bed")
    reg = isoform_expected_outputs("run1", "chr1_0_10", True)
    assert reg["isoforms_bed"] == Path("chr1_0_10.isoforms.bed")
    assert "isoforms_bed_pattern" in reg


def test_build_flair_cmds_basic(tmp_path):
    pairs = [(Path("in.bed"), "tag")]
    cmds = build_flair_cmds(
        "collapse",
        pairs,
        Path("gen.fa"),
        [Path("reads.fa")],
        "run1",
        ["--foo"],
        regionalized=False,
        use_bed=True,
    )
    assert cmds == [["flair", "collapse", "-g", "gen.fa", "-r", "reads.fa", "-q", "in.bed", "-o", "run1", "--foo"]]


def test_count_reads_handles_fastx(tmp_path):
    from flair_test_suite.stages.stage_utils import count_reads
    # FASTA
    fa = tmp_path / "reads.fa"
    fa.write_text(">r1\nAAA\n>r2\nTTT\n")
    assert count_reads(fa) == 2

    # FASTQ
    fq = tmp_path / "reads.fq"
    fq.write_text("@r1\nAAA\n+\nIII\n@r2\nTTT\n+\nIII\n")
    assert count_reads(fq) == 2

    # gzipped FASTA
    fa_gz = tmp_path / "reads.fa.gz"
    import gzip
    with gzip.open(fa_gz, "wt") as fh:
        fh.write(">r1\nAAA\n>r2\nTTT\n")
    assert count_reads(fa_gz) == 2

    # gzipped FASTQ
    fq_gz = tmp_path / "reads.fastq.gz"
    with gzip.open(fq_gz, "wt") as fh:
        fh.write("@r1\nAAA\n+\nIII\n@r2\nTTT\n+\nIII\n")
    assert count_reads(fq_gz) == 2
