import sys
from pathlib import Path
import types

import pandas as pd
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1] / 'src'))

sys.modules.setdefault('intervaltree', types.SimpleNamespace(IntervalTree=object))
sys.modules.setdefault('pysam', types.SimpleNamespace(AlignedSegment=object))

from flair_test_suite.qc import ted
from flair_test_suite.lib import PathBuilder
from flair_test_suite.config_schema import Config as Cfg


def test_collect_regionalized_resolves_peaks(tmp_path: Path, monkeypatch):
    repo_root = tmp_path
    data_dir = repo_root / 'test_data'
    data_dir.mkdir()
    for name in ['exp5.bed', 'exp3.bed', 'ref5.bed', 'ref3.bed']:
        (data_dir / name).write_text('chr1\t0\t1\t.\t0\t+\n')

    run_id = 'run1'
    stage_uuid = '123abc'
    stage_dir = repo_root / 'outputs' / run_id / 'collapse' / stage_uuid
    stage_dir.mkdir(parents=True)

    tag = 'chr1_0_100'
    iso_bed = stage_dir / f'{tag}.isoforms.bed'
    iso_bed.write_text('chr1\t10\t50\tread_ENSG000001.1\t0\t+\n')
    map_txt = stage_dir / f'{tag}.isoform.read.map.txt'
    map_txt.write_text(f'{iso_bed.stem}\tread1\n')

    corr_dir = repo_root / 'outputs' / run_id / 'correct' / 'abcd'
    corr_dir.mkdir(parents=True)
    (corr_dir / f'{tag}_all_corrected.bed').write_text('chr1\t10\t50\tread_ENSG000001.1\t0\t+\n')

    reg_dir = repo_root / 'outputs' / run_id / 'regionalize' / 'reg1'
    reg_dir.mkdir(parents=True)
    (reg_dir / 'region_metrics.tsv').write_text('region_tag\tgene_count\ttranscript_count\nchr1_0_100\t1\t1\n')
    for base in ['exp5.bed', 'exp3.bed', 'ref5.bed', 'ref3.bed']:
        (reg_dir / f'{tag}_{base}').write_text('chr1\t10\t20\t.\t0\t+\n')

    cfg = {
        'run': {'data_dir': 'test_data'},
        'qc': {
            'collapse': {
                'TED': {
                    'experiment_5_prime_regions_bed_file': 'exp5.bed',
                    'experiment_3_prime_regions_bed_file': 'exp3.bed',
                    'reference_5_prime_regions_bed_file': 'ref5.bed',
                    'reference_3_prime_regions_bed_file': 'ref3.bed',
                    'window': 10,
                }
            }
        },
    }

    captured = {}

    def fake_metrics(iso_bed_path, peaks, window, audit_rows, tag_ctx):
        captured.update(peaks)
        return {
            '5prime_precision': 0.1,
            '5prime_recall': 0.1,
            '5prime_f1': 0.1,
            '3prime_precision': 0.2,
            '3prime_recall': 0.2,
            '3prime_f1': 0.2,
            'ref5prime_precision': 0.3,
            'ref5prime_recall': 0.3,
            'ref5prime_f1': 0.3,
            'ref3prime_precision': 0.4,
            'ref3prime_recall': 0.4,
            'ref3prime_f1': 0.4,
        }

    monkeypatch.setattr(ted, '_tss_tts_metrics_full', fake_metrics)

    reg_pb = PathBuilder(repo_root / 'outputs', run_id, 'regionalize', 'reg1')
    ted.collect(stage_dir, cfg, upstreams={'regionalize': reg_pb})

    assert captured['prime5'] == reg_dir / f'{tag}_exp5.bed'
    assert captured['prime3'] == reg_dir / f'{tag}_exp3.bed'
    assert captured['ref_prime5'] == reg_dir / f'{tag}_ref5.bed'
    assert captured['ref_prime3'] == reg_dir / f'{tag}_ref3.bed'

    df = pd.read_csv(stage_dir / 'TED.tsv', sep='\t')
    assert df.loc[0, '5prime_precision'] == 0.1
    assert df.loc[0, '3prime_precision'] == 0.2
    assert df.loc[0, 'ref5prime_precision'] == 0.3
    assert df.loc[0, 'ref3prime_precision'] == 0.4


def test_collect_regionalized_uses_run_level_inputs(tmp_path: Path, monkeypatch):
    """Run-level fields provide peaks when QC config omits them."""
    repo_root = tmp_path
    data_dir = repo_root / 'test_data'
    data_dir.mkdir()
    for name in ['exp5.bed', 'exp3.bed', 'ref5.bed', 'ref3.bed']:
        (data_dir / name).write_text('chr1\t0\t1\t.\t0\t+\n')

    run_id = 'run1'
    stage_uuid = '123abc'
    stage_dir = repo_root / 'outputs' / run_id / 'collapse' / stage_uuid
    stage_dir.mkdir(parents=True)

    tag = 'chr1_0_100'
    iso_bed = stage_dir / f'{tag}.isoforms.bed'
    iso_bed.write_text('chr1\t10\t50\tread_ENSG000001.1\t0\t+\n')

    reg_dir = repo_root / 'outputs' / run_id / 'regionalize' / 'reg1'
    reg_dir.mkdir(parents=True)
    (reg_dir / 'region_metrics.tsv').write_text('region_tag\tgene_count\ttranscript_count\nchr1_0_100\t1\t1\n')
    for base in ['exp5.bed', 'exp3.bed', 'ref5.bed', 'ref3.bed']:
        (reg_dir / f'{tag}_{base}').write_text('chr1\t10\t20\t.\t0\t+\n')

    cfg_dict = {
        'run': {
            'version': '1',
            'conda_env': 'env',
            'work_dir': 'outputs',
            'data_dir': 'test_data',
            'experiment_5_prime_regions_bed_file': 'exp5.bed',
            'experiment_3_prime_regions_bed_file': 'exp3.bed',
            'reference_5_prime_regions_bed_file': 'ref5.bed',
            'reference_3_prime_regions_bed_file': 'ref3.bed',
            'stages': [{"name": "collapse", "requires": []}],
        },
        'qc': {'collapse': {'TED': {'window': 10}}},
    }
    cfg = Cfg.model_validate(cfg_dict)

    captured = {}

    def fake_metrics(iso_bed_path, peaks, window, audit_rows, tag_ctx):
        captured.update(peaks)
        return {}

    monkeypatch.setattr(ted, '_tss_tts_metrics_full', fake_metrics)

    reg_pb = PathBuilder(repo_root / 'outputs', run_id, 'regionalize', 'reg1')
    ted.collect(stage_dir, cfg, upstreams={'regionalize': reg_pb})

    assert captured['prime5'] == reg_dir / f'{tag}_exp5.bed'
    assert captured['prime3'] == reg_dir / f'{tag}_exp3.bed'
    assert captured['ref_prime5'] == reg_dir / f'{tag}_ref5.bed'
    assert captured['ref_prime3'] == reg_dir / f'{tag}_ref3.bed'


def test_collect_regionalized_uses_stage_flags(tmp_path: Path, monkeypatch):
    """When QC config omits peaks, run.stage flags supply them."""
    repo_root = tmp_path
    data_dir = repo_root / 'test_data'
    data_dir.mkdir()
    for name in ['exp5.bed', 'exp3.bed', 'ref5.bed', 'ref3.bed']:
        (data_dir / name).write_text('chr1\t0\t1\t.\t0\t+\n')

    run_id = 'run1'
    stage_uuid = '123abc'
    stage_dir = repo_root / 'outputs' / run_id / 'collapse' / stage_uuid
    stage_dir.mkdir(parents=True)

    tag = 'chr1_0_100'
    iso_bed = stage_dir / f'{tag}.isoforms.bed'
    iso_bed.write_text('chr1\t10\t50\tread_ENSG000001.1\t0\t+\n')

    reg_dir = repo_root / 'outputs' / run_id / 'regionalize' / 'reg1'
    reg_dir.mkdir(parents=True)
    (reg_dir / 'region_metrics.tsv').write_text('region_tag\tgene_count\ttranscript_count\nchr1_0_100\t1\t1\n')
    for base in ['exp5.bed', 'exp3.bed', 'ref5.bed', 'ref3.bed']:
        (reg_dir / f'{tag}_{base}').write_text('chr1\t10\t20\t.\t0\t+\n')

    cfg = {
        'run': {
            'data_dir': 'test_data',
            'stages': {
                'collapse': {
                    'flags': {
                        'experiment_5_prime_regions_bed_file': 'exp5.bed',
                        'experiment_3_prime_regions_bed_file': 'exp3.bed',
                        'reference_5_prime_regions_bed_file': 'ref5.bed',
                        'reference_3_prime_regions_bed_file': 'ref3.bed',
                    }
                }
            },
        },
        'qc': {'collapse': {'TED': {'window': 10}}},
    }

    captured = {}

    def fake_metrics(iso_bed_path, peaks, window, audit_rows, tag_ctx):
        captured.update(peaks)
        return {}

    monkeypatch.setattr(ted, '_tss_tts_metrics_full', fake_metrics)

    reg_pb = PathBuilder(repo_root / 'outputs', run_id, 'regionalize', 'reg1')
    ted.collect(stage_dir, cfg, upstreams={'regionalize': reg_pb})

    assert captured['prime5'] == reg_dir / f'{tag}_exp5.bed'
    assert captured['prime3'] == reg_dir / f'{tag}_exp3.bed'
    assert captured['ref_prime5'] == reg_dir / f'{tag}_ref5.bed'
    assert captured['ref_prime3'] == reg_dir / f'{tag}_ref3.bed'
