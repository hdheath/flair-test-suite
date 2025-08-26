import sys
from pathlib import Path
import types

import pandas as pd
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1] / 'src'))

from flair_test_suite.qc import sqanti
from flair_test_suite.config_schema import Config as Cfg


def _make_cfg(tmp_path: Path):
    data_dir = tmp_path / 'data'
    data_dir.mkdir()
    for name in ['ref.gtf', 'ref.fa', 'ref5.bed', 'ref3.bed', 'sj.bed']:
        (data_dir / name).write_text('')
    cfg_dict = {
        'run': {
            'version': '3.0',
            'conda_env': 'flair',
            'work_dir': str(tmp_path / 'out'),
            'data_dir': str(data_dir),
            'gtf': 'ref.gtf',
            'genome_fa': 'ref.fa',
            'reference_5_prime_regions_bed_file': 'ref5.bed',
            'reference_3_prime_regions_bed_file': 'ref3.bed',
            'junctions': 'sj.bed',
            'stages': [{"name": "collapse"}],
        },
        'qc': {'collapse': {'SQANTI': {'conda_env': 'sqenv', 'cpus': 1}}},
    }
    return Cfg.model_validate(cfg_dict)


def _write_classification(path: Path):
    header = [
        'isoform','chrom','strand','length','exons','structural_category',
        'associated_gene','associated_transcript','ref_length','ref_exons',
        'diff_to_TSS','diff_to_TTS','diff_to_gene_TSS','diff_to_gene_TTS',
        'subcategory','RTS_stage','all_canonical','min_sample_cov','min_cov',
        'min_cov_pos','sd_cov','FL','n_indels','n_indels_junc','bite',
        'iso_exp','gene_exp','ratio_exp','FSM_class','coding','ORF_length',
        'CDS_length','CDS_start','CDS_end','CDS_genomic_start','CDS_genomic_end',
        'predicted_NMD','perc_A_downstream_TTS','seq_A_downstream_TTS',
        'dist_to_cage_peak','within_cage_peak','pos_cage_peak',
        'dist_to_polya_site','within_polya_site','polyA_motif','polyA_dist',
        'ORF_seq','TSS_genomic_coord','TTS_genomic_coord','experiment_id',
        'entry_id','LRGASP_id'
    ]
    row = [
        'iso1','chr1','+','1000','3','full-splice_match',
        'Gene1','Tx1','1000','3',
        '0','0','0','0',
        'sub','NA','True','1','1',
        '1','0','True','0','0','NA',
        '1','1','1','FSM','True','100',
        '100','1','100','1000','1100',
        'False','0','AAAA',
        '10','True','1',
        '20','True','AATAAA','0',
        'ORF','0','0','exp','id','LRG'
    ]
    with open(path, 'w') as fh:
        fh.write('\t'.join(header) + '\n')
        fh.write('\t'.join(row) + '\n')


def test_sqanti_skips_missing_env(tmp_path: Path, monkeypatch):
    cfg = _make_cfg(tmp_path)
    stage_dir = tmp_path / 'out' / 'run1' / 'collapse' / 'sig'
    stage_dir.mkdir(parents=True)
    (stage_dir / 'sample.isoforms.gtf').write_text('')
    monkeypatch.setattr(sqanti, '_conda_env_exists', lambda env: False)
    out_dir = stage_dir / 'qc' / 'sqanti'
    sqanti.collect(stage_dir, cfg, out_dir=out_dir)
    assert not (out_dir / 'sqanti_results.tsv').exists()


def test_sqanti_summarizes_existing_classification(tmp_path: Path, monkeypatch):
    cfg = _make_cfg(tmp_path)
    stage_dir = tmp_path / 'out' / 'run1' / 'collapse' / 'sig'
    stage_dir.mkdir(parents=True)
    gtf = stage_dir / 'sample.isoforms.gtf'
    gtf.write_text('')
    out_dir = stage_dir / 'qc' / 'sqanti'
    out_dir.mkdir(parents=True)
    class_txt = out_dir / 'sample_classification.txt'
    _write_classification(class_txt)
    monkeypatch.setattr(sqanti, '_conda_env_exists', lambda env: True)
    monkeypatch.setattr(sqanti, '_run_sqanti', lambda *a, **k: (_ for _ in ()).throw(AssertionError('should not run')))
    monkeypatch.setattr(sqanti, 'plot_summary', lambda *a, **k: None)
    sqanti.collect(stage_dir, cfg, out_dir=out_dir)
    tsv = out_dir / 'sqanti_results.tsv'
    assert tsv.exists()
    df = pd.read_csv(tsv, sep='\t')
    assert df.loc[0, 'sample'] == 'sample'
    assert df.loc[0, 'FSM'] == 1
