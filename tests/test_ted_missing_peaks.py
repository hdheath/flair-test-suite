import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1] / 'src'))

from flair_test_suite.qc.ted import _tss_tts_metrics_full


def test_tss_tts_metrics_full_audits_missing(tmp_path):
    iso = tmp_path / 'iso.bed'
    iso.write_text('chr1\t0\t1\tiso1\t0\t+\n')
    audit = []
    peaks = {'prime5': None, 'prime3': None, 'ref_prime5': None, 'ref_prime3': None}
    metrics = _tss_tts_metrics_full(iso, peaks, 50, audit, 'ctx')
    assert metrics['5prime_precision'] is None
    assert len(audit) == 4
    notes = {row['note'] for row in audit}
    assert notes == {'missing'}
    tags = {row['tag'] for row in audit}
    assert tags == {'5', '3', 'ref5', 'ref3'}
