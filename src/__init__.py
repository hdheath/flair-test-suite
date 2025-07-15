"""
flair_automate: per-region FLAIR align → correct → slice → collapse → (opt) SQANTI QC & plotting
"""

from .flair_align_automation    import align_all
from .flair_correct_automation  import correct_all
from .slicer                    import slice_all_regions
from .flair_collapse_automation          import collapse_all_regions
from .parse_region_metrics_from_gtf import parse_all_regions
'''
from .sqanti_runner             import run_one, summarize_one
from .sqanti_plot               import plot_summary
'''

__all__ = [
    "align_all",
    "correct_all",
    "slice_all_regions",
    "collapse_all_regions",
    "parse_all_regions",
    "run_one",
    "summarize_one",
    "plot_summary",
]


