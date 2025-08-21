import types
import sys


# The transcriptome browser module depends on heavy optional libraries such as
# pandas and matplotlib. For the purpose of testing the CLI helper
# ``_parse_region`` we only need to load the module, so we provide lightweight
# stubs for those imports before importing the function.
for name in [
    "pandas",
    "matplotlib",
    "matplotlib.pyplot",
    "numpy",
    "scipy",
]:
    sys.modules.setdefault(name, types.ModuleType(name))

# Provide minimal attributes required during import.
pysam_mod = types.ModuleType("pysam")
class DummyAligned:  # pragma: no cover
    pass
pysam_mod.AlignedSegment = DummyAligned
sys.modules.setdefault("pysam", pysam_mod)
patches = types.ModuleType("matplotlib.patches")
patches.Rectangle = object
sys.modules.setdefault("matplotlib.patches", patches)

collections = types.ModuleType("matplotlib.collections")
collections.PatchCollection = object
collections.LineCollection = object
sys.modules.setdefault("matplotlib.collections", collections)

stats = types.ModuleType("scipy.stats")
stats.gaussian_kde = object
sys.modules.setdefault("scipy.stats", stats)

intervaltree = types.ModuleType("intervaltree")
class DummyTree:  # pragma: no cover - behaviour irrelevant
    pass
intervaltree.IntervalTree = DummyTree
sys.modules.setdefault("intervaltree", intervaltree)

from flair_test_suite.plotting.transcriptome_browser import _parse_region


def test_parse_region_basic():
    assert _parse_region("chr1:100-200") == (("chr1", 100, 200), False)


def test_parse_region_with_dash_in_chrom():
    assert _parse_region("chr1-2_random:100-200") == (("chr1-2_random", 100, 200), False)


def test_parse_region_too_long():
    assert _parse_region("chr1:100-25000") == (("chr1", 100, 25000), True)
