from .align import AlignStage
from .correct import CorrectStage
from .region_test import RegionTestStage
from .collapse import CollapseStage
from .combine import CombineStage
from .transcriptome import TranscriptomeStage
from .quantify import QuantifyStage


STAGE_REGISTRY = {
    AlignStage.name: AlignStage,
    CorrectStage.name: CorrectStage,
    RegionTestStage.name: RegionTestStage,
    CollapseStage.name: CollapseStage,
    CombineStage.name: CombineStage,
    TranscriptomeStage.name: TranscriptomeStage,
    QuantifyStage.name: QuantifyStage,
}
