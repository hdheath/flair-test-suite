from .align import AlignStage
from .correct import CorrectStage
from .regionalize import RegionalizeStage
from .collapse import CollapseStage
from .combine import CombineStage
from .transcriptome import TranscriptomeStage


STAGE_REGISTRY = {
    AlignStage.name: AlignStage,
    CorrectStage.name: CorrectStage,
    RegionalizeStage.name: RegionalizeStage,
    CollapseStage.name: CollapseStage,
    CombineStage.name: CombineStage,
    TranscriptomeStage.name: TranscriptomeStage,
}
