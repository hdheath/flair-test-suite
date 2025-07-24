from .align import AlignStage
from .correct import CorrectStage
from .slice import SliceStage
from .collapse import CollapseStage

STAGE_REGISTRY = {
    AlignStage.name: AlignStage,
    CorrectStage.name: CorrectStage,
    SliceStage.name: SliceStage,
    CollapseStage.name: CollapseStage,
}
