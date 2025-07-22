from .align import AlignStage
from .correct import CorrectStage
from .slice import SliceStage

STAGE_REGISTRY = {
    AlignStage.name: AlignStage,
    CorrectStage.name: CorrectStage,
    SliceStage.name: SliceStage,
}
