from .align import AlignStage
STAGE_REGISTRY = {AlignStage.name: AlignStage}

from .correct import CorrectStage
STAGE_REGISTRY["correct"] = CorrectStage

