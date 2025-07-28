from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any

class StageConfig(BaseModel):
    name: str
    requires: List[str] = Field(default_factory=list)  # <-- always a list
    flags: Optional[Dict[str, Any]] = Field(default_factory=dict)

class RunConfig(BaseModel):
    version: str
    conda_env: str
    work_dir: str
    input_root: str
    data_dir: str
    reads_file: Optional[str] = None
    genome_fa: Optional[str] = None
    stages: List[StageConfig]

class Config(BaseModel):
    run_id: Optional[str] = None
    run: RunConfig