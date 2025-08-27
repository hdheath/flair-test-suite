from pydantic import BaseModel, Field, field_validator
from typing import List, Optional, Dict, Any
from pathlib import Path

def _split_commas_preserve_escapes(s: str) -> List[str]:
    """
    Split on commas, allowing escaped commas '\,' inside paths.
    Examples:
      "a.fq,b.fq"          -> ["a.fq", "b.fq"]
      "dir\\,with\\,comma" -> ["dir,with,comma"]
    """
    parts, buf, esc = [], [], False
    for ch in s:
        if esc:
            buf.append(ch)         # take char literally after backslash
            esc = False
        elif ch == '\\':
            esc = True
        elif ch == ',':
            parts.append(''.join(buf))
            buf = []
        else:
            buf.append(ch)
    parts.append(''.join(buf))
    return parts

class StageConfig(BaseModel):
    name: str
    requires: List[str] = Field(default_factory=list)
    flags: Dict[str, Any] = Field(default_factory=dict)

class RunConfig(BaseModel):
    version: str
    conda_env: str
    work_dir: str
    data_dir: str

    # Accept str or list[str]; normalize to list[str]
    reads_file: Optional[List[str]] = None

    @field_validator("reads_file", mode="before")
    def _coerce_reads_file(cls, v):
        if v is None:
            return None

        def normalize_one(item: Any) -> List[str]:
            # Accept Path or str; split by commas (with \, escape), strip spaces
            if isinstance(item, Path):
                return [str(item)]
            if isinstance(item, str):
                # Handle comma-separated single string or multiline TOML strings
                raw_parts = _split_commas_preserve_escapes(item)
                return [p.strip() for p in raw_parts if p.strip()]
            raise TypeError("reads_file entries must be strings or paths")

        out: List[str] = []
        if isinstance(v, (str, Path)):
            out.extend(normalize_one(v))
        elif isinstance(v, list):
            for item in v:
                out.extend(normalize_one(item))
        else:
            raise TypeError("reads_file must be a string, list of strings, or paths")

        # De-duplicate while preserving order
        seen = set()
        deduped: List[str] = []
        for p in out:
            if p not in seen:
                seen.add(p)
                deduped.append(p)

        if not deduped:
            raise ValueError("reads_file resolved to an empty list")
        return deduped

    genome_fa: Optional[str] = None
    # Shared optional inputs available to multiple stages
    gtf: Optional[str] = None
    regions_tsv: Optional[str] = None
    junctions: Optional[str] = None
    experiment_5_prime_regions_bed_file: Optional[str] = None
    experiment_3_prime_regions_bed_file: Optional[str] = None
    reference_5_prime_regions_bed_file: Optional[str] = None
    reference_3_prime_regions_bed_file: Optional[str] = None
    stages: List[StageConfig]

class Config(BaseModel):
    run_id: Optional[str] = None
    case_name: Optional[str] = None
    run: RunConfig
    qc: Dict[str, Any] = Field(default_factory=dict)


class SQANTIQCConfig(BaseModel):
    """Configuration block for optional SQANTI QC."""
    conda_env: str
    cpus: int = 1
