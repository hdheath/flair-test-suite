import os

THIS_DIR = os.path.dirname(__file__)
DATA_ROOT = os.path.join(THIS_DIR, "..", "tests", "data")

CONFIG = {
    "data_dirs": {
        "human_GENCODEv48_WTC11_PacBio_CAGE_QuantSeq": {
            "reads_fasta": "WTC11.PacBio.3reps.fasta",
            "genome_fasta": "GRCh38.primary_assembly.genome.fa",
            "annotation_gtf": "gencode.v48.annotation.gtf",
            "junctions": "WTC11_all.SJ.out.tab",
            "ref_tss_bed": "human_ref_TSS.bed",
            "ref_tts_bed": "human_ref_TTS.bed",
            "experimental_tss_bed": "CAGE_TSS_human.bed",
            "experimental_tts_bed": "WTC11_all_polyApeaks_fixed.bed",
        },
    },
    "runs": [
        {
            "version": "2.0.0",
            "conda_env": "flair-v2",
            "regions": [
                {"chr20": ["3218000:3250000"]}
            ],
            "steps": {
                "align": {
                    "options": [
                        {
                            "name": "default",
                            "flags": {"--threads": 12}
                        }
                    ]
                },
                "correct": {
                    "options": [
                        {
                            "name": "default",
                            "flags": {
                                "--threads": 12,
                                "--junc": os.path.join(
                                    DATA_ROOT,
                                    "human_GENCODEv48_WTC11_PacBio_CAGE_QuantSeq",
                                    "WTC11_all.SJ.out.tab"
                                )
                            }
                        },
                        {
                            "name": "with_ss_window",
                            "flags": {
                                "--threads": 12,
                                "--junc": os.path.join(
                                    DATA_ROOT,
                                    "human_GENCODEv48_WTC11_PacBio_CAGE_QuantSeq",
                                    "WTC11_all.SJ.out.tab"
                                ),
                                "--ss_window": 5
                            }
                        }
                    ]
                },
                "collapse": {
                    "options": [
                        {
                            "name": "default",
                            "flags": {
                                "--support": 1,
                                "--check_splice": True,
                                "--filter": "stringent",
                                "--threads": 8
                            }
                        }
                    ]
                },
                "quantify": {
                    "options": [
                        {
                            "name": "default",
                            "flags": {
                                "--support": 1,
                                "--check_splice": True,
                                "--filter": "stringent",
                                "--threads": 8
                            }
                        }
                    ]
                }
            }
        },
        {
            "version": "3.0.0",
            "conda_env": "flair-v3",
            "regions": [
                {"chr20": ["3218000:3250000"]}
            ],
            "steps": {
                "align": {
                    "options": [
                        {
                            "name": "default",
                            "flags": {"--threads": 12}
                        },
                        {
                            "name": "with_nvrna",
                            "flags": {"--threads": 12, "--nvrna": True}
                        }
                    ]
                },
                "correct": {
                    "options": [
                        {
                            "name": "default",
                            "flags": {
                                "--threads": 12,
                                "--junc": os.path.join(
                                    DATA_ROOT,
                                    "human_GENCODEv48_WTC11_PacBio_CAGE_QuantSeq",
                                    "WTC11_all.SJ.out.tab"
                                )
                            }
                        },
                        {
                            "name": "ss_and_window",
                            "flags": {
                                "--threads": 12,
                                "--junc": os.path.join(
                                    DATA_ROOT,
                                    "human_GENCODEv48_WTC11_PacBio_CAGE_QuantSeq",
                                    "WTC11_all.SJ.out.tab"
                                ),
                                "--ss_window": 5
                            }
                        }
                    ]
                },
                "collapse": {
                    "options": [
                        {
                            "name": "default",
                            "flags": {
                                "--support": 3,
                                "--filter": "best-only",
                                "--threads": 8
                            }
                        }
                    ]
                },
                "quantify": {
                    "options": [
                        {
                            "name": "default",
                            "flags": {
                                "--support": 3,
                                "--filter": "best-only",
                                "--threads": 8
                            }
                        }
                    ]
                }
            }
        }
    ]
}

