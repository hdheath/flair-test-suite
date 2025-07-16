"""
src/flair_test_suite/print_config.py

Pretty-print the loaded TOML config, now including sample_name and output paths.
"""
from pathlib import Path

def print_config_summary(cfg):
    """
    Nicely print out:
    - sample_name and output base directory
    - top-level values
    - data_dir details
    - whole_sample flag
    - regions block
    - stages and their flags
    """
    # Header
    print("\n🔧  FLAIR Test Suite Run Configuration")
    print("=" * 60)

    # Sample name and output directory
    sample = getattr(cfg, 'sample_name', None)
    if sample:
        out_base = Path("outputs") / sample
        print(f"Sample Name  : {sample}")
        print(f"Output Base  : {out_base}")
        print("-" * 60)

    # Run ID, Version, Conda env
    if hasattr(cfg, "run_id"):
        print(f"Run ID       : {cfg.run_id}")
    if hasattr(cfg, "run"):
        print(f"Version      : {cfg.run.version}")
        print(f"Conda Env    : {cfg.run.conda_env}")

    # Input root
    print(f"Input Root   : {cfg.input_root}")

    # Data Dir breakdown
    dd = cfg.data_dir
    if isinstance(dd, str):
        print(f"Data Dir     : {dd}")
    else:
        print(f"Data Dir     : {getattr(dd, 'name', dd)}")
        for key, val in vars(dd).items():
            if key == 'name':
                continue
            if hasattr(val, '__dict__'):
                print(f"  {key:15}:")
                for subk, subv in vars(val).items():
                    print(f"      {subk:13} = {subv}")
            else:
                print(f"  {key:15} = {val}")

    # Whole sample mode
    print(f"Whole Sample : {cfg.whole_sample}")
    print("-" * 60)

    # Regions
    regions = getattr(cfg.run, 'regions', None)
    if regions:
        print("Regions to process:")
        for r in regions:
            chrom = getattr(r, 'chr', '<no chr>')
            spans = getattr(r, 'ranges', [])
            print(f"  • {chrom:8} {', '.join(spans)}")
        print("-" * 60)

    # Stages & Flags
    steps = getattr(cfg.run, 'steps', None)
    if steps:
        print("Enabled stages and flags:")
        for stage in ('align','correct','collapse','quantify'):
            stage_obj = getattr(steps, stage, None)
            if stage_obj is not None:
                flags = getattr(stage_obj, 'flags', None)
                print(f"\n▶ {stage.upper()}:')")
                if flags:
                    for k,v in vars(flags).items():
                        print(f"    • {k:15} = {v}")
                else:
                    print("    (no flags set)")
    print("\n" + "=" * 60 + "\n")
