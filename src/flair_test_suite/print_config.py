def print_config_summary(cfg):
    """
    Nicely print out the top‐level values, the regions block,
    and which stages have flags configured.
    """
    print("\n🔧  FLAIR Test Suite Configuration")
    print("=" * 50)

    # Top‐level
    if hasattr(cfg, "run_id"):
        print(f"Run ID       : {cfg.run_id}")
    if hasattr(cfg, "run"):
        print(f"Version      : {cfg.run.version}")
        print(f"Conda Env    : {cfg.run.conda_env}")
    print(f"Input Root   : {cfg.input_root}")
    print(f"Data Dir     : {cfg.data_dir}")
    print(f"Whole Sample : {cfg.whole_sample}")
    print("-" * 50)

    # Regions
    regions = getattr(cfg.run, "regions", None)
    if regions:
        print("Regions to process:")
        for r in regions:
            chr_    = getattr(r, "chr", "<no chr>")
            ranges  = getattr(r, "ranges", [])
            print(f"  • {chr_:8} {', '.join(ranges)}")
        print("-" * 50)

    # Stages & flags
    steps = getattr(cfg.run, "steps", None)
    if steps:
        print("Enabled stages and flags:")
        for stage in ("align", "correct", "collapse", "quantify"):
            stage_obj = getattr(steps, stage, None)
            if stage_obj is not None:
                flags = getattr(stage_obj, "flags", None)
                print(f"\n▶ {stage.upper()}:")
                if flags:
                    for k, v in vars(flags).items():
                        print(f"    • {k:15} = {v}")
                else:
                    print("    (no flags set)")
    print("\n" + "=" * 50 + "\n")
