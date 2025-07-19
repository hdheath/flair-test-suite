# src/flair_test_suite/util/dag.py


def topological_sort(stage_cfg_list):
    """Return stages sorted by their 'requires' fields."""
    graph = {s.name: set(getattr(s, "requires", [])) for s in stage_cfg_list}
    ordered = []
    while graph:
        ready = [k for k, v in graph.items() if not v]
        if not ready:
            raise ValueError("Stage dependency cycle detected")
        ordered.extend(ready)
        for r in ready:
            graph.pop(r)
            for deps in graph.values():
                deps.discard(r)
    name2cfg = {s.name: s for s in stage_cfg_list}
    return [name2cfg[n] for n in ordered]
