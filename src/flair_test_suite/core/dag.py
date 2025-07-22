# src/flair_test_suite/util/dag.py
# ---------------------------------
# Utility for performing a topological sort on a list of stage configurations,
# ensuring that each stage runs after its dependencies (the 'requires' field).

from typing import List


def topological_sort(stage_cfg_list: List) -> List:
    """
    Given a list of stage configuration objects, each with:
      - .name: unique identifier for the stage
      - .requires: iterable of names of stages that must run before it
    Return a new list sorted so that dependencies come before dependents.

    Raises:
      - ValueError if there is a cycle in the dependency graph.
    """
    # 1) Build a dependency graph:
    #    { stage_name: set(of required stage names) }
    graph = {
        s.name: set(getattr(s, "requires", []))
        for s in stage_cfg_list
    }

    ordered: List[str] = []  # list of stage names in sorted order

    # 2) Kahn's algorithm: repeatedly pick stages with no remaining requirements
    while graph:
        # Find all stages whose requirement set is empty => ready to run
        ready = [name for name, deps in graph.items() if not deps]
        if not ready:
            # No ready nodes found, but graph not empty => cycle exists
            raise ValueError("Stage dependency cycle detected")

        # Add all ready stages to the output order
        ordered.extend(ready)

        # Remove each ready stage from the graph, and from other stages' dep sets
        for r in ready:
            graph.pop(r)
            for deps in graph.values():
                deps.discard(r)

    # 3) Map sorted names back to the original config objects
    name2cfg = {s.name: s for s in stage_cfg_list}
    return [name2cfg[name] for name in ordered]

