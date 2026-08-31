"""Independent validation of the notebook's bipartite matching stage.

The production implementation uses Hopcroft--Karp followed by Konig cover
recovery.  This validator deliberately uses a simpler augmenting-path
implementation and, on small graphs, a brute-force minimum vertex-cover
oracle.  The separate code path reduces the risk that the validation merely
repeats an implementation defect in the production routine.
"""

from __future__ import annotations

import argparse
import csv
import json
from itertools import combinations
from pathlib import Path

import numpy as np


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
DEFAULT_OUTPUT = SCRIPT_DIRECTORY / "chapter5_validation_results"
DEFAULT_SEED = 20260822


def resolve_notebook(requested_path=None):
    """Resolve an explicitly supplied notebook or a standard neighbouring name."""
    if requested_path is not None:
        notebook = Path(requested_path).expanduser().resolve()
        if not notebook.is_file():
            raise FileNotFoundError(f"Notebook not found: {notebook}")
        return notebook

    candidates = [
        SCRIPT_DIRECTORY / "orthogonal_line_center_complete_experiments.ipynb",
        SCRIPT_DIRECTORY / "orthogonal_line_center_complete_experiments(2).ipynb",
    ]
    existing = [path for path in candidates if path.is_file()]
    if len(existing) == 1:
        return existing[0]
    if len(existing) > 1:
        raise RuntimeError(
            "Two candidate notebooks were found. Rename the final notebook or "
            "supply --notebook explicitly."
        )
    raise FileNotFoundError(
        "No experiment notebook was found beside this script. Place the final "
        "notebook in the same folder or supply --notebook PATH."
    )


def load_notebook_definitions(notebook):
    """Execute import/definition cells only, never the experiment call cells."""
    notebook_data = json.loads(notebook.read_text(encoding="utf-8"))
    namespace = {"__name__": "chapter5_notebook_definitions"}
    for index, cell in enumerate(notebook_data["cells"]):
        if cell.get("cell_type") == "code" and index <= 30:
            source = "".join(cell.get("source", []))
            exec(compile(source, f"{notebook.name}:cell-{index}", "exec"), namespace)
    return namespace


def independent_augmenting_path_cover(graph):
    """Minimum cover via elementary augmenting paths, independent of HK code."""
    left = list(graph["vertical_nodes"])
    right = list(graph["horizontal_nodes"])
    adjacency = {u: list(graph["adjacency"].get(u, [])) for u in left}
    matched_right = {v: None for v in right}

    def augment(u, visited_right):
        for v in adjacency[u]:
            if v in visited_right:
                continue
            visited_right.add(v)
            previous = matched_right[v]
            if previous is None or augment(previous, visited_right):
                matched_right[v] = u
                return True
        return False

    matching_size = 0
    for u in left:
        if augment(u, set()):
            matching_size += 1

    matched_left = {u: None for u in left}
    for v, u in matched_right.items():
        if u is not None:
            matched_left[u] = v

    # Alternating reachability from unmatched left vertices.  Starting from
    # left, traverse unmatched L->R edges and matched R->L edges.
    reachable_left = {u for u in left if matched_left[u] is None}
    reachable_right = set()
    stack = list(reachable_left)
    while stack:
        u = stack.pop()
        for v in adjacency[u]:
            if matched_left[u] == v:
                continue
            if v not in reachable_right:
                reachable_right.add(v)
                previous = matched_right[v]
                if previous is not None and previous not in reachable_left:
                    reachable_left.add(previous)
                    stack.append(previous)

    cover = (set(left) - reachable_left) | reachable_right
    return cover, matching_size


def brute_force_cover_size(graph):
    nodes = list(graph["vertical_nodes"]) + list(graph["horizontal_nodes"])
    for size in range(len(nodes) + 1):
        for subset in combinations(nodes, size):
            chosen = set(subset)
            if all(u in chosen or v in chosen for u, v in graph["edges"]):
                return size
    raise AssertionError("No vertex cover exists")


def cover_valid(graph, cover):
    return all(u in cover or v in cover for u, v in graph["edges"])


def scalar_nearest_line_distances(points, vertical_lines, horizontal_lines):
    """Direct scalar implementation, intentionally not using the notebook helper."""
    distances = []
    for x, y in np.asarray(points, dtype=float):
        vertical_distance = min(
            (abs(float(x) - float(line)) for line in vertical_lines),
            default=float("inf"),
        )
        horizontal_distance = min(
            (abs(float(y) - float(line)) for line in horizontal_lines),
            default=float("inf"),
        )
        distances.append(min(vertical_distance, horizontal_distance))
    return distances


def independent_geometry_check(result, object_type, atol=1.0e-8):
    if object_type == "point":
        points = result["points"]
        direct = scalar_nearest_line_distances(
            points, result["selected_v_lines"], result["selected_h_lines"]
        )
        coverage = all(distance <= result["r"] + atol for distance in direct)
        maximum_consistent = abs(
            max(direct) - result["max_distance_to_selected"]
        ) <= atol
        reduction_equivalence = True
    else:
        points = result["disk_centres"]
        centre_direct = scalar_nearest_line_distances(
            points, result["selected_v_lines"], result["selected_h_lines"]
        )
        direct = [max(0.0, value - result["disk_radius"]) for value in centre_direct]
        coverage = all(distance <= result["r"] + atol for distance in direct)
        maximum_consistent = abs(
            max(direct) - result["max_disk_distance_to_selected"]
        ) <= atol
        reduction_equivalence = all(
            (disk_distance <= result["r"] + atol)
            == (centre_distance <= result["effective_radius"] + atol)
            for disk_distance, centre_distance in zip(direct, centre_direct)
        )

    expected_count = len(result["selected_v_lines"]) + len(result["selected_h_lines"])
    checks = {
        "scalar_coverage": coverage,
        "maximum_consistent": maximum_consistent,
        "reduction_equivalence": reduction_equivalence,
        "line_count_consistent": result["selected_line_count"] == expected_count,
        "candidate_vertical_bound": result["candidate_vertical_count"] <= result["n"],
        "candidate_horizontal_bound": result["candidate_horizontal_count"] <= result["n"],
        "edge_bound": result["graph_edge_count"] <= result["n"],
        "konig_equality": result["selected_line_count"] == result["matching_size"],
        "status_consistent": (
            result["status"] == "WITHIN_2K"
        ) == (result["selected_line_count"] <= 2 * result["k"]),
    }
    return checks, max(direct)


def validate_graph(graph, production_solver, source, trial, brute_force=False):
    production_cover, production_matching = production_solver(graph)
    independent_cover, independent_matching = independent_augmenting_path_cover(graph)
    brute_size = brute_force_cover_size(graph) if brute_force else None
    checks = {
        "production_cover_valid": cover_valid(graph, production_cover),
        "independent_cover_valid": cover_valid(graph, independent_cover),
        "production_konig_equality": len(production_cover) == production_matching,
        "independent_konig_equality": len(independent_cover) == independent_matching,
        "backend_size_agreement": len(production_cover) == len(independent_cover),
        "bruteforce_agreement": (
            brute_size is None
            or len(production_cover) == len(independent_cover) == brute_size
        ),
    }
    return {
        "source": source,
        "trial": trial,
        "left_vertices": len(graph["vertical_nodes"]),
        "right_vertices": len(graph["horizontal_nodes"]),
        "edges": len(graph["edges"]),
        "production_cover_size": len(production_cover),
        "independent_cover_size": len(independent_cover),
        "bruteforce_cover_size": "" if brute_size is None else brute_size,
        **checks,
        "passed": all(checks.values()),
    }


def random_bipartite_graph(rng, left_size, right_size, probability):
    left = [f"V{i}" for i in range(left_size)]
    right = [f"H{j}" for j in range(right_size)]
    edges = [
        (u, v)
        for u in left
        for v in right
        if rng.random() < probability
    ]
    if not edges:
        edges = [(left[0], right[0])]
    edges = sorted(set(edges))
    return {
        "vertical_nodes": left,
        "horizontal_nodes": right,
        "adjacency": {u: sorted(v for x, v in edges if x == u) for u in left},
        "edges": edges,
        "edge_points": {},
    }


def main(notebook, output, seed):
    namespace = load_notebook_definitions(notebook)
    production_solver = namespace["builtin_hopcroft_karp_minimum_vertex_cover"]
    build_assignments = namespace["build_sparse_grid_assignments"]
    build_graph = namespace["build_bipartite_graph"]
    generate_points = namespace["generate_points"]
    solve_points = namespace["two_factor_kolcr"]
    solve_disks = namespace["two_factor_disk_kolcr"]

    rng = np.random.default_rng(seed)
    rows = []

    # Small random graphs are checked against a genuine exhaustive oracle.
    for trial in range(200):
        left_size = int(rng.integers(1, 6))
        right_size = int(rng.integers(1, 6))
        probability = float(rng.uniform(0.1, 0.9))
        graph = random_bipartite_graph(rng, left_size, right_size, probability)
        rows.append(validate_graph(
            graph, production_solver, "small_random_graph", trial, brute_force=True
        ))

    # Larger random graphs exercise backend agreement without exponential search.
    for trial in range(300):
        left_size = int(rng.integers(2, 31))
        right_size = int(rng.integers(2, 31))
        probability = float(rng.uniform(0.03, 0.5))
        graph = random_bipartite_graph(rng, left_size, right_size, probability)
        rows.append(validate_graph(
            graph, production_solver, "larger_random_graph", trial, brute_force=False
        ))

    distributions = ["uniform", "clustered", "grid", "near_horizontal", "repeated"]
    # Graphs induced by actual point and equal-radius-disk instances are also checked.
    for object_type in ("point", "disk"):
        for trial in range(250):
            n = int(rng.integers(2, 51))
            distribution = distributions[trial % len(distributions)]
            points = generate_points(
                n, distribution, seed=seed + 1000 * (object_type == "disk") + trial
            )
            r = float(rng.uniform(0.2, 3.0))
            rho = float(rng.uniform(0.1, 1.5)) if object_type == "disk" else 0.0
            _, _, assignments, _ = build_assignments(points, r + rho)
            graph = build_graph(assignments)
            rows.append(validate_graph(
                graph, production_solver, f"{object_type}_induced_graph", trial
            ))

    geometry_rows = []
    for object_type in ("point", "disk"):
        for trial in range(250):
            n = int(rng.integers(2, 51))
            k = int(rng.integers(1, min(8, n) + 1))
            distribution = distributions[trial % len(distributions)]
            instance_seed = seed + 5000 + 1000 * (object_type == "disk") + trial
            points = generate_points(n, distribution, seed=instance_seed, scale=40.0)
            r = float(rng.uniform(0.2, 3.0))
            rho = float(rng.uniform(0.1, 1.5))
            result = (
                solve_points(points, k, r, matching_backend="builtin")
                if object_type == "point"
                else solve_disks(points, k, r, rho, matching_backend="builtin")
            )
            checks, maximum_distance = independent_geometry_check(result, object_type)
            geometry_rows.append({
                "object_type": object_type,
                "trial": trial,
                "distribution": distribution,
                "instance_seed": instance_seed,
                "n": n,
                "k": k,
                "r": r,
                "rho": "" if object_type == "point" else rho,
                "selected_lines": result["selected_line_count"],
                "maximum_direct_distance": maximum_distance,
                **checks,
                "passed": all(checks.values()),
            })

    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "independent_backend_validation.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    geometry_csv = output / "independent_geometry_validation.csv"
    with geometry_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(geometry_rows[0]))
        writer.writeheader()
        writer.writerows(geometry_rows)

    failed = [row for row in rows if not row["passed"]]
    by_source = {}
    for source in sorted({row["source"] for row in rows}):
        selected = [row for row in rows if row["source"] == source]
        by_source[source] = {
            "cases": len(selected),
            "passed": sum(bool(row["passed"]) for row in selected),
            "maximum_vertices": max(
                row["left_vertices"] + row["right_vertices"] for row in selected
            ),
            "maximum_edges": max(row["edges"] for row in selected),
        }
    summary = {
        "seed": seed,
        "notebook": str(notebook),
        "production_backend": "builtin Hopcroft-Karp plus Konig recovery",
        "independent_backend": "elementary DFS augmenting paths plus alternating recovery",
        "total_cases": len(rows),
        "passed": len(rows) - len(failed),
        "failed": len(failed),
        "small_graphs_checked_by_bruteforce": sum(
            row["source"] == "small_random_graph" for row in rows
        ),
        "by_source": by_source,
        "independent_scalar_geometry": {
            "total_cases": len(geometry_rows),
            "passed": sum(bool(row["passed"]) for row in geometry_rows),
            "failed": sum(not bool(row["passed"]) for row in geometry_rows),
            "point_cases": sum(row["object_type"] == "point" for row in geometry_rows),
            "disk_cases": sum(row["object_type"] == "disk" for row in geometry_rows),
        },
    }
    (output / "independent_backend_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    if failed:
        raise AssertionError(f"Independent backend validation failed in {len(failed)} case(s)")
    geometry_failed = [row for row in geometry_rows if not row["passed"]]
    if geometry_failed:
        raise AssertionError(
            f"Independent geometry validation failed in {len(geometry_failed)} case(s)"
        )


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run the dissertation's independent Chapter 5 validation."
    )
    parser.add_argument(
        "--notebook",
        type=Path,
        default=None,
        help="Path to the final experiment notebook. By default, search beside this script.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Directory for the validation CSV and JSON files.",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_arguments()
    selected_notebook = resolve_notebook(arguments.notebook)
    selected_output = arguments.output.expanduser().resolve()
    print("Notebook:", selected_notebook)
    print("Output directory:", selected_output)
    main(selected_notebook, selected_output, arguments.seed)
