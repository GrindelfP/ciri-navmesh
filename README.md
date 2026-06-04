# ciri-navmesh: A NavMesh Triangulation Library 

> A C++ library for navigation mesh construction via triangulation, with shortest-path planning using A\* and the Funnel Algorithm. Compares **Delaunay**, **Greedy MWT**, and **Quasi-Greedy MWT** triangulations against one another on weight, build time, and path quality metrics.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Architecture](#architecture)
- [Algorithms](#algorithms)
  - [Triangulation](#triangulation)
  - [Pathfinding](#pathfinding)
- [Project Structure](#project-structure)
- [Building](#building)
- [Usage](#usage)
  - [Benchmark CLI (`navmesh_bench`)](#benchmark-cli-navmesh_bench)
  - [Visualiser CLI (`navmesh_vis`)](#visualiser-cli-navmesh_vis)
  - [Input Format](#input-format)
- [Output](#output)
  - [Benchmark Table](#benchmark-table)
  - [CSV Export](#csv-export)
  - [SVG Visualisation](#svg-visualisation)
- [Data Structures](#data-structures)
  - [DCEL (Half-Edge Mesh)](#dcel-half-edge-mesh)
  - [NavMesh Graph](#navmesh-graph)
- [Metrics](#metrics)
- [Complexity Reference](#complexity-reference)
- [Sample Files](#sample-files)
- [References](#references)

---

## Overview

Minimum Weight Triangulation (MWT) is an NP-hard problem: finding the triangulation of a planar point set that minimises the total edge length. This project implements and benchmarks three practical approaches:

| Algorithm | Strategy | Provable guarantee |
|---|---|---|
| **Delaunay (Bowyer–Watson)** | Maximise minimum angle | Not MWT; O(n log n) avg |
| **Greedy MWT** | Accept shortest non-crossing edge | Heuristic; O(n² log n) |
| **Quasi-Greedy MWT** | Lune-skeleton + pocket fill + flips | Constant-factor approx; O(n²) |

Once a triangulation is built, it is converted into a **navigation mesh** (NavMesh) — a weighted adjacency graph of triangles — on which A\* finds a triangle-corridor path and the **Funnel Algorithm** computes the geometrically shortest Euclidean path through that corridor.

---

## Features

- Three triangulation algorithms behind a clean `ITriangulator` interface
- Full **DCEL** (Doubly Connected Edge List / half-edge mesh) implementation with:
  - Stable integer indices (no pointer invalidation)
  - Soft deletion via `dead` flags
  - Edge flip, face split, edge split, twin stitching
  - Debug-mode invariant validation (`validate()`)
- **NavMesh** graph built directly from the DCEL in O(n)
- A\* shortest-path search with centroid-distance heuristic (admissible)
- **Funnel Algorithm** (Simple Stupid Funnel, Mononen 2010) for Euclidean path smoothing
- Benchmark CLI with formatted table output and optional CSV export
- SVG visualiser producing per-algorithm images with triangle mesh, corridor highlight, and smoothed path
- Robust geometric predicates: `orientation`, `inCircle`, `segmentIntersect`, `lineIntersectionPoint`, `isLocallyDelaunay`
- Pure C++ standard library — no external algorithmic dependencies

---

## Architecture

```
┌─────────────────────────────────────────────────────┐
│                     Application                     │
│           navmesh_bench          navmesh_vis        │
└──────────────┬──────────────────────────┬───────────┘
               │                          │
┌──────────────▼──────────────┐  ┌────────▼────────────┐
│       Triangulation         │  │     SVG Writer      │
│  ┌──────────────────────┐   │  └─────────────────────┘
│  │  ITriangulator (ABC) │   │
│  └──────────────────────┘   │
│   Delaunay  Greedy  Quasi   │
└──────────────┬──────────────┘
               │  fills
┌──────────────▼──────────────┐
│         DCEL                │
│  Vertex │ HalfEdge │ Face   │
└──────────────┬──────────────┘
               │  feeds
┌──────────────▼──────────────┐
│           NavMesh           │
│   Node (triangle) + Arc     │
└──────┬───────────────┬──────┘
       │               │
┌──────▼──────┐  ┌─────▼───────┐
│    A*       │  │   Funnel    │
│ (corridor)  │  │  (smooth)   │
└─────────────┘  └─────────────┘
```

---

## Algorithms

### Triangulation

#### Delaunay — Bowyer–Watson

Incremental point-insertion algorithm. A super-triangle enclosing all points is inserted first; then for each input point:

1. Find all triangles whose circumcircle contains the new point (**cavity**) via DFS flood-fill.
2. Extract the polygonal **boundary** of the cavity.
3. Kill cavity triangles and re-connect the new point to each boundary edge.
4. Remove the super-triangle and its incident faces.

The result satisfies the **Delaunay condition**: no input point lies strictly inside the circumcircle of any triangle.

#### Greedy MWT

A classical baseline:

1. Generate all C(n, 2) candidate edges, sorted by Euclidean length.
2. Accept each edge if it does not **properly intersect** any already-accepted edge.
3. Build the DCEL from the accepted edge set.

Guaranteed to produce a valid triangulation of the convex hull (Chazelle et al.), but not optimal.

#### Quasi-Greedy MWT (Levcopoulos & Krznaric 1998)

The best-known polynomial constant-factor MWT heuristic:

1. **Phase 1 — Delaunay:** Run Bowyer–Watson to obtain MWT candidate edges (only Delaunay edges can be in the MWT).
2. **Phase 2 — Lune test:** For each candidate edge (p, q), the open disk with diameter pq must contain no other input point (`dot(r−p, r−q) < 0` test). Passing edges are "safe" / forced.
3. **Phase 3 — Pocket fill:** Safe edges form a planar graph; polygonal faces with > 3 vertices are **pockets** triangulated greedily by diagonal length.
4. **Phase 3.5 — Fallback:** Any remaining un-triangulated gaps are closed with shortest non-crossing edges.
5. **Phase 4 — Weight-reducing flips:** Repeatedly scan interior edges; flip edge u–v → w–x when `dist(w, x) < dist(u, v)`. Converges in < 4 passes typically.

### Pathfinding

#### A\* on NavMesh

Standard A\* on the triangle adjacency graph:

- **g(n)** = accumulated centroid-to-centroid cost from start
- **h(n)** = Euclidean distance from centroid(n) to centroid(goal) ← admissible
- **f(n)** = g(n) + h(n)

Min-heap with lazy deletion. Returns a **triangle corridor** (sequence of NodeIdx).

#### Funnel Algorithm

Converts the triangle corridor into the shortest Euclidean path (Lee & Preparata 1984; Simple Stupid Funnel, Mononen 2010):

- Maintain two deques (`leftDeq`, `rightDeq`) sharing the current **apex** at their front.
- For each portal (left/right vertex of the shared edge between adjacent triangles):
  - If the new vertex tightens the funnel on its side → push it onto the chain.
  - If it crosses to the other side → pop vertices from the opposite chain, emitting each as a **waypoint**, until the funnel is valid; reset apex.
- Close the funnel to the goal point.

Result: minimum waypoints needed for the shortest path, in O(k) time.

---

## Project Structure

```
project/
├── CMakeLists.txt
├── lib/
│   ├── geometry/
│   │   ├── primitives.hpp / primitives.cpp   # Point2D, Segment, Triangle
│   │   ├── predicates.hpp / predicates.cpp   # orientation, inCircle, segmentIntersect …
│   │   └── dcel.hh / dcel.cc                 # DCEL half-edge mesh
│   ├── triangulation/
│   │   ├── i_triangulator.hh                 # Abstract interface + TriangulationResult
│   │   ├── delaunay.hh / delaunay.cc         # Bowyer–Watson
│   │   ├── greedy.hh / greedy.cc             # Greedy MWT
│   │   └── quasi_greedy.hh / quasi_greedy.cc # Quasi-Greedy MWT
│   └── navmesh/
│       ├── navmesh.hh / navmesh.cc           # NavMesh graph
│       ├── astar.hh / astar.cc               # A* search
│       └── funnel.hh / funnel.cc             # Funnel Algorithm
├── app/
│   ├── main.cc                               # navmesh_bench CLI
│   └── vis_main.cc                           # navmesh_vis SVG visualiser
├── data/
│   ├── points.txt                            # Small test set (8 pts)
│   ├── tavern.txt                            # Room with obstacles (~25 pts)
│   ├── labyrinth.txt                         # Labyrinth scene (~70 pts)
│   └── labyrinth_big.txt                     # Large labyrinth (~500 pts)
└── tests/
    └── test_geometry.cpp                     # Unit tests
```

File extension convention: `.hh` / `.cc` for all files, except for `primitives` and `predicates` retain `.hpp` / `.cpp`.

---

## Building

### Prerequisites

- CMake ≥ 3.16
- C++20-capable compiler (GCC 10+, Clang 12+, MSVC 2019+)

### Build steps

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

This produces:

| Target | Binary | Description |
|---|---|---|
| `geometry` | `libgeometry.a` | Geometry static library |
| `navmesh_bench` | `navmesh_bench` | Benchmark CLI |
| `test_geometry` | `test_geometry` | Unit test runner |

Run tests:

```bash
ctest --test-dir build --output-on-failure
```

---

## Usage

### Benchmark CLI (`navmesh_bench`)

```
navmesh_bench [OPTIONS] <input_file>

Options:
  -a, --algo <name>    delaunay | greedy | quasi  (default: all three)
  -s, --start <x,y>   Start point for pathfinding
  -g, --goal  <x,y>   Goal  point for pathfinding
  -n, --no-flip        Disable edge flips in quasi-greedy
  -q, --quiet          Print summary table only
  -o, --output <file>  Write CSV results to file
  -h, --help           Print this help
```

**Examples:**

```bash
# Run all three algorithms on the labyrinth, auto-select start/goal
./navmesh_bench data/labyrinth.txt

# Run Delaunay only with explicit start and goal
./navmesh_bench -a delaunay -s 0.1,0.1 -g 19.9,19.9 data/labyrinth.txt

# Benchmark all algorithms, save results to CSV
./navmesh_bench -o results.csv data/labyrinth_big.txt

# Quasi-greedy without flips, quiet mode
./navmesh_bench -a quasi -n -q data/tavern.txt
```

### Visualiser CLI (`navmesh_vis`)

```
navmesh_vis [OPTIONS] <input_file>

Options:
  -s, --start <x,y>    Start point
  -g, --goal  <x,y>    Goal  point
  -n, --no-flip         Disable edge flips in quasi-greedy
  -o, --output <pfx>    Output filename prefix (default: out)
  -W, --width  <px>     SVG width  (default: 900)
  -H, --height <px>     SVG height (default: 900)
  -h, --help            Print this help
```

Produces three SVG files: `<prefix>_delaunay.svg`, `<prefix>_greedy.svg`, `<prefix>_quasi.svg`.

```bash
./navmesh_vis -o viz/labyrinth data/labyrinth.txt
# → viz/labyrinth_delaunay.svg
# → viz/labyrinth_greedy.svg
# → viz/labyrinth_quasi.svg
```

### Input Format

Plain-text file, one point per line, space-separated doubles. Lines starting with `#` are comments. Duplicate points are silently removed.

```
# Outer boundary
0.0  0.0
20.0 0.0
20.0 20.0
0.0  20.0

# Obstacle 1
5.0 5.0
7.0 5.0
6.0 7.0
```

---

## Output

### Benchmark Table

```
 Algorithm                               | Triangles   | Weight        | Build (ms)   | PathLen   | Turns   |
--------------------------------------------------------------------------------------------------------------
 Delaunay (Bowyer–Watson)                | 1590        | 5677.42       | 107.1        | 85.05     | 29      |
 Greedy MWT                              | 1590        | 5475.05       | 20578.1      | 81.55     | 29      |
 Quasi-Greedy MWT (Levcopoulos–Krznaric) | 1590        | 5492.22       | 17854.1      | 87.46     | 37      |```
```
### CSV Export

When `-o results.csv` is specified, a machine-readable CSV is written:

```
algorithm,triangles,weight,build_ms,path_length,turns,pathfind_ms
Delaunay (Bowyer–Watson),198,4521.720000,12.300000,87.400000,3,0.450000
Greedy MWT,198,4318.550000,43.100000,83.100000,2,0.210000
Quasi-Greedy MWT (Levcopoulos–Krznaric),198,4289.110000,15.800000,81.900000,2,0.230000
```

### SVG Visualisation

Each SVG contains (rendered in layers):

1. **Triangle mesh** — lightly tinted faces + grey edge strokes
2. **A\* corridor** — yellow semi-transparent overlay
3. **Input points** — dark blue circles
4. **Smoothed path** — bold coloured polyline with waypoint dots
5. **Start (S) / Goal (G)** — green / red circle markers
6. **Legend box** — algorithm name, triangle count, weight, build time, path length, turns

Colour theme per algorithm:

| Algorithm | Colour |
|---|---|
| Delaunay | Blue `#1a6bbf` |
| Greedy | Green `#217a3c` |
| Quasi-Greedy | Orange `#c46200` |

---

## Data Structures

### DCEL (Half-Edge Mesh)

The core data structure is a **Doubly Connected Edge List** with three pools of elements:

```
Vertex    pos: Point2D,  edge: HalfEdgeIdx,  dead: bool
HalfEdge  origin, twin, next, prev: HalfEdgeIdx,  face: FaceIdx,  dead: bool
Face      edge: HalfEdgeIdx,  dead: bool
```

Every undirected edge (u, v) is stored as two directed half-edges. Face 0 is always the **outer (unbounded) face**. All cross-references use stable integer indices — appending never invalidates existing indices; deletion sets the `dead` flag.

**Maintained invariants (I1–I7):**

| # | Invariant |
|---|---|
| I1 | `h→twin→twin == h` |
| I2 | `h→next→prev == h` |
| I3 | `h→next→origin == h→twin→origin` |
| I4 | `h→face == h→next→face` |
| I5 | `v→edge→origin == v` |
| I6 | `f→edge→face == f` |
| I7 | Outer face always present at index 0 |

Call `dcel.validate()` in debug builds to check all invariants in O(n).

### NavMesh Graph

Built from the DCEL in two O(n) passes:

- **Pass 1:** One `Node` per live bounded face — stores `faceIdx`, pre-computed `centroid`, and three `VertexIdx`.
- **Pass 2:** For each node, walk its three half-edges; each twin that points to another live face creates an `Arc` with `cost = centroid–centroid distance` and `portal = {left_vertex, right_vertex}`.

The portal orientation follows the CCW winding of the **source** face, so that the Funnel Algorithm receives consistently oriented left/right pairs.

---

## Metrics

| Metric | Description |
|---|---|
| **Total weight** | Sum of all undirected edge lengths in the triangulation |
| **Build time** | Wall-clock time of the triangulation phase (`std::chrono`) |
| **Path length** | Total Euclidean length of the smoothed funnel path |
| **Turn count** | Interior waypoints where the direction change exceeds π/12 (15°) |
| **Pathfind time** | Wall-clock time of A\* search alone |

---

## Complexity Reference

| Component | Time complexity | Space |
|---|---|---|
| Bowyer–Watson | O(n²) worst / O(n log n) avg | O(n) |
| Greedy MWT | O(n² log n) | O(n²) |
| Quasi-Greedy MWT | O(n²) (O(n log n) with spatial index) | O(n) |
| Edge flip pass | O(n) per pass | O(1) |
| NavMesh build | O(n) | O(n) |
| A\* search | O((V+E) log V), V≈2n, E≈6n | O(n) |
| Funnel Algorithm | O(k), k = corridor length | O(k) |
| SVG render | O(n) | O(n) |

---

## Sample Files

| File | Points | Description |
|---|---|---|
| `points.txt` | 8 | Minimal test — small grid |
| `tavern.txt` | ~25 | Square room with 4 triangle obstacles and scattered interior points |
| `labyrinth.txt` | ~70 | 20×20 labyrinth with 18 irregular barrier polygons |
| `labyrinth_big.txt` | ~500 | 50×50 labyrinth with 174 barriers — stress test |

---

## References

1. **Bowyer, A.** (1981). Computing Dirichlet tessellations. *The Computer Journal*, 24(2), 162–166.
2. **Watson, D.F.** (1981). Computing the n-dimensional Delaunay tessellation with application to Voronoi polytopes. *The Computer Journal*, 24(2), 167–172.
3. **Levcopoulos, C. & Krznaric, D.** (1998). Quasi-greedy triangulations approximating the minimum weight triangulation. *Journal of Algorithms*, 27(2), 303–338.
4. **Mulzer, W. & Rote, G.** (2008). Minimum-weight triangulation is NP-hard. *Journal of the ACM*, 55(2).
5. **Lee, D.T. & Preparata, F.P.** (1984). Euclidean shortest paths in the presence of rectilinear barriers. *Networks*, 14(3), 393–410.
6. **Mononen, M.** (2010). Simple Stupid Funnel Algorithm. *Digesting Duck Blog*. https://digestingduck.blogspot.com/2010/03/simple-stupid-funnel-algorithm.html
7. **de Berg, M. et al.** (2008). *Computational Geometry: Algorithms and Applications* (3rd ed.). Springer.
