#pragma once

/**
 * @file astar.hh
 * @brief A* shortest-path search on a NavMesh graph.
 *
 * ## Algorithm overview
 *
 * Standard A* on the navmesh node graph (triangle adjacency):
 *
 *   - **g(n)** = accumulated centroid-to-centroid arc cost from start to n.
 *   - **h(n)** = Euclidean distance from centroid(n) to centroid(goal).
 *                Admissible because arc costs ARE centroid-to-centroid distances,
 *                so h never overestimates.
 *   - **f(n)** = g(n) + h(n)
 *
 * The open set is a min-heap (std::priority_queue with a reverse comparator)
 * keyed by f.  Lazy deletion is used: a node is "expanded" only if its
 * recorded g matches the best known g (stale heap entries are skipped).
 *
 * ## Output
 *
 * Returns a sequence of NodeIdx from start to goal (inclusive).
 * This "triangle corridor" is then passed to the Funnel Algorithm for
 * geometric path smoothing.
 *
 * ## Complexity
 *
 *   O((V + E) log V) where V = triangle count, E = adjacency arcs.
 *   For a triangulation of n points: V ≈ 2n, E ≈ 6n → O(n log n).
 *
 * ## Design notes
 *
 * - AStar is stateless between calls; all working memory is allocated
 *   on the stack/local heap per call.
 * - The NavMesh is passed by const reference — A* does not modify it.
 * - `PathResult` carries both the node sequence (for Funnel) and scalar
 *   metrics (for the HUD and benchmark output).
 */

#include "../navmesh/navmesh.hh"
#include "../geometry/primitives.hpp"

#include <optional>
#include <vector>

namespace pathfinding {

// ─────────────────────────────────────────────────────────────────────────────
//  PathResult
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Result of an A* query.
 *
 * If no path exists (start/goal outside the mesh, or disconnected), the
 * `nodes` vector is empty and `found` is false.
 */
struct PathResult {
    /// Sequence of NavMesh node indices from start to goal (both inclusive).
    /// Passed directly to the Funnel Algorithm.
    std::vector<navmesh::NodeIdx> nodes;

    /// Total cost along the returned path (sum of arc costs = centroid distances).
    double cost{0.0};

    /// Number of nodes expanded during the search (for benchmarking).
    std::size_t nodesExpanded{0};

    /// True if a path was found.
    bool found{false};
};

// ─────────────────────────────────────────────────────────────────────────────
//  AStar
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief A* pathfinder on a NavMesh.
 *
 * ### Usage
 * @code
 *   pathfinding::AStar astar;
 *   auto result = astar.search(nm, startNode, goalNode);
 *   if (result.found) {
 *       // pass result.nodes to FunnelAlgorithm
 *   }
 * @endcode
 *
 * The object is stateless and may be reused across multiple queries.
 */
class AStar {
public:
    AStar() = default;

    /**
     * @brief Runs A* from @p start to @p goal on the given NavMesh.
     *
     * ### Steps
     *
     * 1. Initialize: push (f=h(start), g=0, node=start) onto the open heap.
     * 2. Pop the node with lowest f.
     *    - If it is the goal: reconstruct path and return.
     *    - If g recorded for this node is stale (lazy deletion): skip.
     * 3. For each arc (neighbor, cost) of the current node:
     *    - Compute tentative g' = g_current + arc.cost.
     *    - If g' < best known g for neighbor: update, push to heap.
     * 4. If heap is empty without reaching goal: return empty PathResult.
     *
     * @param nm     The NavMesh to search.
     * @param start  Starting node index (must be valid, i.e. != kInvalidNode).
     * @param goal   Goal node index (must be valid).
     * @return       PathResult with the node corridor and metrics.
     *
     * @note If start == goal, returns a single-node path immediately.
     */
    [[nodiscard]] PathResult search(const navmesh::NavMesh& nm,
                                    navmesh::NodeIdx        start,
                                    navmesh::NodeIdx        goal) const;

private:
    // ── Internal types ────────────────────────────────────────────────────────

    /**
     * @brief An entry in the A* open-set priority queue.
     *
     * Ordered by f ascending (min-heap via std::greater in priority_queue).
     */
    struct HeapEntry {
        double          f;    ///< f = g + h
        double          g;    ///< Accumulated cost from start
        navmesh::NodeIdx node; ///< Node index

        /// Min-heap comparison: smaller f has higher priority.
        bool operator>(const HeapEntry& o) const noexcept { return f > o.f; }
    };

    // ── Helpers ───────────────────────────────────────────────────────────────

    /**
     * @brief Euclidean heuristic: distance between centroids of @p a and @p b.
     *
     * Admissible because arc weights are centroid-to-centroid distances and
     * the straight-line distance is a lower bound.
     *
     * @param nm    The NavMesh (for centroid lookup).
     * @param a     Source node.
     * @param b     Goal node.
     * @return      Euclidean distance between their centroids.
     */
    [[nodiscard]] static double heuristic(const navmesh::NavMesh& nm,
                                          navmesh::NodeIdx        a,
                                          navmesh::NodeIdx        b) noexcept;

    /**
     * @brief Reconstructs the path from the @p parent table.
     *
     * Walks backwards from @p goal to @p start following parent pointers,
     * then reverses the result.
     *
     * @param parent  parent[i] = the node we came from to reach node i
     *                (kInvalidNode for the start).
     * @param start   Start node (used as termination condition).
     * @param goal    Goal node (reconstruction begins here).
     * @return        Ordered node sequence [start, ..., goal].
     */
    [[nodiscard]] static std::vector<navmesh::NodeIdx>
    reconstructPath(const std::vector<navmesh::NodeIdx>& parent,
                    navmesh::NodeIdx                     start,
                    navmesh::NodeIdx                     goal);
};

} // namespace pathfinding
