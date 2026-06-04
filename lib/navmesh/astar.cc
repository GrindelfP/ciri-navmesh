/**
 * @file astar.cc
 * @brief A* shortest-path search on a NavMesh graph.
 */

#include "astar.hh"

#include <cassert>
#include <limits>
#include <queue>
#include <vector>

namespace pathfinding {
    using navmesh::NavMesh;
    using navmesh::NodeIdx;
    using navmesh::kInvalidNode;

    // ─────────────────────────────────────────────────────────────────────────────
    //  search
    // ─────────────────────────────────────────────────────────────────────────────

    PathResult AStar::search(const NavMesh &nm,
                             NodeIdx start,
                             NodeIdx goal) const {
        assert(start != kInvalidNode && "AStar::search — invalid start node");
        assert(goal != kInvalidNode && "AStar::search — invalid goal node");

        PathResult result;

        // Trivial case: start == goal.
        if (start == goal) {
            result.nodes = {start};
            result.cost = 0.0;
            result.found = true;
            return result;
        }

        const std::size_t N = nm.nodeCount();

        // Per-node best g-score. Initialised to infinity.
        constexpr double kInf = std::numeric_limits<double>::infinity();
        std::vector<double> gScore(N, kInf);
        std::vector<NodeIdx> parent(N, kInvalidNode);

        gScore[start] = 0.0;

        // Min-heap ordered by f = g + h.
        std::priority_queue<HeapEntry,
            std::vector<HeapEntry>,
            std::greater<HeapEntry> > open;

        open.push({heuristic(nm, start, goal), 0.0, start});

        while (!open.empty()) {
            const auto [f, g, cur] = open.top();
            open.pop();

            ++result.nodesExpanded;

            // Lazy deletion: skip stale entries.
            if (g > gScore[cur] + 1e-12) continue;

            // Goal reached.
            if (cur == goal) {
                result.found = true;
                result.cost = gScore[goal];
                result.nodes = reconstructPath(parent, start, goal);
                return result;
            }

            // Expand neighbours.
            for (const auto &arc: nm.node(cur).neighbors) {
                const double tentative = gScore[cur] + arc.cost;
                if (tentative < gScore[arc.to]) {
                    gScore[arc.to] = tentative;
                    parent[arc.to] = cur;
                    const double h = heuristic(nm, arc.to, goal);
                    open.push({tentative + h, tentative, arc.to});
                }
            }
        }

        // No path found — return empty result (found == false).
        return result;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  heuristic
    // ─────────────────────────────────────────────────────────────────────────────

    double AStar::heuristic(const NavMesh &nm,
                            NodeIdx a,
                            NodeIdx b) noexcept {
        return nm.centroid(a).dist(nm.centroid(b));
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  reconstructPath
    // ─────────────────────────────────────────────────────────────────────────────

    std::vector<NodeIdx> AStar::reconstructPath(
        const std::vector<NodeIdx> &parent,
        NodeIdx start,
        NodeIdx goal) {
        std::vector<NodeIdx> path;
        for (NodeIdx cur = goal; cur != kInvalidNode; cur = parent[cur]) {
            path.push_back(cur);
            if (cur == start) break;
        }
        std::reverse(path.begin(), path.end());
        return path;
    }
} // namespace pathfinding
