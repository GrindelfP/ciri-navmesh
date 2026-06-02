#pragma once

/**
 * @file navmesh.hh
 * @brief Navigation mesh: a weighted adjacency graph built from a DCEL
 *        triangulation, used as input to A* and the Funnel Algorithm.
 *
 * ## Conceptual model
 *
 * After triangulation, each bounded face (triangle) becomes a **node** in the
 * navmesh graph.  Two nodes are connected by an **arc** if the corresponding
 * triangles share an edge (i.e. are face-adjacent in the DCEL).
 *
 * Arc weight = Euclidean distance between the centroids of the two triangles.
 * This gives A* an admissible heuristic (straight-line distance to goal) and
 * makes path cost meaningful for comparison metrics.
 *
 * ## Data stored per node
 *
 *   - `faceIdx`   : original FaceIdx in the source DCEL (stable reference)
 *   - `centroid`  : pre-computed centroid of the triangle
 *   - `vertices`  : the three VertexIdx of the triangle (for Funnel Algorithm)
 *   - `neighbors` : up to 3 arcs to adjacent nodes
 *
 * ## Data stored per arc
 *
 *   - `to`        : destination node index (NodeIdx)
 *   - `cost`      : centroid-to-centroid Euclidean distance
 *   - `portal`    : the shared edge (two VertexIdx), used by Funnel Algorithm
 *                   to determine the passable "gate" between triangles
 *
 * ## Coordinate system
 *
 * NavMesh works entirely in the same 2-D coordinate space as the DCEL.
 * No transformation is applied.
 *
 * ## Rebuild vs reuse
 *
 * `NavMesh::build()` rebuilds from scratch each time.  Since a navmesh is
 * typically built once per triangulation switch (D/G/Q keys in the demo),
 * this is acceptable.  The build is O(n) in the number of triangles.
 *
 * ## Point location
 *
 * `NavMesh::findNode()` locates the triangle containing a query point.
 * It delegates to `DCEL::findFace()` and then looks up the NodeIdx via the
 * internal `faceToNode_` map.  O(n) — a spatial index can be added later.
 */

#include "../geometry/dcel.hh"
#include "../geometry/primitives.hpp"

#include <array>
#include <optional>
#include <vector>

namespace navmesh {

// ─────────────────────────────────────────────────────────────────────────────
//  Index type
// ─────────────────────────────────────────────────────────────────────────────

/// @brief Index into NavMesh::nodes_.  Distinct from geometry::FaceIdx.
using NodeIdx = std::size_t;

/// @brief Sentinel: no node.
inline constexpr NodeIdx kInvalidNode = std::numeric_limits<NodeIdx>::max();

// ─────────────────────────────────────────────────────────────────────────────
//  Arc  (directed adjacency entry)
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief A directed arc between two navmesh nodes.
 *
 * Because the graph is undirected, each undirected edge is represented by
 * two Arc entries (one in each node's neighbor list), symmetrically.
 */
struct Arc {
    NodeIdx to;     ///< Destination node index

    /// @brief Euclidean distance between the centroids of the two triangles.
    double cost;

    /**
     * @brief The "portal" edge shared by the two triangles.
     *
     * Stored as the two vertex indices (from the source DCEL) of the shared
     * half-edge, in the order they appear in the SOURCE triangle's face loop
     * (i.e. CCW around the source face).
     *
     * The Funnel Algorithm uses portals to constrain path smoothing:
     * the path must pass through the portal segment when moving from
     * this node to `to`.
     */
    std::array<geometry::VertexIdx, 2> portal;
};

// ─────────────────────────────────────────────────────────────────────────────
//  Node
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief A node in the navmesh graph, corresponding to one triangle.
 */
struct Node {
    geometry::FaceIdx faceIdx;   ///< Corresponding face in the source DCEL

    geometry::Point2D centroid;  ///< Pre-computed centroid (a+b+c)/3

    /// @brief The three vertex indices of the triangle in CCW order.
    std::array<geometry::VertexIdx, 3> vertices;

    /// @brief Up to 3 arcs to adjacent triangles (hull edges have no arc).
    std::vector<Arc> neighbors;
};

// ─────────────────────────────────────────────────────────────────────────────
//  NavMesh
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Navigation mesh graph built from a DCEL triangulation.
 *
 * ### Typical usage
 * @code
 *   navmesh::NavMesh nm;
 *   nm.build(dcel);
 *
 *   NodeIdx start = nm.findNode(startPt);
 *   NodeIdx goal  = nm.findNode(goalPt);
 *   // ... pass nm to AStar ...
 * @endcode
 */
class NavMesh {
public:
    NavMesh() = default;

    // ── Construction ──────────────────────────────────────────────────────────

    /**
     * @brief Builds the navmesh from a triangulated DCEL.
     *
     * Iterates over all live bounded faces, creates one Node per face,
     * then links adjacent faces via Arc entries.
     *
     * Previous navmesh data is discarded.
     *
     * @param dcel  A fully-triangulated DCEL (all bounded faces are triangles).
     *              The DCEL must remain valid for the lifetime of this NavMesh
     *              if `vertexPos()` is called; it is NOT stored by reference
     *              internally — only indices are kept.
     *
     * @note O(n) in the number of triangles.
     */
    void build(const geometry::DCEL& dcel);

    /// @brief Discards all nodes and resets to empty state.
    void clear() noexcept;

    // ── Accessors ──────────────────────────────────────────────────────────────

    /// @brief Number of nodes (triangles) in the navmesh.
    [[nodiscard]] std::size_t nodeCount() const noexcept;

    /// @brief Access node by index.  Asserts valid index.
    [[nodiscard]] const Node& node(NodeIdx i) const;

    /**
     * @brief Finds the navmesh node whose triangle contains point @p p.
     *
     * Delegates to DCEL::findFace() and maps the result through faceToNode_.
     *
     * @param p     Query point.
     * @param dcel  The same DCEL used in build().
     * @return      NodeIdx, or kInvalidNode if p is outside the triangulation.
     */
    [[nodiscard]] NodeIdx findNode(geometry::Point2D     p,
                                   const geometry::DCEL& dcel) const;

    /**
     * @brief Returns the centroid of node @p i.
     *
     * Convenience wrapper around node(i).centroid.
     */
    [[nodiscard]] geometry::Point2D centroid(NodeIdx i) const;

    /**
     * @brief Returns the portal (shared edge) on the arc from @p from
     *        to @p to, as a pair of world-space Point2D.
     *
     * Looks up the Arc in from's neighbor list and retrieves vertex positions
     * from the stored vertex coordinates.
     *
     * @param from  Source node.
     * @param to    Destination node (must be adjacent to from).
     * @param dcel  Source DCEL (needed to dereference VertexIdx → Point2D).
     * @return      {left, right} portal points in CCW order around `from`.
     *              Returns std::nullopt if `to` is not adjacent to `from`.
     */
    [[nodiscard]] std::optional<std::array<geometry::Point2D, 2>>
    portalPoints(NodeIdx               from,
                 NodeIdx               to,
                 const geometry::DCEL& dcel) const;

    // ── Statistics ────────────────────────────────────────────────────────────

    /**
     * @brief Returns the total number of arcs (counting both directions).
     */
    [[nodiscard]] std::size_t arcCount() const noexcept;

    /**
     * @brief Average node degree (arcs per node).
     */
    [[nodiscard]] double averageDegree() const noexcept;

private:
    std::vector<Node> nodes_;

    /**
     * @brief Maps geometry::FaceIdx → NodeIdx.
     *
     * Built during build() and used by findNode().
     * Size = dcel.faceCount() (includes dead faces, mapped to kInvalidNode).
     */
    std::vector<NodeIdx> faceToNode_;
};

} // namespace navmesh
