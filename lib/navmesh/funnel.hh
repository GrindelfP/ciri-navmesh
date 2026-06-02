#pragma once

/**
 * @file funnel.hh
 * @brief Funnel Algorithm for geometric path smoothing over a triangle corridor.
 *
 * ## Background
 *
 * A* on the navmesh returns a sequence of triangles (a "corridor").  The path
 * through centroids is valid but not optimal — it needlessly bends at every
 * triangle boundary.  The Funnel Algorithm computes the shortest Euclidean
 * path through the corridor in O(k) time, where k is the corridor length.
 *
 * ## Algorithm overview (Lee & Preparata 1984, Chazelle 1982)
 *
 * Conceptually, maintain a "funnel" — a wedge-shaped region rooted at the
 * last confirmed waypoint, bounded by two chains of vertices (left and right
 * "tangent chains").  As each portal (shared edge between adjacent triangles)
 * is processed:
 *
 *   1. If the new portal vertex tightens the funnel on the same side → push it.
 *   2. If it crosses to the other side → pop vertices from the opposite chain
 *      until the funnel is valid again, emitting waypoints for each popped
 *      vertex that "wraps" around an apex corner.
 *   3. When the goal point is reached, close the funnel to it.
 *
 * The result is the shortest path from start to goal through the corridor,
 * touching only the vertices of the portals where turns are required.
 *
 * ## Implementation variant
 *
 * We implement the "simple funnel" variant (Mononen 2010, redblobgames):
 *
 *   - The funnel is represented as two deques (left chain, right chain).
 *   - The apex (last confirmed turn vertex) is tracked separately.
 *   - Cross-product sign tests replace angle comparisons (faster, no trig).
 *
 * ## Input
 *
 * - A start point and a goal point (arbitrary positions inside their
 *   respective triangles — not necessarily centroids).
 * - The portal sequence: for each step from triangle i to triangle i+1,
 *   two vertices (left, right) forming the shared edge, oriented CCW
 *   around the left face.
 *
 * ## Output
 *
 * - A `SmoothPath` containing the waypoint sequence (including start and goal).
 * - Scalar metrics: total length, turn count (inflection points with angle
 *   above a configurable threshold).
 *
 * ## Complexity
 *
 *   O(k) time and space, where k = number of portals = corridor length − 1.
 *
 * ## References
 *
 *   Lee, D.T. & Preparata, F.P. (1984).
 *   "Euclidean shortest paths in the presence of rectilinear barriers."
 *   Networks, 14(3), 393–410.
 *
 *   Mononen, M. (2010). "Simple Stupid Funnel Algorithm."
 *   https://digestingduck.blogspot.com/2010/03/simple-stupid-funnel-algorithm.html
 */

#include "../navmesh/navmesh.hh"
#include "../geometry/primitives.hpp"

#include <vector>

namespace pathfinding {

// ─────────────────────────────────────────────────────────────────────────────
//  Portal
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief A portal between two adjacent triangles in the A* corridor.
 *
 * Stores the two world-space endpoints of the shared edge, oriented so that:
 *   - `left`  is the vertex to the LEFT  of the direction of travel
 *   - `right` is the vertex to the RIGHT of the direction of travel
 *
 * "Direction of travel" = from the current triangle toward the next triangle.
 *
 * The FunnelAlgorithm builds the portal list from the navmesh arc data.
 */
struct Portal {
    geometry::Point2D left;   ///< Left portal vertex (CCW side)
    geometry::Point2D right;  ///< Right portal vertex (CW side)
};

// ─────────────────────────────────────────────────────────────────────────────
//  SmoothPath
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief The result of the Funnel Algorithm.
 */
struct SmoothPath {
    /// Ordered waypoints from start to goal (both included).
    /// These are the minimum number of points required to define the
    /// shortest path through the corridor.
    std::vector<geometry::Point2D> waypoints;

    /// Total Euclidean length of the path (sum of distances between consecutive
    /// waypoints).
    double totalLength{0.0};

    /**
     * @brief Number of "turn" waypoints — points where the path changes
     *        direction by more than @p turnThreshold radians.
     *
     * Start and goal are not counted as turns.
     * Filled by FunnelAlgorithm::smooth() using the threshold passed in.
     */
    std::size_t turnCount{0};
};

// ─────────────────────────────────────────────────────────────────────────────
//  FunnelAlgorithm
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Computes the shortest Euclidean path through a triangle corridor
 *        using the Funnel (Simple Stupid Funnel) Algorithm.
 *
 * ### Usage
 * @code
 *   pathfinding::FunnelAlgorithm funnel;
 *   auto path = funnel.smooth(start, goal, portals);
 *   std::cout << path.totalLength << "\n";
 * @endcode
 *
 * The object is stateless and may be reused.
 */
class FunnelAlgorithm {
public:
    FunnelAlgorithm() = default;

    /**
     * @brief Runs the Funnel Algorithm and returns the smooth path.
     *
     * ### Preconditions
     *
     *   - `start` is inside the first triangle of the corridor.
     *   - `goal`  is inside the last  triangle of the corridor.
     *   - `portals` is ordered from start-triangle to goal-triangle.
     *   - Each portal's (left, right) is oriented consistently CCW around
     *     the source triangle (use NavMesh::portalPoints()).
     *
     * ### Degenerate cases
     *
     *   - Empty portal list (start == goal triangle): returns [start, goal].
     *   - Single portal: straightforward funnel close.
     *   - start or goal coincides with a portal vertex: handled correctly.
     *
     * @param start    World-space start position.
     * @param goal     World-space goal  position.
     * @param portals  Sequence of portals between consecutive corridor triangles.
     *                 Length = corridor.size() − 1.
     * @param turnThreshold  Angle (radians) above which a waypoint is counted
     *                       as a "turn" for the metric. Default: π/12 (15°).
     * @return SmoothPath with waypoints and metrics.
     */
    [[nodiscard]] SmoothPath smooth(
        geometry::Point2D            start,
        geometry::Point2D            goal,
        const std::vector<Portal>&   portals,
        double                       turnThreshold = 3.14159265358979323846 / 12.0
    ) const;

    /**
     * @brief Convenience overload: builds portals from the A* node corridor
     *        and the NavMesh, then runs smooth().
     *
     * This is the typical entry point in production code — the caller
     * doesn't need to manually build the portal list.
     *
     * @param start    World-space start position.
     * @param goal     World-space goal  position.
     * @param corridor Node sequence from AStar::search().nodes.
     * @param nm       The NavMesh (for portal and centroid data).
     * @param dcel     The source DCEL (for vertex position lookup).
     * @param turnThreshold  Turn-counting threshold in radians.
     * @return SmoothPath.
     */
    [[nodiscard]] SmoothPath smooth(
        geometry::Point2D                        start,
        geometry::Point2D                        goal,
        const std::vector<navmesh::NodeIdx>&     corridor,
        const navmesh::NavMesh&                  nm,
        const geometry::DCEL&                    dcel,
        double                                   turnThreshold = 3.14159265358979323846 / 12.0
    ) const;

private:
    // ── Internal funnel state ─────────────────────────────────────────────────

    /**
     * @brief The funnel state maintained during portal processing.
     *
     * The funnel is a triangular region:
     *   apex → left chain tip  (left tangent)
     *   apex → right chain tip (right tangent)
     *
     * `leftChain` and `rightChain` are stacks (front = apex side).
     */
    struct Funnel {
        geometry::Point2D apex;                 ///< Current apex (last confirmed turn)
        std::vector<geometry::Point2D> left;    ///< Left tangent chain (apex at back)
        std::vector<geometry::Point2D> right;   ///< Right tangent chain (apex at back)
    };

    // ── Helpers ───────────────────────────────────────────────────────────────

    /**
     * @brief Builds the Portal list from a NavMesh node corridor.
     *
     * For each consecutive pair (corridor[i], corridor[i+1]), looks up the
     * arc from i to i+1 and retrieves the portal world-space points via
     * NavMesh::portalPoints().
     *
     * The left/right orientation is determined by the CCW face order of the
     * source triangle: the portal vertex that is "to the left" when walking
     * from centroid[i] toward centroid[i+1] is assigned to Portal::left.
     *
     * @param corridor  Node indices from A*.
     * @param nm        NavMesh for arc/portal lookup.
     * @param dcel      DCEL for vertex positions.
     * @return          Ordered portal sequence (size = corridor.size() − 1).
     */
    [[nodiscard]] static std::vector<Portal>
    buildPortals(const std::vector<navmesh::NodeIdx>& corridor,
                 const navmesh::NavMesh&              nm,
                 const geometry::DCEL&               dcel);

    /**
     * @brief Computes the total length and turn count of the given waypoint list.
     *
     * @param waypoints       Ordered list of path points.
     * @param turnThreshold   Minimum angle (radians) to count a turn.
     * @param[out] length     Receives the total Euclidean path length.
     * @param[out] turns      Receives the number of turns.
     */
    static void computeMetrics(const std::vector<geometry::Point2D>& waypoints,
                               double                                 turnThreshold,
                               double&                                length,
                               std::size_t&                           turns) noexcept;

    /**
     * @brief Returns the 2D cross product of vectors (apex→a) and (apex→b).
     *
     * Sign convention:
     *   > 0  → b is to the LEFT  of (apex→a)  (CCW)
     *   < 0  → b is to the RIGHT of (apex→a)  (CW)
     *   = 0  → collinear
     *
     * Used to decide whether a new portal vertex widens or tightens the funnel.
     *
     * @param apex  Funnel apex.
     * @param a     Existing tangent chain tip.
     * @param b     New portal vertex to test.
     */
    [[nodiscard]] static double cross2D(const geometry::Point2D& apex,
                                        const geometry::Point2D& a,
                                        const geometry::Point2D& b) noexcept;
};

} // namespace pathfinding
