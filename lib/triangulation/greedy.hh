#pragma once

/**
 * @file greedy.hh
 * @brief Greedy minimum-weight triangulation heuristic.
 *
 * ## Algorithm overview
 *
 * The greedy triangulation is a classical baseline for MWT approximation:
 *
 *   1. Generate all C(n,2) candidate edges between input points.
 *   2. Sort them by Euclidean length (ascending).
 *   3. Greedily accept each edge if:
 *      a) It does not properly intersect any already-accepted edge, AND
 *      b) Adding it does not violate the triangulation topology
 *         (i.e. it lies inside the convex hull and does not over-triangulate
 *         an already-fully-triangulated region).
 *   4. The accepted edge set is guaranteed to form a triangulation of the
 *      convex hull of the input (by a result of Chazelle et al.).
 *   5. Build the DCEL from the accepted edge set.
 *
 * ## Complexity
 *
 *   - Candidate generation: O(n²)
 *   - Sort:                 O(n² log n)
 *   - Intersection checks:  O(n² · k) where k = accepted edges so far
 *     In practice k = O(n), so this is O(n³) worst case.
 *     For the project's n ≤ ~1000 this is acceptable.
 *   - DCEL build:           O(n log n)
 *   - Overall:              O(n² log n) dominated by sort (practical n ≤ 1000)
 *
 * ## Notes
 *
 * - "Proper intersection" means the segments cross at an interior point
 *   (IntersectionType::Proper). Shared endpoints (Endpoint) are fine.
 * - Collinear input points produce degenerate triangles; caller should
 *   pre-filter or accept that some triangles may be near-degenerate.
 * - The result is NOT necessarily Delaunay — it typically has lower total
 *   weight than Delaunay for uniform random point sets, but is not optimal.
 */

#include "i_triangulator.hh"
#include "../geometry/primitives.hpp"

#include <vector>

namespace triangulation {

/**
 * @brief Greedy minimum-weight triangulation.
 *
 * Generates all candidate edges, sorts by length, and greedily adds
 * non-crossing edges until a full triangulation is formed.
 *
 * ### Usage
 * @code
 *   GreedyTriangulator gt;
 *   geometry::DCEL dcel;
 *   auto result = gt.triangulate(points, dcel);
 * @endcode
 */
class GreedyTriangulator : public ITriangulator {
public:
    GreedyTriangulator() = default;

    /**
     * @copydoc ITriangulator::triangulate()
     *
     * ### Implementation steps
     *
     * 1. Build candidate list: all pairs (i, j), i < j, stored as index pairs
     *    with precomputed squared length for cheap sorting.
     * 2. Sort by squared length (avoids sqrt).
     * 3. For each candidate (u, v):
     *    - Check proper intersection against all accepted edges.
     *    - If clear: accept the edge.
     * 4. From the accepted edge set, extract triangular faces using a
     *    half-edge construction pass.
     * 5. Wrap result into the DCEL.
     *
     * @note The DCEL is built in a separate pass after all edges are accepted,
     *       rather than incrementally, because the greedy algorithm works on
     *       an abstract edge set and does not need DCEL topology during
     *       the acceptance phase.
     */
    [[nodiscard]] TriangulationResult
    triangulate(const std::vector<geometry::Point2D>& points,
                geometry::DCEL&                       dcel) override;

    [[nodiscard]] std::string_view name() const noexcept override {
        return "Greedy MWT";
    }

private:
    // ── Internal types ────────────────────────────────────────────────────────

    /**
     * @brief A candidate or accepted edge, stored as point indices.
     */
    struct Edge {
        std::size_t u;      ///< Index into the input point array (u < v)
        std::size_t v;      ///< Index into the input point array
        double      lenSq;  ///< Squared Euclidean length (for sorting)
    };

    // ── Internal helpers ──────────────────────────────────────────────────────

    /**
     * @brief Generates and returns all C(n,2) candidate edges, sorted by
     *        squared length ascending.
     *
     * @param points  Input point set.
     * @return        Sorted candidate edge list.
     */
    std::vector<Edge>
    buildCandidates(const std::vector<geometry::Point2D>& points) const;



    /**
     * @brief Tests whether segment (points[u], points[v]) properly intersects
     *        any edge in @p accepted.
     *
     * "Proper" means IntersectionType::Proper only — shared endpoints are
     * allowed (IntersectionType::Endpoint is NOT a conflict).
     *
     * @param u        First endpoint index of the candidate edge.
     * @param v        Second endpoint index of the candidate edge.
     * @param points   Input point set.
     * @param accepted Already-accepted edges.
     * @return         true if a proper intersection exists (reject candidate).
     */
    template<class EdgeVec>
    bool properlyIntersectsAny(std::size_t u, std::size_t v, const std::vector<geometry::Point2D> &points,
                               const EdgeVec &accepted) noexcept;

    /**
     * @brief Builds a DCEL from the accepted triangulation edge set.
     *
     * Identifies all triangular faces from the edge adjacency, then
     * populates the DCEL using addVertex / addTriangle / stitchTwins.
     *
     * @param points   Input point set.
     * @param edges    Accepted edges that form a valid triangulation.
     * @param dcel     Target DCEL (must be cleared before call).
     */
    void buildDCEL(const std::vector<geometry::Point2D>& points,
                   const std::vector<Edge>&              edges,
                   geometry::DCEL&                       dcel) const;
};

} // namespace triangulation
