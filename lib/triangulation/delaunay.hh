#pragma once

/**
 * @file delaunay.hh
 * @brief Delaunay triangulation via the Bowyer–Watson algorithm.
 *
 * ## Algorithm overview
 *
 * Bowyer–Watson is an incremental point-insertion algorithm:
 *
 *   1. Create a "super-triangle" that contains all input points.
 *   2. For each input point p:
 *      a. Find all triangles whose circumcircle contains p  → "bad triangles".
 *      b. Find the boundary of the cavity formed by bad triangles
 *         (the polygonal hole).
 *      c. Remove bad triangles from the DCEL.
 *      d. Re-triangulate the cavity by connecting p to every boundary edge.
 *   3. Remove all triangles that share a vertex with the super-triangle.
 *
 * The result satisfies the Delaunay condition: no point lies strictly inside
 * the circumcircle of any triangle.
 *
 * ## Complexity
 *
 *   - Worst case:  O(n²)  — degenerate point sets (e.g. grid, circle)
 *   - Average case: O(n log n) — random point sets
 *   - Space: O(n)
 *
 * ## Implementation notes
 *
 * The super-triangle uses a bounding box expanded by a safety margin.  Its
 * vertices are removed via DCEL::removeVertex() at step 3, which also kills
 * all incident faces.
 *
 * The "bad triangle" flood-fill uses a simple stack-based DFS from the
 * triangle containing p.  This is cache-friendly and avoids the O(n) linear
 * scan of the naive approach.
 *
 * For the Delaunay property we rely on predicates::inCircle() which uses
 * the standard 3×3 determinant formulation.
 */

#include "i_triangulator.hh"

namespace triangulation {

/**
 * @brief Delaunay triangulation using the Bowyer–Watson incremental algorithm.
 *
 * ### Usage
 * @code
 *   DelaunayTriangulator dt;
 *   geometry::DCEL dcel;
 *   auto result = dt.triangulate(points, dcel);
 *   std::cout << result.totalWeight << "\n";
 * @endcode
 *
 * The object is stateless between calls and may be reused freely.
 */
class DelaunayTriangulator : public ITriangulator {
public:
    DelaunayTriangulator() = default;

    /**
     * @copydoc ITriangulator::triangulate()
     *
     * ### Steps (detailed)
     *
     * **Pre-sort**: input points are sorted by x (tie-break y) to improve
     * circumcircle-search locality.  The original order is not preserved.
     *
     * **Super-triangle**: vertices at indices 0, 1, 2 of the DCEL are the
     * super-triangle.  All real vertices start at index 3.
     *
     * **Cavity flood-fill**: starting from the face found by DCEL::findFace(),
     * a DFS/BFS expands across shared edges while the inCircle predicate holds.
     * This is O(k) per insertion where k is the size of the cavity.
     *
     * **Cleanup**: DCEL::removeVertex() is called for the three
     * super-triangle vertices, which cascades to kill all incident faces.
     *
     * @note Points that lie exactly on a circumcircle (co-circular)
     *       are treated as "outside" (OnCircle → do not remove the triangle).
     *       This may produce non-unique triangulations for co-circular inputs
     *       but is numerically safe.
     */
    TriangulationResult
    triangulate(const std::vector<geometry::Point2D> &points,
                geometry::DCEL &dcel) override;

    [[nodiscard]] std::string_view name() const noexcept override {
        return "Delaunay (Bowyer–Watson)";
    }

private:
    // ── Internal helpers ──────────────────────────────────────────────────────

    /**
     * @brief Computes a bounding box of the point set and returns a
     *        super-triangle that strictly contains all points.
     *
     * The super-triangle is added to @p dcel as face 1 (bounded face index 1),
     * with vertices at indices 0, 1, 2.
     *
     * @param points  Input point set (must be non-empty).
     * @param dcel    Target DCEL (must be freshly cleared).
     * @return        Array of the three super-triangle vertex indices {0, 1, 2}.
     */
    std::array<geometry::VertexIdx, 3>
    buildSuperTriangle(std::vector<geometry::Point2D> &points,
                       geometry::DCEL &dcel);

    /**
     * @brief Inserts a single point into the current triangulation.
     *
     * Performs the full Bowyer–Watson step for point @p p:
     *   1. Find the cavity (set of faces whose circumcircle contains p).
     *   2. Find the cavity boundary edges (edges shared with non-bad faces).
     *   3. Kill all cavity faces.
     *   4. Connect p to each boundary edge, creating new triangles.
     *
     * @param p       The point to insert.
     * @param dcel    The DCEL to modify.
     * @return        VertexIdx of the newly inserted vertex.
     */
    geometry::VertexIdx insertPoint(geometry::Point2D p, geometry::DCEL& dcel);

    /**
     * @brief Flood-fills the "bad triangle" cavity around point @p p.
     *
     * Starting from @p startFace (which must already contain p in its
     * circumcircle), expands via adjacency.
     *
     * @param p          The query point.
     * @param startFace  Seed face — must satisfy inCircle(p).
     * @param dcel       The DCEL (read-only during this step).
     * @return           Set of face indices that form the cavity.
     */
    std::vector<geometry::FaceIdx>
    findCavity(geometry::Point2D   p,
               geometry::FaceIdx   startFace,
               const geometry::DCEL& dcel) const;

    /**
     * @brief Extracts the ordered boundary half-edges of the cavity.
     *
     * A boundary half-edge is one whose twin belongs to a non-bad face
     * (i.e. a live face not in the cavity set).
     *
     * The returned half-edges are the twins (belonging to non-cavity faces),
     * ordered to form a CCW polygon around the cavity.  Connecting p to the
     * origin/destination of each twin re-triangulates the hole correctly.
     *
     * @param cavity  Output of findCavity() — set of bad face indices.
     * @param dcel    The DCEL.
     * @return        Ordered list of boundary half-edge indices
     *                (from non-cavity faces, pointing INTO the cavity).
     */
    std::vector<geometry::HalfEdgeIdx>
    extractBoundary(const std::vector<geometry::FaceIdx>& cavity,
                    const geometry::DCEL&                 dcel) const;
};

} // namespace triangulation
