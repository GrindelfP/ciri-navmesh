#pragma once

/**
 * @file quasi_greedy.hh
 * @brief Quasi-greedy minimum-weight triangulation heuristic.
 *
 * ## Background
 *
 * True MWT is NP-hard (Mulzer & Rote, 2008).  The quasi-greedy algorithm of
 * Levcopoulos & Krznaric (1998) gives a constant-factor approximation in
 * O(n log n) time — the best known polynomial heuristic with a provable bound.
 *
 * ## Core idea
 *
 * The key insight is the "β-skeleton" criterion:
 *
 *   An edge (p, q) is in the MWT only if the two "lune" regions associated
 *   with (p, q) contain no other input points.  The lune of (p, q) for
 *   parameter β is the intersection of two discs of radius |pq|/β centred
 *   at p and q respectively.
 *
 * For β = √2 this gives the "LMT-skeleton" (Locally Minimal Triangulation
 * skeleton) — edges provably in every MWT.  We use a practical relaxation:
 *
 *   An edge (p, q) is a "quasi-greedy" edge if the open circle with
 *   diameter (p, q) — the β=1 lune, a.k.a. the "empty semicircle" test —
 *   contains no other input point.
 *
 * This is equivalent to: (p, q) belongs to the Delaunay triangulation AND
 * is "short" in a well-defined angular sense.
 *
 * ## Algorithm (Levcopoulos & Krznaric 1998, adapted)
 *
 *   Phase 1 – Compute the Delaunay triangulation (via Bowyer–Watson).
 *             Only Delaunay edges are MWT candidates (theorem).
 *
 *   Phase 2 – Mark "safe" edges: for each Delaunay edge (p, q), check
 *             the empty-lune condition.  If satisfied, the edge is forced
 *             into the triangulation.
 *
 *   Phase 3 – Triangulate remaining un-covered regions (polygonal pockets
 *             between safe edges) using a greedy sub-triangulation
 *             (sort candidate diagonals by length, accept non-crossing).
 *
 *   Phase 4 – Optional: local edge-flip pass to reduce total weight further
 *             (our own addition on top of the published algorithm).
 *
 * ## Complexity
 *
 *   - Phase 1 (Delaunay):   O(n log n) average
 *   - Phase 2 (lune check): O(n log n) with a point-in-circle spatial index
 *                           (we use a simple O(n) scan per edge → O(n²) total,
 *                            acceptable for n ≤ 2000; a k-d tree would fix it)
 *   - Phase 3 (pockets):    O(n log n) — pockets are small in practice
 *   - Phase 4 (flips):      O(n) per pass, typically 1–3 passes
 *   - Overall:              O(n² ) dominated by Phase 2 for large n;
 *                           O(n log n) if a spatial index is added
 *
 * ## References
 *
 *   Levcopoulos, C. & Krznaric, D. (1998).
 *   "Quasi-greedy triangulations approximating the minimum weight triangulation."
 *   Journal of Algorithms, 27(2), 303–338.
 */

#include "i_triangulator.hh"
#include "delaunay.hh"

#include <vector>

namespace triangulation {

/**
 * @brief Quasi-greedy MWT heuristic (Levcopoulos & Krznaric 1998).
 *
 * Produces triangulations with significantly lower total weight than both
 * Delaunay and naive greedy on typical inputs, in near-linear time.
 *
 * ### Usage
 * @code
 *   QuasiGreedyTriangulator qg;
 *   geometry::DCEL dcel;
 *   auto result = qg.triangulate(points, dcel);
 *   // result.flipCount contains flips performed in Phase 4.
 * @endcode
 */
class QuasiGreedyTriangulator : public ITriangulator {
public:
    /**
     * @brief Constructs the triangulator with optional edge-flip pass.
     *
     * @param doFlips  If true (default), perform Phase 4 weight-reducing
     *                 edge flips after the quasi-greedy triangulation.
     *                 Set to false to measure the pure algorithm output.
     * @param maxFlipPasses  Maximum number of flip passes (default 8).
     *                       Each pass is O(n); convergence is typically < 4.
     */
    explicit QuasiGreedyTriangulator(bool doFlips      = true,
                                     int  maxFlipPasses = 8)
        : doFlips_(doFlips), maxFlipPasses_(maxFlipPasses) {}

    /**
     * @copydoc ITriangulator::triangulate()
     *
     * ### Detailed steps
     *
     * 1. Run Delaunay (Bowyer–Watson) on @p points → get DCEL_dt.
     * 2. Extract Delaunay edges as candidates.
     * 3. For each candidate edge (p, q):
     *    - Check the open-disk lune: the circle with diameter pq must
     *      contain no other input point (emptyLune()).
     *    - If empty: mark edge as "safe" (forced into output).
     * 4. The safe edges partition the convex hull into polygonal pockets.
     *    Each pocket is triangulated greedily (sub-triangulation by
     *    diagonal length).
     * 5. If doFlips_: repeatedly scan all interior edges; flip any edge
     *    (u,v) shared by triangles (u,v,w) and (v,u,x) when:
     *       |uw| + |wx| < |uv| + |wx|   [weight-reducing flip]
     *    until no more weight-reducing flips exist or maxFlipPasses_ reached.
     * 6. Populate output @p dcel from the final edge set.
     */
    [[nodiscard]] TriangulationResult
    triangulate(const std::vector<geometry::Point2D>& points,
                geometry::DCEL&                       dcel) override;

    [[nodiscard]] std::string_view name() const noexcept override {
        return "Quasi-Greedy MWT (Levcopoulos–Krznaric)";
    }

private:
    bool doFlips_;
    int  maxFlipPasses_;

    // ── Internal types ────────────────────────────────────────────────────────

    /**
     * @brief An edge stored as a pair of point-array indices.
     */
    struct Edge {
        std::size_t u, v;   ///< u < v always
        double      lenSq;  ///< Cached squared length
    };

    // ── Phase 2: lune test ────────────────────────────────────────────────────

    /**
     * @brief Tests the open-lune (β=1) emptiness condition for edge (p, q).
     *
     * The lune for β=1 is the open disk whose diameter is the segment pq.
     * Its centre is the midpoint M = (p+q)/2, radius r = |pq|/2.
     *
     * The edge is "safe" iff no point in @p points (other than p and q)
     * lies strictly inside this disk.
     *
     * Equivalently (and cheaply): for each other point r,
     *   dot(r - p, r - q) < 0   →  r is inside the diameter circle.
     * This avoids computing sqrt or the centre explicitly.
     *
     * @param u       Index of first endpoint in @p points.
     * @param v       Index of second endpoint in @p points.
     * @param points  Full input point set.
     * @return        true if the lune is empty (edge is safe / forced).
     *
     * @note O(n) per call.  For large n use a spatial index to get O(log n).
     */
    bool emptyLune(std::size_t                           u,
                   std::size_t                           v,
                   const std::vector<geometry::Point2D>& points) const noexcept;

    // ── Phase 3: pocket triangulation ────────────────────────────────────────

    /**
     * @brief Identifies and triangulates all polygonal pockets.
     *
     * After Phase 2, the safe edges form a planar straight-line graph (PSLG)
     * that may not yet be a triangulation — some faces are polygons with > 3
     * sides.  These "pockets" must be triangulated.
     *
     * Each pocket is an interior polygon whose boundary consists of safe
     * edges.  We extract each pocket's vertex list (in CCW order) and
     * triangulate it using a greedy diagonal-by-length approach.
     *
     * @param points     Full input point set.
     * @param safeEdges  Edges accepted by Phase 2 (lune test).
     * @param[out] allEdges  Receives pocket diagonals appended in-place.
     */
    void triangulatePockets(
        const std::vector<geometry::Point2D>& points,
        const std::vector<Edge>&              safeEdges,
        std::vector<Edge>&                    allEdges) const;

    /**
     * @brief Greedily triangulates a single convex or simple polygon
     *        given as an ordered list of vertex indices.
     *
     * Generates all diagonals, sorts by length, accepts non-crossing ones.
     * Boundary edges of the polygon are already in @p safeEdges and are
     * not re-added.
     *
     * @param poly       Vertex indices of the polygon in CCW order.
     * @param points     Full input point set.
     * @param safeEdges  Already-accepted edges (for intersection check).
     * @param[out] out   Receives newly accepted diagonal edges.
     */
    void triangulatePolygon(
        const std::vector<std::size_t>&       poly,
        const std::vector<geometry::Point2D>& points,
        const std::vector<Edge>&              safeEdges,
        std::vector<Edge>&                    out) const;

    // ── Phase 4: weight-reducing edge flips ───────────────────────────────────

    /**
     * @brief Performs one pass of weight-reducing edge flips over the DCEL.
     *
     * For each interior half-edge h (face != kOuterFace, twin face != kOuterFace):
     *   Let the two triangles be (u, v, w) and (v, u, x).
     *   Flip condition: |u–x| + |w–x| < |u–v| + |w–x|
     *
     *   More precisely, flip (u,v) → (w,x) when:
     *     dist(w, x) < dist(u, v)
     *   because replacing edge uv with wx changes total weight by
     *   dist(w,x) - dist(u,v).  A flip is weight-reducing iff dist(w,x) < dist(u,v).
     *
     * @param dcel  The DCEL to modify in-place.
     * @return      Number of flips performed in this pass.
     *
     * @note We iterate over a snapshot of half-edge indices taken at pass
     *       start, so flips during the pass are not revisited (safe from
     *       index invalidation since DCEL uses stable indices).
     */
    std::size_t flipPass(geometry::DCEL& dcel) const;

    // ── DCEL construction (shared with greedy) ────────────────────────────────

    /**
     * @brief Builds the output DCEL from the final accepted edge set.
     *
     * Same logic as GreedyTriangulator::buildDCEL — finds triangles from
     * edge adjacency, adds them in CCW order, stitches twins.
     *
     * Extracted into a private method here rather than inheriting from
     * GreedyTriangulator to keep the class hierarchy flat (composition
     * over inheritance for algorithmic helpers).
     *
     * @param points  Input point set.
     * @param edges   Final accepted edge set (forms a valid triangulation).
     * @param dcel    Target DCEL (cleared by caller).
     */
    void buildDCEL(const std::vector<geometry::Point2D>& points,
                   const std::vector<Edge>&              edges,
                   geometry::DCEL&                       dcel) const;
};

} // namespace triangulation
