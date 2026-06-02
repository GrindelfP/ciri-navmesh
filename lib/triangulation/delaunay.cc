/**
 * @file delaunay.cc
 * @brief Bowyer–Watson Delaunay triangulation implementation.
 */

#include "delaunay.hh"
#include <unordered_map>

#include "../geometry/predicates.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <unordered_set>
#include <stack>

namespace triangulation {

using geometry::DCEL;
using geometry::FaceIdx;
using geometry::HalfEdgeIdx;
using geometry::Point2D;
using geometry::VertexIdx;
using geometry::kInvalidIdx;
using geometry::kOuterFace;
namespace pred = geometry::predicates;

// ─────────────────────────────────────────────────────────────────────────────
//  Public entry point
// ─────────────────────────────────────────────────────────────────────────────

TriangulationResult
DelaunayTriangulator::triangulate(const std::vector<Point2D>& points, DCEL& dcel)
{
    if (points.size() < 3) {
        throw std::invalid_argument(
            "DelaunayTriangulator: need at least 3 points");
    }

    dcel.clear();

    // Reserve capacity: Euler formula gives ~2n triangles, ~6n half-edges.
    dcel.reserve(points.size() + 3, 2 * points.size() + 1);

    // Sort input by x (tie-break y) for better insertion locality.
    // We work on a local sorted copy; caller's order is unaffected.
    std::vector<Point2D> sorted(points.begin(), points.end());
    std::sort(sorted.begin(), sorted.end());
    // Remove consecutive exact duplicates (robustness).
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());

    if (sorted.size() < 3) {
        throw std::invalid_argument(
            "DelaunayTriangulator: fewer than 3 distinct points after dedup");
    }

    auto t0 = std::chrono::steady_clock::now();

    // Step 1: super-triangle.
    auto superVerts = buildSuperTriangle(sorted, dcel);

    // Step 2: incremental insertion.
    for (const Point2D& p : sorted) {
        insertPoint(p, dcel);
    }

    // Step 3: remove super-triangle vertices (kills all incident faces).
    // Remove in reverse order to avoid dangling references during removal.
    for (int i = 2; i >= 0; --i) {
        dcel.removeVertex(superVerts[static_cast<std::size_t>(i)]);
    }

    auto t1 = std::chrono::steady_clock::now();

    TriangulationResult res;
    res.elapsed       = t1 - t0;
    res.totalWeight   = dcel.totalWeight();
    res.triangleCount = dcel.liveTriangleCount();
    res.flipCount     = 0; // BW does no explicit flips; Delaunay is maintained
                           // structurally via cavity re-triangulation.
    return res;
}

// ─────────────────────────────────────────────────────────────────────────────
//  buildSuperTriangle
// ─────────────────────────────────────────────────────────────────────────────

std::array<VertexIdx, 3>
DelaunayTriangulator::buildSuperTriangle(std::span<const Point2D> points,
                                         DCEL&                    dcel)
{
    // Compute axis-aligned bounding box.
    double minX = points[0].x, maxX = points[0].x;
    double minY = points[0].y, maxY = points[0].y;
    for (const auto& p : points) {
        minX = std::min(minX, p.x);  maxX = std::max(maxX, p.x);
        minY = std::min(minY, p.y);  maxY = std::max(maxY, p.y);
    }

    // Expand bounding box by a generous margin so that the super-triangle
    // circumcircle test is never numerically ambiguous.
    const double dx     = maxX - minX;
    const double dy     = maxY - minY;
    const double margin = std::max({dx, dy, 1.0}) * 10.0;

    // Super-triangle: equilateral-ish, CCW order.
    //   v0 = bottom-left  (well below and left)
    //   v1 = bottom-right (well below and right)
    //   v2 = top-center   (well above)
    const double cx = (minX + maxX) * 0.5;
    const double cy = (minY + maxY) * 0.5;

    const Point2D sv0{cx - margin,     cy - margin};
    const Point2D sv1{cx + margin * 2, cy - margin};
    const Point2D sv2{cx,              cy + margin * 2};

    // Add vertices.
    VertexIdx vi0 = dcel.addVertex(sv0);
    VertexIdx vi1 = dcel.addVertex(sv1);
    VertexIdx vi2 = dcel.addVertex(sv2);

    // Add the single bounded face (triangle 0,1,2) and link all half-edges.
    dcel.addTriangle(vi0, vi1, vi2);

    return {vi0, vi1, vi2};
}

// ─────────────────────────────────────────────────────────────────────────────
//  insertPoint
// ─────────────────────────────────────────────────────────────────────────────

VertexIdx
DelaunayTriangulator::insertPoint(Point2D p, DCEL& dcel)
{
    // ── 1. Locate the face containing p ──────────────────────────────────────
    FaceIdx startFace = dcel.findFace(p);
    // findFace() returns kInvalidIdx if p is outside the triangulation.
    // This can happen if the super-triangle is too small (shouldn't occur with
    // our 10× margin) or if p is on a boundary edge of the convex hull.
    // In both cases we fall back to a linear scan of all live faces.
    if (startFace == kInvalidIdx) {
        // Fallback: find any face whose circumcircle contains p.
        for (FaceIdx fi = 1; fi < dcel.faceCount(); ++fi) {
            if (dcel.face(fi).dead) continue;
            auto [va, vb, vc] = dcel.faceVertices(fi);
            const Point2D& a = dcel.vertex(va).pos;
            const Point2D& b = dcel.vertex(vb).pos;
            const Point2D& c = dcel.vertex(vc).pos;
            if (pred::isInCircle(a, b, c, p)) {
                startFace = fi;
                break;
            }
        }
    }

    // If still not found, p might be exactly on a convex hull edge or an
    // existing vertex.  Insert it anyway into the closest face.
    if (startFace == kInvalidIdx) {
        // Last resort: just pick the first live non-outer face.
        for (FaceIdx fi = 1; fi < dcel.faceCount(); ++fi) {
            if (!dcel.face(fi).dead) { startFace = fi; break; }
        }
    }

    assert(startFace != kInvalidIdx && "No live face found for insertion");

    // ── 2. Find the cavity (all faces whose circumcircle contains p) ──────────
    std::vector<FaceIdx> cavity = findCavity(p, startFace, dcel);
    assert(!cavity.empty());

    // ── 3. Extract the ordered boundary of the cavity ─────────────────────────
    std::vector<HalfEdgeIdx> boundary = extractBoundary(cavity, dcel);

    // ── 4. Kill all cavity faces (marks them dead, frees their face slots) ────
    for (FaceIdx fi : cavity) {
        dcel.killFace(fi);
    }

    // ── 5. Add new vertex ─────────────────────────────────────────────────────
    VertexIdx newVi = dcel.addVertex(p);

    // ── 6. Re-triangulate the cavity ─────────────────────────────────────────
    //
    // Each boundary half-edge h belongs to a NON-cavity face, directed
    // h.origin → h.twin.origin (going along the cavity boundary CCW).
    //
    // We create a new triangle: (newVi, h.origin, h.twin.origin).
    // The edge h.origin → h.twin.origin already has a half-edge in the DCEL
    // (the twin of h, which was a cavity-face half-edge, now dead).  We reuse
    // that dead half-edge slot by re-linking, but it's simpler and safer to
    // just call addTriangle() and then stitch the outer edge to h.
    //
    // For each boundary edge we:
    //   a) Add a triangle (p, u, v) where u→v is the boundary edge direction.
    //   b) Stitch the new triangle's edge u→v to h (which is v→u direction
    //      from the non-cavity side).

    // Collect new face indices so we can stitch twins among them afterward.
    struct BoundaryTriangle {
        FaceIdx      faceIdx;
        HalfEdgeIdx  outerHE;   // half-edge on the boundary (twin of h)
        HalfEdgeIdx  boundaryH; // the original boundary half-edge h
    };
    std::vector<BoundaryTriangle> newTris;
    newTris.reserve(boundary.size());

    for (HalfEdgeIdx h : boundary) {
        const HalfEdge& he = dcel.halfEdge(h);
        VertexIdx u = he.origin;
        VertexIdx v = dcel.halfEdge(he.twin).origin;

        // addTriangle creates (newVi → u → v) in CCW order, returns face idx.
        FaceIdx fi = dcel.addTriangle(newVi, u, v);

        // The half-edge from u→v in the new triangle needs to be stitched
        // with h (which goes from v→u in the non-cavity adjacent face).
        // faceHalfEdges(fi) = {newVi→u, u→v, v→newVi}
        auto [h0, h1, h2] = dcel.faceHalfEdges(fi);
        // h1 is u→v; its twin should be h (v→u, the boundary edge).
        dcel.stitchTwins(h1, h);

        newTris.push_back({fi, h1, h});
    }

    // Stitch the "inner" edges of the new fan (the edges connecting p to the
    // cavity boundary vertices).  Adjacent new triangles share an edge:
    //   triangle[i] has edge  v→newVi  (i.e. h2 of triangle i)
    //   triangle[i+1] has edge newVi→v (i.e. h0 of triangle i+1)
    // where v is the shared vertex on the boundary.
    const std::size_t n = newTris.size();
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t j = (i + 1) % n;

        FaceIdx fi = newTris[i].faceIdx;
        FaceIdx fj = newTris[j].faceIdx;

        auto [fi_h0, fi_h1, fi_h2] = dcel.faceHalfEdges(fi); // newVi→u, u→v, v→newVi
        auto [fj_h0, fj_h1, fj_h2] = dcel.faceHalfEdges(fj); // newVi→u', u'→v', v'→newVi

        // fi_h2: v → newVi
        // fj_h0: newVi → v  (same v, shared boundary vertex)
        // These two are twins.
        dcel.stitchTwins(fi_h2, fj_h0);
    }

    return newVi;
}

// ─────────────────────────────────────────────────────────────────────────────
//  findCavity
// ─────────────────────────────────────────────────────────────────────────────

std::vector<FaceIdx>
DelaunayTriangulator::findCavity(Point2D            p,
                                  FaceIdx            startFace,
                                  const DCEL&        dcel) const
{
    std::vector<FaceIdx>              cavity;
    std::unordered_set<FaceIdx>       visited;
    std::stack<FaceIdx>               stack;

    // Check whether a face is "bad" (circumcircle contains p).
    auto isBad = [&](FaceIdx fi) -> bool {
        if (fi == kInvalidIdx || fi == kOuterFace) return false;
        if (dcel.face(fi).dead) return false;
        auto [va, vb, vc] = dcel.faceVertices(fi);
        const Point2D& a = dcel.vertex(va).pos;
        const Point2D& b = dcel.vertex(vb).pos;
        const Point2D& c = dcel.vertex(vc).pos;
        // isInCircle expects CCW order; DCEL guarantees CCW for bounded faces.
        return pred::isInCircle(a, b, c, p);
    };

    if (!isBad(startFace)) {
        // startFace doesn't contain p in its circumcircle — still need at
        // least one bad face.  The caller already confirmed p is inside
        // startFace geometrically, so it must be bad.  If not, it means p is
        // exactly on the circumcircle (OnCircle).  Treat as not bad; BW still
        // works correctly — we just produce a locally non-Delaunay edge for
        // co-circular inputs, which is acceptable.
        // Expand to neighbors anyway to find truly bad faces.
        visited.insert(startFace);
        auto neighbors = dcel.faceNeighbors(startFace);
        for (FaceIdx nb : neighbors) {
            if (nb != kInvalidIdx && isBad(nb)) {
                stack.push(nb);
                visited.insert(nb);
                cavity.push_back(nb);
            }
        }
    } else {
        stack.push(startFace);
        visited.insert(startFace);
        cavity.push_back(startFace);
    }

    while (!stack.empty()) {
        FaceIdx fi = stack.top(); stack.pop();
        auto neighbors = dcel.faceNeighbors(fi);
        for (FaceIdx nb : neighbors) {
            if (nb == kInvalidIdx) continue;
            if (visited.count(nb)) continue;
            visited.insert(nb);
            if (isBad(nb)) {
                cavity.push_back(nb);
                stack.push(nb);
            }
        }
    }

    return cavity;
}

// ─────────────────────────────────────────────────────────────────────────────
//  extractBoundary
// ─────────────────────────────────────────────────────────────────────────────

std::vector<HalfEdgeIdx>
DelaunayTriangulator::extractBoundary(const std::vector<FaceIdx>& cavity,
                                       const DCEL&                 dcel) const
{
    // Build a set of cavity face indices for O(1) lookup.
    std::unordered_set<FaceIdx> cavitySet(cavity.begin(), cavity.end());

    // A boundary half-edge h is in a cavity face, but h.twin is NOT in a
    // cavity face (it's in a live non-cavity face or on the outer boundary).
    // We collect the TWINS (they are the half-edges from non-cavity faces
    // pointing into the cavity).
    std::vector<HalfEdgeIdx> boundary;
    boundary.reserve(cavity.size() + 2);

    for (FaceIdx fi : cavity) {
        auto halfEdges = dcel.faceHalfEdges(fi);
        for (HalfEdgeIdx h : halfEdges) {
            const HalfEdge& he   = dcel.halfEdge(h);
            FaceIdx         twinFace = dcel.halfEdge(he.twin).face;
            if (!cavitySet.count(twinFace)) {
                // h is a boundary half-edge; collect its twin.
                boundary.push_back(he.twin);
            }
        }
    }

    // Order the boundary edges into a consistent CCW polygon.
    // Build a vertex → boundary-he map: for each boundary half-edge b,
    // b.origin is the "start" vertex.
    std::unordered_map<VertexIdx, HalfEdgeIdx> startMap;
    for (HalfEdgeIdx b : boundary) {
        startMap[dcel.halfEdge(b).origin] = b;
    }

    // Walk the chain: start from boundary[0], follow origin of each he's twin.
    std::vector<HalfEdgeIdx> ordered;
    ordered.reserve(boundary.size());
    HalfEdgeIdx cur = boundary[0];
    const std::size_t maxSteps = boundary.size() + 1; // guard against bad DCEL
    for (std::size_t step = 0; step < maxSteps && ordered.size() < boundary.size(); ++step) {
        ordered.push_back(cur);
        // The destination of cur (= origin of cur.twin) is the next start vertex.
        VertexIdx nextOrigin = dcel.halfEdge(dcel.halfEdge(cur).twin).origin;
        auto it = startMap.find(nextOrigin);
        if (it == startMap.end()) break; // shouldn't happen in a valid DCEL
        cur = it->second;
        if (cur == boundary[0]) break; // completed the loop
    }

    assert(ordered.size() == boundary.size() && "Cavity boundary is not a simple polygon");
    return ordered;
}

} // namespace triangulation
