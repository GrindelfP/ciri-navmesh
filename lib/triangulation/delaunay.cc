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
    using geometry::HalfEdge;
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
    DelaunayTriangulator::triangulate(const std::vector<Point2D> &points, DCEL &dcel) {
        if (points.size() < 3) {
            throw std::invalid_argument(
                "DelaunayTriangulator: need at least 3 points");
        }

        dcel.clear();

        // Reserve capacity: Euler formula gives ~2n triangles, ~6n half-edges.
        dcel.reserve(points.size() + 3, 2 * points.size() + 1);

        // Sort input by x (tie-break y) for better insertion locality.
        std::vector<Point2D> sorted(points.begin(), points.end());
        std::sort(sorted.begin(), sorted.end());
        // Remove consecutive exact duplicates (robustness).
        sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());

        if (sorted.size() < 3) {
            throw std::invalid_argument(
                "DelaunayTriangulator: fewer than 3 distinct points after dedup");
        }

        auto t0 = std::chrono::steady_clock::now();

        // Step 1: super-triangle
        auto superVerts = buildSuperTriangle(sorted, dcel);

        // Step 2: incremental insertion.
        for (const Point2D &p: sorted) {
            insertPoint(p, dcel);
        }

        // Step 3: Extract valid triangles and rebuild DCEL cleanly
        struct TriIdx { VertexIdx a, b, c; };
        std::vector<TriIdx> validTris;
        validTris.reserve(dcel.faceCount());

        std::vector<VertexIdx> sortedToOrig(sorted.size());
        for (std::size_t i = 0; i < sorted.size(); ++i) {
            for (std::size_t j = 0; j < points.size(); ++j) {
                if (sorted[i].nearlyEqual(points[j])) {
                    sortedToOrig[i] = static_cast<VertexIdx>(j);
                    break;
                }
            }
        }

        for (FaceIdx fi = 1; fi < dcel.faceCount(); ++fi) {
            if (dcel.face(fi).dead) continue;
            auto [va, vb, vc] = dcel.faceVertices(fi);

            if (va == superVerts[0] || va == superVerts[1] || va == superVerts[2] ||
                vb == superVerts[0] || vb == superVerts[1] || vb == superVerts[2] ||
                vc == superVerts[0] || vc == superVerts[1] || vc == superVerts[2]) {
                continue;
            }

            validTris.push_back({sortedToOrig[va - 3],
                                 sortedToOrig[vb - 3],
                                 sortedToOrig[vc - 3]});
        }

        dcel.clear();
        dcel.reserve(points.size(), validTris.size() * 3);

        for (const auto& p : points) {
            dcel.addVertex(p);
        }

        struct EdgeKey {
            std::size_t lo, hi;
            bool operator==(const EdgeKey& o) const noexcept { return lo == o.lo && hi == o.hi; }
        };
        struct EdgeKeyHash {
            std::size_t operator()(const EdgeKey& k) const noexcept {
                return std::hash<std::size_t>{}(k.lo) ^ (std::hash<std::size_t>{}(k.hi) * 0x9e3779b9u);
            }
        };

        std::unordered_map<EdgeKey, HalfEdgeIdx, EdgeKeyHash> heMap;
        heMap.reserve(validTris.size() * 3);

        for (const auto& tri : validTris) {
            FaceIdx fi = dcel.addTriangle(tri.a, tri.b, tri.c);
            auto [h0, h1, h2] = dcel.faceHalfEdges(fi);

            auto storeHE = [&](HalfEdgeIdx h, VertexIdx src, VertexIdx dst) {
                heMap[{static_cast<std::size_t>(src), static_cast<std::size_t>(dst)}] = h;
            };
            storeHE(h0, tri.a, tri.b);
            storeHE(h1, tri.b, tri.c);
            storeHE(h2, tri.c, tri.a);
        }

        std::unordered_set<HalfEdgeIdx> stitched;
        stitched.reserve(heMap.size());

        for (auto& [key, h] : heMap) {
            if (stitched.count(h)) continue;
            EdgeKey revKey{key.hi, key.lo};
            auto it = heMap.find(revKey);
            if (it != heMap.end()) {
                HalfEdgeIdx twinH = it->second;
                if (!stitched.count(twinH)) {
                    dcel.stitchTwins(h, twinH);
                    stitched.insert(h);
                    stitched.insert(twinH);
                }
            }
        }

        auto t1 = std::chrono::steady_clock::now();

        TriangulationResult res;
        res.elapsed = t1 - t0;
        res.totalWeight = dcel.totalWeight();
        res.triangleCount = dcel.liveTriangleCount();
        res.flipCount = 0;
        return res;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  buildSuperTriangle
    // ─────────────────────────────────────────────────────────────────────────────

    std::array<VertexIdx, 3>
    DelaunayTriangulator::buildSuperTriangle(std::vector<geometry::Point2D> &points,
                                             DCEL &dcel) {
        // Compute axis-aligned bounding box.
        double minX = points[0].x, maxX = points[0].x;
        double minY = points[0].y, maxY = points[0].y;
        for (const auto &p: points) {
            minX = std::min(minX, p.x);
            maxX = std::max(maxX, p.x);
            minY = std::min(minY, p.y);
            maxY = std::max(maxY, p.y);
        }

        // Expand bounding box by a generous margin so that the super-triangle
        // circumcircle test is never numerically ambiguous.
        const double dx = maxX - minX;
        const double dy = maxY - minY;
        const double margin = std::max({dx, dy, 1.0}) * 10.0;

        // Super-triangle: equilateral-ish, CCW order.
        //   v0 = bottom-left  (well below and left)
        //   v1 = bottom-right (well below and right)
        //   v2 = top-center   (well above)
        const double cx = (minX + maxX) * 0.5;
        const double cy = (minY + maxY) * 0.5;

        const Point2D sv0{cx - margin, cy - margin};
        const Point2D sv1{cx + margin * 2, cy - margin};
        const Point2D sv2{cx, cy + margin * 2};

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
    DelaunayTriangulator::insertPoint(Point2D p, DCEL &dcel) {
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
                const Point2D &a = dcel.vertex(va).pos;
                const Point2D &b = dcel.vertex(vb).pos;
                const Point2D &c = dcel.vertex(vc).pos;
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
                if (!dcel.face(fi).dead) {
                    startFace = fi;
                    break;
                }
            }
        }

        assert(startFace != kInvalidIdx && "No live face found for insertion");

        // ── 2. Find the cavity (all faces whose circumcircle contains p) ──────────
        std::vector<FaceIdx> cavity = findCavity(p, startFace, dcel);
        assert(!cavity.empty());

        // ── 3. Extract the ordered boundary of the cavity ─────────────────────────
        std::vector<HalfEdgeIdx> boundary = extractBoundary(cavity, dcel);

        struct BoundaryEdgeData {
            HalfEdgeIdx h;
            VertexIdx u;
            VertexIdx v;
        };
        std::vector<BoundaryEdgeData> bEdges;
        bEdges.reserve(boundary.size());
        for (HalfEdgeIdx h: boundary) {
            const HalfEdge &he = dcel.halfEdge(h);
            bEdges.push_back({
                h,
                he.origin,
                dcel.halfEdge(he.twin).origin
            });
        }

        // ── 4. Kill all cavity faces (marks them dead, frees their face slots) ────
        for (FaceIdx fi: cavity) {
            dcel.killFace(fi);
        }

        // ── 5. Add new vertex ─────────────────────────────────────────────────────
        VertexIdx newVi = dcel.addVertex(p);

        // ── 6. Re-triangulate the cavity ─────────────────────────────────────────
        struct BoundaryTriangle {
            FaceIdx faceIdx;
            HalfEdgeIdx outerHE;
            HalfEdgeIdx boundaryH;
        };
        std::vector<BoundaryTriangle> newTris;
        newTris.reserve(boundary.size());

        for (const auto &edgeData: bEdges) {
            VertexIdx u = edgeData.u;
            VertexIdx v = edgeData.v;
            HalfEdgeIdx h = edgeData.h;

            // (newVi, v, u) создает строго CCW-ориентированный треугольник!
            FaceIdx fi = dcel.addTriangle(newVi, v, u);

            auto [h0, h1, h2] = dcel.faceHalfEdges(fi);
            // h1 направлено как v->u, а h как u->v. Сшиваем их как честных близнецов.
            dcel.stitchTwins(h1, h);

            newTris.push_back({fi, h1, h});
        }

        const std::size_t n = newTris.size();
        std::size_t compStart = 0;

        for (std::size_t i = 0; i < n; ++i) {
            std::size_t j = i + 1;

            if (j == n || bEdges[j].u != bEdges[i].v) {
                j = compStart;
                compStart = i + 1;
            }

            FaceIdx fi = newTris[i].faceIdx;
            FaceIdx fj = newTris[j].faceIdx;

            auto [fi_h0, fi_h1, fi_h2] = dcel.faceHalfEdges(fi);
            auto [fj_h0, fj_h1, fj_h2] = dcel.faceHalfEdges(fj);

            dcel.stitchTwins(fi_h0, fj_h2);
        }

        return newVi;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  findCavity
    // ─────────────────────────────────────────────────────────────────────────────

    std::vector<FaceIdx>
    DelaunayTriangulator::findCavity(Point2D p,
                                     FaceIdx startFace,
                                     const DCEL &dcel) const {
        std::vector<FaceIdx> cavity;
        std::unordered_set<FaceIdx> visited;
        std::stack<FaceIdx> stack;

        // Check whether a face is "bad" (circumcircle contains p).
        auto isBad = [&](FaceIdx fi) -> bool {
            if (fi == kInvalidIdx || fi == kOuterFace) return false;
            if (dcel.face(fi).dead) return false;
            auto [va, vb, vc] = dcel.faceVertices(fi);
            const Point2D &a = dcel.vertex(va).pos;
            const Point2D &b = dcel.vertex(vb).pos;
            const Point2D &c = dcel.vertex(vc).pos;
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
            for (FaceIdx nb: neighbors) {
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
            FaceIdx fi = stack.top();
            stack.pop();
            auto neighbors = dcel.faceNeighbors(fi);
            for (FaceIdx nb: neighbors) {
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
    DelaunayTriangulator::extractBoundary(const std::vector<FaceIdx> &cavity,
                                          const DCEL &dcel) const {
        std::unordered_set<FaceIdx> cavitySet(cavity.begin(), cavity.end());
        std::vector<HalfEdgeIdx> boundary;
        boundary.reserve(cavity.size() + 2);

        for (FaceIdx fi: cavity) {
            auto halfEdges = dcel.faceHalfEdges(fi);
            for (HalfEdgeIdx h: halfEdges) {
                const HalfEdge &he = dcel.halfEdge(h);
                FaceIdx twinFace = dcel.halfEdge(he.twin).face;
                if (!cavitySet.count(twinFace)) {
                    boundary.push_back(he.twin);
                }
            }
        }

        // Build map: origin vertex -> boundary half-edge starting there.
        // In a valid simple-polygon boundary each vertex appears exactly once.
        std::unordered_map<VertexIdx, HalfEdgeIdx> startMap;
        startMap.reserve(boundary.size());
        for (HalfEdgeIdx b: boundary) {
            VertexIdx orig = dcel.halfEdge(b).origin;
            // Duplicate origin means the cavity boundary is not a simple polygon
            // (can happen with degenerate/collinear configurations). We still
            // insert the first occurrence so the walk can make partial progress,
            // but we detect and handle this below.
            startMap.emplace(orig, b);
        }

        std::vector<HalfEdgeIdx> ordered;
        ordered.reserve(boundary.size());

        // Walk ALL connected components of the boundary.
        // A well-formed cavity has exactly one component, but collinear /
        // cocircular inputs can produce a pinch-point boundary with two
        // components sharing a single vertex.  We collect all of them so
        // that the re-triangulation step covers the full cavity.
        std::unordered_set<HalfEdgeIdx> visited;
        visited.reserve(boundary.size());

        for (HalfEdgeIdx seed : boundary) {
            if (visited.count(seed)) continue;

            const VertexIdx seedOrigin = dcel.halfEdge(seed).origin;
            HalfEdgeIdx cur = seed;

            for (std::size_t step = 0; step < boundary.size(); ++step) {
                if (visited.count(cur)) break;
                visited.insert(cur);
                ordered.push_back(cur);

                VertexIdx nextOrigin =
                    dcel.halfEdge(dcel.halfEdge(cur).twin).origin;
                if (nextOrigin == seedOrigin) break; // closed this component

                auto it = startMap.find(nextOrigin);
                if (it == startMap.end()) break;    // open chain (degenerate)
                cur = it->second;
            }
        }

        // ordered must cover the full boundary; anything else is a geometry bug.
        assert(ordered.size() == boundary.size() &&
               "Cavity boundary is not a simple polygon");
        return ordered;
    }
} // namespace triangulation
