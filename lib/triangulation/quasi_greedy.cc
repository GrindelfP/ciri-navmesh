/**
 * @file quasi_greedy.cc
 * @brief Quasi-greedy MWT triangulation implementation.
 */

#include "quasi_greedy.hh"

#include "../geometry/predicates.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

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
    //  Local helpers (file-scope, not exposed in header)
    // ─────────────────────────────────────────────────────────────────────────────

    namespace {
        /// Undirected edge key: always stores (lo, hi) with lo < hi.
        struct EdgeKey {
            std::size_t lo, hi;

            bool operator==(const EdgeKey &o) const noexcept {
                return lo == o.lo && hi == o.hi;
            }
        };

        struct EdgeKeyHash {
            std::size_t operator()(const EdgeKey &k) const noexcept {
                return std::hash<std::size_t>{}(k.lo) ^
                       (std::hash<std::size_t>{}(k.hi) * 0x9e3779b9u);
            }
        };

        inline EdgeKey makeKey(std::size_t a, std::size_t b) noexcept {
            return {std::min(a, b), std::max(a, b)};
        }

        /// Returns true if segment (points[u], points[v]) properly intersects
        /// any edge in `accepted` that does not share an endpoint with (u, v).
        template<typename EdgeVec>
        bool properlyIntersectsAny(
            std::size_t u,
            std::size_t v,
            const std::vector<Point2D> &points,
            const EdgeVec &accepted) noexcept {
            const Point2D &pu = points[u];
            const Point2D &pv = points[v];
            for (const auto &e: accepted) {
                if (e.u == u || e.u == v || e.v == u || e.v == v) continue;
                auto t = pred::segmentIntersect(pu, pv, points[e.u], points[e.v]);
                if (t == pred::IntersectionType::Proper ||
                    t == pred::IntersectionType::Overlap)
                    return true;
            }
            return false;
        }
    } // namespace

    // ─────────────────────────────────────────────────────────────────────────────
    //  Public entry point
    // ─────────────────────────────────────────────────────────────────────────────

    TriangulationResult
    QuasiGreedyTriangulator::triangulate(const std::vector<Point2D> &points,
                                         DCEL &dcel) {
        if (points.size() < 3) {
            throw std::invalid_argument(
                "QuasiGreedyTriangulator: need at least 3 points");
        }

        dcel.clear();

        auto t0 = std::chrono::steady_clock::now();

        // ── Phase 1: Delaunay triangulation ──────────────────────────────────────
        //
        // We use a temporary DCEL to hold the Delaunay result, then extract
        // its edges as candidates.  The output `dcel` is populated in Phase 4/5.

        DCEL delaunayDcel;
        DelaunayTriangulator dt;
        dt.triangulate(points, delaunayDcel);

        // Extract Delaunay edges as (u, v, lenSq) triples.
        // We iterate over every other half-edge (h < twin) to visit each
        // undirected edge exactly once.
        //
        // Important: vertex indices in delaunayDcel correspond to the
        // insertion order used by Bowyer–Watson, which pre-sorts points by x.
        // We need to map DCEL vertex positions back to indices in `points`.

        // Build position → index map for the original `points` array.
        std::unordered_map<std::size_t, // hash of packed (x,y)
                    std::size_t> // index in points[]
                posToIdx;
        posToIdx.reserve(points.size() * 2);

        // Use a simple spatial map: key = index into points[], value confirmed
        // by coordinate comparison.  Since we pre-sorted in BW, we search linearly.
        // For correctness we do a full position-match lookup.
        auto findIdx = [&](const Point2D &pos) -> std::size_t {
            for (std::size_t i = 0; i < points.size(); ++i) {
                if (points[i].nearlyEqual(pos)) return i;
            }
            return kInvalidIdx; // should not happen
        };

        std::vector<Edge> candidates;
        candidates.reserve(3 * points.size());

        for (HalfEdgeIdx hi = 0; hi < delaunayDcel.halfEdgeCount(); ++hi) {
            if (!delaunayDcel.isHalfEdgeLive(hi)) continue;
            const HalfEdge &he = delaunayDcel.halfEdge(hi);
            if (he.dead) continue;
            if (hi > he.twin) continue; // visit each undirected edge once

            // Skip boundary edges (one side is outer face).
            if (he.face == kOuterFace) continue;
            if (delaunayDcel.halfEdge(he.twin).face == kOuterFace) continue;

            const Point2D &pa = delaunayDcel.vertex(he.origin).pos;
            const Point2D &pb = delaunayDcel.vertex(
                delaunayDcel.halfEdge(he.twin).origin).pos;

            std::size_t ia = findIdx(pa);
            std::size_t ib = findIdx(pb);
            if (ia == kInvalidIdx || ib == kInvalidIdx) continue;

            if (ia > ib) std::swap(ia, ib);
            candidates.push_back({ia, ib, points[ia].distSq(points[ib])});
        }

        // Also include convex-hull edges (they are always in any triangulation).
        // Hull edges are those where twin->face == kOuterFace.
        for (HalfEdgeIdx hi = 0; hi < delaunayDcel.halfEdgeCount(); ++hi) {
            if (!delaunayDcel.isHalfEdgeLive(hi)) continue;
            const HalfEdge &he = delaunayDcel.halfEdge(hi);
            if (he.dead) continue;
            if (delaunayDcel.halfEdge(he.twin).face != kOuterFace) continue;

            const Point2D &pa = delaunayDcel.vertex(he.origin).pos;
            const Point2D &pb = delaunayDcel.vertex(
                delaunayDcel.halfEdge(he.twin).origin).pos;

            std::size_t ia = findIdx(pa);
            std::size_t ib = findIdx(pb);
            if (ia == kInvalidIdx || ib == kInvalidIdx) continue;
            if (ia > ib) std::swap(ia, ib);
            candidates.push_back({ia, ib, points[ia].distSq(points[ib])});
        }

        // De-duplicate candidates (BW may produce the same edge twice via hull).
        std::sort(candidates.begin(), candidates.end(),
                  [](const Edge &a, const Edge &b) {
                      if (a.u != b.u) return a.u < b.u;
                      return a.v < b.v;
                  });
        candidates.erase(
            std::unique(candidates.begin(), candidates.end(),
                        [](const Edge &a, const Edge &b) {
                            return a.u == b.u && a.v == b.v;
                        }),
            candidates.end());

        // ── Phase 2: lune test — mark safe (forced) edges ────────────────────────
        //
        // Sort candidates by length so safe edges are processed short-first.
        std::sort(candidates.begin(), candidates.end(),
                  [](const Edge &a, const Edge &b) { return a.lenSq < b.lenSq; });

        std::vector<Edge> safeEdges;
        safeEdges.reserve(candidates.size());

        for (const Edge &e: candidates) {
            if (emptyLune(e.u, e.v, points)) {
                safeEdges.push_back(e);
            }
        }

        // ── Phase 3: triangulate pockets ─────────────────────────────────────────
        std::vector<Edge> allEdges = safeEdges;
        triangulatePockets(points, safeEdges, allEdges);

        // ── Phase 3.5: Fallback для гарантии полной триангуляции ─────────────────
        // Если safeEdges несвязны, карманы могут не закрыть все пустоты.
        // Жадно добиваем оставшиеся ребра для формирования макс. планарного графа.
        std::vector<Edge> remainingCands;
        const std::size_t n = points.size();
        remainingCands.reserve(n * (n - 1) / 2);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                remainingCands.push_back({i, j, points[i].distSq(points[j])});
            }
        }
        std::sort(remainingCands.begin(), remainingCands.end(),
                  [](const Edge &a, const Edge &b) { return a.lenSq < b.lenSq; });

        std::unordered_set<EdgeKey, EdgeKeyHash> acceptedSet;
        for (const auto &e: allEdges) acceptedSet.insert(makeKey(e.u, e.v));

        for (const Edge &e: remainingCands) {
            if (acceptedSet.count(makeKey(e.u, e.v))) continue;
            if (!properlyIntersectsAny(e.u, e.v, points, allEdges)) {
                allEdges.push_back(e);
                acceptedSet.insert(makeKey(e.u, e.v));
            }
        }

        // ── Phase 4: build DCEL ───────────────────────────────────────────────────
        buildDCEL(points, allEdges, dcel);

        // ── Phase 5: weight-reducing edge flips ───────────────────────────────────

        std::size_t totalFlips = 0;
        if (doFlips_) {
            for (int pass = 0; pass < maxFlipPasses_; ++pass) {
                std::size_t flips = flipPass(dcel);
                totalFlips += flips;
                if (flips == 0) break; // converged
            }
        }

        auto t1 = std::chrono::steady_clock::now();

        TriangulationResult res;
        res.elapsed = t1 - t0;
        res.totalWeight = dcel.totalWeight();
        res.triangleCount = dcel.liveTriangleCount();
        res.flipCount = totalFlips;
        return res;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  Phase 2: emptyLune
    // ─────────────────────────────────────────────────────────────────────────────

    bool QuasiGreedyTriangulator::emptyLune(
        std::size_t u,
        std::size_t v,
        const std::vector<Point2D> &points) const noexcept {
        // The open disk with diameter (p, q):
        //   centre = (p + q) / 2,  radius = |pq| / 2.
        //
        // A point r is inside iff |r - centre| < radius
        //   ⟺  |r - centre|² < |pq|² / 4
        //   ⟺  dot(r - p, r - q) < 0          [classic inscribed-angle identity]
        //
        // The dot-product form is cheaper and numerically identical.

        const Point2D &p = points[u];
        const Point2D &q = points[v];

        for (std::size_t i = 0; i < points.size(); ++i) {
            if (i == u || i == v) continue;
            const Point2D &r = points[i];
            // dot(r - p, r - q) < 0  →  r is strictly inside the diameter circle.
            const double dot = (r.x - p.x) * (r.x - q.x)
                               + (r.y - p.y) * (r.y - q.y);
            if (dot < -geometry::kEps) return false; // lune is not empty
        }
        return true;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  Phase 3: triangulatePockets
    // ─────────────────────────────────────────────────────────────────────────────

    void QuasiGreedyTriangulator::triangulatePockets(
        const std::vector<Point2D> &points,
        const std::vector<Edge> &safeEdges,
        std::vector<Edge> &allEdges) const {
        const std::size_t n = points.size();

        // Build adjacency from safe edges.
        std::vector<std::vector<std::size_t> > adj(n);
        std::unordered_set<EdgeKey, EdgeKeyHash> safeSet;
        safeSet.reserve(safeEdges.size() * 2);

        for (const Edge &e: safeEdges) {
            adj[e.u].push_back(e.v);
            adj[e.v].push_back(e.u);
            safeSet.insert(makeKey(e.u, e.v));
        }

        // A "pocket" is a face of the PSLG formed by safe edges that has more
        // than 3 vertices.  We find pockets by looking for faces in the safe-edge
        // graph.
        //
        // Strategy: for each safe edge (u, v), the face to the left of u→v is
        // determined by walking: starting from v, find the neighbour of v that
        // makes the smallest CCW angle with the direction v←u, then continue.
        // This traces the boundary of the face.  Faces with > 3 vertices are pockets.
        //
        // Implementation: sort each adjacency list by angle around each vertex,
        // then walk faces using the "next half-edge" rule (next CCW neighbour).

        // For each vertex, sort neighbours by angle.
        std::vector<std::vector<std::size_t> > sortedAdj(n);
        for (std::size_t vi = 0; vi < n; ++vi) {
            sortedAdj[vi] = adj[vi];
            std::sort(sortedAdj[vi].begin(), sortedAdj[vi].end(),
                      [&](std::size_t a, std::size_t b) {
                          // Angle of edge vi→a vs vi→b.
                          double angA = std::atan2(points[a].y - points[vi].y,
                                                   points[a].x - points[vi].x);
                          double angB = std::atan2(points[b].y - points[vi].y,
                                                   points[b].x - points[vi].x);
                          return angA < angB;
                      });
        }

        // "Next" in the planar half-edge sense: given directed edge u→v,
        // the next CCW edge around v's face is v→w where w is the neighbour of v
        // that comes immediately AFTER u in v's CW-sorted adjacency
        // (i.e. the first neighbour clockwise from the direction v←u).
        //
        // Equivalently in our CCW-sorted list: the neighbour just BEFORE u
        // in the CCW order around v.
        auto nextHalfEdge = [&](std::size_t u, std::size_t v) -> std::size_t {
            const auto &nb = sortedAdj[v];
            auto it = std::find(nb.begin(), nb.end(), u);
            if (it == nb.end()) return kInvalidIdx;
            // Previous in CCW order = next in CW order = the element before it,
            // wrapping around.
            if (it == nb.begin()) return nb.back();
            return *std::prev(it);
        };

        // Trace all faces of the safe-edge PSLG.
        // Use a visited set on directed edges (u→v) to avoid re-tracing.
        std::unordered_set<std::size_t> visitedHE; // encode as u*n + v
        visitedHE.reserve(safeEdges.size() * 4);

        for (const Edge &startEdge: safeEdges) {
            // Try both directed edges (u→v) and (v→u).
            for (int dir = 0; dir < 2; ++dir) {
                std::size_t u0 = (dir == 0) ? startEdge.u : startEdge.v;
                std::size_t v0 = (dir == 0) ? startEdge.v : startEdge.u;

                std::size_t heKey = u0 * n + v0;
                if (visitedHE.count(heKey)) continue;

                // Trace the face starting from u0→v0.
                std::vector<std::size_t> face;
                face.reserve(8);
                std::size_t cu = u0, cv = v0;

                for (;;) {
                    std::size_t faceKey = cu * n + cv;
                    if (visitedHE.count(faceKey)) break; // cycle completed
                    visitedHE.insert(faceKey);
                    face.push_back(cu);

                    std::size_t nv = nextHalfEdge(cu, cv);
                    if (nv == kInvalidIdx || nv == u0) break;
                    cu = cv;
                    cv = nv;
                }

                // The outer face is traced CW (negative signed area) — skip it.
                // Compute signed area of the traced polygon.
                double signedArea = 0.0;
                const std::size_t m = face.size();
                for (std::size_t i = 0; i < m; ++i) {
                    const Point2D &pi = points[face[i]];
                    const Point2D &pj = points[face[(i + 1) % m]];
                    signedArea += pi.x * pj.y - pj.x * pi.y;
                }
                // CCW → positive signed area → interior face.
                if (signedArea <= 0.0) continue;
                // Already a triangle → no pocket.
                if (m <= 3) continue;

                // Pocket found: triangulate it.
                triangulatePolygon(face, points, allEdges, allEdges);
            }
        }
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  Phase 3: triangulatePolygon
    // ─────────────────────────────────────────────────────────────────────────────

    void QuasiGreedyTriangulator::triangulatePolygon(
        const std::vector<std::size_t> &poly,
        const std::vector<Point2D> &points,
        const std::vector<Edge> &existingEdges,
        std::vector<Edge> &out) const {
        const std::size_t m = poly.size();
        if (m < 3) return;
        if (m == 3) return; // triangle — already done, boundary edges are in safeEdges

        // Build set of already-existing edges (boundary + previously accepted).
        std::unordered_set<EdgeKey, EdgeKeyHash> existingSet;
        existingSet.reserve(existingEdges.size() * 2);
        for (const Edge &e: existingEdges) {
            existingSet.insert(makeKey(e.u, e.v));
        }

        // Generate all diagonals of this polygon (non-boundary edges).
        // A diagonal connects poly[i] and poly[j] where |i-j| > 1 (mod m).
        struct Diag {
            std::size_t i, j; // indices into poly[]
            double lenSq;
        };
        std::vector<Diag> diags;
        diags.reserve(m * (m - 3) / 2);

        for (std::size_t i = 0; i < m; ++i) {
            for (std::size_t j = i + 2; j < m; ++j) {
                if (i == 0 && j == m - 1) continue; // boundary edge
                std::size_t pi = poly[i], pj = poly[j];
                EdgeKey key = makeKey(pi, pj);
                if (existingSet.count(key)) continue; // already accepted
                diags.push_back({i, j, points[pi].distSq(points[pj])});
            }
        }

        std::sort(diags.begin(), diags.end(),
                  [](const Diag &a, const Diag &b) { return a.lenSq < b.lenSq; });

        // Greedily accept diagonals that don't cross existing edges.
        for (const Diag &d: diags) {
            std::size_t pi = poly[d.i], pj = poly[d.j];
            if (!properlyIntersectsAny(pi, pj, points, out)) {
                out.push_back({
                    std::min(pi, pj), std::max(pi, pj),
                    points[pi].distSq(points[pj])
                });
                existingSet.insert(makeKey(pi, pj));
            }
        }
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  Phase 4: flipPass
    // ─────────────────────────────────────────────────────────────────────────────

    std::size_t QuasiGreedyTriangulator::flipPass(DCEL &dcel) const {
        // Snapshot current half-edge count (stable indices — new edges appended).
        const std::size_t heCount = dcel.halfEdgeCount();
        std::size_t flips = 0;

        for (HalfEdgeIdx hi = 0; hi < heCount; ++hi) {
            if (!dcel.isHalfEdgeLive(hi)) continue;
            const HalfEdge &he = dcel.halfEdge(hi);
            // if (he.dead) continue;
            // Only process each undirected edge once.
            if (hi > he.twin) continue;
            // Skip boundary edges.
            if (he.face == kOuterFace) continue;
            if (dcel.halfEdge(he.twin).face == kOuterFace) continue;

            // Vertices of the quadrilateral formed by the two adjacent triangles:
            //   Triangle 1: u → v → w   (he: u→v, he->next->twin->origin = w)
            //   Triangle 2: v → u → x   (twin: v→u, twin->next->twin->origin = x)
            //
            // Current shared edge: u–v  (length = dist(u,v))
            // Candidate new edge:  w–x  (length = dist(w,x))
            // Flip is weight-reducing iff dist(w,x) < dist(u,v).

            VertexIdx vu = he.origin;
            VertexIdx vv = dcel.halfEdge(he.twin).origin;

            // ── w: third vertex of triangle 1 ──
            HalfEdgeIdx he_vw = dcel.halfEdge(hi).next; // v→w
            VertexIdx w_idx = dcel.halfEdge(
                dcel.halfEdge(he_vw).next).origin; // w→u, origin = w

            // ── x: third vertex of triangle 2 ──
            HalfEdgeIdx he_ux = dcel.halfEdge(he.twin).next; // u→x  (next after v→u)
            VertexIdx x_idx = dcel.halfEdge(
                dcel.halfEdge(he_ux).next).origin; // x→v, origin = x

            const Point2D &pu = dcel.vertex(vu).pos;
            const Point2D &pv = dcel.vertex(vv).pos;
            const Point2D &pw = dcel.vertex(w_idx).pos;
            const Point2D &px = dcel.vertex(x_idx).pos;

            // Weight-reducing condition: dist(w, x) < dist(u, v).
            // Compare squared distances to avoid sqrt.
            const double uvSq = pu.distSq(pv);
            const double wxSq = pw.distSq(px);

            if (wxSq < uvSq - geometry::kEps) {
                // Sanity check: flipping must produce a valid (non-degenerate)
                // quadrilateral.  The new edge w-x must lie inside the quad,
                // i.e. w and x must be on opposite sides of line u-v.
                auto owx = pred::orientation(pu, pv, pw);
                auto oxv = pred::orientation(pu, pv, px);
                if (owx == pred::Orientation::Collinear ||
                    oxv == pred::Orientation::Collinear)
                    continue;
                if (owx == oxv) continue; // same side — non-convex quad, skip

                dcel.flipEdge(hi);
                ++flips;
            }
        }
        return flips;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  buildDCEL  (mirrors GreedyTriangulator::buildDCEL)
    // ─────────────────────────────────────────────────────────────────────────────

    void QuasiGreedyTriangulator::buildDCEL(const std::vector<Point2D> &points,
                                            const std::vector<Edge> &edges,
                                            DCEL &dcel) const {
        const std::size_t n = points.size();
        dcel.reserve(n, 2 * n);

        for (const Point2D &p: points) {
            dcel.addVertex(p);
        }

        // Build adjacency list.
        std::vector<std::vector<std::size_t> > adj(n);
        std::unordered_set<EdgeKey, EdgeKeyHash> edgeSet;
        edgeSet.reserve(edges.size() * 2);

        for (const Edge &e: edges) {
            adj[e.u].push_back(e.v);
            adj[e.v].push_back(e.u);
            edgeSet.insert(makeKey(e.u, e.v));
        }

        auto hasEdge = [&](std::size_t a, std::size_t b) -> bool {
            return edgeSet.count(makeKey(a, b)) > 0;
        };

        auto isTriangleEmpty = [&](std::size_t a, std::size_t b, std::size_t c) {
            auto pa = points[a];
            auto pb = points[b];
            auto pc = points[c];

            if (!pred::isCCW(pa, pb, pc)) {
                std::swap(pb, pc);
            }

            auto cross = [](const Point2D& p1, const Point2D& p2, const Point2D& p3) {
                return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
            };

            constexpr double eps = -1e-9;

            for (std::size_t i = 0; i < n; ++i) {
                if (i == a || i == b || i == c) continue;

                if (cross(pa, pb, points[i]) >= eps &&
                    cross(pb, pc, points[i]) >= eps &&
                    cross(pc, pa, points[i]) >= eps) {
                    return false;
                    }
            }
            return true;
        };

        // Find all triangles (u < v < w with all three edges present).
        struct TriIdx {
            std::size_t a, b, c;
        };
        std::vector<TriIdx> triangles;
        triangles.reserve(2 * n);

        for (std::size_t u = 0; u < n; ++u) {
            for (std::size_t v: adj[u]) {
                if (v <= u) continue;
                for (std::size_t w: adj[u]) {
                    if (w <= v) continue;
                    if (!hasEdge(v, w)) continue;

                    if (isTriangleEmpty(u, v, w)) {
                        triangles.push_back({u, v, w});
                    }
                }
            }
        }

        // Directed edge → HalfEdgeIdx.  Key = src * n + dst.
        std::unordered_map<std::size_t, HalfEdgeIdx> dirHeMap;
        dirHeMap.reserve(edges.size() * 4);

        for (const TriIdx &tri: triangles) {
            const Point2D &pa = points[tri.a];
            const Point2D &pb = points[tri.b];
            const Point2D &pc = points[tri.c];

            VertexIdx va = static_cast<VertexIdx>(tri.a);
            VertexIdx vb = static_cast<VertexIdx>(tri.b);
            VertexIdx vc = static_cast<VertexIdx>(tri.c);

            if (!pred::isCCW(pa, pb, pc)) std::swap(vb, vc);

            FaceIdx fi = dcel.addTriangle(va, vb, vc);
            auto [h0, h1, h2] = dcel.faceHalfEdges(fi);
            // h0: va→vb, h1: vb→vc, h2: vc→va

            dirHeMap[static_cast<std::size_t>(va) * n +
                     static_cast<std::size_t>(vb)] = h0;
            dirHeMap[static_cast<std::size_t>(vb) * n +
                     static_cast<std::size_t>(vc)] = h1;
            dirHeMap[static_cast<std::size_t>(vc) * n +
                     static_cast<std::size_t>(va)] = h2;
        }

        // Stitch twins: for directed edge (u→v), twin is (v→u).
        std::unordered_set<HalfEdgeIdx> stitched;
        stitched.reserve(dirHeMap.size());

        for (auto &[key, h]: dirHeMap) {
            if (stitched.count(h)) continue;
            std::size_t src = key / n;
            std::size_t dst = key % n;
            std::size_t revKey = dst * n + src;
            auto it = dirHeMap.find(revKey);
            if (it != dirHeMap.end() && !stitched.count(it->second)) {
                dcel.stitchTwins(h, it->second);
                stitched.insert(h);
                stitched.insert(it->second);
            }
        }
    }
} // namespace triangulation
