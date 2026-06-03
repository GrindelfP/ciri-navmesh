/**
 * @file greedy.cc
 * @brief Greedy minimum-weight triangulation implementation.
 */

#include "greedy.hh"

#include "../geometry/predicates.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace triangulation {

using geometry::DCEL;
using geometry::FaceIdx;
using geometry::HalfEdgeIdx;
using geometry::Point2D;
using geometry::VertexIdx;
using geometry::kInvalidIdx;
namespace pred = geometry::predicates;

// ─────────────────────────────────────────────────────────────────────────────
//  Public entry point
// ─────────────────────────────────────────────────────────────────────────────

TriangulationResult
GreedyTriangulator::triangulate(const std::vector<Point2D>& points,
                                 DCEL&                       dcel)
{
    if (points.size() < 3) {
        throw std::invalid_argument(
            "GreedyTriangulator: need at least 3 points");
    }

    dcel.clear();

    auto t0 = std::chrono::steady_clock::now();

    // Phase 1: generate and sort all candidate edges.
    std::vector<Edge> candidates = buildCandidates(points);

    // Phase 2: greedy acceptance.
    std::vector<Edge> accepted;
    accepted.reserve(3 * points.size()); // ~3n edges in a triangulation

    for (const Edge& e : candidates) {
        if (!properlyIntersectsAny(e.u, e.v, points, accepted)) {
            accepted.push_back(e);
        }
    }

    // Phase 3: build DCEL from accepted edges.
    buildDCEL(points, accepted, dcel);

    auto t1 = std::chrono::steady_clock::now();

    TriangulationResult res;
    res.elapsed       = t1 - t0;
    res.totalWeight   = dcel.totalWeight();
    res.triangleCount = dcel.liveTriangleCount();
    res.flipCount     = 0;
    return res;
}

// ─────────────────────────────────────────────────────────────────────────────
//  buildCandidates
// ─────────────────────────────────────────────────────────────────────────────

std::vector<GreedyTriangulator::Edge>
GreedyTriangulator::buildCandidates(const std::vector<Point2D>& points) const
{
    const std::size_t n = points.size();
    std::vector<Edge> cands;
    cands.reserve(n * (n - 1) / 2);

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            cands.push_back({i, j, points[i].distSq(points[j])});
        }
    }

    std::sort(cands.begin(), cands.end(),
              [](const Edge& a, const Edge& b) { return a.lenSq < b.lenSq; });

    return cands;
}

// ─────────────────────────────────────────────────────────────────────────────
//  properlyIntersectsAny
// ─────────────────────────────────────────────────────────────────────────────

bool GreedyTriangulator::properlyIntersectsAny(
    std::size_t                 u,
    std::size_t                 v,
    const std::vector<Point2D>& points,
    const std::vector<Edge>&    accepted) const
{
    const Point2D& pu = points[u];
    const Point2D& pv = points[v];

    for (const Edge& e : accepted) {
        // Shared endpoint → not a proper intersection; skip.
        if (e.u == u || e.u == v || e.v == u || e.v == v) continue;

        auto type = pred::segmentIntersect(pu, pv, points[e.u], points[e.v]);
        if (type == pred::IntersectionType::Proper) return true;
        // Overlap: two collinear edges — also reject.
        if (type == pred::IntersectionType::Overlap) return true;
    }
    return false;
}

// ─────────────────────────────────────────────────────────────────────────────
//  buildDCEL
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Key type for an undirected edge used in adjacency lookup.
 *
 * Stores (min, max) so that (u,v) and (v,u) map to the same key.
 */
struct EdgeKey {
    std::size_t lo, hi;
    bool operator==(const EdgeKey& o) const noexcept {
        return lo == o.lo && hi == o.hi;
    }
};

struct EdgeKeyHash {
    std::size_t operator()(const EdgeKey& k) const noexcept {
        return std::hash<std::size_t>{}(k.lo) ^
               (std::hash<std::size_t>{}(k.hi) * 0x9e3779b9u);
    }
};

void GreedyTriangulator::buildDCEL(const std::vector<Point2D>& points,
                                    const std::vector<Edge>&    edges,
                                    DCEL&                       dcel) const
{
    const std::size_t n = points.size();
    dcel.reserve(n, 2 * n);

    // ── Step 1: add all vertices ──────────────────────────────────────────────
    // Vertex index in DCEL == index in points[] (guaranteed by sequential add).
    for (const Point2D& p : points) {
        dcel.addVertex(p);
    }

    // ── Step 2: build adjacency list from accepted edges ─────────────────────
    // adj[i] = sorted list of vertices adjacent to i.
    std::vector<std::vector<std::size_t>> adj(n);
    for (const Edge& e : edges) {
        adj[e.u].push_back(e.v);
        adj[e.v].push_back(e.u);
    }

    // ── Step 3: find all triangular faces ────────────────────────────────────
    //
    // A triangle (i, j, k) exists iff all three edges (i,j), (j,k), (i,k)
    // are in the accepted set, AND the triangle has CCW orientation.
    //
    // Strategy: for each edge (u, v), find all w that are adjacent to both
    // u and v — these form candidate triangles.  Use a hash set for O(1) lookup.

    // Build edge set for O(1) membership test.
    std::unordered_set<EdgeKey, EdgeKeyHash> edgeSet;
    edgeSet.reserve(edges.size() * 2);
    for (const Edge& e : edges) {
        edgeSet.insert({e.u, e.v});
    }

    auto hasEdge = [&](std::size_t a, std::size_t b) -> bool {
        return edgeSet.count({std::min(a, b), std::max(a, b)}) > 0;
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

    // Collect triangles, avoiding duplicates.
    // A triangle {i,j,k} with i<j<k is stored once.
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

    // ── Step 4: add triangles to DCEL in CCW order ────────────────────────────
    //
    // We call addTriangle(va, vb, vc) which expects CCW vertex order.
    // Check orientation and swap if needed.

    // Map from EdgeKey → HalfEdgeIdx for twin-stitching after all adds.
    // addTriangle returns a FaceIdx; we retrieve half-edges via faceHalfEdges.
    std::unordered_map<EdgeKey, HalfEdgeIdx, EdgeKeyHash> heMap;
    heMap.reserve(edges.size() * 4);

    for (const TriIdx& tri : triangles) {
        const Point2D& pa = points[tri.a];
        const Point2D& pb = points[tri.b];
        const Point2D& pc = points[tri.c];

        // Ensure CCW order.
        VertexIdx va = static_cast<VertexIdx>(tri.a);
        VertexIdx vb = static_cast<VertexIdx>(tri.b);
        VertexIdx vc = static_cast<VertexIdx>(tri.c);

        if (!pred::isCCW(pa, pb, pc)) {
            std::swap(vb, vc); // flip to CCW
        }

        FaceIdx fi = dcel.addTriangle(va, vb, vc);

        // Record the half-edges for twin-stitching.
        // faceHalfEdges returns {h0: va→vb, h1: vb→vc, h2: vc→va}.
        auto [h0, h1, h2] = dcel.faceHalfEdges(fi);

        // Store directed edge → halfEdge index.
        // The twin of (a→b) has key stored as directed (b→a), so we key by
        // the directed pair and look up the reverse direction.
        // Convention: heMap[{src, dst}] = halfEdge going src→dst.
        auto storeHE = [&](HalfEdgeIdx h, VertexIdx src, VertexIdx dst) {
            heMap[{static_cast<std::size_t>(src),
                   static_cast<std::size_t>(dst)}] = h;
        };
        storeHE(h0, va, vb);
        storeHE(h1, vb, vc);
        storeHE(h2, vc, va);
    }

    // ── Step 5: stitch twins ──────────────────────────────────────────────────
    //
    // For each directed half-edge (u→v) in heMap, its twin is (v→u).
    // If (v→u) exists in heMap, stitch them; otherwise the edge is a hull
    // boundary (its twin stays pointing to the outer face, which addTriangle
    // already set up).

    // Track which pairs have been stitched to avoid double-calling.
    std::unordered_set<HalfEdgeIdx> stitched;
    stitched.reserve(heMap.size());

    for (auto& [key, h] : heMap) {
        if (stitched.count(h)) continue;
        EdgeKey revKey{key.hi, key.lo}; // reversed direction
        auto it = heMap.find(revKey);
        if (it != heMap.end()) {
            HalfEdgeIdx twinH = it->second;
            if (!stitched.count(twinH)) {
                dcel.stitchTwins(h, twinH);
                stitched.insert(h);
                stitched.insert(twinH);
            }
        }
        // If no twin found: boundary edge — leave as is (outer face twin).
    }
}

} // namespace triangulation
