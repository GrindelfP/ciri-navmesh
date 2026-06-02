#pragma once

/**
 * @file dcel.hh
 * @brief Doubly Connected Edge List (DCEL) — half-edge mesh data structure.
 *
 * ## Structure overview
 *
 * A DCEL stores a planar subdivision as three interlinked pools:
 *
 *   Vertex    — a point in the plane + one outgoing HalfEdge
 *   HalfEdge  — directed half-edge with: origin, twin, next, prev, face
 *   Face      — a bounded or unbounded region + one incident HalfEdge
 *
 * Every undirected edge (u, v) is represented as TWO half-edges:
 *   h  : u → v   (belongs to the left face)
 *   h* : v → u   (belongs to the right face, twin of h)
 *
 * ## Invariants (must hold after every mutating operation)
 *
 *   I1.  h->twin->twin == h
 *   I2.  h->next->prev == h
 *   I3.  h->next->origin == h->twin->origin   (consistent vertex pointers)
 *   I4.  h->face == h->next->face             (all edges of a face agree)
 *   I5.  v->edge->origin == v                 (vertex outgoing edge is correct)
 *   I6.  f->edge->face   == f                 (face edge points back to face)
 *   I7.  The unbounded face is always present (index 0 by convention).
 *
 * ## Index stability
 *
 * Elements are stored in std::vector.  Indices (VertexIdx, HalfEdgeIdx,
 * FaceIdx) are stable as long as no element is REMOVED — appending never
 * invalidates existing indices.  Removal sets a "dead" flag; dead slots are
 * NOT reused automatically (keeps indices stable for the navmesh graph that
 * references them).
 *
 * ## Coordinate system
 *
 * CCW face orientation is assumed throughout (consistent with predicates.hpp).
 * The unbounded outer face wraps CW half-edges on the convex hull boundary.
 */

#include "primitives.hpp"

#include <cassert>
#include <cstddef>
#include <limits>
#include <vector>

namespace geometry {

// ─────────────────────────────────────────────────────────────────────────────
//  Index types  (strong typedefs to avoid mixing up vertex/edge/face indices)
// ─────────────────────────────────────────────────────────────────────────────

/// @brief Sentinel value meaning "no element" / null pointer in index form.
inline constexpr std::size_t kInvalidIdx = std::numeric_limits<std::size_t>::max();

/// @brief Index into DCEL::vertices_.
using VertexIdx   = std::size_t;
/// @brief Index into DCEL::halfEdges_.
using HalfEdgeIdx = std::size_t;
/// @brief Index into DCEL::faces_.
using FaceIdx     = std::size_t;

/// @brief Conventional index of the unbounded (outer) face.
inline constexpr FaceIdx kOuterFace = 0;

// ─────────────────────────────────────────────────────────────────────────────
//  Forward declaration
// ─────────────────────────────────────────────────────────────────────────────

class DCEL;

// ─────────────────────────────────────────────────────────────────────────────
//  Vertex
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief A vertex in the DCEL.
 *
 * Stores its 2-D position and the index of ONE outgoing half-edge.
 * All other incident half-edges are reachable by walking
 * edge->twin->next repeatedly (the "vertex star" traversal).
 */
struct Vertex {
    Point2D     pos;                    ///< 2-D coordinates
    HalfEdgeIdx edge{kInvalidIdx};      ///< One outgoing half-edge from this vertex
    bool        dead{false};            ///< True if logically deleted

    explicit Vertex(Point2D p, HalfEdgeIdx e = kInvalidIdx)
        : pos(p), edge(e) {}
};

// ─────────────────────────────────────────────────────────────────────────────
//  HalfEdge
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief A directed half-edge in the DCEL.
 *
 * Together with its twin, represents one undirected edge of the mesh.
 *
 * Traversal patterns:
 *   - Face loop:        h, h->next, h->next->next, ... until back to h
 *   - Vertex star CCW:  h, h->twin->next, (h->twin->next)->twin->next, ...
 *   - Destination:      h->twin->origin
 */
struct HalfEdge {
    VertexIdx   origin{kInvalidIdx};    ///< Start vertex of this half-edge
    HalfEdgeIdx twin  {kInvalidIdx};    ///< Opposite half-edge (same edge, opposite dir)
    HalfEdgeIdx next  {kInvalidIdx};    ///< Next half-edge around the same face (CCW)
    HalfEdgeIdx prev  {kInvalidIdx};    ///< Previous half-edge around the same face
    FaceIdx     face  {kInvalidIdx};    ///< Face to the left of this half-edge
    bool        dead  {false};          ///< True if logically deleted

    HalfEdge() = default;
    HalfEdge(VertexIdx o, HalfEdgeIdx tw, HalfEdgeIdx nx,
             HalfEdgeIdx pv, FaceIdx f)
        : origin(o), twin(tw), next(nx), prev(pv), face(f) {}
};

// ─────────────────────────────────────────────────────────────────────────────
//  Face
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief A face (polygon) in the DCEL.
 *
 * After triangulation every bounded face is a triangle.
 * Face 0 is always the outer (unbounded) face.
 */
struct Face {
    HalfEdgeIdx edge{kInvalidIdx};  ///< One half-edge on the boundary of this face
    bool        dead{false};        ///< True if logically deleted

    explicit Face(HalfEdgeIdx e = kInvalidIdx) : edge(e) {}
};

// ─────────────────────────────────────────────────────────────────────────────
//  DCEL
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Doubly Connected Edge List: the central data structure of the library.
 *
 * Owns three pools (vertices, half-edges, faces) stored as contiguous vectors.
 * All cross-references use integer indices (not pointers) for stability and
 * serializability.
 *
 * ### Typical usage
 * @code
 *   DCEL dcel;
 *   // Triangulator fills it via addVertex / addTriangle / flipEdge ...
 *   for (FaceIdx fi = 1; fi < dcel.faceCount(); ++fi) {
 *       if (dcel.face(fi).dead) continue;
 *       auto [a, b, c] = dcel.faceVertices(fi);
 *       // use triangle (a, b, c)
 *   }
 * @endcode
 */
class DCEL {
public:
    // ── Construction ─────────────────────────────────────────────────────────

    /**
     * @brief Default constructor.
     *
     * Creates the outer (unbounded) face at index 0.
     * All triangulators call this before inserting any geometry.
     */
    DCEL();

    /// @brief Clears all geometry and resets to initial state (outer face only).
    void clear();

    /**
     * @brief Reserves memory for expected element counts.
     *
     * Call before bulk-inserting n points to avoid repeated reallocations.
     * For n input points a Delaunay triangulation has ~2n triangles,
     * ~6n half-edges, ~3n edges.
     *
     * @param nVertices   Expected number of vertices
     * @param nFaces      Expected number of bounded faces (triangles)
     */
    void reserve(std::size_t nVertices, std::size_t nFaces);

    // ── Element accessors (bounds-checked in debug, unchecked in release) ────

    /// @brief Access vertex by index. Asserts index is valid.
    [[nodiscard]] Vertex&       vertex(VertexIdx i);
    [[nodiscard]] const Vertex& vertex(VertexIdx i) const;

    /// @brief Access half-edge by index.
    [[nodiscard]] HalfEdge&       halfEdge(HalfEdgeIdx i);
    [[nodiscard]] const HalfEdge& halfEdge(HalfEdgeIdx i) const;

    /// @brief Access face by index.
    [[nodiscard]] Face&       face(FaceIdx i);
    [[nodiscard]] const Face& face(FaceIdx i) const;

    // ── Size queries ──────────────────────────────────────────────────────────

    /// @brief Total allocated vertices (including dead ones).
    [[nodiscard]] std::size_t vertexCount()   const noexcept;
    /// @brief Total allocated half-edges (including dead ones).
    [[nodiscard]] std::size_t halfEdgeCount() const noexcept;
    /// @brief Total allocated faces (including dead ones, including outer face).
    [[nodiscard]] std::size_t faceCount()     const noexcept;

    /// @brief Number of live (non-dead) bounded faces (triangles).
    [[nodiscard]] std::size_t liveTriangleCount() const noexcept;

    // ── Low-level allocation ──────────────────────────────────────────────────

    /**
     * @brief Allocates a new vertex and returns its index.
     * @param pos  Position of the vertex.
     * @return VertexIdx of the newly created vertex.
     */
    VertexIdx addVertex(Point2D pos);

    /**
     * @brief Allocates a pair of twin half-edges and returns the index of
     *        the first one (the second is first + 1 by construction).
     *
     * @note  Fields (origin, next, prev, face) are left as kInvalidIdx —
     *        the caller must fill them in before any traversal.
     * @return HalfEdgeIdx of the first half-edge (twin = returned + 1).
     */
    [[nodiscard]] HalfEdgeIdx addHalfEdgePair();

    /**
     * @brief Allocates a new face and returns its index.
     * @return FaceIdx of the newly created face.
     */
    [[nodiscard]] FaceIdx addFace();

    // ── High-level construction ───────────────────────────────────────────────

    /**
     * @brief Inserts an isolated triangle (a, b, c) in CCW order.
     *
     * Creates:
     *   - 3 half-edges for the inner face (CCW)
     *   - 3 half-edges for the outer boundary (CW twins)
     *   - 1 bounded face
     *
     * The outer twins are linked to kOuterFace but NOT wired to each other's
     * next/prev — that is done by the triangulator when adjacent triangles
     * are stitched together.
     *
     * @pre  a, b, c are valid, live VertexIdx values.
     * @pre  The triangle (pos[a], pos[b], pos[c]) is CCW (not checked here).
     * @return FaceIdx of the newly created triangle face.
     */
    FaceIdx addTriangle(VertexIdx a, VertexIdx b, VertexIdx c);

    /**
     * @brief Stitches two boundary half-edges that should be twins.
     *
     * Called when two adjacent triangles are inserted and share an edge.
     * Sets h1->twin = h2 and h2->twin = h1, and updates the origin pointers
     * so that h1->origin == h2->twin->origin.
     *
     * @pre  h1 and h2 are boundary half-edges (face == kOuterFace or
     *       newly created, not yet twinned to a real face).
     * @pre  h1->twin->origin == h2->origin  (they represent opposite
     *       directions of the same geometric edge).
     */
    void stitchTwins(HalfEdgeIdx h1, HalfEdgeIdx h2);

    // ── Core mutation: edge flip ──────────────────────────────────────────────

    /**
     * @brief Flips the shared edge between two adjacent triangles.
     *
     * Given half-edge @p h (pointing from vertex u to vertex v), and its
     * twin pointing from v to u, the two triangles sharing this edge are:
     *   Left  triangle:  u → v → w  (w = h->next->twin->origin)
     *   Right triangle:  v → u → x  (x = h->twin->next->twin->origin)
     *
     * After the flip, the shared edge becomes w → x (and x → w).
     *
     * @verbatim
     *   Before:          After:
     *     w                w
     *    /|\              / \
     *   / | \            /   \
     *  u  |  v   ==>   u-----v
     *   \ | /            \   /
     *    \|/              \ /
     *     x                x
     * @endverbatim
     *
     * All half-edge, face, and vertex pointers are updated to maintain
     * invariants I1–I7.
     *
     * @param h  Any half-edge of the edge to flip.
     *           Must not be a boundary edge (face or twin->face == kOuterFace).
     * @return   HalfEdgeIdx of the new edge (h is reused, pointing w → x).
     *
     * @pre  Both faces adjacent to h are triangles (face loop has length 3).
     * @pre  h is not a boundary edge.
     */
    HalfEdgeIdx flipEdge(HalfEdgeIdx h);

    // ── Vertex insertion ──────────────────────────────────────────────────────

    /**
     * @brief Splits a triangular face by inserting a point inside it.
     *
     * The face @p fi (triangle a→b→c) is replaced by three new triangles:
     *   (p, a, b),  (p, b, c),  (p, c, a)
     *
     * The original face fi is reused for (p, a, b); two new faces are created.
     * Six new half-edges are created (three pairs connecting p to a, b, c).
     *
     * @param fi   Index of the triangular face to split.
     * @param p    Position of the new interior point.
     * @return     VertexIdx of the inserted point.
     *
     * @pre  fi is a live, non-outer, triangular face.
     * @pre  p is strictly inside the triangle (not on an edge or outside).
     *       Caller is responsible for this precondition.
     */
    VertexIdx splitFace(FaceIdx fi, Point2D p);

    /**
     * @brief Splits the edge represented by half-edge @p h by inserting
     *        a point at position @p p (typically the midpoint).
     *
     * The edge h (u→v) and its twin (v→u) are replaced by four half-edges:
     *   u→p, p→v  and their twins  v→p, p→u.
     * The two triangles sharing the edge are each split into two triangles.
     *
     * @param h   Half-edge of the edge to split.
     * @param p   Position of the new point (must lie on the edge u-v).
     * @return    VertexIdx of the inserted point.
     *
     * @pre  h is not a boundary half-edge.
     * @pre  p lies on segment (vertex(h->origin).pos, vertex(h->twin->origin).pos).
     */
    VertexIdx splitEdge(HalfEdgeIdx h, Point2D p);

    // ── Deletion ──────────────────────────────────────────────────────────────

    /**
     * @brief Marks a face as dead (logically deleted).
     *
     * Used by Bowyer–Watson to remove triangles whose circumcircle contains
     * the newly inserted point.  Half-edges of the face are NOT freed —
     * the triangulator re-links them when rebuilding the cavity.
     *
     * @param fi  Face index to kill. Must not be kOuterFace.
     */
    void killFace(FaceIdx fi);

    /**
     * @brief Marks a vertex and all its incident half-edges/faces as dead.
     *
     * Used at the end of Bowyer–Watson to remove the super-triangle vertices.
     * After this call, any face that had one of these vertices is also dead.
     *
     * @param vi  Vertex index to remove.
     */
    void removeVertex(VertexIdx vi);

    // ── Traversal helpers ─────────────────────────────────────────────────────

    /**
     * @brief Returns the three vertex indices of a triangular face.
     *
     * Walks the face loop exactly 3 steps and collects origin indices.
     *
     * @param fi  A live, bounded, triangular face.
     * @return    Array {a, b, c} in CCW order.
     *
     * @pre  The face loop has exactly 3 half-edges.
     */
    [[nodiscard]] std::array<VertexIdx, 3> faceVertices(FaceIdx fi) const;

    /**
     * @brief Returns the three half-edge indices that bound a triangular face.
     *
     * @param fi  A live, bounded, triangular face.
     * @return    Array {h0, h1, h2} where h0 == face(fi).edge.
     */
    [[nodiscard]] std::array<HalfEdgeIdx, 3> faceHalfEdges(FaceIdx fi) const;

    /**
     * @brief Returns up to 3 face indices adjacent to face @p fi.
     *
     * An adjacent face is the face on the other side of each boundary edge.
     * Boundary edges (twin->face == kOuterFace) contribute kInvalidIdx.
     *
     * @param fi  A live, bounded, triangular face.
     * @return    Array {n0, n1, n2}; entry is kInvalidIdx for hull edges.
     */
    [[nodiscard]] std::array<FaceIdx, 3> faceNeighbors(FaceIdx fi) const;

    /**
     * @brief Collects all half-edges outgoing from vertex @p vi.
     *
     * Walks the vertex star: starting from vi.edge, repeatedly follows
     * twin->next until returning to the start.
     *
     * @param vi  A live vertex.
     * @return    Vector of HalfEdgeIdx, all with origin == vi.
     *
     * @note  O(degree) time.
     */
    [[nodiscard]] std::vector<HalfEdgeIdx> vertexStar(VertexIdx vi) const;

    /**
     * @brief Finds the half-edge from vertex @p u to vertex @p v, if any.
     *
     * Walks the star of u searching for an edge whose twin->origin == v.
     *
     * @return HalfEdgeIdx if found, kInvalidIdx otherwise.
     * @note  O(degree(u)) time.
     */
    [[nodiscard]] HalfEdgeIdx findHalfEdge(VertexIdx u, VertexIdx v) const;

    /**
     * @brief Finds the triangular face containing point @p p.
     *
     * Simple linear walk: checks every live face.  O(n).
     * For the triangulators this is acceptable; a spatial index can be
     * added later for the navmesh query path.
     *
     * @param p   Query point.
     * @return    FaceIdx of a face whose triangle contains p,
     *            or kInvalidIdx if p is outside the triangulation.
     */
    [[nodiscard]] FaceIdx findFace(Point2D p) const;

    // ── Validation (debug builds) ─────────────────────────────────────────────

    /**
     * @brief Checks all DCEL invariants (I1–I7) for every live element.
     *
     * Runs in O(n).  Intended for use in unit tests and debug assertions.
     *
     * @return true if all invariants hold; false (and logs the first
     *         violation to stderr) otherwise.
     */
    [[nodiscard]] bool validate() const;

    // ── Total edge weight (used for triangulation comparison metric) ──────────

    /**
     * @brief Computes the total weight of the triangulation.
     *
     * Sums the Euclidean length of every undirected edge exactly once
     * (iterates over even-indexed half-edges whose index < twin index).
     *
     * @return Sum of all edge lengths.
     */
    [[nodiscard]] double totalWeight() const noexcept;

private:
    std::vector<Vertex>   vertices_;
    std::vector<HalfEdge> halfEdges_;
    std::vector<Face>     faces_;

    // ── Private helpers ───────────────────────────────────────────────────────

    /**
     * @brief Asserts that index @p i is a valid, live vertex.
     *        No-op in release builds.
     */
    void assertVertex(VertexIdx i) const noexcept;

    /**
     * @brief Asserts that index @p i is a valid, live half-edge.
     */
    void assertHalfEdge(HalfEdgeIdx i) const noexcept;

    /**
     * @brief Asserts that index @p i is a valid face (may be dead or outer).
     */
    void assertFace(FaceIdx i) const noexcept;
};

} // namespace geometry
