#include "dcel.hh"
#include "predicates.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace geometry {
    // ─────────────────────────────────────────────────────────────────────────────
    //  Private assert helpers
    // ─────────────────────────────────────────────────────────────────────────────

    void DCEL::assertVertex(VertexIdx i) const noexcept {
#ifndef NDEBUG
        assert(i < vertices_.size() && !vertices_[i].dead);
#else
        (void) i;
#endif
    }

    void DCEL::assertHalfEdge(HalfEdgeIdx i) const noexcept {
#ifndef NDEBUG
        assert(i < halfEdges_.size() && !halfEdges_[i].dead);
#else
        (void) i;
#endif
    }

    void DCEL::assertFace(FaceIdx i) const noexcept {
#ifndef NDEBUG
        assert(i < faces_.size());
#else
        (void) i;
#endif
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  Construction / reset
    // ─────────────────────────────────────────────────────────────────────────────

    DCEL::DCEL() {
        // Reserve a small amount to avoid first-insertion reallocations.
        vertices_.reserve(8);
        halfEdges_.reserve(24);
        faces_.reserve(8);

        // Create the outer (unbounded) face at index 0.
        faces_.emplace_back(kInvalidIdx); // edge = kInvalidIdx — OK for outer face
    }

    void DCEL::clear() {
        vertices_.clear();
        halfEdges_.clear();
        faces_.clear();
        // Re-create the outer face.
        faces_.emplace_back(kInvalidIdx);
    }

    void DCEL::reserve(std::size_t nVertices, std::size_t nFaces) {
        vertices_.reserve(nVertices);
        faces_.reserve(nFaces + 1); // +1 for outer face
        halfEdges_.reserve(nFaces * 3 * 2); // ~6 half-edges per triangle
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  Element accessors
    // ─────────────────────────────────────────────────────────────────────────────

    Vertex &DCEL::vertex(VertexIdx i) {
        assertVertex(i);
        return vertices_[i];
    }

    const Vertex &DCEL::vertex(VertexIdx i) const {
        assertVertex(i);
        return vertices_[i];
    }

    HalfEdge &DCEL::halfEdge(HalfEdgeIdx i) {
        assertHalfEdge(i);
        return halfEdges_[i];
    }

    const HalfEdge &DCEL::halfEdge(HalfEdgeIdx i) const {
        assertHalfEdge(i);
        return halfEdges_[i];
    }

    Face &DCEL::face(FaceIdx i) {
        assertFace(i);
        return faces_[i];
    }

    const Face &DCEL::face(FaceIdx i) const {
        assertFace(i);
        return faces_[i];
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  Size queries
    // ─────────────────────────────────────────────────────────────────────────────

    std::size_t DCEL::vertexCount() const noexcept { return vertices_.size(); }
    std::size_t DCEL::halfEdgeCount() const noexcept { return halfEdges_.size(); }
    std::size_t DCEL::faceCount() const noexcept { return faces_.size(); }

    std::size_t DCEL::liveTriangleCount() const noexcept {
        std::size_t count = 0;
        for (std::size_t i = 1; i < faces_.size(); ++i) // skip outer face (0)
            if (!faces_[i].dead) ++count;
        return count;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  Low-level allocation
    // ─────────────────────────────────────────────────────────────────────────────

    VertexIdx DCEL::addVertex(Point2D pos) {
        const VertexIdx idx = vertices_.size();
        vertices_.emplace_back(pos);
        return idx;
    }

    HalfEdgeIdx DCEL::addHalfEdgePair() {
        const HalfEdgeIdx idx = halfEdges_.size();
        halfEdges_.emplace_back(); // h
        halfEdges_.emplace_back(); // h->twin  (idx + 1)

        // Wire the twin relationship immediately.
        halfEdges_[idx].twin = idx + 1;
        halfEdges_[idx + 1].twin = idx;

        return idx;
    }

    FaceIdx DCEL::addFace() {
        const FaceIdx idx = faces_.size();
        faces_.emplace_back();
        return idx;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  High-level construction: addTriangle
    // ─────────────────────────────────────────────────────────────────────────────

    FaceIdx DCEL::addTriangle(VertexIdx a, VertexIdx b, VertexIdx c) {
        assertVertex(a);
        assertVertex(b);
        assertVertex(c);

        // Allocate the new face.
        const FaceIdx fi = addFace();

        // Allocate 3 half-edge pairs.
        // Inner half-edges: h0 (a→b), h1 (b→c), h2 (c→a)
        // Outer half-edges: h0t (b→a), h1t (c→b), h2t (a→c)  — the twins
        const HalfEdgeIdx h0 = addHalfEdgePair(); // h0 = idx, h0t = idx+1
        const HalfEdgeIdx h1 = addHalfEdgePair();
        const HalfEdgeIdx h2 = addHalfEdgePair();

        const HalfEdgeIdx h0t = h0 + 1;
        const HalfEdgeIdx h1t = h1 + 1;
        const HalfEdgeIdx h2t = h2 + 1;

        // ── Inner half-edges (face fi, CCW: a→b→c→a) ──────────────────────────
        {
            HalfEdge &e0 = halfEdges_[h0];
            e0.origin = a;
            e0.next = h1;
            e0.prev = h2;
            e0.face = fi;
            // twin already set by addHalfEdgePair()
        }
        {
            HalfEdge &e1 = halfEdges_[h1];
            e1.origin = b;
            e1.next = h2;
            e1.prev = h0;
            e1.face = fi;
        }
        {
            HalfEdge &e2 = halfEdges_[h2];
            e2.origin = c;
            e2.next = h0;
            e2.prev = h1;
            e2.face = fi;
        }

        // ── Outer (twin) half-edges (face kOuterFace) ─────────────────────────
        // next/prev are NOT wired here — the triangulator must stitch them.
        // This is intentional: when adding the first isolated triangle we don't
        // yet know what the neighboring triangles will be.
        {
            HalfEdge &e0t = halfEdges_[h0t];
            e0t.origin = b; // twin of a→b starts at b
            e0t.face = kOuterFace;
            // next, prev = kInvalidIdx until stitched
        }
        {
            HalfEdge &e1t = halfEdges_[h1t];
            e1t.origin = c;
            e1t.face = kOuterFace;
        }
        {
            HalfEdge &e2t = halfEdges_[h2t];
            e2t.origin = a;
            e2t.face = kOuterFace;
        }

        // ── Face pointer ──────────────────────────────────────────────────────
        faces_[fi].edge = h0;

        // ── Vertex outgoing edge (only set if not yet assigned) ───────────────
        // We prefer the inner half-edge so the vertex points into the mesh.
        if (vertices_[a].edge == kInvalidIdx) vertices_[a].edge = h0;
        if (vertices_[b].edge == kInvalidIdx) vertices_[b].edge = h1;
        if (vertices_[c].edge == kInvalidIdx) vertices_[c].edge = h2;

        return fi;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  stitchTwins
    // ─────────────────────────────────────────────────────────────────────────────

    void DCEL::stitchTwins(HalfEdgeIdx h1, HalfEdgeIdx h2) {
        assert(h1 < halfEdges_.size());
        assert(h2 < halfEdges_.size());

        halfEdges_[h1].twin = h2;
        halfEdges_[h2].twin = h1;

        // After stitching the twin relation the dead flag on those outer half-edges
        // should be cleared (they may have been freshly created).
        halfEdges_[h1].dead = false;
        halfEdges_[h2].dead = false;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  flipEdge
    // ─────────────────────────────────────────────────────────────────────────────

    HalfEdgeIdx DCEL::flipEdge(HalfEdgeIdx h) {
        assertHalfEdge(h);

        //
        // Label the two triangles sharing edge h:
        //
        //          w (= dest of h->next)
        //         /|\
        //        / | \
        //      c2  |  c0         c0 = h          (u→v)
        //      /   |   \         c1 = h->next    (v→w)
        //     /    |    \        c2 = h->prev    (w→u)
        //    x  h* | h   v
        //     \    |    /        d0 = h*         (v→u)
        //      \   |   /         d1 = h*->next   (u→x)
        //      d1  |  d2         d2 = h*->prev   (x→v)
        //        \ | /
        //         \|/
        //          u (= h->origin)
        //
        // After flip: edge u↔v becomes edge w↔x.
        //

        const HalfEdgeIdx hIdx = h;
        const HalfEdgeIdx hStar = halfEdges_[h].twin;

        assert(halfEdges_[hIdx ].face != kOuterFace && "Cannot flip a boundary edge");
        assert(halfEdges_[hStar].face != kOuterFace && "Cannot flip a boundary edge");

        // Collect all six half-edges
        const HalfEdgeIdx c0 = hIdx;
        const HalfEdgeIdx c1 = halfEdges_[c0].next;
        const HalfEdgeIdx c2 = halfEdges_[c0].prev;

        const HalfEdgeIdx d0 = hStar;
        const HalfEdgeIdx d1 = halfEdges_[d0].next;
        const HalfEdgeIdx d2 = halfEdges_[d0].prev;

        // Collect the four vertices
        const VertexIdx u = halfEdges_[c0].origin; // c0: u→v
        const VertexIdx v = halfEdges_[d0].origin; // d0: v→u
        // const VertexIdx w = halfEdges_[c1].twin;    // we need dest(c1) = origin(twin(c1))
        // // dest(c1) = origin of twin(c1)
        // const VertexIdx wIdx = halfEdges_[halfEdges_[c1].twin].origin;
        // const VertexIdx xIdx = halfEdges_[halfEdges_[d1].twin].origin;
        // (void)w; (void)wIdx; // just use directly below

        const VertexIdx W = halfEdges_[halfEdges_[c1].twin].origin; // w = dest(c1)
        const VertexIdx X = halfEdges_[halfEdges_[d1].twin].origin; // x = dest(d1)

        // The two faces (reuse existing face objects)
        const FaceIdx fLeft = halfEdges_[c0].face; // triangle u,v,W
        const FaceIdx fRight = halfEdges_[d0].face; // triangle v,u,X

        //
        // Rewire: after the flip
        //   c0 becomes w→x  (was u→v)
        //   d0 becomes x→w  (was v→u)
        //   New triangles:
        //     fLeft  = w, x, u   (half-edges: c0=w→x, d1=x→u... wait, see below)
        //
        // Standard flip rewiring (see de Berg et al. §3.2):
        //
        //   c0: w→x,   next=d1,  prev=c2,  face=fLeft
        //   d0: x→w,   next=c1,  prev=d2,  face=fRight
        //   c1: v→w,   next=c2... NO — c1 stays the same direction but changes face
        //
        // Let me re-derive carefully:
        //   fLeft  gets half-edges: c0 (w→x), d1 (x→u... wait)
        //
        // The cleanest way: enumerate the two new triangles explicitly.
        //
        //   New left  triangle (w, x, u):  c0 (w→x),  d1 (x→u),  c2 (u→w)
        //   New right triangle (x, w, v):  d0 (x→w),  c1 (w→v),  d2 (v→x)
        //
        // Note: dest of c1 = w (we read from halfEdges_[c1].twin.origin = W).
        //       dest of d1 = x (halfEdges_[d1].twin.origin = X).
        //       c2: origin=W (was u→w before? No — c2=h->prev, whose origin = W,
        //           pointing toward u).
        //
        // Let me re-check: c0=h (u→v), c1=h->next (v→W), c2=h->prev (W→u).
        //   So: c2.origin = W, c2 points W→u. ✓
        //       d1.origin = u, d1 points u→X.
        //       d2.origin = X, d2 points X→v.
        //

        // ── Update origins of the flipped half-edges ──────────────────────────
        halfEdges_[c0].origin = W; // c0: w→x
        halfEdges_[d0].origin = X; // d0: x→w

        // ── Re-link left triangle: (W→x→u→W) using c0, d1, c2 ───────────────
        halfEdges_[c0].next = d1;
        halfEdges_[c0].prev = c2;
        halfEdges_[c0].face = fLeft;

        halfEdges_[d1].next = c2;
        halfEdges_[d1].prev = c0;
        halfEdges_[d1].face = fLeft;

        halfEdges_[c2].next = c0;
        halfEdges_[c2].prev = d1;
        halfEdges_[c2].face = fLeft;

        // ── Re-link right triangle: (X→w→v→X) using d0, c1, d2 ──────────────
        halfEdges_[d0].next = c1;
        halfEdges_[d0].prev = d2;
        halfEdges_[d0].face = fRight;

        halfEdges_[c1].next = d2;
        halfEdges_[c1].prev = d0;
        halfEdges_[c1].face = fRight;

        halfEdges_[d2].next = d0;
        halfEdges_[d2].prev = c1;
        halfEdges_[d2].face = fRight;

        // ── Update face->edge pointers (ensure they point to live half-edges) ─
        faces_[fLeft].edge = c0;
        faces_[fRight].edge = d0;

        // ── Update vertex outgoing edges ──────────────────────────────────────
        // u and v no longer have c0/d0 as valid outgoing edges — fix if needed.
        // We reassign to a half-edge we KNOW is still correct.
        // d1.origin == u, c1.origin == v (unchanged)
        vertices_[u].edge = d1; // u→x is still a valid outgoing edge from u
        vertices_[v].edge = d2; // v→x via d2... wait: d2.origin = X. Fix below.

        // Re-derive: after flip, what edges originate from u, v, W, X?
        //   u: d1 (u→X) — correct, d1.origin = u (unchanged)
        //   v: c1 (v→W) — correct, c1.origin = v (unchanged), still valid
        //   W: c0 (W→X) — c0.origin = W now
        //   X: d0 (X→W) — d0.origin = X now
        //
        // Only u and v's stored .edge pointers might be stale if they were
        // pointing to c0 or d0, which changed origin.  Safest: reassign all four.
        vertices_[u].edge = d1; // u → X
        vertices_[v].edge = c1; // v → W  (c1 origin is still v, unchanged)
        vertices_[W].edge = c0; // W → X
        vertices_[X].edge = d0; // X → W

        return c0; // The flipped half-edge (now W→X)
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  splitFace
    // ─────────────────────────────────────────────────────────────────────────────

    VertexIdx DCEL::splitFace(FaceIdx fi, Point2D p) {
        assertFace(fi);
        assert(!faces_[fi].dead && fi != kOuterFace);

        // Collect the three boundary half-edges of the face.
        const HalfEdgeIdx e0 = faces_[fi].edge;
        const HalfEdgeIdx e1 = halfEdges_[e0].next;
        const HalfEdgeIdx e2 = halfEdges_[e1].next;

        // Verify it is indeed a triangle.
        assert(halfEdges_[e2].next == e0 && "splitFace: face is not a triangle");

        const VertexIdx vA = halfEdges_[e0].origin; // a
        const VertexIdx vB = halfEdges_[e1].origin; // b
        const VertexIdx vC = halfEdges_[e2].origin; // c

        // Insert the new vertex p.
        const VertexIdx vP = addVertex(p);

        // Allocate 3 new half-edge pairs connecting p to each corner.
        // Pair i: inner half-edge from p toward vertex[i],
        //         twin half-edge from vertex[i] back toward p.
        const HalfEdgeIdx pa = addHalfEdgePair(); // p→a (pa) and a→p (pa+1)
        const HalfEdgeIdx pb = addHalfEdgePair(); // p→b (pb) and b→p (pb+1)
        const HalfEdgeIdx pc = addHalfEdgePair(); // p→c (pc) and c→p (pc+1)

        const HalfEdgeIdx paInner = pa, apInner = pa + 1;
        const HalfEdgeIdx pbInner = pb, bpInner = pb + 1;
        const HalfEdgeIdx pcInner = pc, cpInner = pc + 1;

        // Allocate two new faces (reuse fi for the first triangle).
        const FaceIdx f1 = fi; // triangle (p, a, b): e0, bpInner, paInner?
        // Let's be explicit:
        const FaceIdx f2 = addFace(); // triangle (p, b, c)
        const FaceIdx f3 = addFace(); // triangle (p, c, a)

        //
        // New triangles:
        //   f1 = (p → a → b → p):  pbInner (p→b... wait
        //
        // Let me be explicit about each triangle in CCW order,
        // using the original half-edges e0 (a→b), e1 (b→c), e2 (c→a).
        //
        //   f1 (p, a, b): p→a=paInner? No...
        //
        // The original face had half-edges e0:a→b, e1:b→c, e2:c→a.
        // After split we want:
        //   f1 (a, b, p): e0 (a→b),  bpInner (b→p),  paInner... wait, pa inner = p→a
        //
        // Let me name the new inner edges clearly:
        //   ap = apInner: origin a, going toward p  (twin of pa which is p→a)
        //   bp = bpInner: origin b, going toward p
        //   cp = cpInner: origin c, going toward p
        //   pa = paInner: origin p, going toward a
        //   pb = pbInner: origin p, going toward b
        //   pc = pcInner: origin p, going toward c
        //
        // Three triangles in CCW order:
        //   f1: a → b → p → a   uses: e0(a→b), bpInner(b→p), paInner(p→a)
        //   f2: b → c → p → b   uses: e1(b→c), cpInner(c→p), pbInner(p→b)
        //   f3: c → a → p → c   uses: e2(c→a), apInner(a→p), pcInner(p→c)
        //

        // Set origins
        halfEdges_[paInner].origin = vP;
        halfEdges_[apInner].origin = vA;
        halfEdges_[pbInner].origin = vP;
        halfEdges_[bpInner].origin = vB;
        halfEdges_[pcInner].origin = vP;
        halfEdges_[cpInner].origin = vC;

        // ── f1: a → b → p ────────────────────────────────────────────────────
        halfEdges_[e0].next = bpInner;
        halfEdges_[e0].prev = paInner;
        halfEdges_[e0].face = f1;
        halfEdges_[bpInner].next = paInner;
        halfEdges_[bpInner].prev = e0;
        halfEdges_[bpInner].face = f1;
        halfEdges_[paInner].next = e0;
        halfEdges_[paInner].prev = bpInner;
        halfEdges_[paInner].face = f1;
        faces_[f1].edge = e0;
        faces_[f1].dead = false;

        // ── f2: b → c → p ────────────────────────────────────────────────────
        halfEdges_[e1].next = cpInner;
        halfEdges_[e1].prev = pbInner;
        halfEdges_[e1].face = f2;
        halfEdges_[cpInner].next = pbInner;
        halfEdges_[cpInner].prev = e1;
        halfEdges_[cpInner].face = f2;
        halfEdges_[pbInner].next = e1;
        halfEdges_[pbInner].prev = cpInner;
        halfEdges_[pbInner].face = f2;
        faces_[f2].edge = e1;

        // ── f3: c → a → p ────────────────────────────────────────────────────
        halfEdges_[e2].next = apInner;
        halfEdges_[e2].prev = pcInner;
        halfEdges_[e2].face = f3;
        halfEdges_[apInner].next = pcInner;
        halfEdges_[apInner].prev = e2;
        halfEdges_[apInner].face = f3;
        halfEdges_[pcInner].next = e2;
        halfEdges_[pcInner].prev = apInner;
        halfEdges_[pcInner].face = f3;
        faces_[f3].edge = e2;

        // ── Internal twin pairs ───────────────────────────────────────────────
        // pa/ap: f1 has paInner (p→a), f3 has apInner (a→p) — these are twins.
        halfEdges_[paInner].twin = apInner;
        halfEdges_[apInner].twin = paInner;
        // pb/bp: f2 has pbInner (p→b), f1 has bpInner (b→p)
        halfEdges_[pbInner].twin = bpInner;
        halfEdges_[bpInner].twin = pbInner;
        // pc/cp: f3 has pcInner (p→c), f2 has cpInner (c→p)
        halfEdges_[pcInner].twin = cpInner;
        halfEdges_[cpInner].twin = pcInner;

        // ── Update vertex outgoing edges ──────────────────────────────────────
        vertices_[vP].edge = paInner; // p → a
        // a, b, c already have valid outgoing edges (e0, e1, e2); only update
        // if their stored edge was kInvalidIdx (shouldn't happen, but defensive).
        if (vertices_[vA].edge == kInvalidIdx) vertices_[vA].edge = e0;
        if (vertices_[vB].edge == kInvalidIdx) vertices_[vB].edge = e1;
        if (vertices_[vC].edge == kInvalidIdx) vertices_[vC].edge = e2;

        return vP;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  splitEdge
    // ─────────────────────────────────────────────────────────────────────────────

    VertexIdx DCEL::splitEdge(HalfEdgeIdx h, Point2D p) {
        assertHalfEdge(h);
        assert(halfEdges_[h].face != kOuterFace && "splitEdge: boundary edge not supported");

        const HalfEdgeIdx ht = halfEdges_[h].twin;

        // Collect surrounding half-edges of the two triangles.
        //
        //   Left  triangle (h  side): h  (u→v),  hn (v→w), hp (w→u)
        //   Right triangle (ht side): ht (v→u), htn (u→x), htp (x→v)
        //
        const HalfEdgeIdx hn = halfEdges_[h].next;
        const HalfEdgeIdx hp = halfEdges_[h].prev;
        const HalfEdgeIdx htn = halfEdges_[ht].next;
        const HalfEdgeIdx htp = halfEdges_[ht].prev;

        const VertexIdx u = halfEdges_[h].origin;
        const VertexIdx v = halfEdges_[ht].origin;
        const VertexIdx w = halfEdges_[hn].twin; // dest(hn) — wrong: twin gives the he not vertex
        // dest(hn) = origin(twin(hn))
        const VertexIdx W = halfEdges_[halfEdges_[hn].twin].origin;
        const VertexIdx X = halfEdges_[halfEdges_[htn].twin].origin;

        const FaceIdx fLeft = halfEdges_[h].face;
        const FaceIdx fRight = halfEdges_[ht].face;
        (void) w;

        // Insert the new vertex.
        const VertexIdx vP = addVertex(p);

        // Allocate new half-edge pairs:
        //   pair1: p→v (inner for left triangle), v→p twin
        //   pair2: p→u (inner for right), u→p twin
        //   (We REUSE h for u→p and ht for v→p, so we only need 2 new pairs.)
        //
        // Actually, the cleanest approach:
        //   Keep  h as u→p  (just change its twin-half and next)
        //   Keep ht as v→p  (analogously)
        //   Add new: pv (p→v) and pu_r (p→u for right triangle)
        //
        const HalfEdgeIdx pv = addHalfEdgePair(); // p→v (left),   v→p twin (pvt = pv+1)
        const HalfEdgeIdx pu = addHalfEdgePair(); // p→u (right),  u→p twin (put = pu+1)

        const HalfEdgeIdx pvInner = pv, vpInner = pv + 1;
        const HalfEdgeIdx puInner = pu, upInner = pu + 1;

        // Two new faces (reuse fLeft and fRight for the "u-side" halves,
        // create new faces for the "v-side" halves).
        const FaceIdx fLeftV = addFace(); // new triangle (p, v, W)
        const FaceIdx fRightV = addFace(); // new triangle (p, X, v) — actually (p, v, X... let's derive)

        // ── Set origins ───────────────────────────────────────────────────────
        // Reuse h as u→p: origin stays u, but twin will be upInner
        halfEdges_[h].twin = upInner;
        halfEdges_[h].origin = u;
        halfEdges_[ht].twin = vpInner;
        halfEdges_[ht].origin = v;

        halfEdges_[pvInner].origin = vP;
        halfEdges_[vpInner].origin = v;
        halfEdges_[puInner].origin = vP;
        halfEdges_[upInner].origin = u;

        // ── Left side ─────────────────────────────────────────────────────────
        // fLeft  (u → p → W → u):   h(u→p),  a new edge p→W?
        //
        // Wait — hn is already v→W.  After the split, left triangle has p instead of v.
        // So the left side becomes TWO triangles:
        //   fLeft  (u, p, W):  h(u→p),  new(p→W),  hp(W→u)
        //   fLeftV (p, v, W):  pvInner(p→v), hn(v→W), new(W→p)
        //
        // This requires a NEW half-edge pair p↔W.
        const HalfEdgeIdx pw = addHalfEdgePair(); // p→W (inner left), W→p twin
        const HalfEdgeIdx pwInner = pw, wpInner = pw + 1;
        halfEdges_[pwInner].origin = vP;
        halfEdges_[wpInner].origin = W;

        // fLeft: u → p → W → u
        halfEdges_[h].next = pwInner;
        halfEdges_[h].prev = hp;
        halfEdges_[h].face = fLeft;
        halfEdges_[pwInner].next = hp;
        halfEdges_[pwInner].prev = h;
        halfEdges_[pwInner].face = fLeft;
        halfEdges_[hp].next = h;
        halfEdges_[hp].prev = pwInner;
        halfEdges_[hp].face = fLeft;
        faces_[fLeft].edge = h;
        faces_[fLeft].dead = false;

        // fLeftV: p → v → W → p
        halfEdges_[pvInner].next = hn;
        halfEdges_[pvInner].prev = wpInner;
        halfEdges_[pvInner].face = fLeftV;
        halfEdges_[hn].next = wpInner;
        halfEdges_[hn].prev = pvInner;
        halfEdges_[hn].face = fLeftV;
        halfEdges_[wpInner].next = pvInner;
        halfEdges_[wpInner].prev = hn;
        halfEdges_[wpInner].face = fLeftV;
        faces_[fLeftV].edge = pvInner;

        halfEdges_[pvInner].twin = vpInner; // ht is now v→p
        halfEdges_[vpInner].twin = pvInner;
        halfEdges_[pwInner].twin = wpInner;
        halfEdges_[wpInner].twin = pwInner;

        // ── Right side ────────────────────────────────────────────────────────
        // Similarly, add p↔X pair.
        const HalfEdgeIdx px = addHalfEdgePair(); // p→X, X→p twin
        const HalfEdgeIdx pxInner = px, xpInner = px + 1;
        halfEdges_[pxInner].origin = vP;
        halfEdges_[xpInner].origin = X;

        // fRight: v → p → X → v
        halfEdges_[ht].next = pxInner;
        halfEdges_[ht].prev = htp;
        halfEdges_[ht].face = fRight;
        halfEdges_[pxInner].next = htp;
        halfEdges_[pxInner].prev = ht;
        halfEdges_[pxInner].face = fRight;
        halfEdges_[htp].next = ht;
        halfEdges_[htp].prev = pxInner;
        halfEdges_[htp].face = fRight;
        faces_[fRight].edge = ht;
        faces_[fRight].dead = false;

        // fRightV: p → u → X → p  (htn is u→X)
        halfEdges_[puInner].next = htn;
        halfEdges_[puInner].prev = xpInner;
        halfEdges_[puInner].face = fRightV;
        halfEdges_[htn].next = xpInner;
        halfEdges_[htn].prev = puInner;
        halfEdges_[htn].face = fRightV;
        halfEdges_[xpInner].next = puInner;
        halfEdges_[xpInner].prev = htn;
        halfEdges_[xpInner].face = fRightV;
        faces_[fRightV].edge = puInner;

        halfEdges_[puInner].twin = upInner;
        halfEdges_[upInner].twin = puInner;
        halfEdges_[pxInner].twin = xpInner;
        halfEdges_[xpInner].twin = pxInner;

        // h and ht are now u→p and v→p; their twins are upInner and vpInner.
        halfEdges_[upInner].twin = h;
        halfEdges_[vpInner].twin = ht;

        // ── Vertex outgoing edges ─────────────────────────────────────────────
        vertices_[vP].edge = pvInner; // p → v  (arbitrary valid choice)
        vertices_[u].edge = h;
        vertices_[v].edge = ht;

        return vP;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  Deletion
    // ─────────────────────────────────────────────────────────────────────────────

    void DCEL::killFace(FaceIdx fi) {
        assert(fi != kOuterFace && "Cannot kill the outer face");
        assertFace(fi);
        faces_[fi].dead = true;
        // Mark the face's half-edges as dead too so traversal skips them.
        // TODO: ERROR FIX (DLN, QGR assertion fail here)
        HalfEdgeIdx h = faces_[fi].edge;
        if (h == kInvalidIdx) return;
        for (int i = 0; i < 3; ++i) {
            halfEdges_[h].dead = true;
            h = halfEdges_[h].next;
        }
    }

    void DCEL::removeVertex(VertexIdx vi) {
        assertVertex(vi);

        std::vector<FaceIdx> facesToKill;
        for (std::size_t i = 0; i < halfEdges_.size(); ++i) {
            if (!halfEdges_[i].dead && halfEdges_[i].origin == vi) {
                FaceIdx fi = halfEdges_[i].face;
                if (fi != kOuterFace && !faces_[fi].dead) {
                    facesToKill.push_back(fi);
                }
            }
        }

        for (FaceIdx fi: facesToKill) {
            if (!faces_[fi].dead) {
                killFace(fi);
            }
        }

        for (std::size_t i = 0; i < halfEdges_.size(); ++i) {
            if (halfEdges_[i].dead) continue;

            bool connectedToVi = (halfEdges_[i].origin == vi);
            if (!connectedToVi && halfEdges_[i].twin != kInvalidIdx) {
                if (halfEdges_[halfEdges_[i].twin].origin == vi) {
                    connectedToVi = true;
                }
            }

            if (connectedToVi) {
                halfEdges_[i].dead = true;
            }
        }

        vertices_[vi].dead = true;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  Traversal helpers
    // ─────────────────────────────────────────────────────────────────────────────

    std::array<VertexIdx, 3> DCEL::faceVertices(FaceIdx fi) const {
        assertFace(fi);
        const HalfEdgeIdx h0 = faces_[fi].edge;
        const HalfEdgeIdx h1 = halfEdges_[h0].next;
        const HalfEdgeIdx h2 = halfEdges_[h1].next;
        return {
            halfEdges_[h0].origin,
            halfEdges_[h1].origin,
            halfEdges_[h2].origin
        };
    }

    std::array<HalfEdgeIdx, 3> DCEL::faceHalfEdges(FaceIdx fi) const {
        assertFace(fi);
        const HalfEdgeIdx h0 = faces_[fi].edge;
        const HalfEdgeIdx h1 = halfEdges_[h0].next;
        const HalfEdgeIdx h2 = halfEdges_[h1].next;
        return {h0, h1, h2};
    }

    std::array<FaceIdx, 3> DCEL::faceNeighbors(FaceIdx fi) const {
        assertFace(fi);
        const auto [h0, h1, h2] = faceHalfEdges(fi);
        const auto neighbor = [&](HalfEdgeIdx h) -> FaceIdx {
            const FaceIdx f = halfEdges_[halfEdges_[h].twin].face;
            return (f == kOuterFace || faces_[f].dead) ? kInvalidIdx : f;
        };
        return {neighbor(h0), neighbor(h1), neighbor(h2)};
    }

    std::vector<HalfEdgeIdx> DCEL::vertexStar(VertexIdx vi) const {
        assertVertex(vi);
        std::vector<HalfEdgeIdx> result;
        const HalfEdgeIdx start = vertices_[vi].edge;
        if (start == kInvalidIdx) return result;

        HalfEdgeIdx cur = start;
        // Safety limit to break infinite loops in a corrupt DCEL.
        const std::size_t maxIter = halfEdges_.size() + 1;
        std::size_t iters = 0;
        do {
            result.push_back(cur);
            cur = halfEdges_[halfEdges_[cur].twin].next;
            if (++iters > maxIter) break; // corrupt DCEL guard
        } while (cur != start);

        return result;
    }

    HalfEdgeIdx DCEL::findHalfEdge(VertexIdx u, VertexIdx v) const {
        for (HalfEdgeIdx he: vertexStar(u)) {
            if (halfEdges_[halfEdges_[he].twin].origin == v)
                return he;
        }
        return kInvalidIdx;
    }

    FaceIdx DCEL::findFace(Point2D p) const {
        for (FaceIdx fi = 1; fi < faces_.size(); ++fi) {
            if (faces_[fi].dead) continue;
            const auto [a, b, c] = faceVertices(fi);
            const Triangle tri{vertices_[a].pos, vertices_[b].pos, vertices_[c].pos};
            if (tri.contains(p)) return fi;
        }
        return kInvalidIdx;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  validate
    // ─────────────────────────────────────────────────────────────────────────────

    bool DCEL::validate() const {
        bool ok = true;

        auto fail = [&](const char *msg, std::size_t idx) {
            std::cerr << "[DCEL::validate] " << msg << " at index " << idx << '\n';
            ok = false;
        };

        // Check half-edge invariants.
        for (std::size_t i = 0; i < halfEdges_.size(); ++i) {
            if (halfEdges_[i].dead) continue;
            const HalfEdge &h = halfEdges_[i];

            // I1: twin->twin == self
            if (h.twin == kInvalidIdx || halfEdges_[h.twin].twin != i)
                fail("I1 violated: twin->twin != self", i);

            // I2: next->prev == self
            if (h.next == kInvalidIdx || halfEdges_[h.next].prev != i)
                fail("I2 violated: next->prev != self", i);

            // I3: next->origin == twin->origin
            if (h.next != kInvalidIdx && h.twin != kInvalidIdx) {
                if (halfEdges_[h.next].origin != halfEdges_[h.twin].origin)
                    fail("I3 violated: next->origin != twin->origin", i);
            }

            // I4: next->face == face
            if (h.next != kInvalidIdx) {
                if (halfEdges_[h.next].face != h.face)
                    fail("I4 violated: next->face != face", i);
            }
        }

        // Check vertex invariants.
        for (std::size_t i = 0; i < vertices_.size(); ++i) {
            if (vertices_[i].dead) continue;
            const Vertex &v = vertices_[i];
            // I5: v->edge->origin == v
            if (v.edge != kInvalidIdx && halfEdges_[v.edge].origin != i)
                fail("I5 violated: vertex->edge->origin != vertex", i);
        }

        // Check face invariants.
        for (std::size_t i = 0; i < faces_.size(); ++i) {
            if (faces_[i].dead) continue;
            const Face &f = faces_[i];
            if (f.edge == kInvalidIdx) continue; // outer face may have no edge
            // I6: f->edge->face == f
            if (halfEdges_[f.edge].face != i)
                fail("I6 violated: face->edge->face != face", i);
        }

        return ok;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  totalWeight
    // ─────────────────────────────────────────────────────────────────────────────

    double DCEL::totalWeight() const noexcept {
        double sum = 0.0;
        for (std::size_t i = 0; i < halfEdges_.size(); ++i) {
            const HalfEdge &h = halfEdges_[i];
            if (h.dead) continue;
            // Count each undirected edge once: only when i < twin index.
            if (h.twin == kInvalidIdx || i > h.twin) continue;
            if (h.face == kOuterFace && halfEdges_[h.twin].face == kOuterFace) continue;
            const Point2D &a = vertices_[h.origin].pos;
            const Point2D &b = vertices_[halfEdges_[h.twin].origin].pos;
            sum += a.dist(b);
        }
        return sum;
    }
} // namespace geometry
