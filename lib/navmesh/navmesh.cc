/**
 * @file navmesh.cc
 * @brief NavMesh implementation: builds an adjacency graph from a DCEL.
 */

#include "navmesh.hh"

#include <cassert>
#include <stdexcept>

namespace navmesh {
    using geometry::DCEL;
    using geometry::FaceIdx;
    using geometry::HalfEdgeIdx;
    using geometry::Point2D;
    using geometry::VertexIdx;
    using geometry::kInvalidIdx;
    using geometry::kOuterFace;

    // ─────────────────────────────────────────────────────────────────────────────
    //  build
    // ─────────────────────────────────────────────────────────────────────────────
    void NavMesh::build(const DCEL &dcel) {
        clear();

        const std::size_t nFaces = dcel.faceCount();

        // faceToNode_[fi] = NodeIdx of that face, or kInvalidNode for dead/outer.
        faceToNode_.assign(nFaces, kInvalidNode);

        // ── Pass 1: create one Node per live bounded face ─────────────────────────

        nodes_.reserve(dcel.liveTriangleCount());

        for (FaceIdx fi = 1; fi < nFaces; ++fi) {
            // fi=0 is outer face
            const auto &f = dcel.face(fi);
            if (f.dead) continue;

            auto [va, vb, vc] = dcel.faceVertices(fi);

            const Point2D &pa = dcel.vertex(va).pos;
            const Point2D &pb = dcel.vertex(vb).pos;
            const Point2D &pc = dcel.vertex(vc).pos;

            Node nd;
            nd.faceIdx = fi;
            nd.centroid = (pa + pb + pc) / 3.0;
            nd.vertices = {va, vb, vc};
            // neighbors filled in Pass 2

            NodeIdx ni = nodes_.size();
            nodes_.push_back(std::move(nd));
            faceToNode_[fi] = ni;
        }

        // ── Pass 2: link adjacent nodes via Arc entries ───────────────────────────
        //
        // For each node, walk the 3 half-edges of its face.
        // For each half-edge h, the twin's face is the adjacent triangle.
        // If the twin's face is live and bounded, create an Arc.
        //
        // We use faceNeighbors() for clarity; it returns up to 3 adjacent FaceIdx
        // (kInvalidIdx for hull edges).
        //
        // Portal convention: the portal for the arc from node A to node B
        // is the shared edge as seen from A's face loop — i.e. the two vertices
        // of the half-edge h going from A to B in A's CCW winding.
        //
        //   h: va → vb  (CCW in face A)
        //   portal = {origin(h), origin(twin(h))}  =  {va, vb}
        //
        // From A's perspective the portal is traversed left-to-right as (va, vb).
        // The Funnel Algorithm uses this to build the sleeve of portals.

        for (NodeIdx ni = 0; ni < nodes_.size(); ++ni) {
            Node &nd = nodes_[ni];
            FaceIdx fi = nd.faceIdx;

            auto halfEdges = dcel.faceHalfEdges(fi);

            for (HalfEdgeIdx h: halfEdges) {
                const auto &he = dcel.halfEdge(h);
                FaceIdx adjF = dcel.halfEdge(he.twin).face;

                if (adjF == kOuterFace) continue; // hull boundary — no arc
                if (dcel.face(adjF).dead) continue; // dead face — skip

                NodeIdx adjNi = faceToNode_[adjF];
                assert(adjNi != kInvalidNode);

                // Portal: origin of h and origin of twin(h), in CCW order of fi.
                VertexIdx portalA = he.origin;
                VertexIdx portalB = dcel.halfEdge(he.twin).origin;

                Arc arc;
                arc.to = adjNi;
                arc.cost = nd.centroid.dist(nodes_[adjNi].centroid);
                arc.portal = {portalA, portalB};

                nd.neighbors.push_back(arc);
            }
        }
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  clear
    // ─────────────────────────────────────────────────────────────────────────────
    void NavMesh::clear() noexcept {
        nodes_.clear();
        faceToNode_.clear();
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  Accessors
    // ─────────────────────────────────────────────────────────────────────────────
    std::size_t NavMesh::nodeCount() const noexcept {
        return nodes_.size();
    }

    const Node &NavMesh::node(NodeIdx i) const {
        assert(i < nodes_.size() && "NavMesh::node — index out of range");
        return nodes_[i];
    }

    NodeIdx NavMesh::findNode(Point2D p, const DCEL &dcel) const {
        FaceIdx fi = dcel.findFace(p);
        if (fi == kInvalidIdx) return kInvalidNode;
        if (fi >= faceToNode_.size()) return kInvalidNode;
        return faceToNode_[fi];
    }

    Point2D NavMesh::centroid(NodeIdx i) const {
        return node(i).centroid;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  portalPoints
    // ─────────────────────────────────────────────────────────────────────────────
    std::optional<std::array<Point2D, 2> >
    NavMesh::portalPoints(NodeIdx from, NodeIdx to, const DCEL &dcel) const {
        if (from >= nodes_.size()) return std::nullopt;

        for (const Arc &arc: nodes_[from].neighbors) {
            if (arc.to != to) continue;

            const Point2D &left = dcel.vertex(arc.portal[0]).pos;
            const Point2D &right = dcel.vertex(arc.portal[1]).pos;
            return std::array<Point2D, 2>{left, right};
        }
        return std::nullopt;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  Statistics
    // ─────────────────────────────────────────────────────────────────────────────
    std::size_t NavMesh::arcCount() const noexcept {
        std::size_t total = 0;
        for (const Node &nd: nodes_) total += nd.neighbors.size();
        return total;
    }

    double NavMesh::averageDegree() const noexcept {
        if (nodes_.empty()) return 0.0;
        return static_cast<double>(arcCount()) /
               static_cast<double>(nodes_.size());
    }
} // namespace navmesh
