/**
 * @file svg_writer.cc
 * @brief SVG renderer implementation for navmesh triangulation visualisation.
 */

#include "svg_writer.hh"

#include <array>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

// ─────────────────────────────────────────────────────────────────────────────
//  Construction
// ─────────────────────────────────────────────────────────────────────────────

SvgWriter::SvgWriter(int width, int height, int margin)
    : width_(width), height_(height), margin_(margin)
{}

// ─────────────────────────────────────────────────────────────────────────────
//  Configuration
// ─────────────────────────────────────────────────────────────────────────────

void SvgWriter::setTitle(std::string_view title) {
    title_ = std::string(title);
}

void SvgWriter::setWorldBounds(double minX, double minY,
                                double maxX, double maxY)
{
    // Guard against degenerate bounds.
    const double kPad = 1e-6;
    worldMinX_ = minX - kPad;
    worldMinY_ = minY - kPad;
    worldMaxX_ = maxX + kPad;
    worldMaxY_ = maxY + kPad;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Coordinate mapping
// ─────────────────────────────────────────────────────────────────────────────

double SvgWriter::toSvgX(double wx) const noexcept {
    const double drawW = static_cast<double>(width_  - 2 * margin_);
    const double span  = worldMaxX_ - worldMinX_;
    return margin_ + (wx - worldMinX_) / span * drawW;
}

double SvgWriter::toSvgY(double wy) const noexcept {
    // SVG Y grows downward; world Y grows upward → flip.
    const double drawH = static_cast<double>(height_ - 2 * margin_);
    const double span  = worldMaxY_ - worldMinY_;
    return (height_ - margin_) - (wy - worldMinY_) / span * drawH;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Colour helpers
// ─────────────────────────────────────────────────────────────────────────────

std::string_view SvgWriter::strokeColour(Palette p) noexcept {
    switch (p) {
        case Palette::Delaunay: return "#1a6bbf";   // blue
        case Palette::Greedy:   return "#217a3c";   // green
        case Palette::Quasi:    return "#c46200";   // orange
    }
    return "#333333";
}

std::string_view SvgWriter::fillColour(Palette p) noexcept {
    switch (p) {
        case Palette::Delaunay: return "#ddeeff";
        case Palette::Greedy:   return "#ddf5e4";
        case Palette::Quasi:    return "#fff0d8";
    }
    return "#f0f0f0";
}

// ─────────────────────────────────────────────────────────────────────────────
//  Internal SVG helpers (lambdas / free functions in anonymous namespace)
// ─────────────────────────────────────────────────────────────────────────────

namespace {

/// Formats a double with fixed precision into a string.
std::string fmt(double v, int prec = 2) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(prec) << v;
    return ss.str();
}

/// Escapes & < > in a string for safe embedding in SVG text.
std::string xmlEscape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        switch (c) {
            case '&': out += "&amp;";  break;
            case '<': out += "&lt;";   break;
            case '>': out += "&gt;";   break;
            default:  out += c;
        }
    }
    return out;
}

} // anonymous namespace

// ─────────────────────────────────────────────────────────────────────────────
//  drawTriangulation
// ─────────────────────────────────────────────────────────────────────────────

void SvgWriter::drawTriangulation(const geometry::DCEL& dcel, Palette palette)
{
    std::ostringstream ss;

    const std::string fCol = std::string(fillColour(palette));
    const std::string eCol = "#8899aa";   // neutral grey for edges

    // ── Filled triangles: draw large first so small ones render on top ────────
    struct TriSvg { double ax,ay,bx,by,cx,cy,area; };
    std::vector<TriSvg> tris;
    tris.reserve(dcel.liveTriangleCount());
    for (geometry::FaceIdx fi = 1; fi < dcel.faceCount(); ++fi) {
        if (dcel.face(fi).dead) continue;
        const auto [va, vb, vc] = dcel.faceVertices(fi);
        const double ax = toSvgX(dcel.vertex(va).pos.x);
        const double ay = toSvgY(dcel.vertex(va).pos.y);
        const double bx = toSvgX(dcel.vertex(vb).pos.x);
        const double by = toSvgY(dcel.vertex(vb).pos.y);
        const double cx = toSvgX(dcel.vertex(vc).pos.x);
        const double cy = toSvgY(dcel.vertex(vc).pos.y);
        const double area = std::abs((bx-ax)*(cy-ay)-(cx-ax)*(by-ay))*0.5;
        tris.push_back({ax,ay,bx,by,cx,cy,area});
    }
    std::sort(tris.begin(), tris.end(),
              [](const TriSvg& a, const TriSvg& b){ return a.area > b.area; });
    ss << "  <!-- triangulation faces -->\n";
    for (const auto& t : tris) {
        ss << "  <polygon points=\"" << fmt(t.ax) << "," << fmt(t.ay)
           << " " << fmt(t.bx) << "," << fmt(t.by)
           << " " << fmt(t.cx) << "," << fmt(t.cy)
           << "\" fill=\"" << fCol << "\" stroke=\"none\"/>\n";
    }

    // ── Edges ─────────────────────────────────────────────────────────────────
    ss << "  <!-- triangulation edges -->\n";
    for (geometry::HalfEdgeIdx hi = 0; hi < dcel.halfEdgeCount(); ++hi) {
        if (!dcel.isHalfEdgeLive(hi)) continue;
        const geometry::HalfEdge& he = dcel.halfEdge(hi);
        if (he.dead) continue;
        if (hi > he.twin) continue;                  // each edge once
        if (he.face == geometry::kOuterFace &&
            dcel.halfEdge(he.twin).face == geometry::kOuterFace) continue;

        const double x1 = toSvgX(dcel.vertex(he.origin).pos.x);
        const double y1 = toSvgY(dcel.vertex(he.origin).pos.y);
        const double x2 = toSvgX(dcel.vertex(dcel.halfEdge(he.twin).origin).pos.x);
        const double y2 = toSvgY(dcel.vertex(dcel.halfEdge(he.twin).origin).pos.y);

        ss << "  <line x1=\"" << fmt(x1) << "\" y1=\"" << fmt(y1)
           << "\" x2=\"" << fmt(x2) << "\" y2=\"" << fmt(y2)
           << "\" stroke=\"" << eCol
           << "\" stroke-width=\"0.7\" opacity=\"0.7\"/>\n";
    }

    body_ += ss.str();
}

// ─────────────────────────────────────────────────────────────────────────────
//  drawPoints
// ─────────────────────────────────────────────────────────────────────────────

void SvgWriter::drawPoints(const std::vector<geometry::Point2D>& pts,
                            double radius)
{
    std::ostringstream ss;
    ss << "  <!-- input points -->\n";
    for (const auto& p : pts) {
        ss << "  <circle cx=\"" << fmt(toSvgX(p.x))
           << "\" cy=\"" << fmt(toSvgY(p.y))
           << "\" r=\"" << fmt(radius)
           << "\" fill=\"#445566\" stroke=\"none\" opacity=\"0.8\"/>\n";
    }
    body_ += ss.str();
}

// ─────────────────────────────────────────────────────────────────────────────
//  drawCorridor
// ─────────────────────────────────────────────────────────────────────────────

void SvgWriter::drawCorridor(const std::vector<navmesh::NodeIdx>& corridor,
                              const navmesh::NavMesh&              nm,
                              const geometry::DCEL&                dcel)
{
    if (corridor.empty()) return;

    std::ostringstream ss;
    ss << "  <!-- A* corridor -->\n";

    for (navmesh::NodeIdx ni : corridor) {
        const navmesh::Node& nd = nm.node(ni);
        const auto [va, vb, vc] = dcel.faceVertices(nd.faceIdx);

        const double ax = toSvgX(dcel.vertex(va).pos.x);
        const double ay = toSvgY(dcel.vertex(va).pos.y);
        const double bx = toSvgX(dcel.vertex(vb).pos.x);
        const double by = toSvgY(dcel.vertex(vb).pos.y);
        const double cx = toSvgX(dcel.vertex(vc).pos.x);
        const double cy = toSvgY(dcel.vertex(vc).pos.y);

        ss << "  <polygon points=\""
           << fmt(ax) << "," << fmt(ay) << " "
           << fmt(bx) << "," << fmt(by) << " "
           << fmt(cx) << "," << fmt(cy)
           << "\" fill=\"#ffee55\" stroke=\"none\" opacity=\"0.35\"/>\n";
    }

    body_ += ss.str();
}

// ─────────────────────────────────────────────────────────────────────────────
//  drawPath
// ─────────────────────────────────────────────────────────────────────────────

void SvgWriter::drawPath(const std::vector<geometry::Point2D>& waypoints,
                          Palette palette)
{
    if (waypoints.size() < 2) return;

    std::ostringstream ss;
    const std::string col = std::string(strokeColour(palette));

    // ── Polyline ──────────────────────────────────────────────────────────────
    ss << "  <!-- smoothed path -->\n"
       << "  <polyline points=\"";
    for (const auto& p : waypoints) {
        ss << fmt(toSvgX(p.x)) << "," << fmt(toSvgY(p.y)) << " ";
    }
    ss << "\" fill=\"none\" stroke=\"" << col
       << "\" stroke-width=\"2.5\" stroke-linejoin=\"round\""
       << " stroke-linecap=\"round\"/>\n";

    // ── Waypoint dots (interior only) ─────────────────────────────────────────
    for (std::size_t i = 1; i + 1 < waypoints.size(); ++i) {
        ss << "  <circle cx=\"" << fmt(toSvgX(waypoints[i].x))
           << "\" cy=\"" << fmt(toSvgY(waypoints[i].y))
           << "\" r=\"3.5\" fill=\"" << col
           << "\" stroke=\"white\" stroke-width=\"1\"/>\n";
    }

    body_ += ss.str();
}

// ─────────────────────────────────────────────────────────────────────────────
//  drawStartGoal
// ─────────────────────────────────────────────────────────────────────────────

void SvgWriter::drawStartGoal(geometry::Point2D start, geometry::Point2D goal)
{
    std::ostringstream ss;
    ss << "  <!-- start / goal markers -->\n";

    auto marker = [&](geometry::Point2D p, const char* colour,
                      const char* label)
    {
        const double cx = toSvgX(p.x);
        const double cy = toSvgY(p.y);
        ss << "  <circle cx=\"" << fmt(cx) << "\" cy=\"" << fmt(cy)
           << "\" r=\"9\" fill=\"" << colour
           << "\" stroke=\"white\" stroke-width=\"2\"/>\n"
           << "  <text x=\"" << fmt(cx) << "\" y=\"" << fmt(cy + 4.5)
           << "\" text-anchor=\"middle\" font-family=\"monospace\""
           << " font-size=\"11\" font-weight=\"bold\" fill=\"white\">"
           << label << "</text>\n";
    };

    marker(start, "#2aa832", "S");
    marker(goal,  "#cc2222", "G");

    body_ += ss.str();
}

// ─────────────────────────────────────────────────────────────────────────────
//  drawLegend
// ─────────────────────────────────────────────────────────────────────────────

void SvgWriter::drawLegend(const LegendMetrics& m)
{
    const int bx = margin_;
    const int by = height_ - margin_ - 120;
    const int bw = 280;
    const int bh = 110;

    std::ostringstream ss;
    ss << "  <!-- legend box -->\n"
       << "  <rect x=\"" << bx << "\" y=\"" << by
       << "\" width=\"" << bw << "\" height=\"" << bh
       << "\" fill=\"white\" fill-opacity=\"0.88\""
       << " stroke=\"#aabbcc\" stroke-width=\"1\" rx=\"6\"/>\n";

    // Helper: one text line
    auto line = [&](int row, const std::string& txt, bool bold = false) {
        const int y = by + 20 + row * 18;
        ss << "  <text x=\"" << (bx + 10) << "\" y=\"" << y
           << "\" font-family=\"monospace\" font-size=\"12\"";
        if (bold) ss << " font-weight=\"bold\"";
        ss << " fill=\"#222\">" << xmlEscape(txt) << "</text>\n";
    };

    line(0, m.algoName, /*bold=*/true);
    line(1, "Triangles : " + std::to_string(m.triangleCount));
    line(2, "Weight    : " + fmt(m.totalWeight, 2));
    line(3, "Build     : " + fmt(m.buildMs, 1) + " ms");
    if (m.pathFound) {
        line(4, "PathLen   : " + fmt(m.pathLength, 2));
        line(5, "Turns     : " + std::to_string(m.turnCount));
    } else {
        line(4, "Path      : N/A");
    }

    body_ += ss.str();
}

// ─────────────────────────────────────────────────────────────────────────────
//  str / save
// ─────────────────────────────────────────────────────────────────────────────

std::string SvgWriter::str() const
{
    std::ostringstream doc;

    // ── SVG header ────────────────────────────────────────────────────────────
    doc << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<svg xmlns=\"http://www.w3.org/2000/svg\""
        << " width=\""  << width_  << "\""
        << " height=\"" << height_ << "\""
        << " viewBox=\"0 0 " << width_ << " " << height_ << "\">\n"

        // Background
        << "  <rect width=\"" << width_ << "\" height=\"" << height_
        << "\" fill=\"#f7f9fb\"/>\n"

        // Title
        << "  <text x=\"" << (width_ / 2) << "\" y=\"" << (margin_ / 2 + 8)
        << "\" text-anchor=\"middle\" font-family=\"sans-serif\""
        << " font-size=\"15\" font-weight=\"bold\" fill=\"#334\">"
        << xmlEscape(title_) << "</text>\n";

    // ── Body (all accumulated elements) ──────────────────────────────────────
    doc << body_;

    // ── Footer ────────────────────────────────────────────────────────────────
    doc << "</svg>\n";

    return doc.str();
}

void SvgWriter::save(const std::string& path) const
{
    std::ofstream f(path);
    if (!f) throw std::runtime_error("SvgWriter: cannot write to " + path);
    f << str();
}
