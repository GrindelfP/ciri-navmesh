/**
 * @file vis_main.cc
 * @brief CLI entry point for the navmesh SVG visualiser.
 *
 * Reads the same input format as navmesh_bench and produces one SVG file
 * per triangulation algorithm:
 *
 *   <prefix>_delaunay.svg
 *   <prefix>_greedy.svg
 *   <prefix>_quasi.svg
 *
 * Each SVG shows:
 *   - All triangle faces (lightly filled) and edges
 *   - Input point cloud
 *   - A* triangle corridor (yellow overlay)
 *   - Smoothed funnel path (coloured bold polyline)
 *   - Start (S) and goal (G) markers
 *   - A legend box with algorithm metrics
 *
 * ## Usage
 *
 * @code
 *   navmesh_vis [OPTIONS] <input_file>
 *
 *   Options:
 *     -s, --start <x,y>   Start point (default: left-most point)
 *     -g, --goal  <x,y>   Goal  point (default: right-most point)
 *     -n, --no-flip        Disable edge flips in quasi-greedy
 *     -o, --output <pfx>   Output filename prefix (default: "utils")
 *     -W, --width  <px>    SVG width  in pixels (default: 900)
 *     -H, --height <px>    SVG height in pixels (default: 900)
 *     -h, --help           Print this help
 * @endcode
 */

#include "svg_writer.hh"

#include "geometry/primitives.hpp"
#include "geometry/dcel.hh"
#include "triangulation/delaunay.hh"
#include "triangulation/greedy.hh"
#include "triangulation/quasi_greedy.hh"
#include "navmesh/navmesh.hh"
#include "navmesh/astar.hh"
#include "navmesh/funnel.hh"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
//  Config
// ─────────────────────────────────────────────────────────────────────────────

struct Config {
    std::string inputFile;
    std::optional<geometry::Point2D> start;
    std::optional<geometry::Point2D> goal;
    std::string outputPrefix{"utils"};
    int svgWidth{900};
    int svgHeight{900};
    bool doFlips{true};
};

static geometry::Point2D parsePoint(const std::string& s) {
    const auto c = s.find(',');
    if (c == std::string::npos)
        throw std::invalid_argument("Expected x,y: " + s);
    return {std::stod(s.substr(0, c)), std::stod(s.substr(c + 1))};
}

static void printHelp(const std::string& prog) {
    std::cout
        << "Usage: " << prog << " [OPTIONS] <input_file>\n\n"
        << "Options:\n"
        << "  -s, --start <x,y>    Start point\n"
        << "  -g, --goal  <x,y>    Goal  point\n"
        << "  -n, --no-flip         Disable edge flips in quasi-greedy\n"
        << "  -o, --output <pfx>    Output prefix (default: utils)\n"
        << "  -W, --width  <px>     SVG width  (default: 900)\n"
        << "  -H, --height <px>     SVG height (default: 900)\n"
        << "  -h, --help            Print this help\n";
}

static std::optional<Config> parseArgs(int argc, char* argv[]) {
    if (argc < 2) {
        printHelp(argv[0]);
        throw std::invalid_argument("No input file.");
    }
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (++i >= argc) throw std::invalid_argument("Missing value for " + a);
            return argv[i];
        };
        if (a == "-h" || a == "--help") { printHelp(argv[0]); return std::nullopt; }
        else if (a == "-s" || a == "--start")  cfg.start  = parsePoint(next());
        else if (a == "-g" || a == "--goal")   cfg.goal   = parsePoint(next());
        else if (a == "-o" || a == "--output") cfg.outputPrefix = next();
        else if (a == "-W" || a == "--width")  cfg.svgWidth  = std::stoi(next());
        else if (a == "-H" || a == "--height") cfg.svgHeight = std::stoi(next());
        else if (a == "-n" || a == "--no-flip") cfg.doFlips = false;
        else if (a[0] == '-') throw std::invalid_argument("Unknown flag: " + a);
        else cfg.inputFile = a;
    }
    if (cfg.inputFile.empty()) throw std::invalid_argument("No input file.");
    return cfg;
}

// ─────────────────────────────────────────────────────────────────────────────
//  readPoints  (same as in main.cc)
// ─────────────────────────────────────────────────────────────────────────────

static std::vector<geometry::Point2D> readPoints(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open: " + path);

    std::vector<geometry::Point2D> pts;
    std::string line;
    int ln = 0;
    while (std::getline(in, line)) {
        ++ln;
        const auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#') continue;
        std::istringstream ss(line);
        double x{}, y{};
        if (!(ss >> x >> y))
            throw std::runtime_error("Parse error line " + std::to_string(ln));
        pts.push_back({x, y});
    }
    std::sort(pts.begin(), pts.end());
    pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
    if (pts.size() < 3)
        throw std::runtime_error("Need >= 3 distinct points.");
    return pts;
}

// ─────────────────────────────────────────────────────────────────────────────
//  defaultStartGoal
// ─────────────────────────────────────────────────────────────────────────────

static std::pair<geometry::Point2D, geometry::Point2D>
defaultStartGoal(const std::vector<geometry::Point2D>& pts) {
    auto minX = *std::min_element(pts.begin(), pts.end(),
                    [](auto& a, auto& b){ return a.x < b.x; });
    auto maxX = *std::max_element(pts.begin(), pts.end(),
                    [](auto& a, auto& b){ return a.x < b.x; });
    constexpr double kIn = 1e-3;
    return {{minX.x + kIn, minX.y}, {maxX.x - kIn, maxX.y}};
}

// ─────────────────────────────────────────────────────────────────────────────
//  worldBounds
// ─────────────────────────────────────────────────────────────────────────────

static void worldBounds(const std::vector<geometry::Point2D>& pts,
                         double& minX, double& minY,
                         double& maxX, double& maxY)
{
    minX = minY = std::numeric_limits<double>::infinity();
    maxX = maxY = -std::numeric_limits<double>::infinity();
    for (const auto& p : pts) {
        minX = std::min(minX, p.x); maxX = std::max(maxX, p.x);
        minY = std::min(minY, p.y); maxY = std::max(maxY, p.y);
    }
    // Ensure square-ish aspect: add a 5 % margin on the shorter axis.
    const double dx = maxX - minX, dy = maxY - minY;
    const double pad = std::max(dx, dy) * 0.05;
    minX -= pad; minY -= pad; maxX += pad; maxY += pad;
}

// ─────────────────────────────────────────────────────────────────────────────
//  renderOne  — triangulate + path-find + write SVG
// ─────────────────────────────────────────────────────────────────────────────

static void renderOne(triangulation::ITriangulator&         tri,
                      SvgWriter::Palette                    palette,
                      const std::vector<geometry::Point2D>& pts,
                      geometry::Point2D                     start,
                      geometry::Point2D                     goal,
                      double minX, double minY,
                      double maxX, double maxY,
                      const Config&                         cfg,
                      const std::string&                    svgPath)
{
    // ── Triangulate ───────────────────────────────────────────────────────────
    geometry::DCEL dcel;
    const auto triRes = tri.triangulate(pts, dcel);

    // ── NavMesh ───────────────────────────────────────────────────────────────
    navmesh::NavMesh nm;
    nm.build(dcel);

    // ── A* + Funnel ───────────────────────────────────────────────────────────
    LegendMetrics legend;
    legend.algoName      = std::string(tri.name());
    legend.triangleCount = triRes.triangleCount;
    legend.totalWeight   = triRes.totalWeight;
    legend.buildMs       = triRes.elapsed.count() * 1000.0;

    std::vector<navmesh::NodeIdx> corridor;

    const navmesh::NodeIdx sNode = nm.findNode(start, dcel);
    const navmesh::NodeIdx gNode = nm.findNode(goal,  dcel);

    std::vector<geometry::Point2D> waypoints;

    if (sNode != navmesh::kInvalidNode && gNode != navmesh::kInvalidNode) {
        pathfinding::AStar astar;
        auto aRes = astar.search(nm, sNode, gNode);

        if (aRes.found) {
            corridor = aRes.nodes;

            pathfinding::FunnelAlgorithm funnel;
            auto smooth = funnel.smooth(start, goal, aRes.nodes, nm, dcel);

            waypoints        = smooth.waypoints;
            legend.pathLength = smooth.totalLength;
            legend.turnCount  = smooth.turnCount;
            legend.pathFound  = true;
        }
    }

    // ── Build SVG ─────────────────────────────────────────────────────────────
    SvgWriter w(cfg.svgWidth, cfg.svgHeight);
    w.setTitle(std::string(tri.name()) + "  —  " + cfg.inputFile);
    w.setWorldBounds(minX, minY, maxX, maxY);

    w.drawTriangulation(dcel, palette);

    if (!corridor.empty())
        w.drawCorridor(corridor, nm, dcel);

    w.drawPoints(pts);

    if (!waypoints.empty())
        w.drawPath(waypoints, palette);

    w.drawStartGoal(start, goal);
    w.drawLegend(legend);

    w.save(svgPath);

    std::cout << "  Written: " << svgPath << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
//  main
// ─────────────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[])
{
    // ── 1. Parse args ─────────────────────────────────────────────────────────
    std::optional<Config> cfgOpt;
    try {
        cfgOpt = parseArgs(argc, argv);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Argument error: " << e.what() << "\n";
        return 1;
    }
    if (!cfgOpt) return 0;
    const Config& cfg = *cfgOpt;

    // ── 2. Read points ────────────────────────────────────────────────────────
    std::vector<geometry::Point2D> pts;
    try {
        pts = readPoints(cfg.inputFile);
    } catch (const std::runtime_error& e) {
        std::cerr << "Input error: " << e.what() << "\n";
        return 2;
    }
    std::cout << "Loaded " << pts.size() << " points from " << cfg.inputFile << "\n";

    // ── 3. Start / goal ───────────────────────────────────────────────────────
    auto [start, goal] = defaultStartGoal(pts);
    if (cfg.start) start = *cfg.start;
    if (cfg.goal)  goal  = *cfg.goal;
    std::cout << "Start: (" << start.x << ", " << start.y << ")"
              << "  Goal: ("  << goal.x  << ", " << goal.y  << ")\n\n";

    // ── 4. World bounds ───────────────────────────────────────────────────────
    double minX, minY, maxX, maxY;
    worldBounds(pts, minX, minY, maxX, maxY);

    // ── 5. Render each algorithm ──────────────────────────────────────────────
    struct AlgoEntry {
        std::unique_ptr<triangulation::ITriangulator> tri;
        SvgWriter::Palette                            palette;
        std::string                                   suffix;
    };

    // We can't copy unique_ptr, so we build and render inline.
    auto run = [&](triangulation::ITriangulator& t,
                   SvgWriter::Palette p,
                   const std::string& suffix)
    {
        const std::string path = cfg.outputPrefix + "_" + suffix + ".svg";
        try {
            renderOne(t, p, pts, start, goal,
                      minX, minY, maxX, maxY, cfg, path);
        } catch (const std::exception& e) {
            std::cerr << "  [" << t.name() << "] FAILED: " << e.what() << "\n";
        }
    };

    {
        triangulation::DelaunayTriangulator dt;
        run(dt, SvgWriter::Palette::Delaunay, "delaunay");
    }
    {
        triangulation::GreedyTriangulator gt;
        run(gt, SvgWriter::Palette::Greedy, "greedy");
    }
    {
        triangulation::QuasiGreedyTriangulator qg(cfg.doFlips);
        run(qg, SvgWriter::Palette::Quasi, "quasi");
    }

    std::cout << "\nDone. Open the .svg files in any browser or SVG viewer.\n";
    return 0;
}
