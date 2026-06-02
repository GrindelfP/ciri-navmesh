/**
 * @file main.cc
 * @brief CLI entry point for the navmesh triangulation benchmark.
 *
 * ## Usage
 *
 * @code
 *   navmesh_bench [OPTIONS] <input_file>
 *
 *   Options:
 *     -a, --algo <name>    Triangulation algorithm: delaunay | greedy | quasi
 *                          (default: all three, runs benchmark comparison)
 *     -s, --start <x,y>   Start point for pathfinding (default: centroid of
 *                          a random triangle)
 *     -g, --goal  <x,y>   Goal  point for pathfinding (default: centroid of
 *                          a random triangle far from start)
 *     -n, --no-flip        Disable edge-flip post-processing in quasi-greedy
 *     -q, --quiet          Suppress per-triangle output; print summary only
 *     -o, --output <file>  Write benchmark results as CSV to <file>
 *     -h, --help           Print this help message
 * @endcode
 *
 * ## Input format
 *
 * A plain-text file where each line contains two space-separated doubles:
 * @code
 *   0.0 0.0
 *   1.0 0.0
 *   0.5 0.866
 *   ...
 * @endcode
 * Lines starting with '#' are treated as comments and ignored.
 * Duplicate points are silently removed.
 *
 * ## Output
 *
 * When run in benchmark mode (no --algo specified), prints a comparison table:
 * @code
 *   Algorithm          | Triangles | Weight      | Build (ms) | PathLen | Turns
 *   -------------------|-----------|-------------|------------|---------|------
 *   Delaunay (BW)      |       198 |    4521.72  |      12.3  |   87.4  |   3
 *   Greedy MWT         |       198 |    4318.55  |      43.1  |   83.1  |   2
 *   Quasi-Greedy MWT   |       198 |    4289.11  |      15.8  |   81.9  |   2
 * @endcode
 *
 * ## Exit codes
 *
 *   0 — Success
 *   1 — Invalid arguments or missing input file
 *   2 — Input file parse error (e.g. malformed coordinates)
 *   3 — Triangulation failed (e.g. all points collinear)
 *   4 — Pathfinding failed (start or goal outside triangulation)
 */

#include "../lib/geometry/primitives.hpp"
#include "../lib/geometry/dcel.hh"
#include "../lib/triangulation/delaunay.hh"
#include "../lib/triangulation/greedy.hh"
#include "../lib/triangulation/quasi_greedy.hh"
#include "../lib/navmesh/navmesh.hh"
#include "astar.hh"
#include "funnel.hh"

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
//  CLI argument parsing
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Parsed command-line configuration.
 */
struct Config {
    std::string              inputFile;
    std::optional<std::string> algo;      ///< "delaunay" | "greedy" | "quasi" | nullopt (all)
    std::optional<geometry::Point2D> start;
    std::optional<geometry::Point2D> goal;
    std::optional<std::string> outputCsv;
    bool doFlips{true};
    bool quiet{false};
};

/**
 * @brief Parses a "x,y" string into a Point2D.
 *
 * @param s   Input string, e.g. "1.5,3.0".
 * @return    Parsed point.
 * @throws    std::invalid_argument on parse failure.
 */
static geometry::Point2D parsePoint(const std::string& s) {
    const auto comma = s.find(',');
    if (comma == std::string::npos)
        throw std::invalid_argument("Expected x,y format: " + s);
    return { std::stod(s.substr(0, comma)),
             std::stod(s.substr(comma + 1)) };
}

/**
 * @brief Prints usage information to @p os.
 *
 * @param prog  Program name (argv[0]).
 * @param os    Output stream (default: std::cout).
 */
static void printHelp(const std::string& prog, std::ostream& os = std::cout) {
    os << "Usage: " << prog << " [OPTIONS] <input_file>\n\n"
       << "Options:\n"
       << "  -a, --algo <name>    delaunay | greedy | quasi  (default: all)\n"
       << "  -s, --start <x,y>   Start point for pathfinding\n"
       << "  -g, --goal  <x,y>   Goal  point for pathfinding\n"
       << "  -n, --no-flip        Disable edge flips in quasi-greedy\n"
       << "  -q, --quiet          Print summary table only\n"
       << "  -o, --output <file>  Write CSV results to file\n"
       << "  -h, --help           Print this help\n";
}

/**
 * @brief Parses argv into a Config.
 *
 * @param argc  Argument count.
 * @param argv  Argument vector.
 * @return      Filled Config, or nullopt if --help was requested.
 * @throws      std::invalid_argument on unknown flags or missing values.
 */
static std::optional<Config> parseArgs(int argc, char* argv[]) {
    if (argc < 2) {
        printHelp(argv[0], std::cerr);
        throw std::invalid_argument("No input file specified.");
    }

    Config cfg;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto nextArg = [&]() -> std::string {
            if (++i >= argc) throw std::invalid_argument("Missing value for " + arg);
            return argv[i];
        };

        if (arg == "-h" || arg == "--help") {
            printHelp(argv[0]);
            return std::nullopt;
        } else if (arg == "-a" || arg == "--algo") {
            cfg.algo = nextArg();
        } else if (arg == "-s" || arg == "--start") {
            cfg.start = parsePoint(nextArg());
        } else if (arg == "-g" || arg == "--goal") {
            cfg.goal  = parsePoint(nextArg());
        } else if (arg == "-o" || arg == "--output") {
            cfg.outputCsv = nextArg();
        } else if (arg == "-n" || arg == "--no-flip") {
            cfg.doFlips = false;
        } else if (arg == "-q" || arg == "--quiet") {
            cfg.quiet = true;
        } else if (arg[0] == '-') {
            throw std::invalid_argument("Unknown flag: " + arg);
        } else {
            cfg.inputFile = arg;
        }
    }

    if (cfg.inputFile.empty())
        throw std::invalid_argument("No input file specified.");
    return cfg;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Input parsing
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Reads a point set from a plain-text file.
 *
 * Each non-comment line must contain exactly two doubles.
 * Lines starting with '#' (after stripping whitespace) are skipped.
 * Duplicate points (exact equality) are removed.
 *
 * @param path  File path.
 * @return      Deduplicated point set.
 * @throws      std::runtime_error on file-open failure or parse error.
 */
static std::vector<geometry::Point2D> readPoints(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open file: " + path);

    std::vector<geometry::Point2D> pts;
    std::string line;
    int lineNum = 0;
    while (std::getline(in, line)) {
        ++lineNum;
        // Strip leading whitespace for comment detection.
        const auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#') continue;

        std::istringstream ss(line);
        double x{}, y{};
        if (!(ss >> x >> y))
            throw std::runtime_error("Parse error at line " +
                                     std::to_string(lineNum) + ": " + line);
        pts.push_back({x, y});
    }

    // Deduplicate (exact equality — duplicates are rare but crash DCEL).
    std::sort(pts.begin(), pts.end());
    pts.erase(std::unique(pts.begin(), pts.end()), pts.end());

    if (pts.size() < 3)
        throw std::runtime_error("Need at least 3 distinct points; got " +
                                 std::to_string(pts.size()));
    return pts;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Single-algorithm run
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Aggregate metrics for one algorithm run.
 */
struct RunMetrics {
    std::string                   algoName;
    triangulation::TriangulationResult tri;  ///< Triangulation timing + weight
    double                        pathLength{0.0};
    std::size_t                   turnCount{0};
    double                        pathfindMs{0.0};
    bool                          pathFound{false};
};

/**
 * @brief Runs one triangulation algorithm plus pathfinding and returns metrics.
 *
 * @param triangulator  The algorithm to use.
 * @param points        Input point set.
 * @param start         Start point for pathfinding.
 * @param goal          Goal  point for pathfinding.
 * @param quiet         If false, prints per-step progress.
 * @return              RunMetrics filled with results.
 */
static RunMetrics runOne(triangulation::ITriangulator&         triangulator,
                         const std::vector<geometry::Point2D>& points,
                         geometry::Point2D                     start,
                         geometry::Point2D                     goal,
                         bool                                  quiet) {
    RunMetrics m;
    m.algoName = std::string(triangulator.name());

    if (!quiet) std::cout << "  [" << m.algoName << "] Triangulating...";
    geometry::DCEL dcel;
    m.tri = triangulator.triangulate(points, dcel);
    if (!quiet) {
        std::cout << " done. Triangles=" << m.tri.triangleCount
                  << " Weight=" << std::fixed << std::setprecision(2)
                  << m.tri.totalWeight << " Build="
                  << std::setprecision(1)
                  << m.tri.elapsed.count() * 1000.0 << "ms\n";
    }

    // Build navmesh
    navmesh::NavMesh nm;
    nm.build(dcel);

    // Locate start/goal nodes
    const auto startNode = nm.findNode(start, dcel);
    const auto goalNode  = nm.findNode(goal,  dcel);
    if (startNode == navmesh::kInvalidNode || goalNode == navmesh::kInvalidNode) {
        if (!quiet) std::cout << "  [" << m.algoName
                              << "] WARNING: start or goal outside mesh.\n";
        return m;
    }

    // A* pathfinding
    auto t0 = std::chrono::steady_clock::now();
    pathfinding::AStar astar;
    auto aResult = astar.search(nm, startNode, goalNode);
    auto t1 = std::chrono::steady_clock::now();
    m.pathfindMs = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (!aResult.found) {
        if (!quiet) std::cout << "  [" << m.algoName << "] No path found.\n";
        return m;
    }
    m.pathFound = true;

    // Funnel smoothing
    pathfinding::FunnelAlgorithm funnel;
    auto smooth = funnel.smooth(start, goal, aResult.nodes, nm, dcel);
    m.pathLength = smooth.totalLength;
    m.turnCount  = smooth.turnCount;

    if (!quiet) {
        std::cout << "  [" << m.algoName << "] Path: len="
                  << std::setprecision(2) << m.pathLength
                  << " turns=" << m.turnCount
                  << " pathfind=" << std::setprecision(2) << m.pathfindMs << "ms\n";
    }
    return m;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Output formatting
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Prints a comparison table to @p os.
 *
 * @param results  One RunMetrics per algorithm.
 * @param os       Output stream.
 */
static void printTable(const std::vector<RunMetrics>& results,
                       std::ostream&                  os = std::cout) {
    // Column widths
    constexpr int W0 = 26, W1 = 11, W2 = 13, W3 = 12, W4 = 9, W5 = 7;
    const std::string sep(W0+W1+W2+W3+W4+W5+5*3, '-');

    auto cell = [&](const std::string& s, int w) {
        os << " " << std::left << std::setw(w) << s << " |";
    };

    os << "\n";
    cell("Algorithm",           W0); cell("Triangles", W1);
    cell("Weight",              W2); cell("Build (ms)", W3);
    cell("PathLen",             W4); cell("Turns",      W5);
    os << "\n" << sep << "\n";

    for (const auto& r : results) {
        cell(r.algoName, W0);
        cell(std::to_string(r.tri.triangleCount), W1);

        std::ostringstream w;
        w << std::fixed << std::setprecision(2) << r.tri.totalWeight;
        cell(w.str(), W2);

        std::ostringstream b;
        b << std::fixed << std::setprecision(1) << r.tri.elapsed.count()*1000.0;
        cell(b.str(), W3);

        if (r.pathFound) {
            std::ostringstream pl;
            pl << std::fixed << std::setprecision(2) << r.pathLength;
            cell(pl.str(), W4);
            cell(std::to_string(r.turnCount), W5);
        } else {
            cell("N/A", W4);
            cell("N/A", W5);
        }
        os << "\n";
    }
    os << "\n";
}

/**
 * @brief Writes results as CSV to @p path.
 *
 * CSV header:
 * @code
 *   algorithm,triangles,weight,build_ms,path_length,turns,pathfind_ms
 * @endcode
 *
 * @param results  One RunMetrics per algorithm.
 * @param path     Output file path.
 * @throws         std::runtime_error on file-open failure.
 */
static void writeCsv(const std::vector<RunMetrics>& results,
                     const std::string&             path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write CSV: " + path);

    f << "algorithm,triangles,weight,build_ms,path_length,turns,pathfind_ms\n";
    for (const auto& r : results) {
        f << r.algoName << ","
          << r.tri.triangleCount << ","
          << std::fixed << std::setprecision(6) << r.tri.totalWeight << ","
          << r.tri.elapsed.count() * 1000.0 << ","
          << r.pathLength << ","
          << r.turnCount << ","
          << r.pathfindMs << "\n";
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Default start/goal selection
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Picks a default start/goal pair from the point set if not specified.
 *
 * Chooses the bottom-left point as start and the top-right point as goal
 * (by bounding-box extremes), which generally gives an interesting path
 * across the triangulation.
 *
 * @param pts  Input point set.
 * @return     {start, goal} pair.
 */
static std::pair<geometry::Point2D, geometry::Point2D>
defaultStartGoal(const std::vector<geometry::Point2D>& pts) {
    auto minX = *std::min_element(pts.begin(), pts.end(),
                    [](auto& a, auto& b){ return a.x < b.x; });
    auto maxX = *std::max_element(pts.begin(), pts.end(),
                    [](auto& a, auto& b){ return a.x < b.x; });
    // Slight offset inward to avoid boundary issues.
    constexpr double kInset = 1e-3;
    auto start = geometry::Point2D{ minX.x + kInset, minX.y };
    auto goal  = geometry::Point2D{ maxX.x - kInset, maxX.y };
    return {start, goal};
}

// ─────────────────────────────────────────────────────────────────────────────
//  main()
// ─────────────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {
    // ── 1. Parse arguments ────────────────────────────────────────────────────
    std::optional<Config> cfgOpt;
    try {
        cfgOpt = parseArgs(argc, argv);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Argument error: " << e.what() << "\n";
        return 1;
    }
    if (!cfgOpt) return 0;  // --help was printed
    const Config& cfg = *cfgOpt;

    // ── 2. Read input ─────────────────────────────────────────────────────────
    std::vector<geometry::Point2D> pts;
    try {
        pts = readPoints(cfg.inputFile);
    } catch (const std::runtime_error& e) {
        std::cerr << "Input error: " << e.what() << "\n";
        return 2;
    }

    if (!cfg.quiet)
        std::cout << "Loaded " << pts.size() << " points from " << cfg.inputFile << "\n";

    // ── 3. Resolve start/goal ─────────────────────────────────────────────────
    auto [start, goal] = defaultStartGoal(pts);
    if (cfg.start) start = *cfg.start;
    if (cfg.goal)  goal  = *cfg.goal;

    if (!cfg.quiet)
        std::cout << "Start: (" << start.x << ", " << start.y << ")  "
                  << "Goal:  (" << goal.x  << ", " << goal.y  << ")\n\n";

    // ── 4. Run selected algorithms ────────────────────────────────────────────
    std::vector<RunMetrics> results;

    auto runAlgo = [&](triangulation::ITriangulator& algo) {
        try {
            results.push_back(runOne(algo, pts, start, goal, cfg.quiet));
        } catch (const std::exception& e) {
            std::cerr << "  [" << algo.name() << "] FAILED: " << e.what() << "\n";
        }
    };

    const bool runAll = !cfg.algo.has_value();
    const std::string algoStr = cfg.algo.value_or("");

    if (runAll || algoStr == "delaunay") {
        triangulation::DelaunayTriangulator dt;
        runAlgo(dt);
    }
    if (runAll || algoStr == "greedy") {
        triangulation::GreedyTriangulator gt;
        runAlgo(gt);
    }
    if (runAll || algoStr == "quasi") {
        triangulation::QuasiGreedyTriangulator qg(cfg.doFlips);
        runAlgo(qg);
    }

    if (results.empty()) {
        std::cerr << "No results produced. Check algorithm name and input.\n";
        return 3;
    }

    // ── 5. Print / save results ───────────────────────────────────────────────
    printTable(results);

    if (cfg.outputCsv) {
        try {
            writeCsv(results, *cfg.outputCsv);
            if (!cfg.quiet)
                std::cout << "CSV written to " << *cfg.outputCsv << "\n";
        } catch (const std::runtime_error& e) {
            std::cerr << "CSV write error: " << e.what() << "\n";
            return 1;
        }
    }

    return 0;
}
