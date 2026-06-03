#pragma once

/**
 * @file svg_writer.hh
 * @brief SVG renderer for navmesh triangulation and pathfinding results.
 *
 * Produces a self-contained SVG file showing:
 *   - All triangulation edges (thin grey lines)
 *   - Triangle faces (lightly filled, colour-coded per algorithm)
 *   - Input point cloud (small circles)
 *   - Start / goal markers
 *   - Raw A* triangle corridor (semi-transparent highlight)
 *   - Smoothed funnel path (bold coloured polyline with waypoint dots)
 *   - A legend and title box
 *
 * No external dependencies — pure C++ standard library.
 *
 * ### Coordinate mapping
 *
 * Input coordinates are in "world space" (arbitrary doubles).
 * The writer maps them into SVG viewport with a configurable margin,
 * flipping Y so that Y-up world coords render correctly in SVG (Y-down).
 *
 * ### Usage
 * @code
 *   SvgWriter w(800, 800);
 *   w.setTitle("Delaunay triangulation — labyrinth.txt");
 *   w.setWorldBounds(minX, minY, maxX, maxY);
 *   w.drawTriangulation(dcel, SvgWriter::Palette::Delaunay);
 *   w.drawPoints(pts);
 *   w.drawCorridor(corridorFaces, dcel, nm);
 *   w.drawPath(waypoints, SvgWriter::Palette::Delaunay);
 *   w.drawStartGoal(start, goal);
 *   w.drawLegend(metrics);
 *   w.save("out_delaunay.svg");
 * @endcode
 */

#include "geometry/dcel.hh"
#include "geometry/primitives.hpp"
#include "navmesh/navmesh.hh"

#include <string>
#include <string_view>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
//  LegendMetrics  — data shown in the info box
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Scalar metrics displayed in the SVG legend box.
 */
struct LegendMetrics {
    std::string algoName;       ///< e.g. "Delaunay (Bowyer–Watson)"
    std::size_t triangleCount{0};
    double      totalWeight{0.0};
    double      buildMs{0.0};
    double      pathLength{0.0};
    std::size_t turnCount{0};
    bool        pathFound{false};
};

// ─────────────────────────────────────────────────────────────────────────────
//  SvgWriter
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Builds and saves an SVG image of a triangulation + path.
 *
 * All draw* calls accumulate SVG elements in an internal string buffer.
 * Call save() (or str()) to obtain the final document.
 *
 * Draw order matters: elements added later appear on top.
 * Recommended call order:
 *   1. drawTriangulation()
 *   2. drawCorridor()      (optional)
 *   3. drawPoints()
 *   4. drawPath()
 *   5. drawStartGoal()
 *   6. drawLegend()
 */
class SvgWriter {
public:
    // ── Visual theme per algorithm ────────────────────────────────────────────

    /**
     * @brief Per-algorithm colour palette.
     *
     * Each algorithm gets a distinct hue used for face fill, path stroke,
     * and legend label.
     */
    enum class Palette {
        Delaunay,   ///< Blue tones
        Greedy,     ///< Green tones
        Quasi       ///< Orange tones
    };

    // ── Construction ──────────────────────────────────────────────────────────

    /**
     * @brief Creates an SVG canvas of @p width × @p height pixels.
     *
     * @param width   SVG viewport width in pixels  (default 900).
     * @param height  SVG viewport height in pixels (default 900).
     * @param margin  Padding around the drawing area in pixels (default 50).
     */
    explicit SvgWriter(int width  = 900,
                       int height = 900,
                       int margin = 60);

    // ── Configuration ─────────────────────────────────────────────────────────

    /**
     * @brief Sets the title shown at the top of the SVG.
     * @param title  UTF-8 string.
     */
    void setTitle(std::string_view title);

    /**
     * @brief Defines the world-space bounding box used for coordinate mapping.
     *
     * Must be called before any draw* method.
     * Typically derived from the input point set bounding box.
     *
     * @param minX, minY  Lower-left corner in world space.
     * @param maxX, maxY  Upper-right corner in world space.
     */
    void setWorldBounds(double minX, double minY,
                        double maxX, double maxY);

    // ── Drawing primitives ────────────────────────────────────────────────────

    /**
     * @brief Draws all live triangles and edges of a DCEL triangulation.
     *
     * Each triangle is filled with a very light tint of the palette colour.
     * Edges are drawn as thin grey strokes.
     *
     * @param dcel     Source triangulation.
     * @param palette  Colour theme (one per algorithm).
     */
    void drawTriangulation(const geometry::DCEL& dcel, Palette palette);

    /**
     * @brief Draws the input point cloud as small filled circles.
     *
     * @param pts    Points to draw.
     * @param radius Circle radius in pixels (default 3).
     */
    void drawPoints(const std::vector<geometry::Point2D>& pts,
                    double radius = 3.0);

    /**
     * @brief Highlights the A* triangle corridor as semi-transparent overlays.
     *
     * Each corridor triangle is filled with a pale yellow tint.
     *
     * @param corridor  Node index sequence from A*.
     * @param nm        NavMesh for face lookup.
     * @param dcel      DCEL for vertex positions.
     */
    void drawCorridor(const std::vector<navmesh::NodeIdx>& corridor,
                      const navmesh::NavMesh&              nm,
                      const geometry::DCEL&                dcel);

    /**
     * @brief Draws the smoothed funnel path as a bold polyline.
     *
     * Waypoints are rendered as small filled circles; straight segments
     * as strokes in the palette colour.
     *
     * @param waypoints  Ordered path points (start … goal).
     * @param palette    Colour theme.
     */
    void drawPath(const std::vector<geometry::Point2D>& waypoints,
                  Palette                               palette);

    /**
     * @brief Draws start (green circle + "S") and goal (red circle + "G") markers.
     *
     * @param start  World-space start position.
     * @param goal   World-space goal  position.
     */
    void drawStartGoal(geometry::Point2D start, geometry::Point2D goal);

    /**
     * @brief Draws a legend box in the lower-left corner with algorithm metrics.
     *
     * @param m  Metrics to display.
     */
    void drawLegend(const LegendMetrics& m);

    // ── Output ────────────────────────────────────────────────────────────────

    /**
     * @brief Returns the complete SVG document as a string.
     */
    [[nodiscard]] std::string str() const;

    /**
     * @brief Writes the SVG document to @p path.
     *
     * @throws std::runtime_error if the file cannot be opened.
     */
    void save(const std::string& path) const;

private:
    int    width_;
    int    height_;
    int    margin_;
    double worldMinX_{0}, worldMinY_{0};
    double worldMaxX_{1}, worldMaxY_{1};
    std::string title_;
    std::string body_;   ///< Accumulated SVG element markup

    // ── Coordinate helpers ────────────────────────────────────────────────────

    /// Maps world X to SVG X (left-to-right).
    [[nodiscard]] double toSvgX(double wx) const noexcept;
    /// Maps world Y to SVG Y (flips vertical axis).
    [[nodiscard]] double toSvgY(double wy) const noexcept;

    // ── Colour helpers ────────────────────────────────────────────────────────

    /// Returns the primary stroke colour for the given palette.
    [[nodiscard]] static std::string_view strokeColour(Palette p) noexcept;
    /// Returns the light face-fill colour for the given palette.
    [[nodiscard]] static std::string_view fillColour(Palette p) noexcept;
};
