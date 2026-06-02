#pragma once

/**
 * @file i_triangulator.hh
 * @brief Abstract interface for triangulation algorithms.
 *
 * Every triangulator takes a set of 2-D points and fills a DCEL with the
 * resulting triangulation.  The interface is intentionally thin: callers
 * deal with the DCEL directly for all subsequent operations (navmesh
 * construction, edge flips, metrics).
 *
 * ### Implementing a new triangulator
 *
 * 1. Inherit from ITriangulator.
 * 2. Override triangulate().
 * 3. Clear the target DCEL at the start (dcel.clear()) or document that
 *    you expect it to be empty.
 * 4. Leave face 0 as the outer face (invariant I7 from dcel.hh).
 *
 * ### Preconditions shared by all implementations
 *
 *   - Input must contain at least 3 non-collinear points.
 *   - Duplicate points produce undefined behaviour (filter upstream).
 *   - Points are in "reasonable" coordinate range (|x|, |y| < 1e6),
 *     consistent with geometry::kEps assumptions.
 */

#include "../geometry/dcel.hh"
#include "../geometry/primitives.hpp"

#include <chrono>
#include <string_view>

namespace triangulation {

// ─────────────────────────────────────────────────────────────────────────────
//  Result metadata returned by triangulate()
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Metadata collected during a single triangulate() call.
 *
 * All fields are filled by the triangulator implementation.
 * Clients use this for the comparison metrics described in README.md.
 */
struct TriangulationResult {
    /// Elapsed wall-clock time of the triangulation algorithm itself.
    std::chrono::duration<double> elapsed{};

    /// Total weight of the triangulation (sum of all edge lengths).
    /// Filled via dcel.totalWeight() after the DCEL is complete.
    double totalWeight{0.0};

    /// Number of live (bounded) triangles in the output DCEL.
    std::size_t triangleCount{0};

    /// Number of edge flips performed (non-zero only for algorithms that flip).
    std::size_t flipCount{0};
};

// ─────────────────────────────────────────────────────────────────────────────
//  ITriangulator
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Pure abstract base class for all triangulation algorithms.
 *
 * Ownership model: the caller owns the DCEL and passes it by reference.
 * This allows reusing the same DCEL object across multiple triangulations
 * (e.g. when switching algorithms interactively in the visualizer).
 */
class ITriangulator {
public:
    virtual ~ITriangulator() = default;

    // Non-copyable, non-movable (stateless algorithms are fine as singletons).
    ITriangulator(const ITriangulator&)            = delete;
    ITriangulator& operator=(const ITriangulator&) = delete;
    ITriangulator(ITriangulator&&)                 = delete;
    ITriangulator& operator=(ITriangulator&&)      = delete;

    /**
     * @brief Triangulates the given point set and fills @p dcel.
     *
     * The DCEL is cleared at the start of this call — any existing geometry
     * is discarded.
     *
     * @param points  View of input points.  Must have >= 3 non-collinear pts.
     *                The span need not remain valid after the call returns.
     * @param dcel    Output DCEL.  Will be in a valid, fully-linked state on
     *                return (validate() returns true).
     * @return        TriangulationResult with timing and metric data.
     *
     * @throws std::invalid_argument if fewer than 3 points are provided.
     * @note   Does NOT throw for collinear-only inputs; the result is an
     *         empty DCEL (0 bounded faces) and totalWeight == 0.
     */
    virtual TriangulationResult
    triangulate(const std::vector<geometry::Point2D>& points,
                geometry::DCEL&                       dcel) = 0;

    /**
     * @brief Human-readable name of the algorithm (e.g. "Delaunay (BW)").
     *
     * Used in HUD display and benchmark output.
     */
    [[nodiscard]] virtual std::string_view name() const noexcept = 0;

protected:
    ITriangulator() = default;
};

} // namespace triangulation
