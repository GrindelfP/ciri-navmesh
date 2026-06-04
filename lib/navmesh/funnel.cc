/**
 * @file funnel.cc
 * @brief Funnel Algorithm — shortest Euclidean path through a triangle corridor.
 *
 * ## Portal convention (from navmesh.cc)
 *
 *   arc.portal[0] → the vertex that is CCW-first in the source face half-edge,
 *                   i.e. the LEFT  vertex when walking src→dst.
 *   arc.portal[1] → RIGHT vertex.
 *
 * portalPoints(from, to, dcel) returns {left, right} in this order.
 *
 * ## Funnel state
 *
 * We maintain two deques that together describe the "funnel" — the region of
 * the plane visible from the current apex:
 *
 *   leftDeque  — chain of vertices on the LEFT  side of the funnel.
 *                leftDeque.front() == apex.
 *   rightDeque — chain of vertices on the RIGHT side of the funnel.
 *                rightDeque.front() == apex.
 *
 * When a new portal (L, R) is processed:
 *
 *   - New L vertex: tighten or extend the left  chain; if it crosses the right
 *     chain, the right chain's apex becomes the new apex (emit waypoints).
 *   - New R vertex: symmetric.
 *
 * ## Cross product sign convention
 *
 *   cross2D(apex, a, b) = (a - apex) × (b - apex)
 *   > 0 → b is CCW (LEFT)  of ray apex→a
 *   < 0 → b is CW  (RIGHT) of ray apex→a
 *   = 0 → collinear
 */

#include "funnel.hh"

#include <cassert>
#include <cmath>
#include <deque>

namespace pathfinding {
    using geometry::Point2D;
    using geometry::DCEL;
    using navmesh::NavMesh;
    using navmesh::NodeIdx;
    using navmesh::kInvalidNode;

    // ─────────────────────────────────────────────────────────────────────────────
    //  cross2D  (static helper, defined before use)
    // ─────────────────────────────────────────────────────────────────────────────
    double FunnelAlgorithm::cross2D(const Point2D &apex,
                                    const Point2D &a,
                                    const Point2D &b) noexcept {
        // (a - apex) × (b - apex)
        const double ax = a.x - apex.x, ay = a.y - apex.y;
        const double bx = b.x - apex.x, by = b.y - apex.y;
        return ax * by - ay * bx;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  buildPortals
    // ─────────────────────────────────────────────────────────────────────────────
    std::vector<Portal>
    FunnelAlgorithm::buildPortals(const std::vector<NodeIdx> &corridor,
                                  const NavMesh &nm,
                                  const DCEL &dcel) {
        std::vector<Portal> portals;
        portals.reserve(corridor.size() > 0 ? corridor.size() - 1 : 0);

        for (std::size_t i = 0; i + 1 < corridor.size(); ++i) {
            auto pts = nm.portalPoints(corridor[i], corridor[i + 1], dcel);
            assert(pts.has_value() && "buildPortals: non-adjacent nodes in corridor");

            // portalPoints returns {left, right} per navmesh.cc convention.
            portals.push_back({(*pts)[0], (*pts)[1]});
        }
        return portals;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  computeMetrics
    // ─────────────────────────────────────────────────────────────────────────────
    void FunnelAlgorithm::computeMetrics(const std::vector<Point2D> &waypoints,
                                         double turnThreshold,
                                         double &length,
                                         std::size_t &turns) noexcept {
        length = 0.0;
        turns = 0;

        const std::size_t n = waypoints.size();
        if (n < 2) return;

        for (std::size_t i = 0; i + 1 < n; ++i)
            length += waypoints[i].dist(waypoints[i + 1]);

        // Count interior waypoints whose turn angle exceeds the threshold.
        // The turn angle at waypoints[i] is the angle between vectors
        // (waypoints[i-1]→waypoints[i]) and (waypoints[i]→waypoints[i+1]).
        for (std::size_t i = 1; i + 1 < n; ++i) {
            const Point2D v1 = waypoints[i] - waypoints[i - 1];
            const Point2D v2 = waypoints[i + 1] - waypoints[i];
            const double d = v1.norm() * v2.norm();
            if (d < 1e-12) continue; // degenerate zero-length segment

            // cos θ = (v1 · v2) / (|v1| |v2|)
            const double cosA = std::clamp(v1.dot(v2) / d, -1.0, 1.0);
            // θ is the turning angle (π - interior angle).
            const double angle = std::acos(cosA);
            if (angle > turnThreshold) ++turns;
        }
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  smooth  (core: takes pre-built portal list)
    // ─────────────────────────────────────────────────────────────────────────────
    SmoothPath FunnelAlgorithm::smooth(Point2D start,
                                       Point2D goal,
                                       const std::vector<Portal> &portals,
                                       double turnThreshold) const {
        SmoothPath result;
        std::vector<Point2D> &waypoints = result.waypoints;
        waypoints.push_back(start);

        // Trivial: no portals (start and goal in same triangle).
        if (portals.empty()) {
            waypoints.push_back(goal);
            computeMetrics(waypoints, turnThreshold, result.totalLength, result.turnCount);
            return result;
        }

        // ── Funnel deques ─────────────────────────────────────────────────────────
        // Both deques share the apex at their front.
        // leftDeque[0] == rightDeque[0] == current apex.

        std::deque<Point2D> leftDeq = {start};
        std::deque<Point2D> rightDeq = {start};

        // ── Process each portal ───────────────────────────────────────────────────

        for (const Portal &portal: portals) {
            const Point2D &newLeft = portal.left;
            const Point2D &newRight = portal.right;

            // ── Update LEFT side ──────────────────────────────────────────────────

            // cross2D(apex, leftDeq.back(), newLeft):
            //   > 0 → newLeft is further CCW  → tighten (stay on left side)
            //   < 0 → newLeft crosses to right side → new apex event

            if (cross2D(leftDeq.front(), leftDeq.back(), newLeft) >= 0.0) {
                // newLeft is within or on the left boundary — update left chain.
                // Pop vertices that are no longer on the convex hull.
                while (leftDeq.size() > 1 &&
                       cross2D(leftDeq.front(),
                               leftDeq[leftDeq.size() - 2],
                               newLeft) >= 0.0) {
                    leftDeq.pop_back();
                }
                leftDeq.push_back(newLeft);
            } else {
                // newLeft crosses to the right of the right chain:
                // emit waypoints by walking the right chain until we find the
                // new apex.
                while (rightDeq.size() > 1 &&
                       cross2D(rightDeq.front(),
                               rightDeq[rightDeq.size() - 2],
                               newLeft) <= 0.0) {
                    // Emit the next right-chain vertex as a waypoint.
                    waypoints.push_back(rightDeq.back());
                    rightDeq.pop_back();
                }
                // The new apex is the remaining tip of the right chain.
                const Point2D newApex = rightDeq.back();

                // Reset both chains from the new apex.
                leftDeq.clear();
                rightDeq.clear();
                leftDeq.push_back(newApex);
                rightDeq.push_back(newApex);
                leftDeq.push_back(newLeft);
            }

            // ── Update RIGHT side ─────────────────────────────────────────────────

            if (cross2D(rightDeq.front(), rightDeq.back(), newRight) <= 0.0) {
                // newRight is within or on the right boundary.
                while (rightDeq.size() > 1 &&
                       cross2D(rightDeq.front(),
                               rightDeq[rightDeq.size() - 2],
                               newRight) <= 0.0) {
                    rightDeq.pop_back();
                }
                rightDeq.push_back(newRight);
            } else {
                // newRight crosses to the left of the left chain.
                while (leftDeq.size() > 1 &&
                       cross2D(leftDeq.front(),
                               leftDeq[leftDeq.size() - 2],
                               newRight) >= 0.0) {
                    waypoints.push_back(leftDeq.back());
                    leftDeq.pop_back();
                }
                const Point2D newApex = leftDeq.back();

                leftDeq.clear();
                rightDeq.clear();
                leftDeq.push_back(newApex);
                rightDeq.push_back(newApex);
                rightDeq.push_back(newRight);
            }
        }

        // ── Close funnel to goal ──────────────────────────────────────────────────
        //
        // Treat goal as a degenerate zero-width portal and apply the same logic.
        // We walk whichever side of the funnel the goal is visible from.

        const Point2D &apex = leftDeq.front(); // == rightDeq.front()

        if (cross2D(apex, leftDeq.back(), goal) >= 0.0 &&
            cross2D(apex, rightDeq.back(), goal) <= 0.0) {
            // Goal is directly visible from apex — straight line.
            waypoints.push_back(goal);
        } else {
            // Determine which side goal falls on and walk that chain.
            if (cross2D(apex, leftDeq.back(), goal) < 0.0) {
                // Goal is to the right of the left chain — walk right chain.
                while (rightDeq.size() > 1 &&
                       cross2D(rightDeq.front(),
                               rightDeq[rightDeq.size() - 2],
                               goal) <= 0.0) {
                    waypoints.push_back(rightDeq.back());
                    rightDeq.pop_back();
                }
            } else {
                // Goal is to the left of the right chain — walk left chain.
                while (leftDeq.size() > 1 &&
                       cross2D(leftDeq.front(),
                               leftDeq[leftDeq.size() - 2],
                               goal) >= 0.0) {
                    waypoints.push_back(leftDeq.back());
                    leftDeq.pop_back();
                }
            }
            waypoints.push_back(goal);
        }

        // ── Deduplicate consecutive identical waypoints ───────────────────────────
        // Can occur when start/goal coincides with a portal vertex.
        waypoints.erase(
            std::unique(waypoints.begin(), waypoints.end(),
                        [](const Point2D &a, const Point2D &b) {
                            return a.distSq(b) < 1e-18;
                        }),
            waypoints.end());

        computeMetrics(waypoints, turnThreshold, result.totalLength, result.turnCount);
        return result;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  smooth  (convenience overload: builds portals from corridor)
    // ─────────────────────────────────────────────────────────────────────────────
    SmoothPath FunnelAlgorithm::smooth(Point2D start,
                                       Point2D goal,
                                       const std::vector<NodeIdx> &corridor,
                                       const NavMesh &nm,
                                       const DCEL &dcel,
                                       double turnThreshold) const {
        if (corridor.empty()) {
            SmoothPath r;
            r.waypoints = {start, goal};
            computeMetrics(r.waypoints, turnThreshold, r.totalLength, r.turnCount);
            return r;
        }

        const std::vector<Portal> portals = buildPortals(corridor, nm, dcel);
        return smooth(start, goal, portals, turnThreshold);
    }
} // namespace pathfinding
