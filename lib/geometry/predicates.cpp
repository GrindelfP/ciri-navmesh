#include "predicates.hpp"

#include <cmath>

namespace geometry {
    namespace predicates {
        // ─────────────────────────────────────────────────────────────────────────────
        //  orientation
        // ─────────────────────────────────────────────────────────────────────────────

        Orientation orientation(const Point2D &a,
                                const Point2D &b,
                                const Point2D &c,
                                double eps) noexcept {
            // det = (b-a) × (c-a)
            const double det = (b - a).cross(c - a);

            const double scale = a.distSq(b) + b.distSq(c) + c.distSq(a);
            const double threshold = eps * std::sqrt(scale + kEps);

            if (det > threshold) return Orientation::CCW;
            if (det < -threshold) return Orientation::CW;
            return Orientation::Collinear;
        }

        // ─────────────────────────────────────────────────────────────────────────────
        //  inCircle
        // ─────────────────────────────────────────────────────────────────────────────

        CirclePosition inCircle(const Point2D &a,
                                const Point2D &b,
                                const Point2D &c,
                                const Point2D &d,
                                double eps) noexcept {
            const double ax = a.x - d.x, ay = a.y - d.y;
            const double bx = b.x - d.x, by = b.y - d.y;
            const double cx = c.x - d.x, cy = c.y - d.y;

            // det = | ax  ay  ax²+ay² |
            //       | bx  by  bx²+by² |
            //       | cx  cy  cx²+cy² |
            const double az = ax * ax + ay * ay;
            const double bz = bx * bx + by * by;
            const double cz = cx * cx + cy * cy;

            // det = ax*(by*cz - bz*cy) - ay*(bx*cz - bz*cx) + az*(bx*cy - by*cx)
            const double det = ax * (by * cz - bz * cy)
                               - ay * (bx * cz - bz * cx)
                               + az * (bx * cy - by * cx);

            const double scale = (std::abs(ax) + std::abs(ay) + std::abs(az))
                                 * (std::abs(bx) + std::abs(by) + std::abs(bz))
                                 * (std::abs(cx) + std::abs(cy) + std::abs(cz));
            const double threshold = eps * scale;

            if (det > threshold) return CirclePosition::Inside;
            if (det < -threshold) return CirclePosition::Outside;
            return CirclePosition::OnCircle;
        }

        // ─────────────────────────────────────────────────────────────────────────────
        //  onSegment
        // ─────────────────────────────────────────────────────────────────────────────

        bool onSegment(const Point2D &a,
                       const Point2D &b,
                       const Point2D &p,
                       double eps) noexcept {
            const double minX = std::min(a.x, b.x) - eps;
            const double maxX = std::max(a.x, b.x) + eps;
            const double minY = std::min(a.y, b.y) - eps;
            const double maxY = std::max(a.y, b.y) + eps;

            return p.x >= minX && p.x <= maxX &&
                   p.y >= minY && p.y <= maxY;
        }

        // ─────────────────────────────────────────────────────────────────────────────
        //  segmentIntersect
        // ─────────────────────────────────────────────────────────────────────────────

        IntersectionType segmentIntersect(const Point2D &p1,
                                          const Point2D &p2,
                                          const Point2D &p3,
                                          const Point2D &p4,
                                          double eps) noexcept {
            const Orientation d1 = orientation(p1, p2, p3, eps);
            const Orientation d2 = orientation(p1, p2, p4, eps);
            const Orientation d3 = orientation(p3, p4, p1, eps);
            const Orientation d4 = orientation(p3, p4, p2, eps);

            // Собственное пересечение: концы отрезков строго по разные стороны.
            if (d1 != d2 && d3 != d4 &&
                d1 != Orientation::Collinear && d2 != Orientation::Collinear &&
                d3 != Orientation::Collinear && d4 != Orientation::Collinear) {
                return IntersectionType::Proper;
            }

            // Граничные случаи: коллинеарность.
            //
            // ВАЖНО: Overlap проверяем ПЕРВЫМ — иначе при перекрывающихся отрезках
            // срабатывает Endpoint (одна из конечных точек лежит на другом отрезке),
            // и мы теряем семантику Overlap.

            const bool allCollinear = d1 == Orientation::Collinear &&
                                      d2 == Orientation::Collinear &&
                                      d3 == Orientation::Collinear &&
                                      d4 == Orientation::Collinear;

            if (allCollinear) {
                // Вырожденный отрезок (точка) не может дать Overlap — только Endpoint.
                const bool seg1Degen = p1.nearlyEqual(p2, eps);
                const bool seg2Degen = p3.nearlyEqual(p4, eps);

                if (seg1Degen || seg2Degen) {
                    if (seg1Degen && onSegment(p3, p4, p1, eps)) return IntersectionType::Endpoint;
                    if (seg2Degen && onSegment(p1, p2, p3, eps)) return IntersectionType::Endpoint;
                    return IntersectionType::None;
                }

                // Оба ненулевые — проверяем перекрытие bbox-ами по X и Y.
                const double ax = std::min(p1.x, p2.x), bx = std::max(p1.x, p2.x);
                const double cx = std::min(p3.x, p4.x), dx = std::max(p3.x, p4.x);
                const double ay = std::min(p1.y, p2.y), by = std::max(p1.y, p2.y);
                const double cy = std::min(p3.y, p4.y), dy = std::max(p3.y, p4.y);

                if (ax <= dx + eps && cx <= bx + eps &&
                    ay <= dy + eps && cy <= by + eps) {
                    return IntersectionType::Overlap;
                }
                return IntersectionType::None;
            }

            // Не все коллинеарны — проверяем граничные точки (Endpoint).

            // Коллинеарность p1-p2-p3 и p3 на [p1,p2]
            if (d1 == Orientation::Collinear && onSegment(p1, p2, p3, eps))
                return IntersectionType::Endpoint;

            // Коллинеарность p1-p2-p4 и p4 на [p1,p2]
            if (d2 == Orientation::Collinear && onSegment(p1, p2, p4, eps))
                return IntersectionType::Endpoint;

            // Коллинеарность p3-p4-p1 и p1 на [p3,p4]
            if (d3 == Orientation::Collinear && onSegment(p3, p4, p1, eps))
                return IntersectionType::Endpoint;

            // Коллинеарность p3-p4-p2 и p2 на [p3,p4]
            if (d4 == Orientation::Collinear && onSegment(p3, p4, p2, eps))
                return IntersectionType::Endpoint;

            return IntersectionType::None;
        }

        // ─────────────────────────────────────────────────────────────────────────────
        //  lineIntersectionPoint
        // ─────────────────────────────────────────────────────────────────────────────

        bool lineIntersectionPoint(const Point2D &p1,
                                   const Point2D &p2,
                                   const Point2D &p3,
                                   const Point2D &p4,
                                   Point2D &out,
                                   double eps) noexcept {
            // Параметрическая форма: P = p1 + t*(p2-p1), Q = p3 + u*(p4-p3).
            // Из условия P = Q:
            //   (p2-p1) × (p4-p3) * t = (p3-p1) × (p4-p3)
            //
            // denom = (p2-p1) × (p4-p3)
            const Point2D r = p2 - p1;
            const Point2D s = p4 - p3;

            const double denom = r.cross(s);

            if (std::abs(denom) < eps) {
                // Параллельные или совпадающие прямые.
                return false;
            }

            const double t = (p3 - p1).cross(s) / denom;

            out = p1 + r * t;
            return true;
        }

        // ─────────────────────────────────────────────────────────────────────────────
        //  signedAngle
        // ─────────────────────────────────────────────────────────────────────────────

        double signedAngle(const Point2D &a,
                           const Point2D &b,
                           const Point2D &c) noexcept {
            const Point2D ab = b - a;
            const Point2D ac = c - a;
            // atan2(cross, dot) — угол со знаком от ab до ac.
            return std::atan2(ab.cross(ac), ab.dot(ac));
        }

        // ─────────────────────────────────────────────────────────────────────────────
        //  isLocallyDelaunay
        // ─────────────────────────────────────────────────────────────────────────────

        bool isLocallyDelaunay(const Point2D &p,
                               const Point2D &q,
                               const Point2D &r,
                               const Point2D &s,
                               double eps) noexcept {
            // Ребро pq разделяет треугольники (p, q, r) и (q, p, s).
            // Ребро локально-Делоне ↔ s не внутри circumcircle(p, q, r).
            //
            // ВАЖНО: (p, q, r) должны быть в CCW-порядке для корректного inCircle.
            // Если нет — нормализуем.
            const Orientation ori = orientation(p, q, r, eps);

            if (ori == Orientation::CCW) {
                return inCircle(p, q, r, s, eps) != CirclePosition::Inside;
            } else if (ori == Orientation::CW) {
                // Меняем порядок чтобы треугольник был CCW: (p, r, q)
                return inCircle(p, r, q, s, eps) != CirclePosition::Inside;
            } else {
                // Вырожденный треугольник — ребро считается локально-Делоне
                // (flip не поможет).
                return true;
            }
        }
    } // namespace predicates
} // namespace geometry
