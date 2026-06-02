#include "primitives.hpp"

#include <cassert>
#include <cmath>

namespace geometry {

// ─────────────────────────────────────────────────────────────────────────────
//  Triangle — реализация не-inlined методов
// ─────────────────────────────────────────────────────────────────────────────

Point2D Triangle::circumcenter() const noexcept {
    // Используем формулу через векторы от вершины a:
    //   D  = 2 * [(b-a) × (c-a)]
    //   ux = [(b-a)·(b-a)*(c-a).y - (c-a)·(c-a)*(b-a).y] / D
    //   uy = [(c-a)·(c-a)*(b-a).x - (b-a)·(b-a)*(c-a).x] / D
    //
    // Вывод: пересечение перпендикулярных биссектрис AB и AC.
    const Point2D ab = b - a;
    const Point2D ac = c - a;

    const double D = 2.0 * ab.cross(ac);
    if (std::abs(D) < kEps) {
        // Вырожденный треугольник: возвращаем бесконечно далёкую точку.
        // Вызывающий код должен проверять isDegenerate() заранее.
        return {std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity()};
    }

    const double ab2 = ab.dot(ab); // |AB|²
    const double ac2 = ac.dot(ac); // |AC|²

    const double ux = (ab2 * ac.y - ac2 * ab.y) / D;
    const double uy = (ac2 * ab.x - ab2 * ac.x) / D;

    return a + Point2D{ux, uy};
}

double Triangle::circumradius() const noexcept {
    const double A = area();
    if (A < kEps) return std::numeric_limits<double>::infinity();

    // R = (|AB| * |BC| * |CA|) / (4 * Area)
    const double ab = a.dist(b);
    const double bc = b.dist(c);
    const double ca = c.dist(a);
    return (ab * bc * ca) / (4.0 * A);
}

Point2D Triangle::incenter() const noexcept {
    // Весa — длины противоположных сторон
    const double wa = b.dist(c); // сторона напротив a
    const double wb = a.dist(c); // сторона напротив b
    const double wc = a.dist(b); // сторона напротив c
    const double total = wa + wb + wc;

    if (total < kEps) return a; // вырожденный случай

    return (a * wa + b * wb + c * wc) / total;
}

bool Triangle::contains(const Point2D& p, double eps) const noexcept {
    // Метод знаков площадей подтреугольников:
    //   p внутри (a,b,c) ↔ ориентации (a,b,p), (b,c,p), (c,a,p) одного знака.
    //
    // Граничный случай: eps-допуск позволяет считать точки на рёбрах внутренними.
    // Используем знаковые площади (без деления на 2 — только знак важен).

    const auto sign = [](double v, double e) -> int {
        if (v > e)  return  1;
        if (v < -e) return -1;
        return 0;
    };

    const double d1 = (p - a).cross(b - a);
    const double d2 = (p - b).cross(c - b);
    const double d3 = (p - c).cross(a - c);

    // Масштабируем eps относительно площади треугольника,
    // чтобы допуск не зависел от размера.
    const double areaScale = std::abs((b - a).cross(c - a));
    const double scaledEps = eps * (areaScale + kEps);

    const int s1 = sign(d1, scaledEps);
    const int s2 = sign(d2, scaledEps);
    const int s3 = sign(d3, scaledEps);

    // Если хотя бы один знак = 0 (точка на ребре) — считаем «внутри».
    // Если все знаки одного знака — точка внутри.
    const bool hasNeg = (s1 < 0) || (s2 < 0) || (s3 < 0);
    const bool hasPos = (s1 > 0) || (s2 > 0) || (s3 > 0);

    return !(hasNeg && hasPos);
}

} // namespace geometry
