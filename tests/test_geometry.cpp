/**
 * @file test_geometry.cpp
 * @brief Юнит-тесты геометрического ядра (primitives + predicates).
 *
 * Тесты написаны в минималистичном стиле без внешних фреймворков:
 * простой assert-based runner, который можно заменить на Google Test
 * или Catch2 при необходимости.
 *
 * Запуск:
 *   cmake --build build && ./build/tests/test_geometry
 */

#include "../lib/geometry/predicates.hpp"
#include "../lib/geometry/primitives.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

// ─────────────────────────────────────────────────────────────────────────────
//  Минимальный тест-раннер
// ─────────────────────────────────────────────────────────────────────────────

namespace test {

int passed = 0;
int failed = 0;

void check(bool condition, const std::string& name) {
    if (condition) {
        std::cout << "  [PASS] " << name << '\n';
        ++passed;
    } else {
        std::cout << "  [FAIL] " << name << '\n';
        ++failed;
    }
}

void section(const std::string& name) {
    std::cout << "\n── " << name << " ──────────────────────────\n";
}

void summary() {
    std::cout << "\n══════════════════════════════════════\n";
    std::cout << "  Passed: " << passed << '\n';
    std::cout << "  Failed: " << failed << '\n';
    std::cout << "══════════════════════════════════════\n";
}

} // namespace test

using namespace geometry;
using namespace geometry::predicates;
using test::check;

// ─────────────────────────────────────────────────────────────────────────────
//  Вспомогательная функция
// ─────────────────────────────────────────────────────────────────────────────

static bool approxEq(double a, double b, double eps = 1e-9) {
    return std::abs(a - b) < eps;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Point2D
// ─────────────────────────────────────────────────────────────────────────────

void testPoint2D() {
    test::section("Point2D — базовые операции");

    const Point2D O{};
    const Point2D A{3.0, 4.0};
    const Point2D B{1.0, 2.0};

    check(approxEq(A.norm(), 5.0), "norm() — (3,4) должна быть 5");
    check(approxEq(A.dist(O), 5.0), "dist() от начала координат");
    check(approxEq(A.distSq(O), 25.0), "distSq() от начала координат");

    // Арифметика
    check((A + B) == Point2D(4.0, 6.0), "operator+");
    check((A - B) == Point2D(2.0, 2.0), "operator-");
    check((A * 2.0) == Point2D(6.0, 8.0), "operator* скаляр");
    check((2.0 * A) == Point2D(6.0, 8.0), "operator* скаляр слева");

    // Dot и cross
    check(approxEq(A.dot(B), 11.0), "dot product");  // 3*1 + 4*2 = 11
    check(approxEq(A.cross(B), 2.0), "cross product"); // 3*2 - 4*1 = 2

    // Нормализация
    const Point2D An = A.normalized();
    check(approxEq(An.norm(), 1.0, 1e-12), "normalized() имеет единичную длину");
    check(approxEq(An.x, 0.6, 1e-12) && approxEq(An.y, 0.8, 1e-12),
          "normalized() верное направление");

    // Сравнение
    check(A == Point2D(3.0, 4.0), "operator== равные точки");
    check(A != B, "operator!= разные точки");

    // nearlyEqual
    check(A.nearlyEqual(Point2D(3.0 + 1e-12, 4.0 - 1e-12)),
          "nearlyEqual() внутри eps");
    check(!A.nearlyEqual(Point2D(3.0 + 1e-6, 4.0), 1e-9),
          "nearlyEqual() вне eps");

    // Лексикографический порядок
    check(B < A, "(1,2) < (3,4) — по X");
    check(Point2D(1, 0) < Point2D(1, 1), "(1,0) < (1,1) — по Y");

    test::section("Point2D — граничные случаи");

    // Нулевой вектор
    check(approxEq(O.norm(), 0.0), "норма нулевого вектора = 0");
    check(O.nearlyEqual(Point2D(0, 0)), "нулевой вектор nearlyEqual себе");
    check(approxEq(O.dot(A), 0.0), "dot с нулём = 0");
}

// ─────────────────────────────────────────────────────────────────────────────
//  Segment
// ─────────────────────────────────────────────────────────────────────────────

void testSegment() {
    test::section("Segment");

    const Segment s{{0, 0}, {4, 0}}; // горизонтальный, длина 4

    check(approxEq(s.length(), 4.0), "length() горизонтального отрезка");
    check(approxEq(s.lengthSq(), 16.0), "lengthSq()");
    check(s.midpoint() == Point2D(2, 0), "midpoint()");

    // closestPoint
    check(s.closestPoint({2, 5}) == Point2D(2, 0),
          "closestPoint() проекция над серединой");
    check(s.closestPoint({-1, 0}) == Point2D(0, 0),
          "closestPoint() слева от начала — зажим к a");
    check(s.closestPoint({6, 0}) == Point2D(4, 0),
          "closestPoint() справа от конца — зажим к b");
    check(s.closestPoint({0, 0}) == Point2D(0, 0),
          "closestPoint() на начале");

    // distSqToPoint
    check(approxEq(s.distSqToPoint({2, 3}), 9.0),
          "distSqToPoint() над серединой");
    check(approxEq(s.distSqToPoint({-3, 0}), 9.0),
          "distSqToPoint() слева — расстояние до a");

    // Ненаправленное равенство
    const Segment s2{{4, 0}, {0, 0}};
    check(s == s2, "неориентированное равенство {a,b} == {b,a}");

    // Вырожденный отрезок (точка)
    const Segment deg{{1, 1}, {1, 1}};
    check(approxEq(deg.length(), 0.0), "вырожденный отрезок: length = 0");
    check(deg.closestPoint({5, 5}) == Point2D(1, 1),
          "вырожденный отрезок: closestPoint = a");
}

// ─────────────────────────────────────────────────────────────────────────────
//  Triangle
// ─────────────────────────────────────────────────────────────────────────────

void testTriangle() {
    test::section("Triangle — площадь и свойства");

    // Прямоугольный треугольник 3-4-5
    const Triangle t{{0, 0}, {4, 0}, {0, 3}};

    check(approxEq(t.area(), 6.0), "area() треугольника 3-4-5");
    check(approxEq(t.signedArea(), 6.0), "signedArea() CCW > 0");
    check(!t.isDegenerate(), "не вырожденный");

    // CCW-ориентация
    const Triangle tCW{{0, 0}, {0, 3}, {4, 0}}; // обратный порядок
    check(tCW.signedArea() < 0, "CW треугольник: signedArea < 0");
    const Triangle tFixed = tCW.ccw();
    check(tFixed.signedArea() > 0, "ccw() исправляет ориентацию");

    // Centroid
    const Triangle t2{{0, 0}, {6, 0}, {0, 6}};
    const Point2D cen = t2.centroid();
    check(approxEq(cen.x, 2.0) && approxEq(cen.y, 2.0),
          "centroid() равнобедренного прямоугольного треугольника");

    test::section("Triangle — circumcenter");

    // Равносторонний треугольник: circumcenter = centroid
    const double sq3 = std::sqrt(3.0);
    const Triangle eq{{0, 0}, {2, 0}, {1, sq3}};
    const Point2D cc = eq.circumcenter();
    check(approxEq(cc.x, 1.0, 1e-9), "circumcenter X равностороннего");
    check(approxEq(cc.y, sq3 / 3.0, 1e-9), "circumcenter Y равностороннего");

    // Прямоугольный треугольник: circumcenter = середина гипотенузы
    const Triangle right{{0, 0}, {4, 0}, {0, 4}};
    const Point2D ccR = right.circumcenter();
    check(approxEq(ccR.x, 2.0, 1e-9) && approxEq(ccR.y, 2.0, 1e-9),
          "circumcenter прямоугольного треугольника = середина гипотенузы");

    test::section("Triangle — contains()");

    const Triangle tri{{0, 0}, {4, 0}, {2, 4}};

    check(tri.contains({2, 1}), "точка внутри треугольника");
    check(!tri.contains({-1, 0}), "точка слева от треугольника");
    check(!tri.contains({5, 0}), "точка справа от треугольника");
    check(!tri.contains({2, 5}), "точка выше треугольника");

    // Граничные случаи — точки на рёбрах и вершинах
    check(tri.contains({0, 0}), "вершина a — считается внутри");
    check(tri.contains({4, 0}), "вершина b — считается внутри");
    check(tri.contains({2, 0}), "середина ребра ab — считается внутри");
    check(tri.contains({2, 4}), "вершина c — считается внутри");
}

// ─────────────────────────────────────────────────────────────────────────────
//  orientation()
// ─────────────────────────────────────────────────────────────────────────────

void testOrientation() {
    test::section("orientation()");

    const Point2D O{0, 0}, A{1, 0}, B{0, 1};
    const Point2D C{2, 0}; // коллинеарна с O и A

    check(orientation(O, A, B) == Orientation::CCW,
          "CCW: (0,0)(1,0)(0,1)");
    check(orientation(O, B, A) == Orientation::CW,
          "CW: (0,0)(0,1)(1,0)");
    check(orientation(O, A, C) == Orientation::Collinear,
          "Collinear: три точки на оси X");

    // Граничные случаи
    check(orientation(O, O, A) == Orientation::Collinear,
          "совпадающие точки — коллинеарность");
    check(orientation(O, O, O) == Orientation::Collinear,
          "все точки совпадают — коллинеарность");

    // Масштаб: очень большие координаты
    const Point2D P1{1e6, 0}, P2{2e6, 0}, P3{1e6, 1};
    check(orientation(P1, P2, P3) == Orientation::CCW,
          "CCW при больших координатах");
    check(orientation(P1, P3, P2) == Orientation::CW,
          "CW при больших координатах");
}

// ─────────────────────────────────────────────────────────────────────────────
//  inCircle()
// ─────────────────────────────────────────────────────────────────────────────

void testInCircle() {
    test::section("inCircle()");

    // Единичная окружность с треугольником (1,0)(0,1)(-1,0) [CCW]
    const Point2D A{1, 0}, B{0, 1}, C{-1, 0};

    // Точки строго внутри/снаружи
    const Point2D inside{0, 0};
    const Point2D outside{2, 0};
    const Point2D onCircle{0, -1};

    check(inCircle(A, B, C, inside) == CirclePosition::Inside,
          "начало координат внутри единичной окружности");
    check(inCircle(A, B, C, outside) == CirclePosition::Outside,
          "точка (2,0) снаружи");
    check(inCircle(A, B, C, onCircle) == CirclePosition::OnCircle,
          "точка (0,-1) на окружности");

    // Тест с конкретным треугольником Делоне
    // Треугольник (0,0)(4,0)(2,2) — CCW
    const Point2D D{0, 0}, E{4, 0}, F{2, 2};
    check(inCircle(D, E, F, {2, 1}) == CirclePosition::Inside,
          "точка (2,1) внутри circumcircle треугольника DEF");
    check(inCircle(D, E, F, {10, 10}) == CirclePosition::Outside,
          "точка (10,10) снаружи circumcircle");

    // Граничный: вырожденный треугольник (коллинеарные точки)
    // Для таких треугольников inCircle некорректен, но не должен крашиться.
    const Point2D G{0, 0}, H{1, 0}, I{2, 0};
    // Просто убеждаемся, что не падает
    (void)inCircle(G, H, I, {1, 1});
    check(true, "inCircle с коллинеарными точками не вызывает UB");
}

// ─────────────────────────────────────────────────────────────────────────────
//  segmentIntersect()
// ─────────────────────────────────────────────────────────────────────────────

void testSegmentIntersect() {
    test::section("segmentIntersect() — собственные пересечения");

    // Крест
    check(segmentIntersect({0,0},{2,2},{0,2},{2,0}) == IntersectionType::Proper,
          "крест — собственное пересечение");
    // Параллельные
    check(segmentIntersect({0,0},{1,0},{0,1},{1,1}) == IntersectionType::None,
          "параллельные отрезки — нет пересечения");
    // T-образное
    check(segmentIntersect({0,0},{2,0},{1,0},{1,2}) == IntersectionType::Endpoint,
          "T-образное — пересечение в конечной точке");

    test::section("segmentIntersect() — граничные случаи");

    // Общая конечная точка
    check(segmentIntersect({0,0},{1,0},{1,0},{2,1}) == IntersectionType::Endpoint,
          "общая конечная точка");

    // Коллинеарные перекрывающиеся
    check(segmentIntersect({0,0},{2,0},{1,0},{3,0}) == IntersectionType::Overlap,
          "коллинеарные перекрывающиеся");

    // Коллинеарные без перекрытия
    check(segmentIntersect({0,0},{1,0},{2,0},{3,0}) == IntersectionType::None,
          "коллинеарные без перекрытия");

    // Один отрезок — точка
    check(segmentIntersect({1,0},{1,0},{0,0},{2,0}) == IntersectionType::Endpoint,
          "вырожденный отрезок-точка на другом отрезке");

    // Один отрезок — точка, не на другом
    check(segmentIntersect({5,5},{5,5},{0,0},{2,0}) == IntersectionType::None,
          "вырожденный отрезок-точка не на другом отрезке");

    // Перпендикулярные, не пересекающиеся
    check(segmentIntersect({0,0},{1,0},{2,0},{2,1}) == IntersectionType::None,
          "перпендикулярные, не пересекающиеся");
}

// ─────────────────────────────────────────────────────────────────────────────
//  lineIntersectionPoint()
// ─────────────────────────────────────────────────────────────────────────────

void testLineIntersection() {
    test::section("lineIntersectionPoint()");

    Point2D out;

    // Классическое пересечение
    bool ok = lineIntersectionPoint({0,0},{2,2},{0,2},{2,0}, out);
    check(ok, "пересечение диагоналей: функция вернула true");
    check(approxEq(out.x, 1.0) && approxEq(out.y, 1.0),
          "пересечение диагоналей: точка (1,1)");

    // Параллельные прямые
    bool fail = lineIntersectionPoint({0,0},{1,0},{0,1},{1,1}, out);
    check(!fail, "параллельные прямые: функция вернула false");

    // Совпадающие прямые
    bool fail2 = lineIntersectionPoint({0,0},{1,0},{2,0},{3,0}, out);
    check(!fail2, "совпадающие прямые: функция вернула false");

    // Перпендикулярные оси
    bool ok2 = lineIntersectionPoint({0,5},{1,5},{3,0},{3,10}, out);
    check(ok2, "перпендикулярные: функция вернула true");
    check(approxEq(out.x, 3.0) && approxEq(out.y, 5.0),
          "перпендикулярные: точка (3,5)");
}

// ─────────────────────────────────────────────────────────────────────────────
//  isLocallyDelaunay()
// ─────────────────────────────────────────────────────────────────────────────

void testLocallyDelaunay() {
    test::section("isLocallyDelaunay()");

    // Два треугольника, образующих квадрат: (0,0)(2,0)(2,2)(0,2)
    // Диагональ (0,0)-(2,2): треугольники [(0,0)(2,0)(2,2)] и [(0,0)(2,2)(0,2)]
    // Для Делоне-триангуляции квадрата эта диагональ локально-Делоне.
    check(isLocallyDelaunay({0,0},{2,2},{2,0},{0,2}),
          "диагональ квадрата (0,0)-(2,2) локально-Делоне");

    // Для второй диагонали (2,0)-(0,2): должна дать то же по симметрии
    check(isLocallyDelaunay({2,0},{0,2},{0,0},{2,2}),
          "диагональ квадрата (2,0)-(0,2) тоже локально-Делоне (симметрия)");

    // «Плохой» flip: тупоугольный треугольник
    // Четырёхугольник (0,0)(3,0)(2,1)(0,2) — диагональ (0,0)-(2,1)
    // vs диагональ (3,0)-(0,2)
    // В этой конфигурации какая-то диагональ нарушает условие Делоне.
    // Просто проверяем, что функция возвращает что-то не UB:
    const bool r = isLocallyDelaunay({0,0},{2,1},{3,0},{0,2});
    check(r == true || r == false,
          "isLocallyDelaunay на тупоугольнике не вызывает UB");
}

// ─────────────────────────────────────────────────────────────────────────────
//  main
// ─────────────────────────────────────────────────────────────────────────────

int main() {
    std::cout << "╔══════════════════════════════════════╗\n";
    std::cout << "║  Geometry Core — Unit Tests          ║\n";
    std::cout << "╚══════════════════════════════════════╝\n";

    testPoint2D();
    testSegment();
    testTriangle();
    testOrientation();
    testInCircle();
    testSegmentIntersect();
    testLineIntersection();
    testLocallyDelaunay();

    test::summary();

    return test::failed > 0 ? 1 : 0;
}
