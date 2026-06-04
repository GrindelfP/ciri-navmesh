#pragma once

/**
 * @file predicates.hpp
 * @brief Геометрические предикаты для алгоритмов триангуляции и навмешей.
 *
 * Все предикаты используют double-арифметику с epsilon-защитой.
 * Для промышленного применения стоит заменить на точные предикаты Шевчука,
 * но для учебного проекта double + аккуратный epsilon вполне корректны.
 *
 * Соглашение об ориентации: вершины треугольника в порядке
 * counter-clockwise (CCW) считаются «правильными».
 */

#include "primitives.hpp"

namespace geometry::predicates {

// ─────────────────────────────────────────────────────────────────────────────
//  Orientation
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Результат теста ориентации трёх точек.
 */
enum class Orientation {
    CCW,        ///< Counter-clockwise (левый поворот)
    CW,         ///< Clockwise (правый поворот)
    Collinear   ///< Точки коллинеарны
};

/**
 * @brief Определяет ориентацию тройки точек (a, b, c).
 *
 * Вычисляет знак знаковой площади треугольника (a, b, c):
 * @code
 *   det = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
 * @endcode
 *
 * @param a  Первая точка
 * @param b  Вторая точка
 * @param c  Третья точка
 * @param eps Порог для различения коллинеарности от «почти коллинеарности».
 *            По умолчанию kEps (1e-9). Для задач с шумом можно увеличить.
 * @return Orientation::CCW  если det > eps
 *         Orientation::CW   если det < -eps
 *         Orientation::Collinear если |det| <= eps
 *
 * @note Граничный случай: при eps = 0 предикат строг, но нестабилен
 *       на вырожденных конфигурациях. Не устанавливайте eps = 0 в продакшне.
 */
[[nodiscard]] Orientation orientation(const Point2D& a,
                                      const Point2D& b,
                                      const Point2D& c,
                                      double eps = kEps) noexcept;

/**
 * @brief Возвращает true, если точки (a, b, c) образуют CCW-поворот.
 * @see orientation()
 */
[[nodiscard]] inline bool isCCW(const Point2D& a,
                                 const Point2D& b,
                                 const Point2D& c,
                                 double eps = kEps) noexcept {
    return orientation(a, b, c, eps) == Orientation::CCW;
}

/**
 * @brief Возвращает true, если точки (a, b, c) коллинеарны.
 * @see orientation()
 */
[[nodiscard]] inline bool isCollinear(const Point2D& a,
                                      const Point2D& b,
                                      const Point2D& c,
                                      double eps = kEps) noexcept {
    return orientation(a, b, c, eps) == Orientation::Collinear;
}

// ─────────────────────────────────────────────────────────────────────────────
//  InCircle
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Результат теста inCircle.
 */
enum class CirclePosition {
    Inside,    ///< Точка строго внутри описанной окружности
    Outside,   ///< Точка строго снаружи
    OnCircle   ///< Точка лежит на окружности (в пределах eps)
};

/**
 * @brief Проверяет, лежит ли точка @p d внутри описанной окружности
 *        треугольника (a, b, c).
 *
 * Предполагает, что (a, b, c) — CCW. Если треугольник CW,
 * Inside и Outside меняются местами.
 *
 * Вычисляется через детерминант 3×3:
 * @code
 *   | ax-dx  ay-dy  (ax-dx)²+(ay-dy)² |
 *   | bx-dx  by-dy  (bx-dx)²+(by-dy)² |
 *   | cx-dx  cy-dy  (cx-dx)²+(cy-dy)² |
 * @endcode
 * det > 0 → d внутри (для CCW-треугольника).
 *
 * @param a  Первая вершина треугольника (CCW-порядок)
 * @param b  Вторая вершина
 * @param c  Третья вершина
 * @param d  Тестируемая точка
 * @param eps Порог для OnCircle. Масштабируется с размером треугольника.
 * @return CirclePosition
 *
 * @warning Результат некорректен для вырожденных треугольников
 *          (det знаменателя ≈ 0). Проверяйте isDegenerate() заранее.
 */
[[nodiscard]] CirclePosition inCircle(const Point2D& a,
                                      const Point2D& b,
                                      const Point2D& c,
                                      const Point2D& d,
                                      double eps = kEps) noexcept;

/**
 * @brief Строгий тест: возвращает true, если @p d строго внутри описанной
 *        окружности треугольника (a, b, c) в CCW-порядке.
 * @see inCircle()
 */
[[nodiscard]] inline bool isInCircle(const Point2D& a,
                                     const Point2D& b,
                                     const Point2D& c,
                                     const Point2D& d,
                                     double eps = kEps) noexcept {
    return inCircle(a, b, c, d, eps) == CirclePosition::Inside;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Intersection
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Результат теста пересечения двух отрезков.
 */
enum class IntersectionType {
    None,       ///< Отрезки не пересекаются
    Proper,     ///< Собственное пересечение (во внутренних точках обоих)
    Endpoint,   ///< Пересечение в конечной точке одного или обоих отрезков
    Overlap     ///< Отрезки коллинеарны и перекрываются
};

/**
 * @brief Проверяет пересечение двух отрезков.
 *
 * Отрезки: p1–p2 и p3–p4.
 *
 * Алгоритм:
 * 1. Проверяем знаки ориентации: если концы одного отрезка по разные
 *    стороны от другого — собственное пересечение.
 * 2. Граничные случаи (конечные точки, коллинеарность) обрабатываются
 *    отдельно.
 *
 * @param p1 Начало первого отрезка
 * @param p2 Конец первого отрезка
 * @param p3 Начало второго отрезка
 * @param p4 Конец второго отрезка
 * @param eps Порог для различения Endpoint/Overlap от None.
 * @return IntersectionType
 *
 * @note Граничный случай: два совпадающих отрезка вернут Overlap.
 *       Нулевой отрезок (p1 == p2) обрабатывается корректно.
 */
[[nodiscard]] IntersectionType segmentIntersect(const Point2D& p1,
                                                 const Point2D& p2,
                                                 const Point2D& p3,
                                                 const Point2D& p4,
                                                 double eps = kEps) noexcept;

/**
 * @brief Перегрузка через структуры Segment.
 * @see segmentIntersect(Point2D, Point2D, Point2D, Point2D)
 */
[[nodiscard]] inline IntersectionType segmentIntersect(const Segment& s1,
                                                        const Segment& s2,
                                                        double eps = kEps) noexcept {
    return segmentIntersect(s1.a, s1.b, s2.a, s2.b, eps);
}

/**
 * @brief Быстрый тест: возвращает true при любом типе пересечения (кроме None).
 * @see segmentIntersect()
 */
[[nodiscard]] inline bool intersects(const Point2D& p1,
                                     const Point2D& p2,
                                     const Point2D& p3,
                                     const Point2D& p4,
                                     double eps = kEps) noexcept {
    return segmentIntersect(p1, p2, p3, p4, eps) != IntersectionType::None;
}

/**
 * @brief Вычисляет точку собственного пересечения двух прямых (не отрезков).
 *
 * Прямые заданы парами точек: (p1, p2) и (p3, p4).
 * Используется в Funnel Algorithm для вычисления точных пересечений.
 *
 * @param[out] out Точка пересечения (заполняется только при возврате true).
 * @return true если прямые пересекаются (не параллельны).
 * @return false если прямые параллельны или совпадают.
 */
[[nodiscard]] bool lineIntersectionPoint(const Point2D& p1,
                                          const Point2D& p2,
                                          const Point2D& p3,
                                          const Point2D& p4,
                                          Point2D& out,
                                          double eps = kEps) noexcept;

// ─────────────────────────────────────────────────────────────────────────────
//  Дополнительные предикаты
// ─────────────────────────────────────────────────────────────────────────────
/**
 * @brief Проверяет, лежит ли точка @p p на отрезке @p a – @p b.
 *
 * Предполагает, что p коллинеарна с a и b (т.е. ориентация уже проверена).
 * Проверяет только принадлежность по координатам (bbox-тест).
 *
 * @param a Начало отрезка
 * @param b Конец отрезка
 * @param p Тестируемая точка
 */
[[nodiscard]] bool onSegment(const Point2D& a,
                              const Point2D& b,
                              const Point2D& p,
                              double eps = kEps) noexcept;

/**
 * @brief Угол в радианах от вектора (a→b) до вектора (a→c), со знаком.
 *
 * Возвращает значение из диапазона (-π, π].
 * > 0 → CCW-поворот, < 0 → CW-поворот.
 *
 * @param a Вершина угла
 * @param b Первое направление
 * @param c Второе направление
 */
[[nodiscard]] double signedAngle(const Point2D& a,
                                  const Point2D& b,
                                  const Point2D& c) noexcept;

/**
 * @brief Проверяет, является ли рёбро (p, q) локально-Делоне относительно
 *        двух треугольников, которые оно разделяет: (p, q, r) и (q, p, s).
 *
 * Ребро локально-Делоне ↔ точки r и s не лежат внутри
 * описанной окружности друг друга относительно этого ребра.
 *
 * Эквивалентно: s не внутри circumcircle(p, q, r) при CCW-ориентации.
 *
 * @param p  Первая вершина общего ребра
 * @param q  Вторая вершина общего ребра
 * @param r  Третья вершина первого треугольника (напротив ребра pq)
 * @param s  Третья вершина второго треугольника (напротив ребра pq)
 * @param eps Порог для теста inCircle.
 * @return true если ребро локально-Делоне (flip не нужен).
 */
[[nodiscard]] bool isLocallyDelaunay(const Point2D& p,
                                      const Point2D& q,
                                      const Point2D& r,
                                      const Point2D& s,
                                      double eps = kEps) noexcept;

} // namespace geometry::predicates

