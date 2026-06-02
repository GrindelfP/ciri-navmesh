#pragma once

#include <cmath>
#include <functional>
#include <ostream>

namespace geometry {

// ─────────────────────────────────────────────────────────────────────────────
//  Epsilon для сравнений с плавающей точкой
// ─────────────────────────────────────────────────────────────────────────────

/// @brief Машинный эпсилон, масштабированный для геометрических предикатов.
///        Все координаты предполагаются в разумных пределах (< 1e6).
inline constexpr double kEps = 1e-9;

// ─────────────────────────────────────────────────────────────────────────────
//  Point2D
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Точка на плоскости с координатами двойной точности.
 *
 * Все геометрические примитивы работают с этим типом.
 * Намеренно простая структура — только данные и базовая арифметика.
 */
struct Point2D {
    double x{0.0}; ///< Координата X
    double y{0.0}; ///< Координата Y

    /// @brief Конструктор по умолчанию — начало координат.
    constexpr Point2D() noexcept = default;

    /// @brief Конструктор из координат.
    constexpr Point2D(double x_, double y_) noexcept : x(x_), y(y_) {}

    // ── Арифметика ──────────────────────────────────────────────────────────

    [[nodiscard]] constexpr Point2D operator+(const Point2D& o) const noexcept {
        return {x + o.x, y + o.y};
    }
    [[nodiscard]] constexpr Point2D operator-(const Point2D& o) const noexcept {
        return {x - o.x, y - o.y};
    }
    [[nodiscard]] constexpr Point2D operator*(double s) const noexcept {
        return {x * s, y * s};
    }
    [[nodiscard]] constexpr Point2D operator/(double s) const noexcept {
        return {x / s, y / s};
    }

    // ── Скалярные операции ───────────────────────────────────────────────────

    /**
     * @brief Скалярное произведение двух векторов (от начала координат).
     * @param o Второй вектор.
     * @return x*o.x + y*o.y
     */
    [[nodiscard]] constexpr double dot(const Point2D& o) const noexcept {
        return x * o.x + y * o.y;
    }

    /**
     * @brief Z-компонента векторного произведения (2D cross product).
     * @param o Второй вектор.
     * @return x*o.y - y*o.x
     *
     * Знак определяет ориентацию: > 0 — CCW, < 0 — CW, = 0 — коллинеарны.
     */
    [[nodiscard]] constexpr double cross(const Point2D& o) const noexcept {
        return x * o.y - y * o.x;
    }

    /**
     * @brief Квадрат евклидова расстояния до точки @p o.
     * @details Предпочтительнее dist() там, где нужно только сравнение —
     *          избегает вызова sqrt.
     */
    [[nodiscard]] constexpr double distSq(const Point2D& o) const noexcept {
        const double dx = x - o.x;
        const double dy = y - o.y;
        return dx * dx + dy * dy;
    }

    /**
     * @brief Евклидово расстояние до точки @p o.
     */
    [[nodiscard]] double dist(const Point2D& o) const noexcept {
        return std::sqrt(distSq(o));
    }

    /**
     * @brief Длина вектора (расстояние от начала координат).
     */
    [[nodiscard]] double norm() const noexcept {
        return std::sqrt(x * x + y * y);
    }

    /**
     * @brief Нормализованный вектор (единичный).
     * @warning Не вызывать для нулевого вектора — результат неопределён.
     */
    [[nodiscard]] Point2D normalized() const noexcept {
        return *this / norm();
    }

    // ── Сравнение ────────────────────────────────────────────────────────────

    /**
     * @brief Точное равенство (побитовое).
     * @note Для нечёткого сравнения используйте nearlyEqual().
     */
    [[nodiscard]] constexpr bool operator==(const Point2D& o) const noexcept {
        return x == o.x && y == o.y;
    }
    [[nodiscard]] constexpr bool operator!=(const Point2D& o) const noexcept {
        return !(*this == o);
    }

    /**
     * @brief Лексикографический порядок (x, затем y) — для использования
     *        в std::set / std::map и алгоритмах сортировки.
     */
    [[nodiscard]] constexpr bool operator<(const Point2D& o) const noexcept {
        if (x != o.x) return x < o.x;
        return y < o.y;
    }

    /**
     * @brief Нечёткое сравнение с допуском @p eps.
     */
    [[nodiscard]] bool nearlyEqual(const Point2D& o,
                                   double eps = kEps) const noexcept {
        return distSq(o) <= eps * eps;
    }

    // ── Вывод ────────────────────────────────────────────────────────────────

    friend std::ostream& operator<<(std::ostream& os, const Point2D& p) {
        return os << '(' << p.x << ", " << p.y << ')';
    }
};

/// @brief Умножение скаляра слева: s * p.
[[nodiscard]] constexpr inline Point2D operator*(double s,
                                                  const Point2D& p) noexcept {
    return p * s;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Segment
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Отрезок, заданный двумя конечными точками.
 *
 * Никакого порядка точек не предполагается (отрезок ненаправленный),
 * если явно не оговорено иное.
 */
struct Segment {
    Point2D a; ///< Первая конечная точка
    Point2D b; ///< Вторая конечная точка

    constexpr Segment() noexcept = default;
    constexpr Segment(Point2D a_, Point2D b_) noexcept : a(a_), b(b_) {}

    /**
     * @brief Длина отрезка.
     */
    [[nodiscard]] double length() const noexcept { return a.dist(b); }

    /**
     * @brief Квадрат длины (дешевле length() — без sqrt).
     */
    [[nodiscard]] constexpr double lengthSq() const noexcept {
        return a.distSq(b);
    }

    /**
     * @brief Середина отрезка.
     */
    [[nodiscard]] constexpr Point2D midpoint() const noexcept {
        return (a + b) * 0.5;
    }

    /**
     * @brief Вектор направления (b - a), ненормализованный.
     */
    [[nodiscard]] constexpr Point2D direction() const noexcept { return b - a; }

    /**
     * @brief Точка на отрезке при параметре t ∈ [0, 1].
     * @param t  0 → a, 1 → b.
     */
    [[nodiscard]] constexpr Point2D at(double t) const noexcept {
        return a + (b - a) * t;
    }

    /**
     * @brief Ближайшая точка на отрезке к заданной точке @p p.
     * @details Проецируем p на прямую ab и зажимаем t в [0, 1].
     */
    [[nodiscard]] Point2D closestPoint(const Point2D& p) const noexcept {
        const Point2D ab = b - a;
        const double lenSq = ab.dot(ab);
        if (lenSq < kEps * kEps) return a; // вырожденный отрезок
        const double t = std::clamp((p - a).dot(ab) / lenSq, 0.0, 1.0);
        return at(t);
    }

    /**
     * @brief Квадрат расстояния от точки @p p до отрезка.
     */
    [[nodiscard]] double distSqToPoint(const Point2D& p) const noexcept {
        return p.distSq(closestPoint(p));
    }

    // ── Сравнение (ненаправленное) ────────────────────────────────────────

    /**
     * @brief Равенство как неориентированных отрезков: {a,b} == {b,a}.
     */
    [[nodiscard]] bool operator==(const Segment& o) const noexcept {
        return (a == o.a && b == o.b) || (a == o.b && b == o.a);
    }
    [[nodiscard]] bool operator!=(const Segment& o) const noexcept {
        return !(*this == o);
    }

    friend std::ostream& operator<<(std::ostream& os, const Segment& s) {
        return os << '[' << s.a << " — " << s.b << ']';
    }
};

// ─────────────────────────────────────────────────────────────────────────────
//  Triangle
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Треугольник, заданный тремя вершинами в порядке CCW.
 *
 * Вершины хранятся в порядке, переданном при создании. Ориентацию
 * (CCW/CW) проверяйте через predicates::orientation() или signedArea().
 */
struct Triangle {
    Point2D a; ///< Первая вершина
    Point2D b; ///< Вторая вершина
    Point2D c; ///< Третья вершина

    constexpr Triangle() noexcept = default;
    constexpr Triangle(Point2D a_, Point2D b_, Point2D c_) noexcept
        : a(a_), b(b_), c(c_) {}

    /**
     * @brief Ориентированная (знаковая) площадь треугольника.
     * @return > 0 — CCW, < 0 — CW, = 0 — вырожденный (коллинеарные точки).
     *
     * Формула: 0.5 * ((b-a) × (c-a))
     */
    [[nodiscard]] constexpr double signedArea() const noexcept {
        return 0.5 * ((b - a).cross(c - a));
    }

    /**
     * @brief Площадь треугольника (всегда неотрицательна).
     */
    [[nodiscard]] double area() const noexcept {
        return std::abs(signedArea());
    }

    /**
     * @brief Проверка вырожденности (площадь меньше eps).
     */
    [[nodiscard]] bool isDegenerate(double eps = kEps) const noexcept {
        return std::abs(signedArea()) < eps;
    }

    /**
     * @brief Периметр (сумма длин сторон).
     */
    [[nodiscard]] double perimeter() const noexcept {
        return a.dist(b) + b.dist(c) + c.dist(a);
    }

    /**
     * @brief Центроид треугольника.
     */
    [[nodiscard]] constexpr Point2D centroid() const noexcept {
        return (a + b + c) / 3.0;
    }

    /**
     * @brief Центр описанной окружности (circumcenter).
     *
     * Вычисляется через перпендикулярные биссектрисы двух сторон.
     * @warning Для вырожденных треугольников результат некорректен.
     */
    [[nodiscard]] Point2D circumcenter() const noexcept;

    /**
     * @brief Радиус описанной окружности.
     * @warning Для вырожденных треугольников — бесконечность.
     */
    [[nodiscard]] double circumradius() const noexcept;

    /**
     * @brief Центр вписанной окружности (incenter).
     *
     * Взвешенное среднее вершин с весами — длинами противоположных сторон.
     */
    [[nodiscard]] Point2D incenter() const noexcept;

    /**
     * @brief Проверяет, лежит ли точка @p p внутри треугольника
     *        (включая границу).
     *
     * Использует знаки площадей подтреугольников. Граничные случаи
     * (точка на ребре или в вершине) считаются «внутри».
     */
    [[nodiscard]] bool contains(const Point2D& p,
                                double eps = kEps) const noexcept;

    /**
     * @brief Возвращает треугольник с гарантированной CCW-ориентацией.
     */
    [[nodiscard]] constexpr Triangle ccw() const noexcept {
        if (signedArea() >= 0.0) return *this;
        return {a, c, b}; // меняем местами b и c
    }

    /**
     * @brief Сторона i-я (0→ab, 1→bc, 2→ca).
     */
    [[nodiscard]] constexpr Segment side(int i) const noexcept {
        switch (i) {
        case 0:  return {a, b};
        case 1:  return {b, c};
        default: return {c, a};
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const Triangle& t) {
        return os << "Tri[" << t.a << ", " << t.b << ", " << t.c << ']';
    }
};

} // namespace geometry

// ─────────────────────────────────────────────────────────────────────────────
//  std::hash специализации для использования в unordered-контейнерах
// ─────────────────────────────────────────────────────────────────────────────

namespace std {

template <>
struct hash<geometry::Point2D> {
    std::size_t operator()(const geometry::Point2D& p) const noexcept {
        // FNV-подобное смешивание двух double
        auto h1 = std::hash<double>{}(p.x);
        auto h2 = std::hash<double>{}(p.y);
        return h1 ^ (h2 * 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
    }
};

} // namespace std
