import math
import random


def generate_labyrinth(target_n, filename="labyrinth.txt", width=100.0, height=100.0):
    if target_n < 7:
        print("Error: At least 7 vertices are required.")
        return

    frame_points = [
        (0.0, 0.0),
        (width, 0.0),
        (width, height),
        (0.0, height)
    ]

    vertices_left = target_n - 4
    poly_sizes = []

    while vertices_left > 0:
        if 3 <= vertices_left <= 6:
            poly_sizes.append(vertices_left)
            vertices_left = 0
        elif vertices_left < 3:
            if poly_sizes:
                poly_sizes[-1] += vertices_left
            else:
                poly_sizes.append(3)
            vertices_left = 0
        else:
            size = random.randint(3, 6)
            poly_sizes.append(size)
            vertices_left -= size

    num_obstacles = len(poly_sizes)

    grid_size = math.ceil(math.sqrt(num_obstacles))
    cells = [(i, j) for i in range(grid_size) for j in range(grid_size)]
    random.shuffle(cells)

    cell_w = width / grid_size
    cell_h = height / grid_size

    obstacles = []

    for idx, size in enumerate(poly_sizes):
        grid_i, grid_j = cells[idx]

        xmin = grid_i * cell_w
        xmax = (grid_i + 1) * cell_w
        ymin = grid_j * cell_h
        ymax = (grid_j + 1) * cell_h

        pad_x = cell_w * 0.15
        pad_y = cell_h * 0.15

        cx = (xmin + xmax) / 2.0
        cy = (ymin + ymax) / 2.0

        max_rx = (cell_w - 2 * pad_x) / 2.0
        max_ry = (cell_h - 2 * pad_y) / 2.0

        rx = random.uniform(max_rx * 0.5, max_rx)
        ry = random.uniform(max_ry * 0.5, max_ry)

        alpha = random.uniform(0, 2 * math.pi)

        angles = sorted([random.uniform(0, 2 * math.pi) for _ in range(size)])

        poly_points = []
        for angle in angles:
            lx = rx * math.cos(angle)
            ly = ry * math.sin(angle)

            x = cx + lx * math.cos(alpha) - ly * math.sin(alpha)
            y = cy + lx * math.sin(alpha) + ly * math.cos(alpha)

            poly_points.append((round(x, 2), round(y, 2)))

        obstacles.append(poly_points)

    with open(filename, "w", encoding="utf-8") as f:
        f.write("# === Outer edge ===\n")
        for p in frame_points:
            f.write(f"{p[0]} {p[1]}\n")
        f.write("\n")

        for obs_idx, obs in enumerate(obstacles):
            f.write(f"# Barier {obs_idx + 1} (vertices: {len(obs)})\n")
            for p in obs:
                f.write(f"{p[0]} {p[1]}\n")
            f.write("\n")

    total_generated = 4 + sum(len(o) for o in obstacles)
    print(f"Successfully generated '{filename}'!")
    print(f"Vertices: {total_generated} (Edge: 4, Bariers: {num_obstacles})")


if __name__ == "__main__":
    generate_labyrinth(target_n=800, filename="labyrinth_big.txt", width=50.0, height=50.0)
