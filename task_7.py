import math

eps = 1e-6
V = 22.0
T = V


def f(x):
    return 4 * V * x ** 4 - 3 * V * T * x ** 3 + 6 * V * x - 2 * V * T


def p(x):
    return x ** 2


def q(x):
    return x


def check_y(x):
    return V * x ** 2 * (x - T)


def tridiagonal_algorithm(a, b, c, d):
    n = len(a)
    p = [0.0] * (n + 1)
    q = [0.0] * (n + 1)

    for i in range(n):
        p[i + 1] = c[i] / (b[i] - a[i] * p[i])
        q[i + 1] = (a[i] * q[i] - d[i]) / (b[i] - a[i] * p[i])

    res = [0.0] * n
    res[n - 1] = q[n]
    for i in range(n - 2, -1, -1):
        res[i] = p[i + 1] * res[i + 1] + q[i + 1]

    return res


def print_column(data, width=20):
    for value in data:
        print(f"{value:<{width}.5f}", end=' ')
    print()


def print_results(vx, vy):
    k = 10  # Interval for printing

    for i in range(0, len(vx), k):
        m = min(i + k, len(vx))

        # Print each row with proper alignment
        print("x:              ", end='')
        print_column(vx[i:m])

        print("y (computed):   ", end='')
        print_column(vy[i:m])

        print("y (check_y):    ", end='')
        print_column([check_y(x) for x in vx[i:m]])

        print("Difference:     ", end='')
        print_column([vy[j] - check_y(vx[j]) for j in range(i, m)])

        print()

    # Find maximum deviation, its corresponding x, and index
    max_deviation = 0
    max_x = 0
    max_index = 0
    for i in range(len(vy)):
        deviation = abs(vy[i] - check_y(vx[i]))
        if deviation > max_deviation:
            max_deviation = deviation
            max_x = vx[i]
            max_index = i

    print(f"Максимальное отклонение: {max_deviation:.5f} в точке x = {max_x:.5f} (индекс {max_index})")


def solve():
    n = 20000
    l, r = 0.0, T
    h = (r - l) / n

    vx = [l + i * h for i in range(n + 1)]

    a = [0.0] * (n - 1)
    b = [0.0] * (n - 1)
    c = [0.0] * (n - 1)
    d = [0.0] * (n - 1)

    for i in range(n - 1):
        x = vx[i + 1]
        a[i] = (1 / h ** 2 - p(x) / (2 * h))
        b[i] = -(-2 / h ** 2 + q(x))
        c[i] = (1 / h ** 2 + p(x) / (2 * h))
        d[i] = f(x)

    a[0] = 0.0
    c[-1] = 0.0

    vy = tridiagonal_algorithm(a, b, c, d)
    vy.insert(0, 0.0)
    vy.append(0.0)

    print_results(vx, vy)
    print()


if __name__ == "__main__":
    solve()
