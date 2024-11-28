eimport math
from itertools import product

V = 22.0
T = V
eps = 1e-6


def f(x):
    return 4 * V * x ** 4 - 3 * V * T * x ** 3 + 6 * V * x - 2 * V * T


def p(x):
    return x ** 2


def q(x):
    return x


def check_y(x):
    return V * x * x * (x - T)


def phi(x, k):
    k += 1
    return x ** (k + 2) * (x - T)


def phi1(x, k):
    k += 1
    return (k + 2) * x ** (k + 1) * (x - T) + x ** (k + 2)


def phi2(x, k):
    k += 1
    # Вторая производная phi(x, k) = x ** (k + 2) * (x - T)
    term1 = (k + 2) * (k + 1) * x ** k * (x - T)
    term2 = 2 * (k + 2) * x ** (k + 1)
    term3 = x ** k * (k + 1) * k * T
    return term1 + term2 - term3


def gauss(sole):
    n = len(sole)
    m = len(sole[0])
    x_coords = list(range(n))

    for j in range(m - 1):
        mx = max(range(j, n), key=lambda i: abs(sole[i][j]))
        if abs(sole[mx][j]) < eps:
            sole[j], sole[mx] = sole[mx], sole[j]
        else:
            mx = j
            for i in range(j, n):
                if abs(sole[j][i]) > abs(sole[j][mx]):
                    mx = i
            x_coords[j], x_coords[mx] = x_coords[mx], x_coords[j]
            for i in range(n):
                sole[i][mx], sole[i][j] = sole[i][j], sole[i][mx]

        if abs(sole[j][mx]) < eps:
            print("Division by zero!")
            return []

        for i in range(j + 1, n):
            d = sole[i][j] / sole[j][j]
            for k in range(m):
                sole[i][k] -= sole[j][k] * d

    for j in range(n - 1, -1, -1):
        for i in range(j - 1, -1, -1):
            d = sole[i][j] / sole[j][j]
            sole[i][j] -= sole[j][j] * d
            sole[i][m - 1] -= sole[j][m - 1] * d

    s = [0] * n
    for i in range(n):
        s[x_coords[i]] = sole[i][m - 1] / sole[i][i]
    return s


def gauss_with_rhs(A, B):
    for i in range(len(A)):
        A[i].append(B[i])
    return gauss(A)


def print_result(vx, vy):
    k = 10
    print(f"{'x':>12} | {'y':>12} | {'check_y(x)':>12} | {'y - check_y(x)':>16}")
    print("-" * 66)
    for i in range(0, len(vx), k):
        m = min(i + k, len(vx))
        for j in range(i, m):
            print(f"{vx[j]:>12.5f} | {vy[j]:>12.5f} | {check_y(vx[j]):>12.5f} | {(vy[j] - check_y(vx[j])):>16.5f}")
        print()


def solve():
    n = 20
    l, r = 0, T
    h = (r - l) / (n + 1)
    vx = [l + h * (i + 1) for i in range(n)]

    A = [[phi2(vx[j], k) + p(vx[j]) * phi1(vx[j], k) + q(vx[j]) * phi(vx[j], k) for k in range(n)] for j in range(n)]
    B = [f(x) for x in vx]
    va = gauss_with_rhs(A, B)

    vx.insert(0, 0)
    vx.append(T)
    vy = [0] * len(vx)

    for i in range(len(vx)):
        for k in range(n):
            vy[i] += va[k] * phi(vx[i], k)

    print_result(vx, vy)


if __name__ == "__main__":
    solve()
