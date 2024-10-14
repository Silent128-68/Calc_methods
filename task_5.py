import numpy as np

eps = 1e-6


def gauss(sole):
    n = len(sole)
    m = len(sole[0])

    x_coords = list(range(n))

    print("Расширенная матрица СЛАУ:")
    for row in sole:
        print("\t".join(f"{x:.10f}" for x in row))

    print("\nНомера неизвестных:")
    print("\t".join(str(x) for x in x_coords))
    print()

    for j in range(m - 1):
        mx = j
        for i in range(j, n):
            if abs(sole[i][j]) > abs(sole[mx][j]):
                mx = i

        if abs(sole[mx][j]) < eps:
            print("Division by zero!")
            return []

        sole[j], sole[mx] = sole[mx], sole[j]

        for i in range(j, n):
            if abs(sole[j][i]) > abs(sole[j][mx]):
                mx = i

        x_coords[j], x_coords[mx] = x_coords[mx], x_coords[j]
        for i in range(n):
            sole[i][j], sole[i][mx] = sole[i][mx], sole[i][j]

        if abs(sole[j][mx]) < eps:
            print("Division by zero!")
            return []

        for i in range(j + 1, n):
            d = sole[i][j] / sole[j][j]
            for k in range(m):
                sole[i][k] -= sole[j][k] * d

    print("Прямой ход метода Гаусса:")
    for row in sole:
        print("\t".join(f"{x:.10f}" for x in row))

    print("\nНомера неизвестных:")
    print("\t".join(str(x) for x in x_coords))
    print()

    for j in range(n - 1, -1, -1):
        for i in range(j - 1, -1, -1):
            d = sole[i][j] / sole[j][j]
            sole[i][j] -= sole[j][j] * d
            sole[i][m - 1] -= sole[j][m - 1] * d

    print("Обратный ход метода Гаусса:")
    for row in sole:
        print("\t".join(f"{x:.10f}" for x in row))
    print()

    s = [0] * n
    for i in range(n):
        s[x_coords[i]] = sole[i][m - 1] / sole[i][i]

    return s


def gauss_with_b(A, B):
    for i in range(len(A)):
        A[i].append(B[i])
    return gauss(A)


def solve():
    n = 5
    v = 3
    A = [[v / 100 for _ in range(n)] for i in range(n)]

    for i in range(n):
        A[i][i] = v + 2 * i

    V = [v + 2 * i for i in range(n)]
    B = [sum(A[i][j] * V[j] for j in range(n)) for i in range(n)]

    result = gauss_with_b(A, B)

    print("Решение СЛАУ:")
    print("\t".join(f"{x:.10f}" for x in result))


if __name__ == "__main__":
    solve()
