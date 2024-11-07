import numpy as np

eps = 1e-6


def gauss(sole):
    n = len(sole)
    m = len(sole[0])
    x_coords = list(range(n))

    for j in range(m - 1):
        mx = j
        for i in range(j, n):
            if abs(sole[i][j]) > abs(sole[mx][j]):
                mx = i
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
            raise ValueError("Division by zero!")

        for i in range(j + 1, n):
            d = sole[i][j] / sole[j][j]
            for k in range(m):
                sole[i][k] -= sole[j][k] * d

    for j in range(n - 1, -1, -1):
        for i in range(j - 1, -1, -1):
            d = sole[i][j] / sole[j][j]
            sole[i][j] -= sole[j][j] * d
            sole[i][m - 1] -= sole[j][m - 1] * d

    s = np.zeros(n)
    for i in range(n):
        s[x_coords[i]] = sole[i][m - 1] / sole[i][i]
    return s


def gauss_solution(A, B):
    for i in range(len(A)):
        A[i].append(B[i])
    return gauss(A)


def determinant(A):
    n = len(A)
    sg = 1
    for j in range(n):
        mx = j
        for i in range(j, n):
            if abs(A[i][j]) > abs(A[mx][j]):
                mx = i
        if abs(A[mx][j]) < eps:
            sg *= 1 - 2 * (mx != j)
            A[j], A[mx] = A[mx], A[j]
        else:
            mx = j
            for i in range(j, n):
                if abs(A[j][i]) > abs(A[j][mx]):
                    mx = i
            for i in range(n):
                A[i][mx], A[i][j] = A[i][j], A[i][mx]
            sg *= 1 - 2 * (mx != j)

        for i in range(j + 1, n):
            d = A[i][j] / A[j][j]
            for k in range(n):
                A[i][k] -= A[j][k] * d

    res = 1
    for i in range(n):
        res *= A[i][i]
    return res * sg


def inverse(A):
    n = len(A)
    E = np.identity(n)
    B = np.array(A, dtype=float)
    for j in range(n):
        res = gauss_solution(A.copy(), E[j].tolist())
        for i in range(n):
            B[i][j] = res[i]
    return B


def multiply(A, B):
    n, m, q = len(A), len(A[0]), len(B[0])
    C = np.zeros((n, q))
    for i in range(n):
        for j in range(q):
            for k in range(m):
                C[i][j] += A[i][k] * B[k][j]
    return C


def solve():
    A = [
        [-3, 1, 3, 4],
        [3, 0, -1, 4],
        [-5, 2, 3, 0],
        [4, -1, 2, -6],
    ]
    print("Матрица A:")
    for row in A:
        print("\t".join(f"{x:.10f}" for x in row))

    det_A = determinant(A.copy())
    print("\nОпределитель:", det_A)

    print("\nОбратная матрица:")
    inv_A = inverse(A.copy())
    for row in inv_A:
        print("\t".join(f"{x:.10f}" for x in row))

    print("\nПроверка. Произведение матрицы A на обратную:")
    C = multiply(A, inv_A)
    for row in C:
        print("\t".join(f"{x:.10f}" for x in row))


if __name__ == "__main__":
    solve()
