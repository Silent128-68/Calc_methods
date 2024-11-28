import numpy as np

def f(x):
    V = 22
    return V * (4 / 3 * x + 1 / 4 * x ** 2 + 1 / 5 * x ** 3)


def check_y(x):
    V = 22
    return V * x

def gauss(A, B):
    n = len(A)
    A = np.array(A, dtype=float)
    B = np.array(B, dtype=float)
    for j in range(n):
        max_row = max(range(j, n), key=lambda i: abs(A[i, j]))
        if abs(A[max_row, j]) < 1e-6:
            raise ValueError("Division by zero!")
        A[[j, max_row]] = A[[max_row, j]]
        B[[j, max_row]] = B[[max_row, j]]
        for i in range(j + 1, n):
            factor = A[i, j] / A[j, j]
            A[i, j:] -= factor * A[j, j:]
            B[i] -= factor * B[j]

    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (B[i] - np.dot(A[i, i + 1:], x[i + 1:])) / A[i, i]
    return x


def print_results(vx, vy):
    k = 5
    for i in range(0, len(vx), k):
        m = min(i + k, len(vx))
        print("\t".join(f"{vx[j]: .6f}" for j in range(i, m)))
        print("\t".join(f"{vy[j]: .6f}" for j in range(i, m)))
        print("\t".join(f"{check_y(vx[j]): .6f}" for j in range(i, m)))
        print("\t".join(f"{vy[j] - check_y(vx[j]): .6f}" for j in range(i, m)))
        print()


def solve():
    n = 3
    A = [[1 / (i + j + 3) for j in range(n)] for i in range(n)]
    print("A:")
    print("\n".join(" ".join(f"{x: .6f}" for x in row) for row in A))
    for i in range(n):
        A[i][i] += 1
    print("B:")
    B = [12.032777777777777, 9.166666666666666, 7.411904761904761]
    print(" ".join(f"{x}" for x in B))
    vq = gauss(A, B)
    m = 3
    l, r = 0, 1
    h = (r - l) / (m - 1)
    vx = [l + i * h for i in range(m)]

    vy = []
    for x in vx:
        y = f(x)
        for k in range(n):
            y -= vq[k] * x ** (k + 1)
        vy.append(y)

    print_results(vx, vy)



solve()
