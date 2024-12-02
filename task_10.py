import numpy as np

def f(x):
    """Функция правой части."""
    V = 22
    return V * (4 / 3 * x + 1 / 4 * x ** 2 + 1 / 5 * x ** 3)

def kernel(x, t):
    """Ядро интегрального уравнения."""
    return 1 / (1 + x + t)

def gauss(A, B):
    """Решение СЛАУ методом Гаусса."""
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

def solve():
    # Параметры квадратуры
    n = 4  # Число квадратурных узлов
    l, r = 0, 1  # Границы интегрирования
    h = (r - l) / n  # Шаг квадратуры
    t = np.array([l + i * h for i in range(n)])  # Узлы (левые концы интервалов)

    # Формирование матрицы A и вектора B
    lambda_ = 1  # Параметр \lambda
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            A[i, j] = -lambda_ * kernel(t[i], t[j]) * h
        A[i, i] += 1  # Добавляем единичную матрицу (диагональный член)

    # Правая часть
    B = np.array([f(ti) for ti in t])

    # Решение СЛАУ методом Гаусса
    y = gauss(A, B)

    # Вывод результатов
    print("Результаты в узлах квадратуры:")
    for i in range(n):
        print(f"x = {t[i]:.6f}, y(x) = {y[i]:.6f}")

    # Проверка аналитической формулы
    print("\nПроверка аналитической формулы:")
    for i in range(n):
        analytical_y = f(t[i]) + lambda_ * sum(h * kernel(t[i], t[k]) * y[k] for k in range(n))
        print(f"x = {t[i]:.6f}, y(x) = {y[i]:.6f}, Analytical y(x) = {analytical_y:.6f}")

solve()
