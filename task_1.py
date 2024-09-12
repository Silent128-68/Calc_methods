import numpy as np

data = np.array([[-1, 0, 2, 3], [0, 1, 9, 28]])

print("Исходный массив:")
for row in data:
    print(' '.join(f"{num:>4}" for num in row))

x_values = data[0]
f_values = data[1]

n = len(x_values)

matrix = np.zeros((n, n + 1), dtype=int)

for i in range(n):
    for j in range(n):
        matrix[i, j] = x_values[i] ** (n - j - 1)
    matrix[i, -1] = f_values[i]

# Решение СЛАУ
A = matrix[:, :-1]  # Матрица коэффициентов
b = matrix[:, -1]   # Вектор свободных членов

coefficients = np.linalg.solve(A, b)
a0, a1, a2, a3 = coefficients

print("\nКоэффициенты a:")
for i, coeff in enumerate(coefficients):
    print(f"a{i} = {coeff:.2f}")

def p3(x):
    return a3 * x**3 + a2 * x**2 + a1 * x + a0

new_matrix = np.zeros((2, 2*n - 1))

for i in range(n):
    new_matrix[0, 2*i] = x_values[i]
    if i < n - 1:
        new_matrix[0, 2*i + 1] = (x_values[i] + x_values[i + 1]) / 2

for i in range(n):
    new_matrix[1, 2*i] = f_values[i]
    if i < n - 1:
        new_matrix[1, 2*i + 1] = p3((x_values[i] + x_values[i + 1]) / 2)

print("\nРезультирующая матрица:")
for row in new_matrix:
    formatted_row = ' '.join(f"{num:>7.3f}" for num in row)
    print(formatted_row)
