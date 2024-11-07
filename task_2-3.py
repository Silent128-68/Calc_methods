import numpy as np

data = np.array([[-1, 0, 2, 3], [0, 1, 9, 28]])

print("Исходный массив:")
for row in data:
    print(' '.join(f"{num:>4}" for num in row))

x_values = data[0]
f_values = data[1]

n = len(x_values)


def lagrange_interpolation(x):
    result = 0.0
    for i in range(n):
        term = f_values[i]
        for j in range(n):
            if i != j:
                term *= (x - x_values[j]) / (x_values[i] - x_values[j])
        result += term
    return result


def newton_interpolation(x):
    result = x + 1 + 3 * (x + 1) * x + 1/2 * (x + 1) * x * (x - 2)
    return result


new_matrix = np.zeros((2, 2 * n - 1))

for i in range(n):
    new_matrix[0, 2 * i] = x_values[i]
    if i < n - 1:
        new_matrix[0, 2 * i + 1] = (x_values[i] + x_values[i + 1]) / 2

for i in range(n):
    new_matrix[1, 2 * i] = f_values[i]
    if i < n - 1:
        new_matrix[1, 2 * i + 1] = lagrange_interpolation((x_values[i] + x_values[i + 1]) / 2)

new_matrix_1 = np.zeros((2, 2 * n - 1))

for i in range(n):
    new_matrix_1[0, 2 * i] = x_values[i]
    if i < n - 1:
        new_matrix_1[0, 2 * i + 1] = (x_values[i] + x_values[i + 1]) / 2

for i in range(n):
    new_matrix_1[1, 2 * i] = f_values[i]
    if i < n - 1:
        new_matrix_1[1, 2 * i + 1] = newton_interpolation((x_values[i] + x_values[i + 1]) / 2)


print("\nРезультирующая матрица для формы Лагранжа:")
for row in new_matrix:
    formatted_row = ' '.join(f"{num:>7.3f}" for num in row)
    print(formatted_row)

print("\nРезультирующая матрица для формы Ньютона:")
for row in new_matrix:
    formatted_row = ' '.join(f"{num:>7.3f}" for num in row)
    print(formatted_row)
