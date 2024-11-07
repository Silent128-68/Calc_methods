import math

# Константы
eps = 1e-6
V = 16


def f(x, y):
    return 2 * V * x + V * x ** 2 - y


def check_y(x):
    return V * x ** 2


def euler(n, h, x0, y0):
    x = [0] * n
    y = [0] * n
    x[0] = x0
    y[0] = y0
    for i in range(n - 1):
        x[i + 1] = x[i] + h
        y[i + 1] = y[i] + h * f(x[i], y[i])
    return y


def ex_euler(n, h, x0, y0):
    x = [0] * (2 * n)
    y = [0] * (2 * n)
    x[0] = x0
    y[0] = y0
    for i in range(2 * n - 1):
        x[i + 1] = x[i] + h / 2
        if i % 2:
            y[i + 1] = y[i - 1] + h * f(x[i], y[i])
        else:
            y[i + 1] = y[i] + (h / 2) * f(x[i], y[i])
    res = [y[2 * i] for i in range(n)]
    return res


def pre_corr(n, h, x0, y0):
    x = [0] * n
    y = [0] * n
    x[0] = x0
    y[0] = y0
    for i in range(n - 1):
        x[i + 1] = x[i] + h
        y[i + 1] = y[i] + h * f(x[i], y[i])
        y[i + 1] = y[i] + (h / 2) * (f(x[i], y[i]) + f(x[i + 1], y[i + 1]))
    return y


def print_results(vx, vy):
    header = f"{'x':>10} {'y':>15} {'check_y(x)':>15} {'Deviation':>15}"
    print(header)
    print("-" * len(header))

    deviations = [vy[i] - check_y(vx[i]) for i in range(len(vy))]
    for i in range(len(vx)):
        print(f"{vx[i]:>10.5f} {vy[i]:>15.5f} {check_y(vx[i]):>15.5f} {deviations[i]:>15.5f}")

    max_deviation = max(abs(dev) for dev in deviations)
    print("-" * len(header))
    print(f"{'Максимальное отклонение:':<35}{max_deviation:.5f}")


def solve():
    n = 9
    h = 1.0
    x0 = 1.0
    y0 = V
    vx = [x0 + i * h for i in range(n)]

    print("Метод Эйлера:")
    vy1 = euler(n, h, x0, y0)
    print_results(vx, vy1)

    print("\nРасширенный метод Эйлера:")
    vy2 = ex_euler(n, h, x0, y0)
    print_results(vx, vy2)

    print("\nМетод предварительного и корректирующего счета:")
    vy3 = pre_corr(n, h, x0, y0)
    print_results(vx, vy3)


if __name__ == "__main__":
    solve()
