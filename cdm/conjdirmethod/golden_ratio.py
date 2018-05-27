import math

from sympy.abc import d
from sympy import sqrt
from sympy.utilities.lambdify import lambdify

EPS = 0.01  # погрешность


def golden_ratio(f):  # Метод золотого сечения
    """
    Принимает:
    f - функция от одной переменной f(s), которая получилась после подстановки аргументов в иходную функцию
    Возвращает:
    answer (число) - средняя точка на интервале (a, b)
    """
    f = lambdify(d, f)

    a = -3  # границы интервала поиска
    b = 3
    PHI = (math.sqrt(5) - 1) / 2
    i = 0
    print('\nМетод золотого сечения:')
    while abs(b - a) > EPS:  # основной цикл
        i += 1
        alpha = PHI * a + (1 - PHI) * b
        beta = (1 - PHI) * a + PHI * b
        print('\tИтерация', str(i) + ':')
        print('\t\ta =', round(a, 4), '\t\tb =', round(b, 4))

        if f(alpha) > f(beta):  # сравниваем значение функции в точках alpha и beta
            a = alpha
        else:
            b = beta

    answer = (a + b) / 2
    return answer


el_a = 3
el_b = 4
f = 4 * (d + el_b * sqrt(1 - d ** 2 / el_a ** 2))
ans = golden_ratio(f)
print('\nLoss function:', f.evalf(7), '\nd =', round(ans, 4), '\nh =', round(2000 / math.pi / ans ** 2, 4))
