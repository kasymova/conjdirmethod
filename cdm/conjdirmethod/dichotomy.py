import math

from sympy.abc import x
from sympy.utilities.lambdify import lambdify
from sympy import sqrt

EPS = 0.001  # погрешность


def dichotomy(f):  # Метод дихотомии
    """
    Принимает:
    f - функция от одной переменной f(s), которая получилась после подстановки аргументов в иходную функцию
    Возвращает:
    answer (число) - средняя точка на интервале (a, b)
    """
    f = lambdify(x, f)

    a = 0  # границы интервала поиска
    b = 3
    delta = EPS / 2
    i = 0
    print('\nМетод дихотомии:')
    while abs(b - a) > EPS:  # основной цикл
        i += 1
        x1 = (a + b - delta) / 2
        x2 = (a + b + delta) / 2
        print('\tИтерация', str(i) + ':')
        print('\t\ta =', round(a, 4), '\t\tb =', round(b, 4))

        if f(x1) > f(x2):  # сравниваем значение функции в точках x1 и x2
            b = x2
        else:
            a = x1

    answer = (a + b) / 2
    return answer


el_a = 3
el_b = 4
f = 4 * (x + el_b * sqrt(1 - x ** 2 / el_a ** 2))
ans = dichotomy(f)
print('\nLoss function:', f.evalf(7), '\nx =', round(ans, 4), '\ny =',
      round(el_b * math.sqrt(1 - ans ** 2 / el_a ** 2), 4))
