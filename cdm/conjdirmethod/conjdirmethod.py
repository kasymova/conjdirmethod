import math
from functools import reduce

from sympy import symbols
from sympy.abc import s
from sympy.matrices import Matrix, hessian
from sympy.parsing.sympy_parser import parse_expr
from sympy.utilities.lambdify import lambdify

ROUNDING_NUMBER = 3  # число знаков, до которого округлять (только при выводе, считает точно)
MAX_ITER = 20  # максимальное число итераций

NUMBER_OF_VARS = 2  # число переменных
x_vars = symbols('x1:%d' % (NUMBER_OF_VARS + 1))  # создаем эти переменные

loss_function = 'x1**2 + 2*x2**2 + x1*x2 - 7*x1 - 7*x2'  # целевая функция
loss_function = parse_expr(loss_function)
f = lambdify(x_vars, loss_function)  # создаем ее

init_vector = [-15, 10]  # начальный вектор
init_vector = Matrix(init_vector)  # относим его к классу матриц

EPS = 0.001  # погрешность

listmerge = lambda s: reduce(lambda d, el: d.extend(el) or d, s, [])  # системная функция для слияния списков


def rounding_vectors(vector):  # Функция для округления векторов (только при выводе)
    vector = list(vector.evalf())
    return [round(el, ROUNDING_NUMBER) for el in vector]


def golden_ratio(f):  # Метод золотого сечения
    """
    Принимает:
    f - функция от одной переменной f(s), которая получилась после подстановки аргументов в иходную функцию
    Возвращает:
    answer (число) - средняя точка на интервале (a, b)
    """
    f = lambdify(s, f)

    a = -50  # границы интервала поиска
    b = 50
    PHI = (math.sqrt(5) - 1) / 2

    while abs(b - a) > EPS:  # основной цикл
        alpha = PHI * a + (1 - PHI) * b
        beta = (1 - PHI) * a + PHI * b

        if f(alpha) > f(beta):  # сравниваем значение функции в точках alpha и beta
            a = alpha
        else:
            b = beta

    answer = (a + b) / 2
    return answer


def dichotomy(f):  # Метод дихотомии
    """
    Принимает:
    f - функция от одной переменной f(s), которая получилась после подстановки аргументов в иходную функцию
    Возвращает:
    answer (число) - средняя точка на интервале (a, b)
    """
    f = lambdify(s, f)

    a = -50  # границы интервала поиска
    b = 50
    delta = EPS / 2

    while abs(b - a) > EPS:  # основной цикл
        x1 = (a + b - delta) / 2
        x2 = (a + b + delta) / 2

        if f(x1) <= f(x2):  # сравниваем значение функции в точках x1 и x2
            b = x2
        else:
            a = x1

    answer = (a + b) / 2
    return answer


def conjdirmethod():  # Метод сопряженных направлений
    """
    Возвращает:
    x (вектор) - точка минимума исходной функции
    """

    H = hessian(loss_function, x_vars)  # находим матрицу Гессе исходной функции
    v = [el.evalf() for el in listmerge([el[2] for el in H.eigenvects()])]  # находим собственные векторы матрицы Гессе

    x_0 = init_vector  # в качестве x_0 берем начальный вектор

    i = 0  # счетчик числа итераций
    # Внешний цикл - отвечает за число итераций, работает до тех пор, пока норма разности x и x0 больше EPS
    while True:
        i += 1
        print('\nИтерация ' + str(i) + ':')
        print('\tШаг 1:')

        # В качестве направления выбираем первый собственный вектор матрицы Гессе, считаем
        # его в точке x_0 (так как у функций, порядок которых выше, чем 2, матрица Гессе - переменная, зависящая
        # в общем случае от всех переменных x1, x2, ... и, следовательно, ее собственные векторы тоже переменные)
        # Если матрица Гессе состоит из констант, то никакой подстановки не произойдет
        d = [v[0].subs({x_vars[k]: x_0[k] for k in range(NUMBER_OF_VARS)})]

        x = x_0 + s * d[0]
        S = dichotomy(f(*x))  # Находим S из решения задачи одномерной минимизации
        x = x.subs(s, S)  # Подставляем S в x
        print('\t\tx = ', rounding_vectors(x))  # Выводим x на экран

        # Второй цикл - отвечает за количество дальнейших шагов (2, 3, ...), которое зависит от числа переменных
        for j in range(1, NUMBER_OF_VARS):
            print('\tШаг ' + str(j + 1) + ':')

            # Находим y0 как x + (j+1)-ый (2-ой, 3-ий) собственный вектор матрицы Гессе,
            # аналогично считаем его в точке x_0
            y0 = x + v[j].subs({x_vars[k]: x_0[k] for k in range(NUMBER_OF_VARS)})

            y = y0  # для дальнейшего суммирования

            # Внутренний цикл для последовательной минимизации y по направлениям d0, d1, ...
            for k in range(0, j):
                y = y + s * d[k]
                S = dichotomy(f(*y))  # Находим S из решения задачи одномерной минимизации
                y = y.subs(s, S)  # Подставляем S в y

            d.append(x - y)  # Находим следующее d как x - y

            x = x + s * d[j]
            S = dichotomy(f(*x))  # Находим S из решения задачи одномерной минимизации
            x = x.subs(s, S)  # Подставляем S в y
            print('\t\tx = ', rounding_vectors(x))

        # Условие выхода из внешнего цикла - если норма разности x и x_0 становится меньше EPS
        if (x - x_0).norm() <= EPS:
            break
        # Условие выхода из внешнего цикла - если привышено максимальное число итерации (вдруг начало расходиться)
        if i > MAX_ITER:
            break
        # Полагаем x как x_0 для следующей итерации
        x_0 = x

    return x


print('Целевая функция: ', loss_function)
print('Начальный вектор: ', rounding_vectors(init_vector))

# Вызываем метод сопряженных направлений
min_x = conjdirmethod()

print('\nОтвет: x_min = ', rounding_vectors(min_x))
