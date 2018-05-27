import math

from sympy import symbols, diff
from sympy.abc import s
from sympy.matrices import Matrix, eye
from sympy.parsing.sympy_parser import parse_expr

ROUNDING_NUMBER = 3  # число знаков, до которого округлять (только при выводе, считает точно)

NUMBER_OF_VARS = 2  # число переменных
x_vars = symbols('x1:%d' % (NUMBER_OF_VARS + 1))  # создаем эти переменные

func = 'x1**2 + 2*x2**2 + x1*x2 - 7*x1 - 7*x2'  # исходная функция
func = parse_expr(func)

init_vector = [-1, 0]  # начальный вектор
init_vector = Matrix(init_vector)

EPS = 0.01  # погрешность


def gradient(f, p):
    return Matrix(
        [diff(f, x_vars[i]).subs({x_vars[i]: p[i] for i in range(NUMBER_OF_VARS)}) for i in range(NUMBER_OF_VARS)])


def rounding_vectors(vector):  # Функция для округления векторов (только при выводе)
    vector = list(vector.evalf())
    return [round(el, ROUNDING_NUMBER) for el in vector]


def f_s(point, f):  # функция для подстановки аргумента в одномерную функцию
    return f.subs(s, point)


def f(point):  # функция для подстановки аргументов в исходную функцию
    f = func
    for i in range(NUMBER_OF_VARS):
        f = f.subs(x_vars[i], point[i])
    return f


def gold_ratio(f):  # Метод золотого сечения
    a = -50  # границы интервала поиска
    b = 50
    PHI = (math.sqrt(5) - 1) / 2
    alpha = PHI * a + (1 - PHI) * b
    beta = (1 - PHI) * a + PHI * b

    while abs(b - a) > EPS:  # основной цикл
        if f_s(alpha, f) > f_s(beta, f):  # сравниваем значение функции в точках alpha и beta
            a = alpha
            alpha = beta
            beta = a + PHI * (b - a)
        else:
            b = beta
            beta = alpha
            alpha = a + (1 - PHI) * (b - a)

    answer = (a + b) / 2
    return answer


def dich(f):  # Метод дихотомии
    a = -50  # границы интервала поиска
    b = 50
    delta = EPS / 2

    while abs(b - a) > EPS:  # основной цикл
        x1 = (a + b - delta) / 2
        x2 = (a + b + delta) / 2
        if f_s(x1, f) <= f_s(x2, f):  # сравниваем значение функции в точках x1 и x2
            b = x2
        else:
            a = x1

    answer = (a + b) / 2
    return answer


def dfp(method):
    def solve(f):  # Вспомогательная функция
        if method == 1:
            return gold_ratio(f)
        elif method == 2:
            return dich(f)

    D = eye(NUMBER_OF_VARS)
    x0 = init_vector
    x = x0 - s * D * gradient(func, x0)
    S = solve(f(list(x)))  # Находим S из решения задачи одномерной минимизации
    x = x.subs(s, S)  # Подставляем S в y

    while (x - x0).norm() > EPS:
        U = x - x0
        V = gradient(func, x) - gradient(func, x0)
        A = U * U.T / (U.T * V)[0]
        B = - D * V * V.T * D / (V.T * D * V)[0]
        D = D + A + B
        x0 = x
        x = x0 - s * D * gradient(func, x0)
        S = solve(f(list(x)))  # Находим S из решения задачи одномерной минимизации
        x = x.subs(s, S)  # Подставляем S в x

    return x


print('Исходная функция: ', func)

# Вызываем метод сопряженных направлений, передаем параметр method, равный
# 1 - Метод золотого сечения
# 2 - Метод дихотомии
min_x = dfp(1)

print('\nОтвет: x_min = ', rounding_vectors(min_x))
