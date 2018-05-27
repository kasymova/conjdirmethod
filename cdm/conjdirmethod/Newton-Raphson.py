import math

from sympy import symbols, diff
from sympy.abc import s
from sympy.matrices import Matrix, hessian
from sympy.parsing.sympy_parser import parse_expr

NUMBER_OF_VARS = 2  # число переменных
x_vars = symbols('x1:%d' % (NUMBER_OF_VARS + 1))  # создаем эти переменные

func = 'x1**2 + 2*x2**2 + x1*x2 - 7*x1 - 7*x2' # исходная функция
#func = 'x1**2 + 2*(x2-1)**2 + (x3-4)**2'  # исходная функция трех переменных
func = parse_expr(func)

init_vector = [-10, 10]  # начальный вектор
init_vector = Matrix(init_vector)

EPS = 0.01  # погрешность


def gradient(f, p):  # Функция для вычисления градиента функции f в точке p
    return Matrix(
        [diff(f, x_vars[i]).subs({x_vars[i]: p[i] for i in range(NUMBER_OF_VARS)}) for i in range(NUMBER_OF_VARS)])


def rounding_vectors(vector):  # Функция для округления векторов (только при выводе)
    ROUNDING_NUMBER = 3  # число знаков, до которого округлять (только при выводе, считает точно)
    vector = list(vector.evalf())
    return [round(el, ROUNDING_NUMBER) for el in vector]


def f_s(point, f):  # функция для подстановки аргумента в одномерную функцию
    return f.subs(s, point)


def f(point):  # функция для подстановки аргументов в исходную функцию
    f = func
    for i in range(NUMBER_OF_VARS):
        f = f.subs(x_vars[i], point[i])
    return f


def zolotoe_sechenie(f):  # Метод золотого сечения
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


def dichotomiya(f):  # Метод дихотомии
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


def positive_definite(matrix):  # Функция для исследования матрицы на положительную определенность
    eigenvals = matrix.eigenvals()
    for value in eigenvals.keys():
        if value <= 0:
            return False
    return True


def newton_raphson():
    x0 = init_vector
    x = 2 * x0  # Положим такой x, чтобы зашло в цикл

    counter = 0 # счетчик числа итераций
    while (x - x0).norm() > EPS:  # Пока норма разности x и x0 больше EPS
        counter += 1
        print('\tИтерация', counter)
        # Находим матрицу Гессе исходной функции в точке x
        H = hessian(func, x_vars)
        H = H.subs({x_vars[i]: init_vector[i] for i in range(NUMBER_OF_VARS)})

        if positive_definite(H.inv()):  # Если обратная матрица Гессе положительно определена
            d = -H.inv() * gradient(func, x)
        else:  # иначе
            d = -gradient(func, x)

        x0 = x
        x = x0 + s * d
        S = zolotoe_sechenie(f(list(x)))  # Находим S из решения задачи одномерной минимизации методом золотого сечения
        # S = dichotomiya(f(list(x))) # или методом дихотомии
        x = x.subs(s, S)  # Подставляем S в x
        print('\t\tx = ', rounding_vectors(x))

    return x


print('Исходная функция: ', func)

# Вызываем метод Ньютона-Рафсона
min_x = newton_raphson()

print('\nОтвет: x_min = ', rounding_vectors(min_x))
