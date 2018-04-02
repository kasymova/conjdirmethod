from django.shortcuts import render
from django.http import HttpResponse
import json
import math
import re
from sympy import Function, symbols, latex, Symbol, diff
from sympy.parsing.sympy_parser import parse_expr
from sympy.matrices import eye, Matrix, hessian
from functools import reduce

ans = ''
ROUNDING_NUMBER = 4

listmerge = lambda s: reduce(lambda d, el: d.extend(el) or d, s, [])


def index(request):
    def f_s(x_var, f):
        return f.subs(s, x_var)

    def f(x_mass):
        f = func
        for i in range(NUMBER_OF_VARS):
            f = f.subs(x_vars[i], x_mass[i])
        return f

    def golden_ratio(f):
        global ans

        a = -100
        b = 100
        PHI = (math.sqrt(5) - 1) / 2
        alpha = PHI * a + (1 - PHI) * b
        beta = (1 - PHI) * a + PHI * b
        while abs(b - a) > EPS:
            if f_s(alpha, f) > f_s(beta, f):
                a = alpha
                alpha = beta
                beta = a + PHI * (b - a)
            else:
                b = beta
                beta = alpha
                alpha = a + (1 - PHI) * (b - a)

        answer = round((a + b) / 2, ROUNDING_NUMBER)
        ans += 'Метод золотого сечения: ' + '\( s = ' + str(answer) + '\)<br>'
        return answer

    def dichotomy(f):
        global ans

        a = -100
        b = 100
        delta = EPS / 2
        while abs(b - a) > EPS:
            x1 = (a + b - delta) / 2
            x2 = (a + b + delta) / 2
            if f_s(x1, f) <= f_s(x2, f):
                b = x2
            else:
                a = x1

        answer = round((a + b) / 2, ROUNDING_NUMBER)
        ans += 'Метод дихотомии: ' + '\( s = ' + str(answer) + '\)<br>'
        return answer

    def conjdirmethod(method):
        global ans

        def solve(f):
            if method == 1:
                return golden_ratio(f)
            elif method == 2:
                return dichotomy(f)

        H = hessian(func, x_vars)
        ans += 'Матрица Гессе: $$ H = ' + latex(H) + ' $$<br>'

        v = listmerge([el[2] for el in H.eigenvects()])
        ans += 'Собственные векторы матрицы Гессе: $$ v_{}'.format(
            'v_'.join([str(i + 1) + '=' + latex(v[i]) for i in range(NUMBER_OF_VARS)])) + ' $$<br>'
        v = [el.evalf(ROUNDING_NUMBER) for el in v]

        x_0 = init_vector
        i = 0
        while True:
            i += 1
            ans += '<b>Итерация №' + str(i) + ':</b><br>'

            ans += 'Шаг 1:<br>'
            d = v[0]
            ans += '$$ d_0 = v_1 = ' + latex(v[0]) + '$$'

            x = x_0 + s * d
            ans += '$$ x_1 = argmin f(x_0 + s \cdot d_0) = argmin f(' + latex(x_0) + '+ s \cdot' + latex(
                d) + ') = argmin f(' + latex(x) + ') $$'

            S = solve(f(list(x)))
            x = (x.subs(s, S)).evalf(ROUNDING_NUMBER)
            ans += '$$ x_1 = ' + latex(x) + '$$'

            for j in range(1, NUMBER_OF_VARS):
                J = str(j)
                ans += 'Шаг ' + str(j + 1) + ':<br>'
                y0 = x + v[j]
                ans += '$$ y_' + str(j - 1) + '= x_' + J + ' + v_' + str(j + 1) + '=' + latex(x) + '+' + latex(
                    v[j]) + '=' + latex(y0) + '$$'

                y1 = y0 + s * d
                ans += '$$ y_' + J + '= argmin f(y_' + str(j - 1) + '+ s \cdot d_' + str(
                    j - 1) + ')= argmin f(' + latex(y0) + '+ s \cdot' + latex(d) + ') = argmin f(' + latex(y1) + ') $$'

                S = solve(f(list(y1)))
                y1 = (y1.subs(s, S)).evalf(ROUNDING_NUMBER)
                ans += '$$ y_' + str(j) + '=' + latex(y1) + '$$'

                d = x - y1
                ans += '$$ d_' + J + '= x_' + J + '- y_' + J + '=' + latex(x) + '-' + latex(y1) + '=' + latex(d) + '$$'
                ans += '$$ x_' + str(j + 1) + '= argmin f(x_' + J + '+ s \cdot d_' + J + ') = argmin f(' + latex(
                    x) + '+ s \cdot' + latex(d) + ') = argmin f('

                x = x + s * d
                ans += latex(x) + ')$$'

                S = solve(f(list(x)))
                x = (x.subs(s, S)).evalf(ROUNDING_NUMBER)
                ans += '$$ x_' + str(j+1) + '=' + latex(x) + ' $$'

            if (x - x_0).norm() <= EPS:
                break

            x_0 = x
        return x

    if request.method == 'GET':
        return render(request, 'index.html')
    else:
        global ans
        ans = ''
        content = json.loads(request.body.decode('utf-8'))

        func = content['func']
        repl = re.findall(r'_\d\d', func)
        for el in repl:
            func = re.sub(el, el[1] + '**' + el[2], func)
        func = re.sub('_', '', func)
        func = re.sub('\^', '**', func)
        func = parse_expr(func)

        init_vector = content['init_vector']
        init_vector = Matrix([float(v) for v in filter(None, re.split("[, ]+", init_vector))])

        NUMBER_OF_VARS = int(content['number_of_vars'])
        x_vars = symbols('x1:%d' % (NUMBER_OF_VARS + 1))

        if len(init_vector) != NUMBER_OF_VARS:
            return HttpResponse('err', content_type='application/json')

        EPS = float(content['eps'])
        s = Symbol('s')

        answer = conjdirmethod(1)

        res = 'Исходная функция: $$ f(x_{}) = '.format('x_'.join([str(i + 1) for i in range(NUMBER_OF_VARS)])) + latex(
            func) + '\\rightarrow min $$' + 'Вектор начального приближения: $$ x_0 = ' + latex(
            init_vector) + ' $$' + ans + '<b>Ответ:</b> $$ x^* = ' + latex(answer) + ' $$'

    return HttpResponse(res, content_type='application/json')
