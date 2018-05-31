import json
import math
import re
from functools import reduce

from django.http import HttpResponse, JsonResponse
from django.shortcuts import render
from sympy import symbols, latex, Symbol, solve
from sympy.matrices import Matrix, hessian
from sympy.parsing.sympy_parser import parse_expr
from sympy.abc import u, v, a, b

from django.views.decorators.csrf import csrf_exempt

ans = ''
ROUNDING_NUMBER = 4
MAX_ITER = 20

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
            d = [(v[0].subs({x_vars[k]: x_0[k] for k in range(NUMBER_OF_VARS)})).evalf(ROUNDING_NUMBER)]
            ans += '$$ d_0 = v_1 = ' + latex(d[0]) + '$$'

            x = x_0 + s * d[0]
            ans += '$$ x_1 = argmin f(x_0 + s \cdot d_0) = argmin f(' + latex(x_0) + '+ s \cdot' + latex(
                d[0]) + ') = argmin f(' + latex(x) + ') $$'
            S = solve(f(list(x)))
            x = (x.subs(s, S)).evalf(ROUNDING_NUMBER)
            ans += '$$ x_1 = ' + latex(x) + '$$'

            for j in range(1, NUMBER_OF_VARS):
                J = str(j)
                ans += 'Шаг ' + str(j + 1) + ':<br>'
                y0 = x + (v[j].subs({x_vars[k]: x_0[k] for k in range(NUMBER_OF_VARS)})).evalf(ROUNDING_NUMBER)
                ans += '$$ y_' + str(j - 1) + '= x_' + J + ' + v_' + str(j + 1) + '=' + latex(x) + '+' + latex(
                    y0 - x) + '=' + latex(y0) + '$$'

                ans += 'Находим \(y_' + J + '\), минимизируя последовательно функцию по направлениям \(d_{}'.format(
                    ',d_'.join([str(k) for k in range(j)])) + '\):<br>'
                y = y0
                for k in range(0, j):
                    y = y + s * d[k]
                    S = solve(f(list(y)))
                    y = (y.subs(s, S)).evalf(ROUNDING_NUMBER)
                ans += '$$ y_' + J + '=' + 'y_' + str(j - 1) + '+ \sum_{i=0}^' + str(j - 1) + 's_i \cdot d_i=' + latex(
                    y) + '$$'

                d.append(x - y)
                ans += '$$ d_' + J + '= x_' + J + '- y_' + J + '=' + latex(x) + '-' + latex(y) + '=' + latex(
                    d[j]) + '$$'
                ans += '$$ x_' + str(j + 1) + '= argmin f(x_' + J + '+ s \cdot d_' + J + ') = argmin f(' + latex(
                    x) + '+ s \cdot' + latex(d[j]) + ') = argmin f('
                x = x + s * d[j]
                ans += latex(x) + ')$$'

                S = solve(f(list(x)))
                x = (x.subs(s, S)).evalf(ROUNDING_NUMBER)
                ans += '$$ x_' + str(j + 1) + '=' + latex(x) + ' $$'

            if (x - x_0).norm() <= EPS:
                break

            if i > MAX_ITER:
                break

            x_0 = x
        ans += '<b>Ответ:</b> $$ x^* = ' + latex(x) + ' $$'
        return 0

    if request.method == 'GET':
        return render(request, 'index.html')
    else:
        global ans
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

        ans = 'Исходная функция: $$f(x_{}) = '.format('x_'.join([str(i + 1) for i in range(NUMBER_OF_VARS)])) + latex(
            func) + '\\rightarrow min $$' + 'Вектор начального приближения: $$x_0 = ' + latex(
            init_vector) + '$$'

        EPS = float(content['eps'])
        s = Symbol('s')

        conjdirmethod(1)

        res = ans

    return HttpResponse(res, content_type='application/json')


@csrf_exempt
def surface(request):
    if request.method == "POST":
        content = json.loads(request.body.decode('utf-8'))
        Xstr = content['X']
        Ystr = content['Y']
        Zstr = content['Z']

        X1str = content['X1']
        Y1str = content['Y1']
        Z1str = content['Z1']

        X = parse_expr(Xstr)
        Y = parse_expr(Ystr)
        Z = parse_expr(Zstr)
        X1 = parse_expr(X1str)
        Y1 = parse_expr(Y1str)
        Z1 = parse_expr(Z1str)

        solution = solve([X - X1, Y - Y1, Z - Z1], [u, v, a])
        answer = {"curve" + str(i): [str(X.subs({u: solution[i][0], v: solution[i][1]})),
                                     str(Y.subs({u: solution[i][0], v: solution[i][1]})),
                                     str(Z.subs({u: solution[i][0], v: solution[i][1]}))]
                  for i in range(len(solution))}
        resp = JsonResponse(answer)

        resp['Access-Control-Allow-Origin'] = '*'
        return resp

    return HttpResponse('It was GET')
