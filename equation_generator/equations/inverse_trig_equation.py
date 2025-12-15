import random
import sympy
from sympy import (sin, cos, tan, cot, asin, acos, atan, acot,
                   Eq, solveset, Reals, symbols, expand,
                   pi, sqrt, Rational, S, Union)

from ..base_class import TrigonometricEquation


class InverseTrigEquation(TrigonometricEquation):

    def _generate(self):
        while True:
            func_type = random.choice(['arcsin', 'arccos', 'arctg', 'arcctg'])

            if func_type in ['arcsin', 'arccos']:
                V_choices = [0, 1, -1, Rational(1, 2), Rational(-1, 2)]
            else:
                V_choices = [0, 1, -1]

            V = random.choice(V_choices)

            if func_type == 'arcsin':
                alpha = asin(V)
                sympy_func = asin
            elif func_type == 'arccos':
                alpha = acos(V)
                sympy_func = acos
            elif func_type == 'arctg':
                alpha = atan(V)
                sympy_func = atan
            elif func_type == 'arcctg':
                alpha = acot(V)
                if V < 0:
                    alpha = alpha + pi
                sympy_func = acot

            k = random.choice([1, 2, 3, 4])

            rhs = k * alpha

            poly_degree = random.choice([1, 2])

            if poly_degree == 1:
                x1 = random.randint(-5, 5)
                a = random.choice([1, -1, 2, -2, 3])
                b = V - a * x1
                P = a * self.x + b

            else:
                x1 = random.randint(-5, 5)
                x2 = random.randint(-5, 5)
                a = random.choice([1, -1, 2])
                P = expand(a * (self.x - x1) * (self.x - x2)) + V

            P = sympy.simplify(P)

            self.equation_obj = Eq(k * sympy_func(P), rhs)

            self.variables = {
                'func_type': func_type,
                'V': V,
                'alpha': alpha,
                'k': k,
                'rhs': rhs,
                'P': P,
                'poly_degree': poly_degree
            }
            break

    def _solve(self):
        P = self.variables['P']
        V = self.variables['V']

        sol = solveset(Eq(P, V), self.x, domain=Reals)
        self.solution_obj = sol

    def _build_solution_steps(self):
        func_type = self.variables['func_type']
        V = self.variables['V']
        alpha = self.variables['alpha']
        k = self.variables['k']
        rhs = self.variables['rhs']
        P = self.variables['P']

        if func_type == 'arctg':
            func_latex = r"\text{arctg}"
        elif func_type == 'arcctg':
            func_latex = r"\text{arcctg}"
        else:
            func_latex = rf"\{func_type}"

        if func_type == 'arcsin':
            direct_func = r"\sin"
        elif func_type == 'arccos':
            direct_func = r"\cos"
        elif func_type == 'arctg':
            direct_func = r"\text{tg}"
        elif func_type == 'arcctg':
            direct_func = r"\text{ctg}"

        P_latex = sympy.latex(P)

        self.steps.append(("text", f"Маємо рівняння: ${sympy.latex(self.equation_obj)}$"))

        current_rhs = rhs
        if k != 1:
            self.steps.append(("text", rf"Поділимо обидві частини рівняння на ${k}$:"))
            current_rhs = rhs / k
            self.steps.append(("math", rf"{func_latex}({P_latex}) = {sympy.latex(current_rhs)}"))

        range_info = ""
        is_valid = True
        if func_type == 'arcsin':
            range_info = r"\left[-\frac{\pi}{2}; \frac{\pi}{2}\right]"
            if not (-pi / 2 <= current_rhs <= pi / 2): is_valid = False
        elif func_type == 'arccos':
            range_info = r"[0; \pi]"
            if not (0 <= current_rhs <= pi): is_valid = False
        elif func_type == 'arctg':
            range_info = r"\left(-\frac{\pi}{2}; \frac{\pi}{2}\right)"
            if not (-pi / 2 < current_rhs < pi / 2): is_valid = False
        elif func_type == 'arcctg':
            range_info = r"(0; \pi)"
            if not (0 < current_rhs < pi): is_valid = False

        self.steps.append(("text",
                           rf"Перевіримо, чи належить права частина області значень функції ${func_latex}$, тобто проміжку ${range_info}$."))

        if not is_valid:
            self.steps.append(("text", r"Значення виходить за межі області значень. Рівняння розв'язків немає."))
            self.solution_obj = S.EmptySet
            return

        self.steps.append(("text", r"Значення належить області визначення."))

        self.steps.append(("text", r"За означенням обернених тригонометричних функцій:"))
        self.steps.append(("math", rf"{P_latex} = {direct_func}\left({sympy.latex(current_rhs)}\right)"))

        self.steps.append(("text",
                           rf"Оскільки ${direct_func}\left({sympy.latex(current_rhs)}\right) = {sympy.latex(V)}$, отримуємо алгебраїчне рівняння:"))

        alg_eq = Eq(P, V)
        self.steps.append(("math", sympy.latex(alg_eq)))

        if V != 0:
            lhs_final = P - V
            self.steps.append(("math", rf"{sympy.latex(lhs_final)} = 0"))

        self.steps.append(("text", "Розв'язуємо отримане рівняння:"))

        final_sol = sympy.latex(self.solution_obj)
        self.steps.append(("math", f"x \in {final_sol}"))

        if func_type in ['arcsin', 'arccos']:
            self.steps.append(("text", r"Перевірка: знайдені корені повинні задовольняти умову $|P(x)| \le 1$."))
            self.steps.append(
                ("text", rf"Оскільки $P(x) = {sympy.latex(V)}$ і $|{sympy.latex(V)}| \le 1$, корені підходять."))

        self.steps.append(("text", "Кінцева відповідь:"))
        self.steps.append(("math", f"x \in {final_sol}"))