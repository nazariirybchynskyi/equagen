import random
import sympy
from sympy import sin, cos, Eq, solveset, Reals, symbols, pi, S, Intersection, Union

from ..base_class import TrigonometricEquation


class BoundedSumEquation(TrigonometricEquation):

    def _generate(self):
        while True:
            num_terms = random.choice([2, 2, 3])
            golden_roots = [0, pi, pi / 2, 3 * pi / 2]
            x0 = random.choice(golden_roots)

            terms_data = []

            for _ in range(num_terms):
                for attempt in range(50):
                    func = random.choice([sin, cos])
                    k = random.choice([1, 1, 2, 3, 4, 5, 6])

                    val = func(k * x0)
                    val = sympy.simplify(val)

                    if abs(val) == 1:
                        sign = int(val)
                        terms_data.append({
                            'func': func,
                            'k': k,
                            'sign': sign
                        })
                        break
                else:
                    break

            if len(terms_data) < num_terms:
                continue

            signatures = set()
            is_duplicate = False
            for t in terms_data:
                sig = (t['func'].__name__, t['k'])
                if sig in signatures:
                    is_duplicate = True
                    break
                signatures.add(sig)

            if is_duplicate:
                continue

            lhs_expr = 0
            for t in terms_data:
                term_expr = t['sign'] * t['func'](t['k'] * self.x)
                lhs_expr += term_expr

            rhs_value = num_terms

            self.equation_obj = Eq(lhs_expr, rhs_value)

            self.variables = {
                'x0': x0,
                'terms_data': terms_data,
                'rhs_value': rhs_value
            }
            break

    def _solve(self):
        x0 = self.variables['x0']
        n = symbols('n', integer=True)
        self.solution_obj = sympy.ImageSet(sympy.Lambda(n, x0 + 2 * pi * n), sympy.Integers)

    def _build_solution_steps(self):
        terms_data = self.variables['terms_data']
        rhs_value = self.variables['rhs_value']
        x0 = self.variables['x0']

        self.steps.append(("text", f"Маємо рівняння: ${sympy.latex(self.equation_obj)}$"))

        lhs_latex_parts = []
        ineq_parts = []

        for t in terms_data:
            func_name = t['func'].__name__
            k = t['k']
            k_str = f"{k}x" if k > 1 else "x"
            sign_str = "" if t['sign'] == 1 else "-"

            term_str = rf"{sign_str}\{func_name}({k_str})"
            lhs_latex_parts.append(term_str)

            ineq_parts.append(rf"{term_str} \le 1")

        self.steps.append(
            ("text", r"Оцінимо ліву частину рівняння. Відомо, що область значень синуса і косинуса $[-1, 1]$."))

        inequalities = r", \quad ".join(ineq_parts)
        self.steps.append(("math", inequalities))

        sum_ineq = " + ".join(["1"] * len(terms_data))
        self.steps.append(("text", f"Отже, сума цих доданків не перевищує: ${sum_ineq} = {rhs_value}$."))

        self.steps.append(
            ("text", f"Рівність досягається тоді й лише тоді, коли кожен доданок досягає свого максимуму ($1$)."))
        self.steps.append(("text", "Рівняння рівносильне системі:"))

        system_lines = []
        for t in terms_data:
            func_name = t['func'].__name__
            k = t['k']
            k_str = f"{k}x" if k > 1 else "x"
            sign_str = "" if t['sign'] == 1 else "-"
            system_lines.append(rf"{sign_str}\{func_name}({k_str}) = 1")

        latex_system = r"\begin{cases} " + r" \\ ".join(system_lines) + r" \end{cases}"
        self.steps.append(("math", latex_system))

        self.steps.append(("text", "Розв'яжемо перше рівняння системи:"))

        first_term = terms_data[0]
        func_name = first_term['func'].__name__
        k = first_term['k']
        sign = first_term['sign']

        rhs_simple = 1 * sign
        k_str = f"{k}x" if k > 1 else "x"

        self.steps.append(("math", rf"\{func_name}({k_str}) = {rhs_simple}"))

        if x0 == 0:
            x0_latex = ""
            plus_sign = ""
        else:
            x0_latex = sympy.latex(x0)
            plus_sign = " + "

        if k == 1:
            root_series = rf"x = {x0_latex}{plus_sign}2\pi n, n \in \mathbb{{Z}}"
        else:
            base_angle = sympy.asin(rhs_simple) if func_name == 'sin' else sympy.acos(rhs_simple)
            base_latex = sympy.latex(base_angle)

            self.steps.append(("math", rf"{k}x = {base_latex} + 2\pi n"))

            term_1 = base_angle / k
            term_2_num = 2
            term_2_den = k

            term_1_lat = sympy.latex(term_1) if term_1 != 0 else ""
            plus = " + " if term_1 != 0 else ""
            term_2_lat = rf"\frac{{{term_2_num}\pi n}}{{{term_2_den}}}"

            root_series = rf"x = {term_1_lat}{plus}{term_2_lat}, n \in \mathbb{{Z}}"

        self.steps.append(("math", root_series))
        self.steps.append(("text", f"Перевіримо, чи задовольняють ці корені інші рівняння системи."))

        self.steps.append(("text",
                           rf"Підставимо значення $x$, які відповідають $n$ кратному {k} (щоб період $2\pi n$ зберігся), наприклад, серію $x = {sympy.latex(x0)} + 2\pi m$:"))

        for i in range(1, len(terms_data)):
            t = terms_data[i]
            func_name = t['func'].__name__
            tk = t['k']
            tsign = t['sign']
            sign_sym = "-" if tsign == -1 else ""

            arg_sub = rf"{tk}\left({sympy.latex(x0)} + 2\pi m\right)"
            arg_open = rf"{sympy.latex(tk * x0)} + {2 * tk}\pi m"

            val_x0 = t['func'](tk * x0)
            res = tsign * val_x0

            check_line = rf"{sign_sym}\{func_name}\left( {arg_sub} \right) = {sign_sym}\{func_name}\left( {arg_open} \right) = {sign_sym}\{func_name}\left( {sympy.latex(tk * x0)} \right) = {sign_sym}({sympy.latex(val_x0)}) = {res}"

            self.steps.append(("math", check_line))

        self.steps.append(("text", "Усі рівняння системи задовольняються."))

        if x0 == 0:
            ans_latex = r"2\pi n"
        else:
            ans_latex = rf"{sympy.latex(x0)} + 2\pi n"

        self.steps.append(("text", "Кінцева відповідь:"))
        self.steps.append(("math", rf"x = {ans_latex}, n \in \mathbb{{Z}}"))