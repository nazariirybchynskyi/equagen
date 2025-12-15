import random
import sympy
from sympy import sin, cos, tan, Eq, solveset, Reals, symbols, Rational, sqrt, Union, S

from ..base_class import TrigonometricEquation


class TanSubstitutionEquation(TrigonometricEquation):
    def _format_polynomial_latex(self, coeffs, var_name):
        degree = len(coeffs) - 1
        terms = []

        for i, c in enumerate(coeffs):
            power = degree - i
            if c == 0: continue

            if c == 1 and power != 0:
                str_c = ""
            elif c == -1 and power != 0:
                str_c = "-"
            else:
                str_c = sympy.latex(c)

            if power == 0:
                str_var = ""
            elif power == 1:
                str_var = var_name
            else:
                str_var = f"{var_name}^{power}"

            term = f"{str_c}{str_var}"

            if not terms:
                terms.append(term)
            else:
                if c > 0:
                    terms.append(f"+ {term}")
                else:
                    terms.append(f"- {abs(c)}{str_var}")

        if not terms: return "0"
        return " ".join(terms)

    def _generate(self):
        while True:
            target_func = random.choice([sin, cos])
            nice_roots = [
                0,
                1, -1,
                sqrt(3), -sqrt(3),
                sqrt(3) / 3, -sqrt(3) / 3
            ]
            t1 = random.choice(nice_roots)

            A = random.choice([1, -1, 2, -2, 3, -3])
            B = random.choice([1, -1, 2, -2, 3, -3])
            denom = 1 + t1 ** 2

            if target_func == sin:
                val_func = 2 * t1 / denom
            else:
                val_func = (1 - t1 ** 2) / denom

            C = A * val_func + B * t1

            C = sympy.sympify(C)
            C = sympy.simplify(C)

            if not (C.is_integer or (hasattr(C, 'q') and C.q < 10) or C.has(sqrt)):
                continue

            self.equation_obj = Eq(A * target_func(2 * self.x) + B * tan(self.x), C)

            if target_func == sin:
                poly_coeffs = [B, -C, 2 * A + B, -C]
            else:
                poly_coeffs = [B, -A - C, B, A - C]

            if poly_coeffs[0] == 0:
                continue

            self.variables = {
                'target_func': target_func,
                'A': A, 'B': B, 'C': C,
                't1': t1,
                'poly_coeffs': poly_coeffs
            }
            break

    def _solve(self):
        t = symbols('t')
        poly_coeffs = self.variables['poly_coeffs']

        poly = 0
        degree = len(poly_coeffs) - 1
        for i, c in enumerate(poly_coeffs):
            poly += c * t ** (degree - i)

        t_roots = sympy.solve(poly, t)

        final_solution = sympy.EmptySet
        for root in t_roots:
            s_root = sympy.sympify(root)
            if s_root.is_real:
                sol = solveset(Eq(tan(self.x), s_root), self.x, domain=Reals)
                final_solution = Union(final_solution, sol)

        self.solution_obj = final_solution

    def _build_solution_steps(self):
        A = self.variables['A']
        B = self.variables['B']
        C = self.variables['C']
        t1 = self.variables['t1']
        poly_coeffs = self.variables['poly_coeffs']
        target_func = self.variables['target_func']
        func_name = target_func.__name__
        t = symbols('t')

        self.steps.append(("text", f"Маємо рівняння: ${sympy.latex(self.equation_obj)}$"))

        self.steps.append(("text",
                           rf"Рівняння містить ${func_name}(2x)$ та $\text{{tg}}(x)$. Виразимо ${func_name}(2x)$ через $\text{{tg}}(x)$ за формулою:"))

        if target_func == sin:
            formula = r"\sin(2x) = \frac{2\text{tg}(x)}{1+\text{tg}^2(x)}"
            sub_expr = r"\frac{2t}{1+t^2}"
        else:
            formula = r"\cos(2x) = \frac{1-\text{tg}^2(x)}{1+\text{tg}^2(x)}"
            sub_expr = r"\frac{1-t^2}{1+t^2}"

        self.steps.append(("math", formula))

        self.steps.append(("text", r"Нехай $t = \text{tg}(x)$. Підставимо у рівняння:"))

        def fmt(val):
            if val == 1: return ""
            if val == -1: return "-"
            return sympy.latex(val)

        eq_sub = rf"{fmt(A)} \cdot {sub_expr}"

        if B > 0:
            eq_sub += f" + {fmt(B)}t"
        elif B < 0:
            eq_sub += f" - {sympy.latex(abs(B))}t"

        eq_sub += f" = {sympy.latex(C)}"
        self.steps.append(("math", eq_sub))

        self.steps.append(
            ("text", r"Помножимо обидві частини на $(1+t^2)$ (оскільки $1+\text{tg}^2 x \neq 0$) та розкриємо дужки:"))

        if target_func == sin:
            step_inter = rf"{sympy.latex(2 * A)}t + {sympy.latex(B)}t(1+t^2) = {sympy.latex(C)}(1+t^2)"
        else:
            term_A = f"{sympy.latex(A)}(1-t^2)"
            B_sign = "+" if B > 0 else "-" if B < 0 else ""
            term_B = f"{B_sign} {sympy.latex(abs(B))}t(1+t^2)" if B != 0 else ""

            step_inter = rf"{term_A} {term_B} = {sympy.latex(C)}(1+t^2)"

        self.steps.append(("math", step_inter))

        poly_latex = self._format_polynomial_latex(poly_coeffs, "t") + " = 0"
        self.steps.append(
            ("text", "Перенесемо все в одну сторону і зведемо подібні доданки. Отримуємо кубічне рівняння:"))
        self.steps.append(("math", poly_latex))

        t = symbols('t')
        poly_sym = 0
        deg = len(poly_coeffs) - 1
        for i, c in enumerate(poly_coeffs): poly_sym += c * t ** (deg - i)
        all_t_roots = sympy.solve(poly_sym, t)

        if poly_coeffs[-1] == 0:
            self.steps.append(("text", "Вільний член дорівнює 0, тому винесемо спільний множник $t$ за дужки:"))

            coeffs_quad = poly_coeffs[:-1]
            quad_latex = self._format_polynomial_latex(coeffs_quad, "t").replace("= 0", "")
            self.steps.append(("math", rf"t({quad_latex}) = 0"))

            self.steps.append(("text", "Звідси перший корінь $t_1 = 0$."))

            other_roots = [r for r in all_t_roots if r != 0]

            if other_roots:
                self.steps.append(("text", "Розв'яжемо квадратне рівняння в дужках:"))
                roots_display = ", ".join([f"t = {sympy.latex(r)}" for r in other_roots])
                self.steps.append(("math", rf"{quad_latex} = 0 \implies {roots_display}"))

        else:
            self.steps.append(("text",
                               f"Спробуємо підібрати корінь серед дільників вільного члена ({sympy.latex(poly_coeffs[-1])})."))
            self.steps.append(("text", rf"Перевіримо $t = {sympy.latex(t1)}$:"))
            self.steps.append(
                ("math", rf"P({sympy.latex(t1)}) = 0 \implies t_1 = {sympy.latex(t1)} \text{{ є коренем.}}"))

            self.steps.append(
                ("text", rf"Розділимо многочлен на $(t - ({sympy.latex(t1)}))$ і отримаємо квадратне рівняння:"))

            quotient, remainder = sympy.div(poly_sym, t - t1)
            self.steps.append(("math", sympy.latex(Eq(quotient, 0))))

            other_roots = sympy.solve(quotient, t)
            if other_roots:
                self.steps.append(("text", "Корені цього квадратного рівняння:"))
                roots_display = ", ".join([f"t = {sympy.latex(r)}" for r in other_roots])
                self.steps.append(("math", roots_display))
            else:
                self.steps.append(("text", "Це квадратне рівняння не має дійсних коренів."))

        self.steps.append(("text", "Повертаємось до заміни $\text{tg}(x) = t$:"))

        final_sets_latex = []

        unique_roots = []
        for r in all_t_roots:
            s_r = sympy.sympify(r)
            if s_r.is_real:
                if s_r not in unique_roots:
                    unique_roots.append(s_r)

        unique_roots.sort(key=lambda x: float(x))

        for i, root in enumerate(unique_roots, 1):
            angle = sympy.atan(root)

            if isinstance(angle, sympy.atan):
                angle_latex = rf"\text{{arctg}}({sympy.latex(root)})"
            else:
                angle_latex = sympy.latex(angle)

            if angle_latex == "0":
                x_expr = r"\pi n"
            else:
                x_expr = rf"{angle_latex} + \pi n"

            prefix = f"{i}) " if len(unique_roots) > 1 else ""

            self.steps.append(
                ("math", rf"{prefix}\text{{tg}}(x) = {sympy.latex(root)} \implies x = {x_expr}, n \in \mathbb{{Z}}"))

            final_sets_latex.append(rf"\left\{{ {x_expr} \; \middle| \; n \in \mathbb{{Z}} \right\}}")

        self.steps.append(("text", "Кінцева відповідь:"))
        if not final_sets_latex:
            self.steps.append(("text", "Розв'язків немає."))
        else:
            full_response = r" \cup ".join(final_sets_latex)
            self.steps.append(("math", rf"x \in {full_response}"))