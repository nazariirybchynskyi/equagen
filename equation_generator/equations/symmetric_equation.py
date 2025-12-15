import random
import sympy
from sympy import sin, cos, Eq, solveset, Reals, symbols, pi, sqrt, Rational, S, Union
from ..base_class import TrigonometricEquation


class SymmetricEquation(TrigonometricEquation):

    def _format_polynomial_latex(self, A, B, C, var_name):
        def fmt_coeff(val):
            if val == 1: return ""
            if val == -1: return "-"
            lat = sympy.latex(abs(val))
            if getattr(val, 'is_Add', False):
                return f"({lat})"
            return lat

        terms = []

        if A != 0:
            coeff = fmt_coeff(A)
            if A == 1:
                terms.append(f"{var_name}^2")
            elif A == -1:
                terms.append(f"-{var_name}^2")
            else:
                if A > 0:
                    terms.append(f"{sympy.latex(A)} {var_name}^2")
                else:
                    terms.append(f"-{sympy.latex(abs(A))} {var_name}^2")

        if B != 0:
            sign = "+" if B > 0 else "-"
            coeff = fmt_coeff(B)
            if coeff in ["", "-"]:
                terms.append(f"{sign} {var_name}")
            else:
                terms.append(f"{sign} {coeff}{var_name}")

        if C != 0:
            sign = "+" if C > 0 else "-"
            terms.append(f"{sign} {sympy.latex(abs(C))}")

        if not terms:
            return "0 = 0"

        latex_str = " ".join(terms)
        if not latex_str.startswith('-'):
            latex_str = latex_str.lstrip('+ ')

        return latex_str + " = 0"

    def _generate(self):
        while True:
            nice_t_roots = [
                0, 1, -1,
                sqrt(2), -sqrt(2),
                sqrt(2) / 2, -sqrt(2) / 2,
                1 / sqrt(2), -1 / sqrt(2)
            ]

            t1 = random.choice(nice_t_roots)
            t2 = random.choice(nice_t_roots)

            if abs(t1) > sqrt(2) and abs(t2) > sqrt(2):
                continue

            sub_type = random.choice(['plus', 'minus'])

            a_quad = random.choice([1, -1, 2])
            b_quad = -a_quad * (t1 + t2)
            c_quad = a_quad * t1 * t2

            coeffs = [a_quad, b_quad, c_quad]
            denoms = [c.q for c in coeffs if hasattr(c, 'q')]
            if denoms:
                lcm = sympy.lcm(denoms)
                a_quad *= lcm
                b_quad *= lcm
                c_quad *= lcm

            a_quad = sympy.simplify(a_quad)
            b_quad = sympy.simplify(b_quad)
            c_quad = sympy.simplify(c_quad)

            if a_quad == 0:
                continue

            if sub_type == 'plus':
                A_eq = a_quad
                B_eq = b_quad
                C_eq = c_quad + a_quad
                term_linear = sin(self.x) + cos(self.x)
            else:
                A_eq = -a_quad
                B_eq = b_quad
                C_eq = c_quad + a_quad
                term_linear = sin(self.x) - cos(self.x)

            if A_eq == 0 and B_eq == 0:
                continue

            self.equation_obj = Eq(A_eq * sin(2 * self.x) + B_eq * term_linear + C_eq, 0)

            self.variables = {
                'sub_type': sub_type,
                't1': t1, 't2': t2,
                'a_quad': a_quad, 'b_quad': b_quad, 'c_quad': c_quad,
                'A_eq': A_eq, 'B_eq': B_eq, 'C_eq': C_eq
            }
            break

    def _solve(self):
        t1 = self.variables['t1']
        t2 = self.variables['t2']
        sub_type = self.variables['sub_type']

        if sub_type == 'plus':
            expr = sin(self.x) + cos(self.x)
        else:
            expr = sin(self.x) - cos(self.x)

        sol1 = solveset(Eq(expr, t1), self.x, domain=Reals)
        sol2 = solveset(Eq(expr, t2), self.x, domain=Reals)

        self.solution_obj = Union(sol1, sol2)

    def _build_solution_steps(self):
        sub_type = self.variables['sub_type']
        t1 = self.variables['t1']
        t2 = self.variables['t2']
        a_quad = self.variables['a_quad']
        b_quad = self.variables['b_quad']
        c_quad = self.variables['c_quad']

        t = symbols('t')

        self.steps.append(("text", f"Маємо рівняння: ${sympy.latex(self.equation_obj)}$"))

        sign_str = "+" if sub_type == 'plus' else "-"

        self.steps.append(("text",
                           rf"Це симетричне рівняння відносно $\sin x$ та $\cos x$. Введемо заміну $t = \sin x {sign_str} \cos x$."))

        if sub_type == 'plus':
            sq_step = r"t^2 = (\sin x + \cos x)^2 = \sin^2 x + 2\sin x \cos x + \cos^2 x = 1 + \sin 2x"
            sin2x_expr = r"\sin 2x = t^2 - 1"
        else:
            sq_step = r"t^2 = (\sin x - \cos x)^2 = \sin^2 x - 2\sin x \cos x + \cos^2 x = 1 - \sin 2x"
            sin2x_expr = r"\sin 2x = 1 - t^2"

        self.steps.append(("text", "Піднесемо заміну до квадрату:"))
        self.steps.append(("math", sq_step))
        self.steps.append(("math", rf"\implies {sin2x_expr}"))

        self.steps.append(("text", rf"Підставимо вираз для $\sin 2x$ та $t$ у початкове рівняння:"))

        A_eq = self.variables['A_eq']
        B_eq = self.variables['B_eq']
        C_eq = self.variables['C_eq']

        term_sin2x = f"(t^2 - 1)" if sub_type == 'plus' else f"(1 - t^2)"

        sub_eq_str = ""
        if A_eq != 0:
            A_lat = sympy.latex(A_eq) if A_eq != 1 and A_eq != -1 else ("-" if A_eq == -1 else "")
            sub_eq_str += f"{A_lat}{term_sin2x}"

        if B_eq != 0:
            sign = "+" if (B_eq > 0 and sub_eq_str) else "-"
            B_val = sympy.latex(abs(B_eq)) if abs(B_eq) != 1 else ""
            sub_eq_str += f" {sign} {B_val}t"

        if C_eq != 0:
            sign = "+" if (C_eq > 0 and sub_eq_str) else "-"
            sub_eq_str += f" {sign} {sympy.latex(abs(C_eq))}"

        sub_eq_str += " = 0"
        self.steps.append(("math", sub_eq_str))

        self.steps.append(("text", "Розкриємо дужки та отримаємо квадратне рівняння:"))
        quad_latex = self._format_polynomial_latex(a_quad, b_quad, c_quad, "t")
        self.steps.append(("math", quad_latex))

        self.steps.append(("text", "Знаходимо корені:"))

        if t1 == t2:
            self.steps.append(("math", f"t = {sympy.latex(t1)}"))
        else:
            self.steps.append(("math", f"t_1 = {sympy.latex(t1)}, \\quad t_2 = {sympy.latex(t2)}"))

        self.steps.append(("text",
                           rf"Зауважимо, що $|t| = |\sin x {sign_str} \cos x| = |\sqrt{{2}}\sin(x \pm \frac{{\pi}}{{4}})| \le \sqrt{{2}} \approx 1.41$."))

        self.steps.append(("text", "Повертаємось до заміни:"))

        unique_roots = [t1]
        if t1 != t2:
            unique_roots.append(t2)

        solution_sets_latex = []

        for i, root in enumerate(unique_roots, 1):
            prefix = f"{i}) " if len(unique_roots) > 1 else ""

            if abs(root) > sqrt(2):
                self.steps.append(("math",
                                   rf"{prefix}\sin x {sign_str} \cos x = {sympy.latex(root)} \implies \text{{розв'язків немає, бо }} |{sympy.latex(root)}| > \sqrt{{2}}"))
            else:
                if len(unique_roots) > 1:
                    self.steps.append(("text", rf"{i}) Розв'язуємо рівняння:"))
                else:
                    self.steps.append(("text", rf"Розв'язуємо рівняння:"))

                self.steps.append(("math", rf"\sin x {sign_str} \cos x = {sympy.latex(root)}"))

                self.steps.append(("text", r"Поділимо на $\sqrt{2}$ та зведемо до синуса суми/різниці:"))
                val_div = root / sqrt(2)
                sign_pi = "+" if sub_type == "plus" else "-"

                self.steps.append(("math",
                                   rf"\frac{{1}}{{\sqrt{{2}}}}\sin x {sign_str} \frac{{1}}{{\sqrt{{2}}}}\cos x = \frac{{{sympy.latex(root)}}}{{\sqrt{{2}}}}"))
                self.steps.append(("math",
                                   rf"\sin x \cos\frac{{\pi}}{{4}} {sign_str} \cos x \sin\frac{{\pi}}{{4}} = {sympy.latex(val_div)}"))
                self.steps.append(("math", rf"\sin(x {sign_pi} \frac{{\pi}}{{4}}) = {sympy.latex(val_div)}"))

                try:
                    angle = sympy.asin(val_div)
                    if isinstance(angle, sympy.asin): raise ValueError

                    is_neg = False
                    if angle < 0:
                        is_neg = True
                        angle = -angle

                    angle_latex = sympy.latex(angle)

                    power_n = "n+1" if is_neg else "n"

                    self.steps.append(("text",
                                       rf"Оскільки $\arcsin({sympy.latex(val_div)}) = {('-' if is_neg else '')}{angle_latex}$, маємо:"))
                    self.steps.append(("math",
                                       rf"x {sign_pi} \frac{{\pi}}{{4}} = (-1)^{{{power_n}}} \cdot {angle_latex} + \pi n, n \in \mathbb{{Z}}"))

                    final_x_expr = f"(-1)^{{{power_n}}} {angle_latex}"
                    if sign_pi == "+":
                        final_x_expr += r" - \frac{\pi}{4} + \pi n"
                    else:
                        final_x_expr += r" + \frac{\pi}{4} + \pi n"

                    self.steps.append(("math", f"x = {final_x_expr}"))

                    solution_sets_latex.append(rf"\left\{{ {final_x_expr} \; \middle| \; n \in \mathbb{{Z}} \right\}}")

                except:
                    sol = solveset(Eq(sin(self.x + (pi / 4 if sub_type == 'plus' else -pi / 4)), val_div), self.x,
                                   domain=Reals)
                    lat = sympy.latex(sol)
                    self.steps.append(("math", rf"x \in {lat}"))
                    solution_sets_latex.append(lat)

        self.steps.append(("text", f"Кінцева відповідь:"))
        if not solution_sets_latex:
            self.steps.append(("text", "Розв'язків немає."))
        else:
            full_latex = r" \cup ".join(solution_sets_latex)
            self.steps.append(("math", rf"x \in {full_latex}"))