import random
import sympy
from sympy import sin, cos, tan, cot, Eq, solveset, Reals, symbols, Rational, sqrt, Union, S

from ..base_class import TrigonometricEquation


class SumTanCotanEquation(TrigonometricEquation):

    def _format_polynomial_latex(self, A, B, C, var_name):
        def fmt_coeff(val):
            if val == 1: return ""
            if val == -1: return "-"
            lat = sympy.latex(abs(val))
            if getattr(val, 'is_Add', False): return f"({lat})"
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
                terms.append(f"{sign} {coeff} {var_name}")

        if C != 0:
            sign = "+" if C > 0 else "-"
            terms.append(f"{sign} {sympy.latex(abs(C))}")

        if not terms: return "0 = 0"
        latex_str = " ".join(terms)
        if not latex_str.startswith('-'): latex_str = latex_str.lstrip('+ ')
        return latex_str + " = 0"

    def _generate(self):
        while True:
            valid_roots = [2, -2, 4 / sqrt(3), -4 / sqrt(3)]

            dummy_roots = [0, 1, -1]

            if random.random() < 0.5:
                t1 = random.choice(valid_roots)
                t2 = random.choice(valid_roots)
            else:
                t1 = random.choice(valid_roots)
                t2 = random.choice(dummy_roots)

            A_sq = 1
            B_sq = -(t1 + t2)
            C_sq = t1 * t2

            coeffs = [A_sq, B_sq, C_sq]
            exprs = [A_sq, B_sq, C_sq]
            mult = S(1)
            for e in exprs:
                e = sympy.simplify(e)
                if e.has(sqrt(3)):
                    numer, denom = e.as_numer_denom()
                    if denom.has(sqrt(3)):
                        mult = sympy.lcm(mult, denom)
                if hasattr(e, 'q'):
                    mult = sympy.lcm(mult, e.q)

            A_sq = sympy.simplify(A_sq * mult)
            B_sq = sympy.simplify(B_sq * mult)
            C_sq = sympy.simplify(C_sq * mult)

            A_eq = A_sq
            B_eq = B_sq
            C_eq = C_sq + 2 * A_sq

            if A_eq == 0: continue

            self.equation_obj = Eq(
                A_eq * (tan(self.x) ** 2 + cot(self.x) ** 2) +
                B_eq * (tan(self.x) + cot(self.x)) +
                C_eq, 0
            )

            self.variables = {
                'A_eq': A_eq, 'B_eq': B_eq, 'C_eq': C_eq,
                'A_sq': A_sq, 'B_sq': B_sq, 'C_sq': C_sq,
                't1': t1, 't2': t2
            }
            break

    def _solve(self):
        t1 = self.variables['t1']
        t2 = self.variables['t2']

        final_solution = sympy.EmptySet

        for t_val in [t1, t2]:
            if abs(t_val) < 2:
                continue

            y = symbols('y')
            quad_y = y ** 2 - t_val * y + 1
            y_roots = sympy.solve(quad_y, y)

            for root in y_roots:
                sol = solveset(Eq(tan(self.x), root), self.x, domain=Reals)
                final_solution = Union(final_solution, sol)

        self.solution_obj = final_solution

    def _build_solution_steps(self):
        A_eq = self.variables['A_eq']
        B_eq = self.variables['B_eq']
        C_eq = self.variables['C_eq']
        A_sq = self.variables['A_sq']
        B_sq = self.variables['B_sq']
        C_sq = self.variables['C_sq']
        t1 = self.variables['t1']
        t2 = self.variables['t2']

        t = symbols('t')

        self.steps.append(("text", f"Маємо рівняння: ${sympy.latex(self.equation_obj)}$"))

        if B_eq != 0:
            self.steps.append(("text", r"Згрупуємо члени рівняння:"))

            grouped_latex = ""
            if A_eq == 1:
                grouped_latex += r"(\text{tg}^2 x + \text{ctg}^2 x)"
            elif A_eq == -1:
                grouped_latex += r"-(\text{tg}^2 x + \text{ctg}^2 x)"
            else:
                grouped_latex += rf"{sympy.latex(A_eq)}(\text{{tg}}^2 x + \text{{ctg}}^2 x)"

            sign_b = "+" if B_eq > 0 else "-"
            val_b = "" if abs(B_eq) == 1 else sympy.latex(abs(B_eq))
            grouped_latex += rf" {sign_b} {val_b}(\text{{tg}} x + \text{{ctg}} x)"

            if C_eq != 0:
                sign_c = "+" if C_eq > 0 else "-"
                grouped_latex += rf" {sign_c} {sympy.latex(abs(C_eq))}"

            grouped_latex += " = 0"
            self.steps.append(("math", grouped_latex))

        self.steps.append(("text", r"Введемо заміну $t = \text{tg } x + \text{ctg } x$."))
        self.steps.append(("text", r"Піднесемо заміну до квадрату:"))
        self.steps.append(("math",
                           r"t^2 = (\text{tg } x + \text{ctg } x)^2 = \text{tg}^2 x + 2\text{tg } x \text{ctg } x + \text{ctg}^2 x"))
        self.steps.append(("text", r"Оскільки $\text{tg } x \cdot \text{ctg } x = 1$, маємо:"))
        self.steps.append(
            ("math", r"t^2 = \text{tg}^2 x + 2 + \text{ctg}^2 x \implies \text{tg}^2 x + \text{ctg}^2 x = t^2 - 2"))

        self.steps.append(("text", "Підставимо вирази через $t$ у рівняння:"))

        sub_eq_latex = ""
        if A_eq == 1:
            sub_eq_latex += "(t^2 - 2)"
        elif A_eq == -1:
            sub_eq_latex += "-(t^2 - 2)"
        else:
            sub_eq_latex += rf"{sympy.latex(A_eq)}(t^2 - 2)"

        if B_eq != 0:
            sign_b = "+" if B_eq > 0 else "-"
            val_b = "" if abs(B_eq) == 1 else sympy.latex(abs(B_eq))
            sub_eq_latex += rf" {sign_b} {val_b}t"

        if C_eq != 0:
            sign_c = "+" if C_eq > 0 else "-"
            sub_eq_latex += rf" {sign_c} {sympy.latex(abs(C_eq))}"

        sub_eq_latex += " = 0"
        self.steps.append(("math", sub_eq_latex))

        self.steps.append(("text", "Розкриємо дужки та зведемо подібні доданки:"))
        quad_latex = self._format_polynomial_latex(A_sq, B_sq, C_sq, "t")
        self.steps.append(("math", quad_latex))

        self.steps.append(("text", "Знаходимо корені квадратного рівняння:"))
        if t1 == t2:
            self.steps.append(("math", f"t = {sympy.latex(t1)}"))
        else:
            self.steps.append(("math", f"t_1 = {sympy.latex(t1)}, \\quad t_2 = {sympy.latex(t2)}"))

        self.steps.append(("text",
                           r"Зауважимо, що $|t| = |\text{tg } x + \text{ctg } x| = |\text{tg } x + \frac{1}{\text{tg } x}| \ge 2$ (за нерівністю Коші)."))

        self.steps.append(("text", "Повертаємось до заміни:"))

        final_sets_latex = []

        unique_roots = [t1]
        if t1 != t2: unique_roots.append(t2)

        for i, root in enumerate(unique_roots, 1):
            prefix = f"{i}) " if len(unique_roots) > 1 else ""
            root_latex = sympy.latex(root)

            if abs(root) < 2:
                self.steps.append(("math",
                                   rf"{prefix} \text{{tg }} x + \text{{ctg }} x = {root_latex} \implies \text{{розв'язків немає, бо }} |{root_latex}| < 2"))
            else:
                self.steps.append(("math", rf"{prefix} \text{{tg }} x + \text{{ctg }} x = {root_latex}"))
                self.steps.append(("math",
                                   rf"\text{{tg }} x + \frac{{1}}{{\text{{tg }} x}} = {root_latex} \implies \text{{tg}}^2 x - ({root_latex})\text{{tg }} x + 1 = 0"))

                y = symbols('y')
                y_roots = sympy.solve(y ** 2 - root * y + 1, y)

                for yr in y_roots:
                    angle = sympy.atan(yr)
                    if isinstance(angle, sympy.atan):
                        angle_latex = rf"\text{{arctg}}({sympy.latex(yr)})"
                    else:
                        angle_latex = sympy.latex(angle)

                    sol_line = rf"x = {angle_latex} + \pi n, n \in \mathbb{{Z}}"
                    self.steps.append(("math", rf"\text{{tg }} x = {sympy.latex(yr)} \implies {sol_line}"))

                    final_sets_latex.append(
                        rf"\left\{{ {angle_latex} + \pi n \; \middle| \; n \in \mathbb{{Z}} \right\}}")

        self.steps.append(("text", "Кінцева відповідь:"))
        if not final_sets_latex:
            self.steps.append(("text", "Розв'язків немає."))
        else:
            full_response = r" \cup ".join(final_sets_latex)
            self.steps.append(("math", rf"x \in {full_response}"))