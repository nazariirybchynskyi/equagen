import random
import sympy
from sympy import sin, cos, Eq, solveset, Reals, symbols, Rational, pi, Mul, Add, expand, Union, sqrt

from ..base_class import TrigonometricEquation

class ReducibleToHomogeneousEquation(TrigonometricEquation):

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
                terms.append(f"{sign} {coeff} {var_name}")

        if C != 0:
            sign = "+" if C > 0 else "-"
            val_latex = sympy.latex(abs(C))
            if getattr(C, 'is_Add', False):
                val_latex = f"({val_latex})"
            terms.append(f"{sign} {val_latex}")

        latex_str = " ".join(terms)
        if not latex_str.startswith('-'):
            latex_str = latex_str.lstrip('+ ')

        return latex_str + " = 0"

    def _format_homogeneous_latex(self, A, B, C):
        terms = []

        if A != 0:
            if A == 1:
                term = "\\sin^2 x"
            elif A == -1:
                term = "-\\sin^2 x"
            else:
                coeff = sympy.latex(A)
                if getattr(A, 'is_Add', False):
                    coeff = f"({coeff})"
                term = f"{coeff}\\sin^2 x"
            terms.append(term)
        else:
            pass

        if B != 0:
            sign = "+" if B > 0 else "-"
            val = abs(B)
            if val == 1:
                term = f"{sign} \\sin x \\cos x"
            else:
                coeff = sympy.latex(val)
                if getattr(val, 'is_Add', False):
                    coeff = f"({coeff})"
                term = f"{sign} {coeff}\\sin x \\cos x"
            terms.append(term)

        if C != 0:
            sign = "+" if C > 0 else "-"
            val = abs(C)
            if val == 1:
                term = f"{sign} \\cos^2 x"
            else:
                coeff = sympy.latex(val)
                if getattr(val, 'is_Add', False):
                    coeff = f"({coeff})"
                term = f"{sign} {coeff}\\cos^2 x"
            terms.append(term)

        homo_str = " ".join(terms)
        if not homo_str.startswith('-'):
            homo_str = homo_str.lstrip('+ ')

        return homo_str + " = 0"

    def _generate(self):
        while True:
            nice_roots = [
                0,
                1, -1,
                sqrt(3), -sqrt(3),
                sqrt(3) / 3, -sqrt(3) / 3
            ]

            t1 = random.choice(nice_roots)
            t2 = random.choice(nice_roots)

            A_kern = random.choice([1, 2, -1, -2])
            B_kern = -A_kern * (t1 + t2)
            C_kern = A_kern * (t1 * t2)

            coeffs_list = [A_kern, B_kern, C_kern]
            denoms = [c.q for c in coeffs_list if hasattr(c, 'q')]
            if denoms:
                lcm = sympy.lcm(denoms)
                A_kern = A_kern * lcm
                B_kern = B_kern * lcm
                C_kern = C_kern * lcm

            if hasattr(A_kern, 'is_Integer') and A_kern.is_Integer: A_kern = int(A_kern)
            if hasattr(B_kern, 'is_Integer') and B_kern.is_Integer: B_kern = int(B_kern)
            if hasattr(C_kern, 'is_Integer') and C_kern.is_Integer: C_kern = int(C_kern)

            D = random.choice([1, -1, 2, -2, 3])

            A_final = A_kern
            C_final = C_kern

            if A_final == 0 or B_kern == 0 or C_final == 0:
                continue

            A_orig = A_final + D
            C_orig = C_final + D
            B_orig = B_kern

            self.equation_obj = Eq(
                A_orig * sin(self.x) ** 2 +
                B_orig * sin(self.x) * cos(self.x) +
                C_orig * cos(self.x) ** 2,
                D
            )

            self.variables = {
                'A_orig': A_orig, 'B_orig': B_orig, 'C_orig': C_orig, 'D': D,
                'A_final': A_final, 'B_final': B_kern, 'C_final': C_final,
                't1': t1, 't2': t2
            }
            break

    def _solve(self):
        t1 = self.variables['t1']
        t2 = self.variables['t2']

        sol1 = solveset(Eq(sympy.tan(self.x), t1), self.x, domain=Reals)
        sol2 = solveset(Eq(sympy.tan(self.x), t2), self.x, domain=Reals)

        self.solution_obj = Union(sol1, sol2)

    def _build_solution_steps(self):
        A_orig = self.variables['A_orig']
        B_orig = self.variables['B_orig']
        C_orig = self.variables['C_orig']
        D = self.variables['D']
        A_final = self.variables['A_final']
        B_final = self.variables['B_final']
        C_final = self.variables['C_final']
        t1 = self.variables['t1']
        t2 = self.variables['t2']

        t = symbols('t')

        self.steps.append(("text", f"Маємо рівняння: ${self.get_equation_latex()}$"))

        self.steps.append(("text",
                           f"Це рівняння, що зводиться до однорідного, оскільки права частина є константою $D={sympy.latex(D)}$."))

        self.steps.append(
            ("text", rf"Помножимо константу $D={sympy.latex(D)}$ на тригонометричну одиницю $\sin^2 x + \cos^2 x$:"))
        self.steps.append(("math", rf"{sympy.latex(D)} = {sympy.latex(D)}(\sin^2 x + \cos^2 x)"))

        original_LHS_latex = sympy.latex(self.equation_obj.lhs)

        D_latex = sympy.latex(D)
        if getattr(D, 'is_Add', False):
            D_latex = f"({D_latex})"
        D_sub_latex = f"{D_latex} (\\sin^2 x + \\cos^2 x)"

        self.steps.append(("text", "Підставляємо $\\sin^2 x + \\cos^2 x = 1$ у праву частину рівняння:"))
        self.steps.append(("math", rf"{original_LHS_latex} = {D_sub_latex}"))

        final_homo_eq = self._format_homogeneous_latex(A_final, B_final, C_final)
        self.steps.append(("text", "Переносимо всі члени вліво та зводимо подібні, отримуючи однорідне рівняння:"))
        self.steps.append(("math", final_homo_eq))

        self.steps.append(("text", rf"Ділимо рівняння на $\cos^2(x) \neq 0$:"))
        tan_eq_latex = self._format_polynomial_latex(A_final, B_final, C_final, r'\tan x')
        self.steps.append(("math", tan_eq_latex))

        self.steps.append(("text", rf"Вводимо заміну $t = \text{{tg}}(x)$:"))
        t_eq_latex = self._format_polynomial_latex(A_final, B_final, C_final, "t")
        self.steps.append(("math", t_eq_latex))

        self.steps.append(("text", "Розв'язуємо квадратне рівняння:"))
        roots_latex = f"t_1 = {sympy.latex(t1)}"
        if t2 != t1:
            roots_latex += rf", \quad t_2 = {sympy.latex(t2)}"
        self.steps.append(("math", roots_latex))

        self.steps.append(("text", rf"Повертаємось до заміни $\text{{tg}}(x) = t$:"))

        sol1_latex = sympy.latex(solveset(Eq(sympy.tan(self.x), t1), self.x, domain=Reals))
        self.steps.append(("math", rf"1) \text{{tg}}(x) = {sympy.latex(t1)} \implies x = {sol1_latex}"))

        if t2 != t1:
            sol2_latex = sympy.latex(solveset(Eq(sympy.tan(self.x), t2), self.x, domain=Reals))
            self.steps.append(("math", rf"2) \text{{tg}}(x) = {sympy.latex(t2)} \implies x = {sol2_latex}"))

        final_sol_latex = sympy.latex(self.solution_obj)
        self.steps.append(("text", "Об'єднуючи розв'язки, отримуємо кінцеву відповідь:"))
        self.steps.append(("math", f"x = {final_sol_latex}"))