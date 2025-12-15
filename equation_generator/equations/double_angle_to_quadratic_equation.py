import random
import sympy
from sympy import sin, cos, tan, Eq, solveset, Reals, symbols, Rational, sqrt
from ..base_class import TrigonometricEquation


class DoubleAngleToQuadraticEquation(TrigonometricEquation):

    def _generate(self):
        while True:
            target_func = random.choice([sin, cos])

            nice_roots = [0, 1, -1, Rational(1, 2), Rational(-1, 2), sqrt(3) / 2, -sqrt(3) / 2, sqrt(2) / 2,
                          -sqrt(2) / 2]
            other_roots = [2, -2, 3, Rational(3, 2)]

            t1 = random.choice(nice_roots)
            t2 = random.choice(nice_roots + other_roots)

            if t1 + t2 == 0:
                continue

            if (abs(t1) > 1) and (abs(t2) > 1):
                continue

            a_quad = 1
            b_quad = -(t1 + t2)
            c_quad = t1 * t2

            coeffs = [a_quad, b_quad, c_quad]
            denoms = [c.q for c in coeffs if hasattr(c, 'q')]
            lcm = sympy.lcm(denoms)

            a = int(a_quad * lcm)
            b = int(b_quad * lcm)
            c = int(c_quad * lcm)

            if a == 0 and b == 0:
                continue

            if target_func == sin:
                A_eq = -a
                B_eq = 2 * b
                C_eq = a + 2 * c

                term_double = cos(2 * self.x)
                term_linear = sin(self.x)
                formula_latex = r"\cos(2x) = 1 - 2\sin^2(x)"
                substitution_latex = rf"{A_eq}(1 - 2\sin^2(x)) + {B_eq}\sin(x) + {C_eq} = 0" if A_eq != 0 else ""

            else:
                A_eq = a
                B_eq = 2 * b
                C_eq = a + 2 * c

                term_double = cos(2 * self.x)
                term_linear = cos(self.x)
                formula_latex = r"\cos(2x) = 2\cos^2(x) - 1"
                substitution_latex = rf"{A_eq}(2\cos^2(x) - 1) + {B_eq}\cos(x) + {C_eq} = 0" if A_eq != 0 else ""

            if A_eq == 0:
                continue

            self.equation_obj = Eq(A_eq * term_double + B_eq * term_linear + C_eq, 0)

            self.variables = {
                'target_func': target_func,
                't1': t1,
                't2': t2,
                'a': a, 'b': b, 'c': c,
                'formula_latex': formula_latex,
                'A_eq': A_eq, 'B_eq': B_eq, 'C_eq': C_eq
            }
            break

    def _solve(self):
        target_func = self.variables['target_func']
        t1 = self.variables['t1']
        t2 = self.variables['t2']

        sol1 = solveset(Eq(target_func(self.x), t1), self.x, domain=Reals)
        sol2 = solveset(Eq(target_func(self.x), t2), self.x, domain=Reals)

        self.solution_obj = sympy.Union(sol1, sol2)

    def _build_solution_steps(self):
        target_func = self.variables['target_func']
        t1 = self.variables['t1']
        t2 = self.variables['t2']
        a = self.variables['a']
        b = self.variables['b']
        c = self.variables['c']
        formula_latex = self.variables['formula_latex']
        A_eq = self.variables['A_eq']
        B_eq = self.variables['B_eq']
        C_eq = self.variables['C_eq']

        t = symbols('t')
        func_name = target_func.__name__

        self.steps.append(("text", f"Маємо рівняння: ${sympy.latex(self.equation_obj)}$"))

        self.steps.append(("text",
                           rf"Рівняння містить $\cos(2x)$ та ${func_name}(x)$. Використаємо формулу косинуса подвійного кута, щоб звести все до ${func_name}(x)$:"))
        self.steps.append(("math", formula_latex))

        self.steps.append(("text", "Підставимо формулу в рівняння:"))

        if A_eq == 1:
            A_str = ""
        elif A_eq == -1:
            A_str = "-"
        else:
            A_str = str(A_eq)

        if target_func == sin:
            sub_step = f"{A_str}(1 - 2\\sin^2 x)"
        else:
            sub_step = f"{A_str}(2\\cos^2 x - 1)"

        B_sign = "+" if B_eq >= 0 else "-"
        C_sign = "+" if C_eq >= 0 else "-"
        sub_step += f" {B_sign} {abs(B_eq)}{func_name} x {C_sign} {abs(C_eq)} = 0"

        self.steps.append(("math", sub_step))

        self.steps.append(("text", "Розкриємо дужки та зведемо подібні доданки:"))

        quad_eq_func = Eq(a * target_func(self.x) ** 2 + b * target_func(self.x) + c, 0)
        self.steps.append(("math", sympy.latex(quad_eq_func)))

        self.steps.append(("text", rf"Введемо заміну $t = {func_name}(x)$ (де $|t| \le 1$):"))
        quad_eq_t = Eq(a * t ** 2 + b * t + c, 0)
        self.steps.append(("math", sympy.latex(quad_eq_t)))

        self.steps.append(("text", "Знаходимо корені квадратного рівняння:"))
        self.steps.append(("math", f"t_1 = {sympy.latex(t1)}, \\quad t_2 = {sympy.latex(t2)}"))

        self.steps.append(("text", "Повертаємось до заміни:"))

        if abs(t1) <= 1:
            self.steps.append(("math",
                               rf"1) {func_name}(x) = {sympy.latex(t1)} \implies x \in {sympy.latex(solveset(Eq(target_func(self.x), t1), self.x))}"))
        else:
            self.steps.append(("math",
                               rf"1) {func_name}(x) = {sympy.latex(t1)} \implies \text{{розв'язків немає, бо }} |{sympy.latex(t1)}| > 1"))

        if t1 != t2:
            if abs(t2) <= 1:
                self.steps.append(("math",
                                   rf"2) {func_name}(x) = {sympy.latex(t2)} \implies x \in {sympy.latex(solveset(Eq(target_func(self.x), t2), self.x))}"))
            else:
                self.steps.append(("math",
                                   rf"2) {func_name}(x) = {sympy.latex(t2)} \implies \text{{розв'язків немає, бо }} |{sympy.latex(t2)}| > 1"))

        final_solution_latex = sympy.latex(self.solution_obj)
        self.steps.append(("text", f"Об'єднуючи розв'язки, отримуємо кінцеву відповідь:"))
        self.steps.append(("math", f"x = {final_solution_latex}"))