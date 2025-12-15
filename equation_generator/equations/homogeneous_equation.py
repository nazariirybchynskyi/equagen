import random
import sympy
from sympy import sin, cos, tan, Eq, solveset, Reals, symbols, Rational
from ..base_class import TrigonometricEquation

class HomogeneousEquation(TrigonometricEquation):

    def _generate(self):
        t1 = random.choice([-3, -2, -1, 1, 2, 3])
        t2 = random.choice([-1, 1, Rational(1, 2), 2])

        A = 1
        B = -(t1 + t2)
        C = t1 * t2

        if B == int(B): B = int(B)
        if C == int(C): C = int(C)

        self.variables = {'A': A, 'B': B, 'C': C, 't1': t1, 't2': t2}

        self.equation_obj = Eq(
            A * sin(self.x) ** 2 + B * sin(self.x) * cos(self.x) + C * cos(self.x) ** 2,
            0
        )

    def _solve(self):
        t_eq_1 = Eq(tan(self.x), self.variables['t1'])
        t_eq_2 = Eq(tan(self.x), self.variables['t2'])

        sol1 = solveset(t_eq_1, self.x, domain=Reals)
        sol2 = solveset(t_eq_2, self.x, domain=Reals)

        self.solution_obj = sympy.Union(sol1, sol2)

    def _build_solution_steps(self):
        A = self.variables['A']
        B = self.variables['B']
        C = self.variables['C']
        t1 = self.variables['t1']
        t2 = self.variables['t2']

        A_str = "" if A == 1 else str(A)
        B_str = f"+ {B}" if B > 0 else f"{B}"
        C_str = f"+ {C}" if C > 0 else f"{C}"

        self.steps = [
            ("text", f"Маємо рівняння: ${sympy.latex(self.equation_obj)}$"),
            ("text", "Це однорідне тригонометричне рівняння другого порядку."),
            ("text", rf"Розділимо обидві частини рівняння на $cos^2(x) \neq 0$:"),
            ("math", f"{A_str}\\frac{{sin^2(x)}}{{cos^2(x)}} {B_str}\\frac{{sin(x)cos(x)}}{{cos^2(x)}} {C_str}\\frac{{cos^2(x)}}{{cos^2(x)}} = 0"),
            ("text", f"Отримаємо: ${A_str}tg^2(x) {B_str}tg(x) {C_str} = 0$"),
            ("text", f"Введемо заміну $t = tg(x)$. Рівняння набуває вигляду:"),
            ("math", f"{A_str}t^2 {B_str}t {C_str} = 0"),
            ("text", f"Коренями цього квадратного рівняння (наприклад, за теоремою Вієта) є $t_1 = {t1}$ та $t_2 = {t2}$."),
            ("text", "Повертаємось до заміни:"),
            ("math", rf"$tg(x) = {t1} \implies x = {sympy.latex(solveset(Eq(tan(self.x), self.variables['t1']), self.x))}$"),
            ("math", rf"$tg(x) = {t2} \implies x = {sympy.latex(solveset(Eq(tan(self.x), self.variables['t2']), self.x))}$"),
            ("text", f"Об'єднуючи розв'язки, отримуємо кінцеву відповідь.")
        ]