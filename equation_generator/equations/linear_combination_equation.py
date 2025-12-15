import random
import sympy
from sympy import sin, cos, Eq, solveset, Reals, symbols, pi, sqrt, Rational, atan2, Integers
from ..base_class import TrigonometricEquation


class LinearCombinationEquation(TrigonometricEquation):

    def __init__(self):
        self.base_angles = {
            pi / 6: (Rational(1, 2), sqrt(3) / 2),
            pi / 4: (sqrt(2) / 2, sqrt(2) / 2),
            pi / 3: (sqrt(3) / 2, Rational(1, 2)),
        }
        self.base_amplitudes = [2, sqrt(2)]
        self.target_values = [0, 1, -1, Rational(1, 2), Rational(-1, 2), sqrt(2) / 2, -sqrt(2) / 2]

        super().__init__()

    def _get_phi_latex(self, phi):
        if phi == pi / 6: return r'\frac{\pi}{6}'
        if phi == pi / 4: return r'\frac{\pi}{4}'
        if phi == pi / 3: return r'\frac{\pi}{3}'
        return sympy.latex(phi)

    def _generate(self):
        while True:
            phi_base, (sin_phi, cos_phi) = random.choice(list(self.base_angles.items()))
            D_base = random.choice(self.base_amplitudes)
            S_target = random.choice(self.target_values)

            reduction_type = random.choice(['sin_sum', 'cos_diff'])

            if reduction_type == 'sin_sum':
                a = D_base * cos_phi
                b = D_base * sin_phi
            else:
                a = D_base * sin_phi
                b = D_base * cos_phi

            c = D_base * S_target
            D = D_base

            if abs(c) > abs(D):
                continue

            if a == 0 or b == 0:
                continue

            phi = atan2(b, a)

            a = sympy.simplify(a)
            b = sympy.simplify(b)
            c = sympy.simplify(c)
            D = sympy.simplify(D)

            self.variables = {
                'a': a, 'b': b, 'c': c, 'D': D, 'phi': phi, 'S': S_target, 'phi_base': phi_base,
                'reduction_type': reduction_type
            }

            self.equation_obj = Eq(a * sin(self.x) + b * cos(self.x), c)
            break

    def _solve(self):
        a = self.variables['a']
        b = self.variables['b']
        c = self.variables['c']

        self.solution_obj = solveset(Eq(a * sin(self.x) + b * cos(self.x), c), self.x, domain=Reals)

    def _build_solution_steps(self):
        a = self.variables['a']
        b = self.variables['b']
        c = self.variables['c']
        D = self.variables['D']
        phi = self.variables['phi']
        phi_base = self.variables['phi_base']
        S = self.variables['S']

        a_latex = sympy.latex(a)
        b_latex = sympy.latex(b)
        c_latex = sympy.latex(c)
        D_latex = sympy.latex(D)

        self.steps.append(
            ("text", rf"Маємо рівняння виду $a \sin x + b \cos x = c$: ${sympy.latex(self.equation_obj)}$"))

        if abs(c) > abs(D):
            self.steps.append(("text",
                               rf"Перевірка умови розв'язності: $c^2 = {sympy.latex(c ** 2)}$ та $a^2 + b^2 = {sympy.latex(a ** 2 + b ** 2)} = {sympy.latex(D ** 2)}$. Оскільки $|c| > \sqrt{{a^2+b^2}}$, розв'язків немає."))
            return

        self.steps.append(("text", rf"Обчислимо допоміжний множник (амплітуду): $D = \sqrt{{a^2 + b^2}}$"))
        self.steps.append(("math", rf"D = \sqrt{{({a_latex})^2 + ({b_latex})^2}} = {D_latex}"))

        self.steps.append(("text", rf"Поділимо обидві частини рівняння на ${D_latex}$:"))

        div_eq_raw = rf"\frac{{{a_latex}}}{{{D_latex}}} \sin x + \frac{{{b_latex}}}{{{D_latex}}} \cos x = \frac{{{c_latex}}}{{{D_latex}}}"
        self.steps.append(("math", div_eq_raw))

        simplified_div_eq = rf"{sympy.latex(a / D)} \sin x + {sympy.latex(b / D)} \cos x = {sympy.latex(S)}"
        self.steps.append(("math", simplified_div_eq))

        cos_phi = a / D
        sin_phi = b / D
        phi_base_latex = self._get_phi_latex(phi_base)

        self.steps.append(("text", rf"Введемо допоміжний кут $\varphi$, такий що:"))
        self.steps.append(("math",
                           rf"\cos \varphi = {sympy.latex(cos_phi)} \quad \text{{та}} \quad \sin \varphi = {sympy.latex(sin_phi)} \implies \varphi = {phi_base_latex}"))

        self.steps.append(("text", "Замінюємо коефіцієнти тригонометричними функціями кута $\varphi$:"))
        new_step_latex = rf"\cos\left({phi_base_latex}\right) \sin x + \sin\left({phi_base_latex}\right) \cos x = {sympy.latex(S)}"
        self.steps.append(("math", new_step_latex))

        self.steps.append(
            ("text", rf"Використовуємо формулу $\sin(x + \varphi) = \sin x \cos \varphi + \cos x \sin \varphi$:"))
        self.steps.append(("math", rf"\sin\left(x + {phi_base_latex}\right) = {sympy.latex(S)}"))

        self.steps.append(("text", "Розв'язуємо найпростіше тригонометричне рівняння:"))

        self.steps.append(("math", rf"x + {phi_base_latex} = \arcsin({sympy.latex(S)}) + 2\pi n, n \in \mathbb{{Z}}"))

        final_sol_latex = sympy.latex(self.solution_obj)
        self.steps.append(("text", "Виражаємо $x$ та отримуємо кінцеву відповідь:"))
        self.steps.append(("math", f"x = {final_sol_latex}"))