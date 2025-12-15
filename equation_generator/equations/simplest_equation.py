import random
import sympy
from sympy import sin, cos, tan, cot, pi, Eq, solveset, Reals, oo, zoo, nan, symbols, Rational, EmptySet, sqrt, \
    ImageSet, Lambda, Integers

from ..base_class import TrigonometricEquation


class SimplestEquation(TrigonometricEquation):

    def __init__(self):
        self.general_formula_map = {
            'sin': r't = (-1)^n \arcsin(val) + \pi n, n \in \mathbb{Z}',
            'cos': r't = \pm \arccos(val) + 2\pi n, n \in \mathbb{Z}',
            'tan': r't = \arctan(val) + \pi n, n \in \mathbb{Z}',
            'cot': r't = \text{arcctg}(val) + \pi n, n \in \mathbb{Z}'
        }

        n = symbols('n', integer=True)
        self.special_case_map = {
            'sin': {
                0: (pi * n, r't = \pi n, n \in \mathbb{Z}'),
                1: (pi / 2 + 2 * pi * n, r't = \frac{\pi}{2} + 2\pi n, n \in \mathbb{Z}'),
                -1: (-pi / 2 + 2 * pi * n, r't = -\frac{\pi}{2} + 2\pi n, n \in \mathbb{Z}')
            },
            'cos': {
                0: (pi / 2 + pi * n, r't = \frac{\pi}{2} + \pi n, n \in \mathbb{Z}'),
                1: (2 * pi * n, r't = 2\pi n, n \in \mathbb{Z}'),
                -1: (pi + 2 * pi * n, r't = \pi + 2\pi n, n \in \mathbb{Z}')
            }
        }

        super().__init__()

    def _generate(self):
        sqrt2, sqrt3 = sympy.sqrt(2), sympy.sqrt(3)

        tabular_rhs_map = {
            'sin': [
                0, 1, -1,
                Rational(1, 2), Rational(-1, 2),
                sqrt2 / 2, -sqrt2 / 2,
                sqrt3 / 2, -sqrt3 / 2
            ],
            'cos': [
                0, 1, -1,
                Rational(1, 2), Rational(-1, 2),
                sqrt2 / 2, -sqrt2 / 2,
                sqrt3 / 2, -sqrt3 / 2
            ],
            'tan': [
                0, 1, -1,
                sqrt3 / 3, -sqrt3 / 3,
                sqrt3, -sqrt3
            ],
            'cot': [
                0, 1, -1,
                sqrt3 / 3, -sqrt3 / 3,
                sqrt3, -sqrt3
            ]
        }

        while True:
            f = random.choice([sin, cos, tan, cot])
            f_name = f.__name__

            k = random.choice([1, 1, 2, 3])
            b = random.choice([0, 0, pi / 6, pi / 4, pi / 3])
            A = random.choice([1, 1, 2, -1])

            a_rhs = random.choice(tabular_rhs_map[f_name])
            a_simplified = A * a_rhs

            self.variables = {'A': A, 'k': k, 'b': b, 'a_rhs': a_rhs, 'f': f}
            self.equation_obj = Eq(A * f(k * self.x + b), a_simplified)

            if self.equation_obj.has(oo, zoo, nan):
                continue

            break

    def _solve(self):
        A = self.variables['A']
        k = self.variables['k']
        b = self.variables['b']
        a_rhs = self.variables['a_rhs']
        f = self.variables['f']
        f_name = f.__name__
        n = symbols('n', integer=True)

        if (f_name == 'sin' or f_name == 'cos') and (a_rhs > 1 or a_rhs < -1):
            self.solution_obj = EmptySet
            return

        arg_expr = k * self.x + b
        special_case_data = self.special_case_map.get(f_name, {}).get(a_rhs)

        if special_case_data:
            formula_expr, formula_latex = special_case_data
            simple_eq = Eq(arg_expr, formula_expr)
            sub_sol_expr = sympy.solve(simple_eq, self.x)[0].expand()
            self.solution_obj = ImageSet(Lambda(n, sub_sol_expr), Integers)
        else:
            self.solution_obj = solveset(self.equation_obj, self.x, domain=Reals)

    def _build_solution_steps(self):
        if self.equation_obj is None:
            return

        A = self.variables['A']
        k = self.variables['k']
        b = self.variables['b']
        a_rhs = self.variables['a_rhs']
        f = self.variables['f']
        f_name = f.__name__
        t = symbols('t')
        n = symbols('n', integer=True)

        self.steps.append(("text", f"Маємо рівняння: ${sympy.latex(self.equation_obj)}$"))

        simple_eq_obj = Eq(f(k * self.x + b), a_rhs)
        if A != 1:
            self.steps.append(("text", f"Розділимо обидві частини на ${A}:"))
            self.steps.append(("math", sympy.latex(simple_eq_obj)))

        if (f_name == 'sin' or f_name == 'cos') and (a_rhs > 1 or a_rhs < -1):
            self.steps.append(("text",
                               rf"Оскільки ${sympy.latex(a_rhs)}$ не належить області значень $\left[-1, 1\right]$ для функції ${f_name}$,"))
            self.steps.append(("text", "рівняння не має дійсних розв'язків."))
            return

        is_simple_arg = (k == 1 and b == 0)

        if is_simple_arg:
            arg_expr_latex = "x"
            formula_variable = "x"
        else:
            arg_expr = k * self.x + b
            arg_expr_latex = sympy.latex(arg_expr)
            t_equation = Eq(f(t), a_rhs)
            formula_variable = "t"
            self.steps.append(
                ("text", rf"Це найпростіше тригонометричне рівняння. Введемо заміну $t = {arg_expr_latex}$:"))
            self.steps.append(("math", sympy.latex(t_equation)))

        special_case_data = self.special_case_map.get(f_name, {}).get(a_rhs)

        if special_case_data:
            formula_expr, formula_latex = special_case_data
            formula_str = formula_latex.replace('t', formula_variable)
            self.steps.append(("text", "Це окремий випадок. Використовуємо спрощену формулу:"))
            self.steps.append(("math", formula_str))
        else:
            self.steps.append(("text", rf"Застосуємо загальну формулу розв'язку для ${f_name}$:"))
            val_latex = sympy.latex(a_rhs)
            raw_formula = self.general_formula_map[f_name]
            formula_with_val = raw_formula.replace('val', val_latex)
            final_formula = formula_with_val.replace('t', formula_variable)
            self.steps.append(("math", final_formula))

        if not is_simple_arg:
            if special_case_data:
                formula_expr, _ = special_case_data
                t_solution = formula_expr
            else:
                t_solution = solveset(Eq(f(t), a_rhs), t, domain=Reals)

            self.steps.append(("text", f"Підставляємо наше значення та розв'язуємо для $t$:"))
            self.steps.append(("math", f"t = {sympy.latex(t_solution)}"))
            self.steps.append(("text", f"Повертаємось до заміни $t = {arg_expr_latex}$:"))
            self.steps.append(("math", f"{arg_expr_latex} = {sympy.latex(t_solution)}"))
            self.steps.append(("text", f"Виражаємо $x$ та отримуємо кінцеву відповідь:"))
        else:
            self.steps.append(("text", "Обчислюємо та отримуємо кінцеву відповідь:"))

        self.steps.append(("math", f"x = {sympy.latex(self.solution_obj)}"))