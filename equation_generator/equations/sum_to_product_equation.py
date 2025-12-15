import random
import sympy
from sympy import sin, cos, tan, cot, pi, Eq, solveset, Reals, symbols, ImageSet, Lambda, Integers, EmptySet

from ..base_class import TrigonometricEquation


class SumToProductEquation(TrigonometricEquation):

    def __init__(self):
        self.formula_map = {
            'sin+sin': r'\sin(\alpha) + \sin(\beta) = 2\sin\left(\frac{\alpha+\beta}{2}\right)\cos\left(\frac{\alpha-\beta}{2}\right)',
            'sin-sin': r'\sin(\alpha) - \sin(\beta) = 2\sin\left(\frac{\alpha-\beta}{2}\right)\cos\left(\frac{\alpha+\beta}{2}\right)',
            'cos+cos': r'\cos(\alpha) + \cos(\beta) = 2\cos\left(\frac{\alpha+\beta}{2}\right)\cos\left(\frac{\alpha-\beta}{2}\right)',
            'cos-cos': r'\cos(\alpha) - \cos(\beta) = -2\sin\left(\frac{\alpha+\beta}{2}\right)\sin\left(\frac{\alpha-\beta}{2}\right)'
        }

        n_int = symbols('n', integer=True)
        self.special_case_map = {
            'sin': {
                0: pi * n_int
            },
            'cos': {
                0: pi / 2 + pi * n_int
            }
        }

        super().__init__()

    def _generate(self):

        f = random.choice([sin, cos])
        f_name = f.__name__

        while True:
            alpha_arg = random.randint(2, 7)
            beta_arg = random.randint(1, 5)

            if alpha_arg == beta_arg:
                continue

            if beta_arg > alpha_arg:
                alpha_arg, beta_arg = beta_arg, alpha_arg

            if (alpha_arg % 2) == (beta_arg % 2):
                break

        op = random.choice(['+', '-'])

        alpha_expr = alpha_arg * self.x
        beta_expr = beta_arg * self.x

        self.variables = {
            'f': f, 'op': op,
            'alpha_arg': alpha_arg,
            'beta_arg': beta_arg
        }

        if op == '+':
            self.equation_obj = Eq(f(alpha_expr) + f(beta_expr), 0)
        else:
            self.equation_obj = Eq(f(alpha_expr) - f(beta_expr), 0)

        alpha_latex = sympy.latex(alpha_expr)
        beta_latex = sympy.latex(beta_expr)

        self.variables[
            'pretty_latex'] = rf"{f_name}{{\left({alpha_latex} \right)}} {op} {f_name}{{\left({beta_latex} \right)}} = 0"

        self._calculate_solutions()

    def get_equation_latex(self) -> str:
        if 'pretty_latex' in self.variables:
            return self.variables['pretty_latex']
        return super().get_equation_latex()

    def _calculate_solutions(self):
        f = self.variables['f']
        op = self.variables['op']
        alpha_expr = self.variables['alpha_arg'] * self.x
        beta_expr = self.variables['beta_arg'] * self.x
        f_name = f.__name__
        n = symbols('n', integer=True)

        if f_name == 'sin' and op == '+':
            lhs = 2 * sin((alpha_expr + beta_expr) / 2) * cos((alpha_expr - beta_expr) / 2)
        elif f_name == 'sin' and op == '-':
            lhs = 2 * sin((alpha_expr - beta_expr) / 2) * cos((alpha_expr + beta_expr) / 2)
        elif f_name == 'cos' and op == '+':
            lhs = 2 * cos((alpha_expr + beta_expr) / 2) * cos((alpha_expr - beta_expr) / 2)
        else:  # cos-cos
            lhs = -2 * sin((alpha_expr + beta_expr) / 2) * sin((alpha_expr - beta_expr) / 2)

        simplified_lhs = sympy.simplify(lhs)

        factors = simplified_lhs.as_ordered_factors()
        self.variables['factors'] = factors

        final_solution_set = EmptySet
        sub_solutions_list = []

        for factor in factors:
            if factor.is_number:
                continue

            sub_eq = Eq(factor, 0)
            f_sub = factor.func
            arg_expr = factor.args[0]
            rhs = 0

            special_formula_expr = self.special_case_map.get(f_sub.__name__, {}).get(rhs)

            if special_formula_expr:
                simple_eq = Eq(arg_expr, special_formula_expr)
                sub_sol_expr = sympy.solve(simple_eq, self.x)[0].expand()
                sub_sol_set = ImageSet(Lambda(n, sub_sol_expr), Integers)
            else:
                sub_sol_set = solveset(sub_eq, self.x, domain=Reals)

            sub_solutions_list.append(sub_sol_set)
            final_solution_set = final_solution_set.union(sub_sol_set)

        self.variables['sub_solutions'] = sub_solutions_list
        self.solution_obj = final_solution_set

    def _solve(self):
        pass

    def _build_solution_steps(self):
        if self.equation_obj is None: return

        f = self.variables['f']
        op = self.variables['op']
        alpha_str = f"{self.variables['alpha_arg']}*x"
        beta_str = f"{self.variables['beta_arg']}*x"

        f_name = f.__name__
        formula_key = f"{f_name}{op}{f_name}"
        formula = self.formula_map.get(formula_key)

        self.steps.append(("text", f"Маємо рівняння: ${self.variables['pretty_latex']}$"))
        self.steps.append(("text", "Це рівняння, що розв'язується перетворенням суми (або різниці) у добуток."))
        self.steps.append(("text", rf"Застосуємо формулу: ${formula}$"))

        alpha_expr = sympy.sympify(alpha_str)
        beta_expr = sympy.sympify(beta_str)

        self.steps.append(
            ("text", rf"Підставляємо $\alpha = {sympy.latex(alpha_expr)}$ та $\beta = {sympy.latex(beta_expr)}$:"))

        factors = self.variables['factors']
        simplified_lhs = 1
        for f_ in factors: simplified_lhs *= f_

        self.steps.append(("math", f"{sympy.latex(simplified_lhs)} = 0"))
        self.steps.append(("text", "Добуток дорівнює нулю, коли хоча б один із множників дорівнює нулю."))

        sub_sols = self.variables['sub_solutions']
        factor_index = 0
        n_latex = symbols('n')

        for i, factor in enumerate(factors):
            if factor.is_number:
                continue

            sub_eq = Eq(factor, 0)

            f_sub = factor.func
            f_sub_name = f_sub.__name__
            arg_expr = factor.args[0]
            rhs = 0

            special_formula_expr = self.special_case_map.get(f_sub_name, {}).get(rhs)

            self.steps.append(("text", f"Розв'язуємо {factor_index + 1}-й множник:"))
            self.steps.append(("math", sympy.latex(sub_eq)))

            is_simple_arg = (arg_expr == self.x)

            if not is_simple_arg:
                if special_formula_expr:
                    formula_latex = sympy.latex(special_formula_expr.subs(symbols('n', integer=True), n_latex))
                    formula_str = f"{sympy.latex(arg_expr)} = {formula_latex}, n " + r"\in \mathbb{Z}"
                    self.steps.append(("math", formula_str))

            sub_sol_set = sub_sols[factor_index]
            sub_sol_latex = sympy.latex(sub_sol_set)

            if not is_simple_arg:
                self.steps.append(("text", "Виражаємо $x$:"))
            else:
                self.steps.append(("text", "Отримуємо розв'язок:"))

            self.steps.append(("math", f"x = {sub_sol_latex}"))

            factor_index += 1

        final_sol_latex = sympy.latex(self.solution_obj)
        self.steps.append(("text", "Об'єднуючи всі розв'язки, отримуємо кінцеву відповідь:"))
        self.steps.append(("math", f"x = {final_sol_latex}"))