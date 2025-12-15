import random
import sympy
from sympy import sin, cos, tan, cot, pi, Eq, solveset, Reals, symbols, ImageSet, Lambda, Integers, EmptySet, Add, gcd

from ..base_class import TrigonometricEquation


class GroupingEquation(TrigonometricEquation):

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

    def _get_transform(self, f, op, alpha_expr, beta_expr):
        f_name = f.__name__
        if f_name == 'sin' and op == '+':
            return 2 * sin((alpha_expr + beta_expr) / 2) * cos((alpha_expr - beta_expr) / 2)
        elif f_name == 'sin' and op == '-':
            return 2 * sin((alpha_expr - beta_expr) / 2) * cos((alpha_expr + beta_expr) / 2)
        elif f_name == 'cos' and op == '+':
            return 2 * cos((alpha_expr + beta_expr) / 2) * cos((alpha_expr - beta_expr) / 2)
        else:  # cos-cos
            return -2 * sin((alpha_expr + beta_expr) / 2) * sin((alpha_expr - beta_expr) / 2)

    def get_equation_latex(self) -> str:
        return sympy.latex(self.equation_obj)

    def _generate(self):

        k_pool = [1, 2, 3, 4]
        m_n_pool = [2, 3, 4, 5, 6, 7, 8, 9]

        while True:
            k = random.choice(k_pool)

            valid_pool = [x for x in m_n_pool if x > k and (x % 2 != k % 2)]

            if len(valid_pool) < 2:
                valid_pool = [x for x in m_n_pool if x > k and (x % 2 == k % 2)]
                if len(valid_pool) < 2:
                    continue

            m, n = random.sample(valid_pool, 2)

            alpha1_arg, beta1_arg = m + k, m - k
            alpha2_arg, beta2_arg = n + k, n - k

            args_set = {alpha1_arg, beta1_arg, alpha2_arg, beta2_arg}

            if 0 in args_set or len(args_set) < 4:
                continue

            break

        f = random.choice([sin, cos])
        op_transform = random.choice(['+', '-'])

        if op_transform == '-':
            op_group = '-'
        else:
            op_group = random.choice(['+', '-'])

        alpha1_expr, beta1_expr = alpha1_arg * self.x, beta1_arg * self.x
        alpha2_expr, beta2_expr = alpha2_arg * self.x, beta2_arg * self.x

        term1 = f(alpha1_expr)
        term2 = f(beta1_expr) if op_transform == '+' else -f(beta1_expr)
        term3 = f(alpha2_expr)
        term4 = f(beta2_expr) if op_transform == '+' else -f(beta2_expr)

        if op_group == '-':
            term3 = -term3
            term4 = -term4

        terms_list = [term1, term2, term3, term4]
        random.shuffle(terms_list)

        self.equation_obj = Eq(Add(*terms_list), 0)

        self.variables = {
            'f': f, 'op_transform': op_transform, 'op_group': op_group,
            'a1_arg': alpha1_arg, 'b1_arg': beta1_arg,
            'a2_arg': alpha2_arg, 'b2_arg': beta2_arg
        }

    def _solve(self):
        self._transform_pairs()
        self._factor_common_term()
        self._factor_remaining_group()
        self._solve_factors()

    def _transform_pairs(self):
        f = self.variables['f']
        op_t = self.variables['op_transform']
        a1_expr = self.variables['a1_arg'] * self.x
        b1_expr = self.variables['b1_arg'] * self.x
        a2_expr = self.variables['a2_arg'] * self.x
        b2_expr = self.variables['b2_arg'] * self.x

        pair1_lhs = self._get_transform(f, op_t, a1_expr, b1_expr)
        pair2_lhs = self._get_transform(f, op_t, a2_expr, b2_expr)

        self.variables['pair1_lhs'] = pair1_lhs
        self.variables['pair2_lhs'] = pair2_lhs

    def _factor_common_term(self):
        pair1_lhs = self.variables['pair1_lhs']
        pair2_lhs = self.variables['pair2_lhs']
        op_g = self.variables['op_group']

        common_factor_expr = gcd(pair1_lhs, pair2_lhs)
        other_term1_expr = sympy.simplify(pair1_lhs / common_factor_expr)
        other_term2_expr = sympy.simplify(pair2_lhs / common_factor_expr)

        self.variables['common_factor'] = common_factor_expr
        self.variables['other_term1'] = other_term1_expr
        self.variables['other_term2'] = other_term2_expr

        if op_g == '+':
            other_group = other_term1_expr + other_term2_expr
        else:
            other_group = other_term1_expr - other_term2_expr

        self.variables['other_group_expr'] = other_group

    def _factor_remaining_group(self):
        op_g = self.variables['op_group']
        other_term1_expr = self.variables['other_term1']
        other_term2_expr = self.variables['other_term2']

        other_f = other_term1_expr.func
        if not hasattr(other_f, '__name__') or not other_f.__name__ in ('sin', 'cos'):
            other_f = other_term2_expr.func

        final_factors_expr = self._get_transform(other_f, op_g, other_term1_expr.args[0], other_term2_expr.args[0])
        self.variables['final_factors_expr'] = final_factors_expr

        all_factors_expr = self.variables['common_factor'] * final_factors_expr
        all_factors = all_factors_expr.as_ordered_factors()
        self.variables['all_factors'] = all_factors

    def _solve_factors(self):
        all_factors = self.variables['all_factors']
        n_int = symbols('n', integer=True)

        final_solution_set = EmptySet
        sub_solutions_list = []

        for factor_obj in all_factors:
            if factor_obj.is_number:
                continue

            sub_eq = Eq(factor_obj, 0)
            f_sub = factor_obj.func

            if not hasattr(f_sub, '__name__') or not f_sub.__name__ in ('sin', 'cos'):
                continue

            arg_expr = factor_obj.args[0]
            rhs = 0

            special_formula_expr = self.special_case_map.get(f_sub.__name__, {}).get(rhs)

            if special_formula_expr:
                simple_eq = Eq(arg_expr, special_formula_expr)
                sub_sol_expr = sympy.solve(simple_eq, self.x)[0].expand()
                sub_sol_set = ImageSet(Lambda(n_int, sub_sol_expr), Integers)
            else:
                sub_sol_set = solveset(sub_eq, self.x, domain=Reals)

            sub_solutions_list.append(sub_sol_set)
            final_solution_set = final_solution_set.union(sub_sol_set)

        self.variables['sub_solutions'] = sub_solutions_list
        self.solution_obj = final_solution_set

    def _build_solution_steps(self):
        if self.equation_obj is None: return

        f = self.variables['f']
        op_t = self.variables['op_transform']
        op_g = self.variables['op_group']
        f_name = f.__name__
        n_latex = symbols('n')

        a1_lat = sympy.latex(self.variables['a1_arg'] * self.x)
        b1_lat = sympy.latex(self.variables['b1_arg'] * self.x)
        a2_lat = sympy.latex(self.variables['a2_arg'] * self.x)
        b2_lat = sympy.latex(self.variables['b2_arg'] * self.x)

        formula_transform_key = f"{f_name}{op_t}{f_name}"
        formula_transform = self.formula_map.get(formula_transform_key)

        self.steps.append(("text", f"Маємо рівняння: ${self.get_equation_latex()}$"))
        self.steps.append(("text", "Згрупуємо доданки в пари, щоб можна було застосувати формули суми/різниці."))

        grouped_latex = (
            rf"\left( {f_name}{{\left({a1_lat} \right)}} {op_t} {f_name}{{\left({b1_lat} \right)}} \right) "
            rf"{op_g} \left( {f_name}{{\left({a2_lat} \right)}} {op_t} {f_name}{{\left({b2_lat} \right)}} \right) = 0")
        self.steps.append(("math", grouped_latex))

        self.steps.append(("text", rf"Використовуємо формулу: ${formula_transform}$"))

        pair1_lat = sympy.latex(self.variables['pair1_lhs'])
        pair2_lat = sympy.latex(self.variables['pair2_lhs'])
        self.steps.append(("text", "Застосовуємо формулу до кожної пари:"))

        pair1_term = self.variables['pair1_lhs']
        pair2_term = self.variables['pair2_lhs']

        if op_g == '+':
            step2_expr = pair1_term + pair2_term
        else:
            step2_expr = pair1_term - pair2_term

        step2_latex = sympy.latex(step2_expr) + " = 0"
        self.steps.append(("math", step2_latex))

        common_factor_lat = sympy.latex(self.variables['common_factor'])
        other_group_lat = sympy.latex(self.variables['other_group_expr'])

        self.steps.append(("text", rf"Винесемо спільний множник ${common_factor_lat}$ за дужки:"))
        self.steps.append(("math", rf"{common_factor_lat} \cdot \left( {other_group_lat} \right) = 0"))

        other_f = self.variables['other_term1'].func
        if not hasattr(other_f, '__name__') or not other_f.__name__ in ('sin', 'cos'):
            other_f = self.variables['other_term2'].func

        other_f_name = other_f.__name__
        formula_group_key = f"{other_f_name}{op_g}{other_f_name}"
        formula_group = self.formula_map.get(formula_group_key)

        if formula_group:
            self.steps.append(("text", rf"Тепер застосуємо формулу (${formula_group}$) до виразу в дужках:"))

            final_factors_expr = self.variables['common_factor'] * self.variables['final_factors_expr']
            full_final_lat = sympy.latex(sympy.simplify(final_factors_expr)) + " = 0"
            self.steps.append(("math", full_final_lat))

        self.steps.append(("text", "Добуток дорівнює нулю, коли хоча б один із множників дорівнює нулю."))

        sub_sols = self.variables['sub_solutions']
        factor_index = 0

        for i, factor_obj in enumerate(self.variables['all_factors']):
            if factor_obj.is_number:
                continue

            sub_eq = Eq(factor_obj, 0)
            f_sub = factor_obj.func

            if not hasattr(f_sub, '__name__') or not f_sub.__name__ in ('sin', 'cos'):
                continue

            arg_expr = factor_obj.args[0]
            rhs = 0

            special_formula_expr = self.special_case_map.get(f_sub.__name__, {}).get(rhs)

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