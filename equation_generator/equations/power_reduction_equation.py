import random
import sympy
from sympy import sin, cos, tan, cot, pi, Eq, S, solveset, Reals, symbols, ImageSet, Lambda, Integers, EmptySet, Add, gcd, \
    Rational

from ..base_class import TrigonometricEquation


class PowerReductionEquation(TrigonometricEquation):

    def __init__(self):
        self.formula_map = {
            'sin': r'\sin^2(\alpha) = \frac{1 - \cos(2\alpha)}{2}',
            'cos': r'\cos^2(\alpha) = \frac{1 + \cos(2\alpha)}{2}',
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
        if f_name == 'cos' and op == '+':
            return 2 * cos((alpha_expr + beta_expr) / 2) * cos((alpha_expr - beta_expr) / 2)
        elif f_name == 'cos' and op == '-':
            return -2 * sin((alpha_expr + beta_expr) / 2) * sin((alpha_expr - beta_expr) / 2)
        else:
            return None

    def _generate(self):

        path_type = random.choice(["3_terms_const", "4_terms"])
        self.variables = {'path_type': path_type}

        self.variables['f_type'] = random.choice([sin, cos])
        f = self.variables['f_type']

        k_pool = [1, 2, 3]
        m_n_pool = [2, 3, 4, 5, 6, 7, 8, 9]

        while True:
            k = random.choice(k_pool)
            m_n_valid_pool = [x for x in m_n_pool if x > k]

            if len(m_n_valid_pool) < 2:
                continue

            m, n = random.sample(m_n_valid_pool, 2)

            if path_type == "4_terms":
                args_set1 = {m + k, m - k}
                args_set2 = {n + k, n - k}
                if not args_set1.isdisjoint(args_set2):
                    continue

            break

        self.variables['k_arg'] = k
        self.variables['m_arg'] = m
        self.variables['n_arg'] = n

        if path_type == "3_terms_const":
            const = Rational(3, 2)

            while True:
                k_3 = random.randint(1, 3)
                m_3 = random.randint(2, 5)
                if k_3 == m_3 or (k_3 % 2 == m_3 % 2):
                    continue

                a, b, c = (m_3 + k_3), (m_3 - k_3), m_3
                if len({a, b, c}) < 3:
                    continue
                break

            final_args = [a, b, c]
            self.variables['final_args'] = final_args

            self.variables['kernel_args'] = [2 * a * self.x, 2 * b * self.x, 2 * c * self.x]
            self.variables['kernel_grouped_terms'] = [cos(2 * a * self.x), cos(2 * b * self.x), cos(2 * c * self.x)]
            self.variables['k_3'] = k_3
            self.variables['m_3'] = m_3

            shuffled_args = final_args[:]
            random.shuffle(shuffled_args)
            terms_list = [f(arg * self.x) ** 2 for arg in shuffled_args]
            self.equation_obj = Eq(Add(*terms_list), const)
            self.variables['const'] = const

        else:
            a1, b1 = (m + k), (m - k)
            a2, b2 = (n + k), (n - k)

            final_args = [Rational(a1, 2), Rational(b1, 2), Rational(a2, 2), Rational(b2, 2)]
            self.variables['final_args'] = final_args

            self.variables['kernel_args'] = [a1, b1, a2, b2]

            terms_list = [
                f(final_args[0] * self.x) ** 2,
                f(final_args[1] * self.x) ** 2,
                -f(final_args[2] * self.x) ** 2,
                -f(final_args[3] * self.x) ** 2
            ]

            random.shuffle(terms_list)
            self.equation_obj = Eq(Add(*terms_list), 0)
            self.variables['const'] = 0

    def _solve(self):
        if self.variables['path_type'] == "3_terms_const":
            self._solve_3_terms()
        else:
            self.variables['op_group'] = '-'
            self._solve_4_terms()

    def _solve_3_terms(self):
        f = self.variables['f_type']
        a, b, c = self.variables['final_args']

        t1, t2, t3 = cos(2 * a * self.x), cos(2 * b * self.x), cos(2 * c * self.x)

        kernel_eq_lhs = t1 + t2 + t3
        self.variables['kernel_eq'] = Eq(kernel_eq_lhs, 0)

        group1_terms = [cos(2 * a * self.x), cos(2 * b * self.x)]
        group2_term = cos(2 * c * self.x)

        self.variables['grouped_kernel_terms'] = [group1_terms[0], group1_terms[1], group2_term]
        self.variables['group2_term'] = group2_term

        group1_trans = self._get_transform(cos, '+', group1_terms[0].args[0], group1_terms[1].args[0])
        self.variables['group1_trans'] = group1_trans

        common_factor = gcd(group1_trans, group2_term)
        other_group = sympy.simplify(group1_trans / common_factor + group2_term / common_factor)

        self.variables['common_factor'] = common_factor
        self.variables['other_group_expr'] = other_group
        self.variables['all_factors'] = [common_factor, other_group]

        self._solve_factors(self.variables['all_factors'])

    def _solve_4_terms(self):
        k, m, n = self.variables['k_arg'], self.variables['m_arg'], self.variables['n_arg']

        T1_e, T2_e = (m + k) * self.x, (m - k) * self.x
        T3_e, T4_e = (n + k) * self.x, (n - k) * self.x

        kernel_eq_lhs = cos(T1_e) + cos(T2_e)
        kernel_eq_rhs = cos(T3_e) + cos(T4_e)
        self.variables['kernel_eq'] = Eq(kernel_eq_lhs, kernel_eq_rhs)

        trans_lhs = self._get_transform(cos, '+', T1_e, T2_e)
        trans_rhs = self._get_transform(cos, '+', T3_e, T4_e)

        self.variables['trans_lhs'] = trans_lhs
        self.variables['trans_rhs'] = trans_rhs

        common_factor = gcd(trans_lhs, trans_rhs)
        other_term1 = sympy.simplify(trans_lhs / common_factor)
        other_term2 = sympy.simplify(trans_rhs / common_factor)
        other_group = other_term1 - other_term2

        self.variables['common_factor'] = common_factor
        self.variables['other_group_expr'] = other_group

        final_factors_expr = self._get_transform(cos, '-', other_term1.args[0], other_term2.args[0])
        self.variables['final_factors_expr'] = final_factors_expr

        if final_factors_expr is None:
            final_factors_expr = S(1)

        all_factors_expr = common_factor * final_factors_expr
        self.variables['all_factors'] = all_factors_expr.as_ordered_factors()

        self._solve_factors(self.variables['all_factors'])

    def _solve_factors(self, all_factors):
        n_int = symbols('n', integer=True)
        final_solution_set = EmptySet
        sub_solutions_list = []

        for factor_obj in all_factors:
            if factor_obj.is_number: continue

            sub_eq = Eq(factor_obj, 0)

            if factor_obj.func == Add:
                sub_sol_set = solveset(sub_eq, self.x, domain=Reals)
                sub_solutions_list.append(sub_sol_set)
                final_solution_set = final_solution_set.union(sub_sol_set)
                continue

            f_sub = factor_obj.func
            if not hasattr(f_sub, '__name__') or not f_sub.__name__ in ('sin', 'cos'):
                sub_sol_set = solveset(sub_eq, self.x, domain=Reals)
                sub_solutions_list.append(sub_sol_set)
                final_solution_set = final_solution_set.union(sub_sol_set)
                continue

            arg_expr = factor_obj.args[0]
            rhs = 0
            special_formula_expr = self.special_case_map.get(f_sub.__name__, {}).get(rhs)

            if special_formula_expr:
                simple_eq = Eq(arg_expr, special_formula_expr)
                sub_sol_expr_list = sympy.solve(simple_eq, self.x)
                if sub_sol_expr_list:
                    sub_sol_expr = sub_sol_expr_list[0].expand()
                    sub_sol_set = ImageSet(Lambda(n_int, sub_sol_expr), Integers)
                else:
                    sub_sol_set = EmptySet
            else:
                sub_sol_set = solveset(sub_eq, self.x, domain=Reals)

            sub_solutions_list.append(sub_sol_set)
            final_solution_set = final_solution_set.union(sub_sol_set)

        self.variables['sub_solutions'] = sub_solutions_list
        self.solution_obj = final_solution_set

    def _build_solution_steps(self):
        if self.variables['path_type'] == "3_terms_const":
            self._build_steps_3_terms()
        else:
            self._build_steps_4_terms()

    def _build_steps_4_terms(self):
        f = self.variables['f_type']
        f_name = f.__name__
        a, b, c, d = self.variables['final_args']

        self.steps.append(("text", f"Маємо рівняння: ${self.get_equation_latex()}$"))
        self.steps.append(("text", "Перенесемо доданки, щоб згрупувати їх:"))

        grouped_eq_latex = rf"{f_name}^2({sympy.latex(a * self.x)}) + {f_name}^2({sympy.latex(b * self.x)}) = {f_name}^2({sympy.latex(c * self.x)}) + {f_name}^2({sympy.latex(d * self.x)})"
        self.steps.append(("math", grouped_eq_latex))

        formula_reduce = self.formula_map.get(f_name)
        self.steps.append(("text", rf"Застосуємо формулу пониження степеня: ${formula_reduce}$"))

        if f_name == 'sin':
            step2_lhs = rf"\frac{{1 - \cos({2 * a}x)}}{{2}} + \frac{{1 - \cos({2 * b}x)}}{{2}}"
            step2_rhs = rf" = \frac{{1 - \cos({2 * c}x)}}{{2}} + \frac{{1 - \cos({2 * d}x)}}{{2}}"
            step3 = rf"\cos({2 * a}x) + \cos({2 * b}x) = \cos({2 * c}x) + \cos({2 * d}x)"
        else:  # cos
            step2_lhs = rf"\frac{{1 + \cos({2 * a}x)}}{{2}} + \frac{{1 + \cos({2 * b}x)}}{{2}}"
            step2_rhs = rf" = \frac{{1 + \cos({2 * c}x)}}{{2}} + \frac{{1 + \cos({2 * d}x)}}{{2}}"
            step3 = rf"\cos({2 * a}x) + \cos({2 * b}x) = \cos({2 * c}x) + \cos({2 * d}x)"

        self.steps.append(("math", step2_lhs + step2_rhs))
        self.steps.append(("text", "Домножимо на 2 та спростимо:"))
        self.steps.append(("math", step3))

        formula_cos_plus = self.formula_map.get('cos+cos')
        self.steps.append(("text", rf"Застосуємо формулу суми косинусів до обох частин: ${formula_cos_plus}$"))

        trans_lhs_lat = sympy.latex(self.variables['trans_lhs'])
        trans_rhs_lat = sympy.latex(self.variables['trans_rhs'])
        self.steps.append(("math", rf"{trans_lhs_lat} = {trans_rhs_lat}"))

        self.steps.append(("text", "Перенесемо всі доданки вліво:"))
        self.steps.append(("math", rf"{trans_lhs_lat} - {trans_rhs_lat} = 0"))

        common_factor_lat = sympy.latex(self.variables['common_factor'])
        other_group_lat = sympy.latex(self.variables['other_group_expr'])

        self.steps.append(("text", rf"Винесемо спільний множник ${common_factor_lat}$ за дужки:"))
        self.steps.append(("math", rf"{common_factor_lat} \cdot \left( {other_group_lat} \right) = 0"))

        formula_cos_minus = self.formula_map.get('cos-cos')
        self.steps.append(
            ("text", rf"Тепер застосуємо формулу різниці косинусів до виразу в дужках: ${formula_cos_minus}$"))

        final_factors_expr = self.variables['common_factor'] * self.variables['final_factors_expr']
        full_final_lat = sympy.latex(sympy.simplify(final_factors_expr)) + " = 0"
        self.steps.append(("math", full_final_lat))

        self._build_steps_solve_factors()

    def _build_steps_3_terms(self):
        f = self.variables['f_type']
        f_name = f.__name__
        a, b, c = self.variables['final_args']
        const = self.variables['const']

        formula_reduce = self.formula_map.get(f_name)
        self.steps.append(("text", f"Маємо рівняння: ${self.get_equation_latex()}$"))
        self.steps.append(("text", rf"Застосуємо формулу пониження степеня: ${formula_reduce}$"))

        t1_lat = rf"\cos({2 * a}x)"
        t2_lat = rf"\cos({2 * b}x)"
        t3_lat = rf"\cos({2 * c}x)"

        if f_name == 'sin':
            step2_lhs = rf"\frac{{1 - {t1_lat}}}{{2}} + \frac{{1 - {t2_lat}}}{{2}} + \frac{{1 - {t3_lat}}}{{2}}"
            step3 = rf"{t1_lat} + {t2_lat} + {t3_lat} = 0"
        else:  # cos
            step2_lhs = rf"\frac{{1 + {t1_lat}}}{{2}} + \frac{{1 + {t2_lat}}}{{2}} + \frac{{1 + {t3_lat}}}{{2}}"
            step3 = rf"{t1_lat} + {t2_lat} + {t3_lat} = 0"

        self.steps.append(("math", rf"{step2_lhs} = {sympy.latex(const)}"))
        self.steps.append(("text", "Домножимо на 2 та спростимо:"))
        self.steps.append(("math", step3))

        g1, g2, g3 = self.variables['grouped_kernel_terms']
        grouped_lat = rf"\left( {sympy.latex(g1)} + {sympy.latex(g2)} \right) + {sympy.latex(g3)} = 0"
        self.steps.append(("text", rf"Згрупуємо доданки:"))
        self.steps.append(("math", grouped_lat))

        formula_cos_plus = self.formula_map.get('cos+cos')
        self.steps.append(("text", rf"Використовуємо формулу: ${formula_cos_plus}$"))

        group1_trans_lat = sympy.latex(self.variables['group1_trans'])
        group2_term_lat = sympy.latex(self.variables['group2_term'])
        self.steps.append(("math", rf"{group1_trans_lat} + {group2_term_lat} = 0"))

        common_factor_lat = sympy.latex(self.variables['common_factor'])
        other_group_lat = sympy.latex(self.variables['other_group_expr'])

        self.steps.append(("text", rf"Винесемо спільний множник ${common_factor_lat}$ за дужки:"))
        self.steps.append(("math", rf"{common_factor_lat} \cdot \left( {other_group_lat} \right) = 0"))

        self._build_steps_solve_factors()

    def _build_steps_solve_factors(self):
        sub_sols = self.variables['sub_solutions']
        factor_index = 0
        n_latex = symbols('n')

        self.steps.append(("text", "Добуток дорівнює нулю, коли хоча б один із множників дорівнює нулю."))

        for i, factor_obj in enumerate(self.variables['all_factors']):
            if factor_obj.is_number: continue

            sub_eq = Eq(factor_obj, 0)

            if factor_obj.func == Add:
                self.steps.append(("text", f"Розв'язуємо {factor_index + 1}-й множник:"))
                self.steps.append(("math", sympy.latex(sub_eq)))

                if factor_index < len(sub_sols):
                    sub_sol_set = sub_sols[factor_index]
                    sub_sol_latex = sympy.latex(sub_sol_set)

                    if "arccos" in sub_sol_latex or "arcsin" in sub_sol_latex:
                        t_eq = Eq(factor_obj.args[1], -factor_obj.args[0])
                        sub_sol_latex = sympy.latex(solveset(t_eq, self.x, domain=Reals))

                    self.steps.append(("text", "Отримуємо розв'язок:"))
                    self.steps.append(("math", rf"x \in {sub_sol_latex}"))

                factor_index += 1
                continue

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
                    formula_str = rf"{sympy.latex(arg_expr)} = {formula_latex}, n \in \mathbb{{Z}}"
                    self.steps.append(("math", formula_str))

            sub_sol_set = sub_sols[factor_index]
            sub_sol_latex = sympy.latex(sub_sol_set)

            if not is_simple_arg:
                self.steps.append(("text", "Виражаємо $x$:"))

            self.steps.append(("math", rf"x \in {sub_sol_latex}"))

            factor_index += 1

        final_sol_latex = sympy.latex(self.solution_obj)
        self.steps.append(("text", "Об'єднуючи всі розв'язки, отримуємо кінцеву відповідь:"))

        self.steps.append(("math", rf"x \in {final_sol_latex}"))
