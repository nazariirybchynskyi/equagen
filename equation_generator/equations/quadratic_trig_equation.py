import random
import sympy
from sympy import sin, cos, tan, cot, pi, Eq, solveset, Reals, symbols, ImageSet, Lambda, Integers, EmptySet, Add, gcd, \
    Rational, sqrt, Mul, expand

from ..base_class import TrigonometricEquation


class QuadraticTrigEquation(TrigonometricEquation):

    def __init__(self):
        n_int = symbols('n', integer=True)
        self.general_formula_map = {
            'sin': r't = (-1)^n \arcsin(val) + \pi n, n \in \mathbb{Z}',
            'cos': r't = \pm \arccos(val) + 2\pi n, n \in \mathbb{Z}',
            'tan': r't = \arctan(val) + \pi n, n \in \mathbb{Z}',
            'cot': r't = \text{arcctg}(val) + \pi n, n \in \mathbb{Z}'
        }

        self.special_case_map = {
            'sin': {
                0: pi * n_int,
                1: pi / 2 + 2 * pi * n_int,
                -1: -pi / 2 + 2 * pi * n_int
            },
            'cos': {
                0: pi / 2 + pi * n_int,
                1: 2 * pi * n_int,
                -1: pi + 2 * pi * n_int
            },
            'tan': {
                0: pi * n_int
            },
            'cot': {
                0: pi / 2 + pi * n_int
            }
        }

        self.identity_map = {
            'sin': (cos(symbols('alpha')) ** 2, r'\cos^2(\alpha) = 1 - \sin^2(\alpha)'),
            'cos': (sin(symbols('alpha')) ** 2, r'\sin^2(\alpha) = 1 - \cos^2(\alpha)')
        }

        super().__init__()

    def _format_polynomial_latex(self, A, B, C, var_name):
        # Допоміжна функція для форматування коефіцієнта
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

        latex_str = " ".join(terms)
        if not latex_str.startswith('-'):
            latex_str = latex_str.lstrip('+ ')

        return latex_str + " = 0"

    def _generate(self):
        path_type = random.choice(["direct", "reducible"])

        if path_type == "direct":
            f_target = random.choice([sin, cos, tan, cot])
        else:  # reducible
            f_target = random.choice([sin, cos])

        f_name = f_target.__name__

        k = random.choice([1, 1, 1, 2, 3])
        b = random.choice([0, 0, 0, pi / 6, pi / 4])
        arg = k * self.x + b

        # Генерація коренів
        while True:
            if f_name in ('sin', 'cos'):
                nice_roots = [
                    1, -1,
                    Rational(1, 2), Rational(-1, 2),
                    sqrt(2) / 2, -sqrt(2) / 2,
                    sqrt(3) / 2, -sqrt(3) / 2
                ]
                # Прибираємо 0 зі списку nice_roots, щоб C != 0

                trap_roots = [2, -2, 3, Rational(3, 2)]

                t1 = random.choice(nice_roots)
                t2 = random.choice(nice_roots + trap_roots)

                if (abs(t1) > 1) and (abs(t2) > 1):
                    continue
            else:  # tan, cot
                nice_roots = [
                    1, -1,
                    sqrt(3), -sqrt(3),
                    1 / sqrt(3), -1 / sqrt(3)
                ]
                # Прибираємо 0

                t1 = random.choice(nice_roots)
                t2 = random.choice(nice_roots)

            # --- НОВА ПЕРЕВІРКА ---
            # 1. Щоб B != 0, сума коренів не має бути 0 (t1 != -t2)
            if t1 + t2 == 0:
                continue

            # 2. Щоб C != 0, жоден корінь не має бути 0
            # (Ми вже виключили 0 зі списків, але про всяк випадок)
            if t1 == 0 or t2 == 0:
                continue
            # ----------------------

            break

        # Створюємо ЯДРО
        a_raw = 1
        b_raw = -(t1 + t2)
        c_raw = t1 * t2

        coeffs_raw = [a_raw, b_raw, c_raw]
        denoms = [c.q for c in coeffs_raw if hasattr(c, 'q')]
        lcm = sympy.lcm(denoms) if denoms else 1

        mult = random.choice([1, -1])

        A_kernel = a_raw * lcm * mult
        B_kernel = b_raw * lcm * mult
        C_kernel = c_raw * lcm * mult

        A_kernel = sympy.simplify(A_kernel)
        B_kernel = sympy.simplify(B_kernel)
        C_kernel = sympy.simplify(C_kernel)

        if hasattr(A_kernel, 'is_Integer') and A_kernel.is_Integer: A_kernel = int(A_kernel)
        if hasattr(B_kernel, 'is_Integer') and B_kernel.is_Integer: B_kernel = int(B_kernel)
        if hasattr(C_kernel, 'is_Integer') and C_kernel.is_Integer: C_kernel = int(C_kernel)

        t_func = f_target(arg)

        self.variables = {
            'A_kernel': A_kernel, 'B_kernel': B_kernel, 'C_kernel': C_kernel,
            't1': t1, 't2': t2,
            'f_target': f_target, 'arg': arg,
            'path_type': path_type
        }

        if path_type == "direct":
            self.equation_obj = Eq(A_kernel * t_func ** 2 + B_kernel * t_func + C_kernel, 0)
        else:
            f_replace_symbol, _ = self.identity_map[f_name]
            f_replace_func = f_replace_symbol.subs(symbols('alpha'), arg)

            # Перевіряємо, чи не зникнуть коефіцієнти після заміни sin^2 -> 1-cos^2
            # A*sin^2 + B*sin + C = 0  --> A(1-cos^2) + ...
            # -A*cos^2 + B*cos + (A+C) = 0
            # Нам треба, щоб нові коефіцієнти (A_orig, B_orig, C_orig) теж не були нулями

            A_orig = -A_kernel  # A_orig != 0, бо A_kernel != 0
            B_orig = B_kernel  # B_orig != 0, бо B_kernel != 0
            C_orig = C_kernel + A_kernel  # C_orig має бути != 0

            # Якщо C_orig стає нулем, то початкове рівняння буде неповним.
            # C_kernel + A_kernel = 0 => t1*t2 + 1 = 0 (після нормування) => t1*t2 = -1
            # Це стається, коли корені взаємно обернені і протилежні за знаком? Ні.
            # Це стається при певних комбінаціях. Якщо сталося - регенеруємо.
            if C_orig == 0:
                return self._generate()

            self.equation_obj = Eq(A_orig * f_replace_func + B_orig * t_func + C_orig, 0)

            replacement_expr = (1 - t_func ** 2)
            term_with_parentheses = Mul(A_orig, replacement_expr, evaluate=False)
            lhs_initial = Add(term_with_parentheses, B_orig * t_func, C_orig, evaluate=False)

            self.variables['initial_replacement'] = Eq(lhs_initial, 0)

            term_expanded = expand(A_orig * replacement_expr)
            terms = []
            if term_expanded.is_Add:
                terms.extend(term_expanded.args)
            else:
                terms.append(term_expanded)

            terms.append(B_orig * t_func)
            terms.append(C_orig)

            lhs_expanded = Add(*terms, evaluate=False)
            self.variables['expanded_replacement'] = Eq(lhs_expanded, 0)

            self.variables['A_orig'] = A_orig
            self.variables['B_orig'] = B_orig
            self.variables['C_orig'] = C_orig

    def get_equation_latex(self) -> str:
        return sympy.latex(self.equation_obj, order='none')

    def _solve(self):
        t1 = self.variables['t1']
        t2 = self.variables['t2']
        f = self.variables['f_target']
        arg = self.variables['arg']
        f_name = f.__name__
        n_int = symbols('n', integer=True)

        roots = [t1]
        if t2 != t1:
            roots.append(t2)

        final_solution_set = EmptySet
        sub_solutions_list = []

        for t_val in roots:
            is_trap = (f_name in ('sin', 'cos') and abs(t_val) > 1)

            if is_trap:
                sub_sol_set = EmptySet
            else:
                sub_sol_set = self._solve_sub_equation(f, arg, t_val, n_int)

            sub_solutions_list.append(sub_sol_set)
            final_solution_set = final_solution_set.union(sub_sol_set)

        self.variables['sub_solutions'] = sub_solutions_list
        self.solution_obj = final_solution_set

    def _solve_sub_equation(self, f, arg, t_val, n_int):
        f_name = f.__name__

        special_formula_expr = self.special_case_map.get(f_name, {}).get(t_val)

        if special_formula_expr is not None:
            simple_eq = Eq(arg, special_formula_expr)
            sub_sol_expr = sympy.solve(simple_eq, self.x)[0].expand()
            return ImageSet(Lambda(n_int, sub_sol_expr), Integers)
        else:
            return solveset(Eq(f(arg), t_val), self.x, domain=Reals)

    def _build_solution_steps(self):
        A_k = self.variables['A_kernel']
        B_k = self.variables['B_kernel']
        C_k = self.variables['C_kernel']
        t1 = self.variables['t1']
        t2 = self.variables['t2']
        f = self.variables['f_target']
        arg = self.variables['arg']
        path_type = self.variables['path_type']
        sub_solutions = self.variables['sub_solutions']

        f_name = f.__name__
        arg_latex = sympy.latex(arg)
        t = symbols('t')
        n_latex = symbols('n')
        is_simple_arg = (arg == self.x)

        self.steps.append(("text", f"Маємо рівняння: ${self.get_equation_latex()}$"))

        f_name_latex = f"{f_name}({arg_latex})"

        if path_type == "reducible":
            _, identity_latex = self.identity_map[f_name]
            identity_latex_applied = identity_latex.replace(r'\alpha', arg_latex)

            self.steps.append(
                ("text", f"Рівняння містить дві різні функції. Зведемо його до однієї, використавши тотожність:"))
            self.steps.append(
                ("math", rf"\sin^2({arg_latex}) + \cos^2({arg_latex}) = 1 \implies {identity_latex_applied}"))

            self.steps.append(("text", "Підставляємо в рівняння:"))
            self.steps.append(("math", sympy.latex(self.variables['initial_replacement'], order='none')))

            self.steps.append(("text", "Розкриваємо дужки:"))
            self.steps.append(("math", sympy.latex(self.variables['expanded_replacement'], order='none')))

            self.steps.append(
                ("text", f"Зводимо подібні доданки та отримуємо квадратне рівняння відносно ${f_name_latex}$:"))

        final_quadratic_latex = self._format_polynomial_latex(A_k, B_k, C_k, f_name_latex)

        if path_type == "direct":
            self.steps.append(("text", f"Це квадратне рівняння відносно ${f_name_latex}$."))

        self.steps.append(("math", final_quadratic_latex))

        self.steps.append(("text", f"Введемо заміну $t = {f_name_latex}$:"))
        t_eq_latex = self._format_polynomial_latex(A_k, B_k, C_k, "t")
        self.steps.append(("math", t_eq_latex))

        self.steps.append(("text", "Розв'язуємо квадратне рівняння:"))

        roots_latex = f"t_1 = {sympy.latex(t1)}"
        if t2 != t1:
            roots_latex += rf", \quad t_2 = {sympy.latex(t2)}"
        self.steps.append(("math", roots_latex))

        self.steps.append(("text", "Повертаємось до заміни:"))

        unique_roots = [t1]
        if t2 != t1:
            unique_roots.append(t2)

        for i, t_val in enumerate(unique_roots):
            sol_set = sub_solutions[i]
            self._build_steps_for_sub_equation(t_val, f_name, arg_latex, is_simple_arg, sol_set)

        if len(sub_solutions) > 0 and self.solution_obj != EmptySet:
            final_sol_latex = sympy.latex(self.solution_obj)
            self.steps.append(("text", "Об'єднуючи всі розв'язки, отримуємо кінцеву відповідь:"))
            self.steps.append(("math", f"x = {final_sol_latex}"))

    def _build_steps_for_sub_equation(self, t_val, f_name, arg_latex, is_simple_arg, sol_set):
        n_latex = symbols('n')

        self.steps.append(("text", f"Розв'язуємо корінь $t = {sympy.latex(t_val)}$:"))

        is_trap = (f_name in ('sin', 'cos') and (t_val > 1 or t_val < -1))

        if is_trap:
            self.steps.append(("text",
                               rf"Корінь $t = {sympy.latex(t_val)}$ не належить області значень $\left[-1, 1\right]$ для функції ${f_name}$. Розв'язків немає."))
            return

        self.steps.append(("math", f"{f_name}({arg_latex}) = {sympy.latex(t_val)}"))

        special_formula_expr = self.special_case_map.get(f_name, {}).get(t_val)

        formula_variable_latex = "x" if is_simple_arg else arg_latex

        if special_formula_expr:
            self.steps.append(("text", "Це окремий випадок. Використовуємо спрощену формулу:"))
            formula_latex = sympy.latex(special_formula_expr.subs(symbols('n', integer=True), n_latex))
            formula_str = f"{formula_variable_latex} = {formula_latex}, n " + r"\in \mathbb{Z}"
            self.steps.append(("math", formula_str))
        else:
            self.steps.append(("text", rf"Застосуємо загальну формулу розв'язку для ${f_name}$:"))
            val_latex = sympy.latex(t_val)
            raw_formula = self.general_formula_map[f_name]
            formula_with_val = raw_formula.replace('val', val_latex)
            final_formula = formula_with_val.replace('t', formula_variable_latex)
            self.steps.append(("math", final_formula))

        if not is_simple_arg:
            self.steps.append(("text", "Виражаємо $x$:"))

        sub_sol_latex = sympy.latex(sol_set)
        self.steps.append(("math", f"x = {sub_sol_latex}"))