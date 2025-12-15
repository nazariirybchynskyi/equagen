from base_class import TrigonometricEquation
import sympy
from sympy import sin, cos, tan, cot, pi, Eq, solveset, Reals, oo, zoo, nan, symbols
import random

# Кількість унікальних комбінацій параметрів для генерації: 24
# (6 варіантів t1 * 4 варіанти t2)
class HomogeneousEquation(TrigonometricEquation):

    def _generate(self):
        t1 = random.choice([-3, -2, -1, 1, 2, 3])
        t2 = random.choice([-1, 1, 1 / 2, 2])

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

        self.steps = [
            ("text", f"Маємо рівняння: ${sympy.latex(self.equation_obj)}$"),
            ("text", "Це однорідне тригонометричне рівняння другого порядку."),
            ("text", rf"Розділимо обидві частини рівняння на $cos^2(x) \neq 0$:"),
            ("math", f"{A}\\frac{{sin^2(x)}}{{cos^2(x)}} + {B}\\frac{{sin(x)cos(x)}}{{cos^2(x)}} + {C}\\frac{{cos^2(x)}}{{cos^2(x)}} = 0"),
            ("text", f"Отримаємо: ${A}tg^2(x) + {B}tg(x) + {C} = 0$"),
            ("text", f"Введемо заміну $t = tg(x)$. Рівняння набуває вигляду:"),
            ("math", f"{A}t^2 + {B}t + {C} = 0"),
            ("text", f"Коренями цього квадратного рівняння є $t_1 = {t1}$ та $t_2 = {t2}$."),
            ("text", "Повертаємось до заміни:"),
            ("math", rf"1) $tg(x) = {t1} \implies x = {sympy.latex(solveset(Eq(tan(self.x), self.variables['t1']), self.x))}$"),
            ("math", rf"2) $tg(x) = {t2} \implies x = {sympy.latex(solveset(Eq(tan(self.x), self.variables['t2']), self.x))}$"),
            ("text", f"Об'єднуючи розв'язки, отримуємо кінцеву відповідь.")
        ]

# Кількість унікальних комбінацій параметрів для генерації: 2160
class SimplestEquation(TrigonometricEquation):

    def __init__(self):
        self.general_formula_map = {
            'sin': r't = (-1)^n \arcsin(val) + \pi n, n \in \mathbb{Z}',
            'cos': r't = \pm \arccos(val) + 2\pi n, n \in \mathbb{Z}',
            'tan': r't = \arctan(val) + \pi n, n \in \mathbb{Z}',
            'cot': r't = \text{arcctg}(val) + \pi n, n \in \mathbb{Z}'
        }
        super().__init__()

    def _generate(self):
        while True:
            f = random.choice([sin, cos, tan, cot])
            k = random.randint(1, 3)
            b = random.choice([0, pi / 6, pi / 4, pi / 3, pi / 2])
            A = random.choice([1, 2, 3, -1])

            nice_angles = [0, pi / 6, pi / 4, pi / 3, pi / 2, 2 * pi / 3, 3 * pi / 4, 5 * pi / 6, pi]
            base_solution = random.choice(nice_angles)

            argument = k * base_solution + b

            a_value = A * f(argument)

            if a_value.has(oo, zoo, nan):
                continue

            a_simplified = sympy.simplify(a_value)

            self.variables = {'A': A, 'k': k, 'b': b, 'a': a_simplified, 'f': f}

            self.equation_obj = Eq(A * f(k * self.x + b), a_simplified)

            break

    def _solve(self):
        if self.equation_obj is not None:
            self.solution_obj = solveset(self.equation_obj, self.x, domain=Reals)

    def _build_solution_steps(self):
        if self.equation_obj is None:
            return

        A = self.variables['A']
        k = self.variables['k']
        b = self.variables['b']
        a = self.variables['a']
        f = self.variables['f']
        f_name = f.__name__

        t = symbols('t')

        self.steps.append(("text", f"Маємо рівняння:"))
        self.steps.append(("math", sympy.latex(self.equation_obj)))

        if A != 1:
            simple_eq_obj = Eq(f(k * self.x + b), a / A)
            self.steps.append(("text", f"Розділимо обидві частини на ${A}:"))
            self.steps.append(("math", sympy.latex(simple_eq_obj)))
        else:
            simple_eq_obj = self.equation_obj

        argument_expression = simple_eq_obj.lhs.args[0]
        t_equation = Eq(f(t), simple_eq_obj.rhs)

        self.steps.append(
            ("text", f"Це найпростіше тригонометричне рівняння. Введемо заміну $t = {sympy.latex(argument_expression)}$:"))
        self.steps.append(("math", sympy.latex(t_equation)))

        value_latex = sympy.latex(t_equation.rhs)
        general_formula = self.general_formula_map[f_name].replace('val', value_latex)

        self.steps.append(("text", rf"Загальний розв'язок для ${f_name}$ має вигляд:"))
        self.steps.append(("math", general_formula))

        t_solution = sympy.solveset(t_equation, t, domain=Reals)
        self.steps.append(("text", f"Підставляємо наше значення, отримуємо: $t = {sympy.latex(t_solution)}$"))

        self.steps.append(("text", f"Повертаємось до заміни $t = {sympy.latex(argument_expression)}$:"))
        self.steps.append(("math", f"{sympy.latex(argument_expression)} = {sympy.latex(t_solution)}"))

        final_solution_latex = sympy.latex(self.solution_obj)
        self.steps.append(("text", f"Виражаємо $x$ та отримуємо кінцеву відповідь:"))
        self.steps.append(("math", f"$x = {final_solution_latex}$"))

# Ваш тестовий код
eq = HomogeneousEquation()
print(f"Рівняння: {eq.get_equation_latex()}")
print(f"Відповідь: {eq.get_solution_latex()}")
print("\n--- Генеруємо PDF (симуляція) ---")
eq.generate_pdf("homogeneous_solution.pdf")