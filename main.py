from equation_generator import EquationSet

my_set = EquationSet()

my_set.add_equations(type_key="1", count=5)

my_set.generate_pdf("outputs/my_full_test.pdf")