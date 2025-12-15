from .base_class import FixedDocument
from pylatex import Document, Section, Math, Package
from pylatex.utils import NoEscape


class EquationSet:

    def __init__(self):
        self.equations = []

    def add_equations(self, type_key: str, count: int = 1):
        from . import EQUATION_REGISTRY
        klass = EQUATION_REGISTRY.get(type_key)

        if klass is None:
            print(f"Попередження: Тип рівняння '{type_key}' не знайдено у EQUATION_REGISTRY.")
            return

        for _ in range(count):
            try:
                self.equations.append(klass())
            except Exception as e:
                print(f"Помилка при генерації класу {klass.__name__}: {e}")

    def clear(self):
        self.equations = []

    def generate_pdf(self, filename: str):
        geometry_options = {"tmargin": "1in", "lmargin": "1in", "rmargin": "1in"}

        doc = FixedDocument(
            documentclass='extarticle',
            document_options='14pt',
            geometry_options=geometry_options,
            fontenc=None,
            inputenc=None
        )

        doc.packages.add(Package('amsmath'))
        doc.packages.add(Package('amssymb'))
        doc.packages.add(Package('inputenc', options=['utf8']))
        doc.packages.add(Package('fontenc', options=['T2A']))
        doc.packages.add(Package('babel', options=['ukrainian']))

        for i, eq in enumerate(self.equations, 1):

            section_title = f"Завдання {i}"
            with doc.create(Section(section_title, numbering=False)):

                doc.append(NoEscape("Розв'яжіть рівняння:"))
                doc.append(Math(data=NoEscape(eq.get_equation_latex()), escape=False))
                doc.append(NoEscape(r"\vspace{20pt}"))

                doc.append(NoEscape(r"\textbf{Хід розв\'язання}"))
                doc.append(NoEscape(r"\\ \vspace{10pt}"))

                for step_type, step_data in eq.steps:
                    if step_type == "text":
                        doc.append(NoEscape(step_data))
                        doc.append(NoEscape(r"\\ \vspace{5pt}"))
                    elif step_type == "math":
                        doc.append(Math(data=NoEscape(step_data), escape=False))
                        doc.append(NoEscape(r"\\ \vspace{5pt}"))

                doc.append(NoEscape(r"\\ \vspace{10pt}"))
                doc.append(NoEscape(r"\textbf{Відповідь:}"))
                doc.append(Math(data=NoEscape(eq.get_solution_latex()), escape=False))

            if i < len(self.equations):
                doc.append(NoEscape(r"\newpage"))

        try:
            if filename.endswith('.pdf'):
                filename = filename[:-4]

            compiler_path = r"D:\apps\latex\texlive\2025\bin\windows\pdflatex.exe"

            doc.generate_pdf(
                filename,
                clean_tex=True,
                compiler=compiler_path
            )

            print(f"PDF-файл '{filename}.pdf' успішно створено.")

        except Exception as e:
            print(f"Помилка при генерації PDF: {e}")