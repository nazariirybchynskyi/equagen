import abc
import sympy
from pylatex import Document, Section, Math, Package
from pylatex.utils import NoEscape
import subprocess
import os
import shutil


class FixedDocument(Document):

    def generate_tex(self, file_name):
        tex_file = file_name + '.tex'
        with open(tex_file, 'w', encoding='utf-8') as f:
            f.write(self.dumps())

    def compile(self, file_name, clean_tex=True, clean=True,
                compiler=None, compiler_args=None, silent=True):

        tex_file = file_name + '.tex'
        pdf_file = file_name + '.pdf'
        log_file = file_name + '.log'
        aux_file = file_name + '.aux'

        if compiler is None:
            if shutil.which("latexmk") is not None:
                compiler = "latexmk"
            elif shutil.which("pdflatex") is not None:
                compiler = "pdflatex"
            else:
                raise (Exception("No LaTex compiler was found"))

        if compiler_args is None:
            compiler_args = []

        if compiler == 'latexmk':
            compiler_args = ['-pdf'] + compiler_args

        args = [compiler] + compiler_args + [tex_file]

        try:
            FNULL = open(os.devnull, 'w')
            stdout = FNULL if silent else None

            p = subprocess.run(args, stdout=stdout, stderr=subprocess.STDOUT,
                               cwd=self.temp_dir)
        finally:
            FNULL.close()

        with open(log_file, 'r', encoding='latin-1') as f:
            log = f.read()

        if p.returncode != 0:
            raise (Exception(
                'This is pdfTeX, Version 3.14159265-2.6-1.40.18 '
                '(TeX Live 2017)\n' + log
            ))

        if clean_tex:
            os.remove(tex_file)

        if clean:
            for ext in [aux_file, log_file, '.out', '.fls',
                        '.fdb_latexmk']:
                try:
                    os.remove(file_name + ext)
                except FileNotFoundError:
                    pass


class TrigonometricEquation(abc.ABC):

    def __init__(self):
        self.x = sympy.symbols('x')
        self.equation_obj = None
        self.solution_obj = None
        self.variables = {}
        self.steps = []

        self._generate()
        self._solve()
        self._build_solution_steps()

    @abc.abstractmethod
    def _generate(self):
        pass

    @abc.abstractmethod
    def _solve(self):
        pass

    @abc.abstractmethod
    def _build_solution_steps(self):
        pass

    def get_equation_latex(self) -> str:
        if self.equation_obj is not None:
            return sympy.latex(self.equation_obj)
        return "Рівняння не згенеровано."

    def get_solution_latex(self) -> str:
        if self.solution_obj is not None:
            return sympy.latex(self.solution_obj)
        return "Розв'язок не знайдено."