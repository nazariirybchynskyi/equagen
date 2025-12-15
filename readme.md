# Trigonometric Equation Generator 📐

A powerful Python-based system for automated generation of trigonometric equations with detailed, step-by-step solutions.

This project utilizes a **"generative reverse engineering"** approach. Instead of randomizing coefficients (which often leads to "ugly" roots), the system first selects didactically appropriate roots and parameters, then reconstructs the equation. The final output is rendered into high-quality PDF documents using LaTeX.

## 📋 Prerequisites

To run this project, you need the following installed on your system:

1.  **Python 3.8+**
2.  **LaTeX Distribution (CRITICAL)**

⚠️ **IMPORTANT:** The system generates PDFs by compiling `.tex` files. You **must** have a LaTeX engine installed, and it **must be added to your system's PATH**.

* **Windows:** Install [MiKTeX](https://miktex.org/download) or [TeX Live](https://www.tug.org/texlive/). During installation, ensure you check the option **"Add to PATH"**.
* **Linux (Ubuntu/Debian):**
    ```bash
    sudo apt-get install texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-lang-cyrillic
    ```
* **macOS:** Install [MacTeX](https://www.tug.org/mactex/).

## 🛠 Installation

1.  Clone the repository:
    ```bash
    git clone [https://github.com/your-username/trig-generator.git](https://github.com/your-username/trig-generator.git)
    cd trig-generator
    ```

2.  Install Python dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## 🚀 Usage & Generation Process

The core of the system is the `EquationSet` class. It acts as an orchestrator that bridges the user request with the specific mathematical logic and the rendering engine.

### How `EquationSet` Works:

1.  **Initialization**: When you create an instance of `EquationSet`, you pass the target `equation_type_id` and the desired `count` (quantity).
2.  **Mapping**: The class dynamically maps the provided ID to a specific generator class (e.g., ID `"8"` $\rightarrow$ `LinearCombinationEquation`).
3.  **Generation Loop**: It runs a generation loop `count` times. Inside this loop, the specific generator class uses the "reverse engineering" algorithm to create a valid equation object with distinct roots and steps.
4.  **Collection**: All generated equation objects (containing LaTeX strings for the problem statement and the solution) are collected into a list.
5.  **Rendering**: The `to_pdf()` method injects these LaTeX strings into a Jinja2 template and calls the system's `pdflatex` (or equivalent) to build the final PDF.

### Code Example:

```python
from equation_generator import EquationSet

my_set = EquationSet()

my_set.add_equations(type_key="1", count=5)

my_set.generate_pdf("outputs/my_full_test.pdf")

print("PDF successfully generated!")
```
## 🧮 Equation Types (ID Reference)

Use the following IDs when using `add_equations` on `EquationSet` to select the specific equation type.

| ID | Class Name | Опис (Метод розв'язання) |
| :--- | :--- | :--- |
| **"1"** | `SimplestEquation` | Найпростіші тригонометричні рівняння ($\sin x = a$ тощо) |
| **"2"** | `HomogeneousEquation` | Однорідні рівняння (зведення до $\tan x$) |
| **"3"** | `SumToProductEquation` | Перетворення суми/різниці в добуток |
| **"4"** | `GroupingEquation` | Метод групування доданків |
| **"5"** | `PowerReductionEquation` | Метод пониження степеня |
| **"6"** | `QuadraticTrigEquation` | Звідні до квадратних (заміна $t = f(x)$) |
| **"7"** | `DoubleAngleToQuadraticEquation` | Використання формул подвійного аргументу |
| **"8"** | `LinearCombinationEquation` | Метод введення допоміжного кута ($a \sin x + b \cos x = c$) |
| **"9"** | `ReducibleToHomogeneousEquation` | Рівняння, що зводяться до однорідних |
| **"10"** | `SymmetricEquation` | Симетричні рівняння (заміна $t = \sin x \pm \cos x$) |
| **"11"** | `TanSubstitutionEquation` | Універсальна тригонометрична підстановка |
| **"12"** | `SumTanCotanEquation` | Рівняння з $\tan x + \cot x$ |
| **"13"** | `BoundedSumEquation` | Метод оцінки (minimax, обмеженість функцій) |
| **"14"** | `InverseTrigEquation` | Рівняння з оберненими тригонометричними функціями |

## 💻 Technologies

* **Python 3**: Core application logic.
* **SymPy**: Computer Algebra System (CAS) used for symbolic math, root validation, and exact arithmetic.
* **LaTeX**: High-quality typesetting system for generating the final documents.
* **Jinja2**: Template engine used to construct dynamic LaTeX source files.

## 📄 License

This project is licensed under the **MIT License**.

```text
MIT License

Copyright (c) 2025 Nazarii Rybchynskyi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.