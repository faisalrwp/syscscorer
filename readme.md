# 🌐 CyReScoF – Compromise Vector (C) Calculator
A Modern GUI Tool for System Dependency Analysis, Sensitivity Evaluation & Visualization
<p align="center"> <img src="https://img.shields.io/badge/Python-3.8%2B-blue?style=flat-square&logo=python" /> <img src="https://img.shields.io/badge/GUI-ttkbootstrap-success?style=flat-square" /> <img src="https://img.shields.io/badge/License-MIT-green?style=flat-square" /> <img src="https://img.shields.io/badge/Powered%20By-Numpy%20%7C%20Pandas-important?style=flat-square&logo=numpy" /> </p>

## 📘 Overview

CyReScoF (Cyber Resilience Score Calculation Framework) is a Python-based graphical tool that computes the Compromise Vector (C) using system dependency matrices.
It supports:

📊 Compute C vector

🔁 Automatic dependency matrix & score alignment

🧪 Full sensitivity analysis for β (propagation factor)

📈 Embedded graph visualization (C vs β)

💾 Export results & graphs

🎨 Theme customization using ttkbootstrap

This tool is designed for researchers, analysts, and engineers working with cyber-physical systems, network resilience modeling, or cascading failure analysis.

## 🚀 Features at a Glance
### 🧮 Compromise Vector (C) Calculator

Solve: (I − βA) C = (1 − S)

Automatic matrix & score validation

Diagnostic outputs:

Spectral radius

Determinant

Condition number

Regularization fallback

### 🔁 Sensitivity Analysis (β-min → β-max)

User-defined:

- β minimum
- β maximum
- β step size

Computes C for each β value

Displays results in a sortable table

Exportable as CSV

### 📈 Sensitivity Analysis Graph

Interactive plotting using matplotlib

Multi-component selection

Visualizes C-value change vs β

Supports:
- PNG
- JPG
- PDF
- SVG

Perfect for presentations & reports

### 🎨 Theme Customization

Based on ttkbootstrap

Live theme switching

Beautiful modern UI styles

### 📁 Template Generation

One-click creation of CSV templates:

- Dependency Matrix (A)
- Component Scores (S)

## 📷 Screenshots
Coming Soon 

## Project Structure

```
📂 Project Structure
├── cyrescof_gui.py
├── settings.json
├── README.md
└── docs/
    └── images/
```

## 🔧 Installation

### ✔ 1. Clone the Repository
```
git clone https://github.com/your-username/CyReScoF.git
cd CyReScoF
```

### ✔ 2. Install Dependencies
```
pip install numpy pandas ttkbootstrap matplotlib
```

### ▶️ Running the Application
```
python cyrescof_gui.py
```

The GUI will start instantly.

## 🧠 How It Works

### 🟦 Input

A: Dependency matrix (n×n)

S: Component scores (n×1)

### 🟥 Process

CyReScoF solves:

(
𝐼
−
𝛽
𝐴
)
⋅
𝐶
=
(
1
−
𝑆
)
(I−βA)⋅C=(1−S)

with:

- Regularization fallback
- Eigen-based diagnostics
- Automatic alignment of S to A

### 🟩 Output

- C values
- Diagnostics
- Exportable CSV
- Graph plots
- Sensitivity matrix

## 📊 Sensitivity Analysis Formula

For β in:

[
𝛽
𝑚
𝑖
𝑛
,
𝛽
𝑚
𝑎
𝑥
]
  
step
  
Δ
𝛽
[β
min
	​

,β
max
	​

]stepΔβ

We compute:

𝐶
(
𝛽
)
=
(
𝐼
−
𝛽
𝐴
)
−
1
(
1
−
𝑆
)
C(β)=(I−βA)
−1
(1−S)

Results are:
- Listed in table
- Exportable
- Visualizable

## 💡 Use Cases
✓ Cyber Resilience Studies
✓ Risk Propagation Modeling
✓ Dependency Network Analysis
✓ Cascading Failure Simulation
✓ Academic Research & Teaching

## 🗂 Export Options

| Feature	           | Format |
|-|-|
| C Results	           | CSV |
| Sensitivity Results  | CSV |
| Graphs	           | PNG, JPG, PDF, SVG |

## 🛠 Tech Stack
| Component	    | Library |
|-|-|
| GUI	        | ttkbootstrap |
| Data Handling	| pandas |
| Math	        | numpy |
| Plotting	    | matplotlib |
| File Export	| csv |

## 🤝 Contributing

Pull requests are welcome!
Please ensure that:
- Existing functionality stays intact
- New features are modular
- UI/UX remains consistent

## 📜 License

This project is licensed under the MIT License.
You are free to use, modify, and distribute with attribution.

## ⭐ Support & Star the Repo

If this tool helped you:

👉 Please give the repository a ⭐ on GitHub!

It motivates further development.

## 📬 Contact

For feature requests or questions, please open an Issue on GitHub.