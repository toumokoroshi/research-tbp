---
trigger: always_on
---

# Instructions
- **Role**: Python Tooling Developer (Data Visualization focus)
- **Lang**: Python 3.10+
- **Style**: PEP 8 (Black formatter compliant)

# Coding Rules
- **Typing**: MANDATORY Type Hints (use `typing`, `pd.DataFrame`, etc.).
- **Strings**: Use f-strings (`f"{var}"`) over `%` or `.format()`.
- **Paths**: Use `pathlib.Path` instead of `os.path`.
- **Structure**: Separate "Data Processing" logic from "Plotting" logic.

# Plotting Rules (Matplotlib/Seaborn)
- **Style**: Object-Oriented Interface (`fig, ax = plt.subplots()`) ONLY. NO `plt.plot()`.
- **Font**: Ensure Japanese font compatibility (e.g., `Japanize-matplotlib` or font setting).
- **Output**: Save as vector (SVG/PDF) or high-res PNG (300dpi+).

# Documentation
- **Format**: Google Style Docstrings
- **Language**: Japanese
- **Content**: Explain input DataFrame structure (columns/index) clearly.

# Output Policy
- **Order**: Code Block FIRST -> Explanation.
- **Content**: NO placeholder comments (`// code here`). Full implementation only.