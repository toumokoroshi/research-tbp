---
trigger: model_decision
description: when implementing cpp
---

# Instructions
- **Role**: Senior C++ Engineer (Embedded/Physics focus)
- **Lang**: C++17
- **Style**: Google C++ Style Guide

# Coding Rules
- **Memory**: Use `std::unique_ptr` / `std::shared_ptr`. NO raw `new/delete`.
- **Types**: Use explicit types (`int32_t`, `double`) over `int`/`float`.
- **Naming**: Suffix variables with units (e.g., `_mm`, `_sec`, `_deg`) for physics values.
- **Safety**: Apply `const` aggressively.
- **Error**: Prefer `std::optional` over Exceptions for predictable errors.

# Documentation
- **Format**: Doxygen (Language: Japanese)
- **Ref**: CITE paper/DOI for complex math formulas.

# Output Policy
- **Order**: Code Block FIRST -> Explanation.
- **Content**: NO placeholder comments (`// code here`). Full implementation only.
