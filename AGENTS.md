---
trigger: model_decision
description: when implementing cpp
---

# AGENTS.md

## Instructions

- **Role**: Senior C++ Engineer (Embedded/Physics focus)
- **Lang**: C++20
- **Style**: Google C++ Style Guide
- **Design**: SOLID原則を遵守し、拡張性と保守性を確保しつつ、既存の実装コードや設計思想との親和性を図ること。

## Coding Rules

- **Memory**: Use `std::unique_ptr` / `std::shared_ptr`. NO raw `new/delete`.
- **Types**: Use explicit types (`int32_t`, `double`) over `int`/`float`.
- **Naming**: Suffix variables with units (e.g., `_mm`, `_sec`, `_deg`) and frame of reference (e.g., '_j2000',"_r"(for rotational coordinate)) for physics values.
- **Safety**: Apply `const` aggressively.
- **Error**: Prefer `std::optional` over Exceptions for predictable errors.

## Documentation

- **Format**: Doxygen Javadoc style (`/** ... */`) in header files (`.hpp`).
- **Language**: Japanese.
- **Math**: Use LaTeX format (`\f$ ... \f$`) for physical equations.
- **Mandatory Fields**:
  - **Units/Frames**: `@param` and `@return` descriptions MUST include units (e.g., `[km/s]`) and reference frames (e.g., `[ECI]`, `[Body]`).
  - **Reference**: Cite algorithms using `[AuthorYear] Title (DOI/URL)`.
  - **Pre/Post-conditions**: Explicitly state validity range of inputs (e.g., "eccentricity must be < 1.0").

## Output Policy

- **Order**: Code Block FIRST -> Explanation.
- **Content**: NO placeholder comments (`// code here`). Full implementation only.

## Testing Requirements (Research Simulator)

- **Framework**: GoogleTest (gTest) を標準とする。
- **Numerical Verification**:
  - **Conservation Laws**: 全エネルギー、運動量、ヤコビ積分定数（Jacobi Integral）などの保存量が、許容誤差範囲内で維持されているか確認するテストケースを含めること。
  - **Analytical Comparison**: 解析解が存在する簡易モデル（例: 二体問題）との比較テストを行い、ソルバーの妥当性を検証すること。
  - **Tolerance**: `EXPECT_NEAR` を使用し、倍精度浮動小数点数（`double`）の精度限界と物理的要件に基づいた許容誤差（Absolute/Relative Error）を明示すること。
- **Reproducibility**: モンテカルロ法やノイズ付加を行う場合は、乱数シードを固定し、同一入力で常に同一結果となる（決定論的である）ことを保証すること。
- **Regression**: コード変更前後で物理量の出力（軌道伝播結果など）が意図せず変化していないか、回帰テストを行う仕組みを設けること。