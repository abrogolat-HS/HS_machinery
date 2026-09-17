# HSmachinery

**A Wolfram Mathematica package for Higher-Spin (HS) generating systems, Sp(2) spinor algebra, and phase-space star products.**

`HSmachinery` is a symbolic computation engine designed for 2-component spinor calculus, specifically tailored for Higher-Spin gauge theory (e.g., Vasiliev equations). It provides robust tools for manipulating Sp(2) invariants, performing multidimensional oscillatory Gaussian integrations over spinor pairs, computing Moyal star products of phase-space symbols, and evaluating shifted homotopy operators.

## ✨ Features

- **Sp(2) Invariant Algebra**: Automatic handling of antisymmetric spinor convolutions (`conv`) and explicit components (`comp`), with intelligent scalar extraction.
- **Oscillatory Gaussian Integration**: Exact integration of polynomials multiplied by $\exp(i \cdot \text{exponent})$ over 2-component spinor pairs, with automatic reduction of the result back to the Sp(2) invariant basis.
- **Higher-Spin Star Products**: Computes the full phase-space (Moyal) star product and the holomorphic star product using integral representations.
- **Shifted Homotopy Operators**: Evaluates the $h_q$ homotopy operator used to resolve differential equations in the auxiliary $Z$-spinor space.
- **Robust Invariant Reduction**: Uses advanced linear algebra to rigorously convert raw spinor component polynomials back into an elegant basis of `conv` invariants.

## 📦 Installation

1. Download or clone the repository.
2. Place `HSmachinery.m` in a directory recognized by Mathematica (e.g., your `$UserBaseDirectory/Applications/`).
3. Load the package in your Mathematica notebook:
   ```mathematica
   << HSmachinery`

## Example

```
(* 1. Declare spinor variables *)
spinors[z, y, zb, yb];

(* 2. Define two phase-space functions as {polynomial, exponent} *)
f1 = {R, conv[y, yb]};
f2 = {1, conv[z, y]};

(* 3. Compute their full phase-space star product *)
starProduct[f1, f2, {z, y}, {zb, yb}]
```
