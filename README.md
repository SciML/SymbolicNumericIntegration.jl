# SymbolicNumericIntegration.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/SymbolicNumericIntegration/stable/)

[![codecov](https://codecov.io/gh/SciML/SymbolicNumericIntegration.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/SymbolicNumericIntegration.jl)
[![Build Status](https://github.com/SciML/SymbolicNumericIntegration.jl/workflows/CI/badge.svg)](https://github.com/SciML/SymbolicNumericIntegration.jl/actions?query=workflow%3ACI)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

**SymbolicNumericIntegration.jl** is a hybrid symbolic/numerical integration package that works on the [Julia Symbolics](https://github.com/JuliaSymbolics/Symbolics.jl) expressions.

## Tutorials and Documentation

For information on using the package,
[see the stable documentation](https://docs.sciml.ai/SymbolicNumericIntegration/stable/). Use the
[in-development documentation](https://docs.sciml.ai/SymbolicNumericIntegration/dev/) for the version of
the documentation, which contains the unreleased features.

## Example

```julia
julia> using SymbolicNumericIntegration
julia> using Symbolics

julia> @variables x a b

julia> integrate(3x^3 + 2x - 5)
(x^2 + (3//4)*(x^4) - (5//1)*x, 0, 0)

julia> integrate(exp(a * x), x; symbolic = true)
(exp(a*x) / a, 0, 0)

julia> integrate(sin(a * x) * cos(b * x), x; symbolic = true, detailed = false)
(-a*cos(a*x)*cos(b*x) - b*sin(a*x)*sin(b*x)) / (a^2 - (b^2))
```

# Citation

If you use SymbolicNumericIntegration.jl in your research, please cite [this paper](https://arxiv.org/abs/2201.12468):

```

@article{Iravanian2022,
author = {Shahriar Iravanian and Carl Julius Martensen and Alessandro Cheli and Shashi Gowda and Anand Jain and Julia Computing and Yingbo Ma and Chris Rackauckas},
doi = {10.48550/arxiv.2201.12468},
month = {1},
title = {Symbolic-Numeric Integration of Univariate Expressions based on Sparse Regression},
url = {https://arxiv.org/abs/2201.12468v2},
year = {2022},
}

```

```
```
