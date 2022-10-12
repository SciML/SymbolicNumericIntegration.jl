# SymbolicNumericIntegration.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](http://symbolicnumericintegration.sciml.ai/stable/)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/dev/modules/SymbolicNumericIntegration/)

[![codecov](https://codecov.io/gh/SciML/SymbolicNumericIntegration.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/SymbolicNumericIntegration.jl)
[![Build Status](https://github.com/SciML/SymbolicNumericIntegration.jl/workflows/CI/badge.svg)](https://github.com/SciML/SymbolicNumericIntegration.jl/actions?query=workflow%3ACI)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

**SymbolicNumericIntegration.jl** is a hybrid symbolic/numerical integration package that works on the [Julia Symbolics](https://github.com/JuliaSymbolics/Symbolics.jl) expressions.

## Tutorials and Documentation

For information on using the package,
[see the stable documentation](https://symbolicnumericintegration.sciml.ai/stable/). Use the
[in-development documentation](https://symbolicnumericintegration.sciml.ai/dev/) for the version of
the documentation, which contains the unreleased features.

## Example

```julia
using Symbolics
using SymbolicNumericIntegration

@variables x

julia> integrate(3x^3 + 2x - 5)
(x^2 + (3//4)*(x^4) - (5x), 0, 0)
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
