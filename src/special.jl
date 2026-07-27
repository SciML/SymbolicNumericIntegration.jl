# Adopted from mpmath (https://mpmath.org/), which is a Python library for real and
# complex floating-point arithmetic with arbitrary precision

E1(z) = Complex(expint(z))

"""
    Ei(z)

Return the exponential integral evaluated at `z` as a `Complex` value.

## Arguments

  - `z`: Real or complex evaluation point.

## Examples

```julia
julia> Ei(1) isa Complex
true
```
"""
function Ei(z)
    if z isa Real
        return Complex(-expint(Complex(-z)))
    end

    v = -expint(-z)

    if imag(z) > 0
        v += Complex(0, π)
    else
        v -= Complex(0, π)
    end

    return v
end

"""
    Ci(z)

Return the cosine integral evaluated at `z` as a `Complex` value.

## Arguments

  - `z`: Real or complex evaluation point.

## Examples

```julia
julia> Ci(1) isa Complex
true
```
"""
function Ci(z)
    if z isa Real
        return Complex(cosint(z))
    end

    v = (Ei(z * im) + Ei(-z * im)) / 2

    if real(z) ≈ 0
        if imag(z) >= 0
            v += Complex(0, π / 2)
        else
            v -= Complex(0, π / 2)
        end
    end

    if real(z) < 0
        if imag(z) >= 0
            v += Complex(0, π)
        else
            v -= Complex(0, π)
        end
    end
    return v
end

"""
    Si(z)

Return the sine integral evaluated at `z` as a `Complex` value.

## Arguments

  - `z`: Real or complex evaluation point.

## Examples

```julia
julia> Si(1) isa Complex
true
```
"""
function Si(z)
    if z isa Real
        return Complex(sinint(z))
    end

    v = (Ei(z * im) - Ei(-z * im)) / 2im

    if real(z) ≈ 0
        return v
    end

    if real(z) > 0
        v -= π / 2
    else
        v += π / 2
    end

    return v
end

"""
    Li(z)

Return the logarithmic integral evaluated at `z` as a `Complex` value.

## Arguments

  - `z`: Real or complex evaluation point.

## Examples

```julia
julia> Li(3) isa Complex
true
```
"""
function Li(z)
    return Ei(log(z)) - Ei(log(2.0))
end

function Erfi(z)
    return erfi(z)
end
