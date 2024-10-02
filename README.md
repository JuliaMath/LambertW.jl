# LambertW
### Lambert W function and associated omega constant

[![Build Status](https://github.com/JuliaMath/LambertW.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/JuliaMath/LambertW.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Codecov](https://codecov.io/gh/JuliaMath/LambertW.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaMath/LambertW.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)

### lambertw

The [Lambert W function](http://en.wikipedia.org/wiki/Lambert_W_function),
also called the omega function or product logarithm.

```julia
lambertw(z,k)   # Lambert W function for argument z and branch index k
lambertw(z)     # the same as lambertw(z,0)
lambertw_check_convergence(z, k=0) # The same as above but throw an error if the computation failed to converge
```

`z` may be Complex or Real. `k` must be an integer. For Real
`z`, `k` must be either `0` or `-1`.

Examples:

```julia
julia> lambertw(10)
1.7455280027406994

julia> lambertw(e)
1

julia> lambertw(1.0)
0.5671432904097838

julia> lambertw(-pi/2 + 0im)  / pi
4.6681174759251105e-18 + 0.5im
```

#### Note on `lambertw_check_convergence`

You can use this for extra safety. But I have been unable to find any input for which the root finding fails to
converge quickly.

### lambertwbp(x,k=0)

Return `1 + W(-1/e + z)`, for `z` satisfying `0 <= abs(z) < 1/e`,
on the branch of index `k`, where `k` must be either `0` or `-1`. This
function is designed to minimize loss of precision near the branch point `z=-1/e`.
`lambertwbp(z,k)` converges to `Float64` precision for `abs(z) < 0.32`.

If `k=-1` and `imag(z) < 0`, the value on the branch `k=1` is returned.

```julia
julia> lambertwbp(1e-3,-1)
-0.07560894118662498

julia> lambertwbp(0)
-0.0
```

`lambertwbp` uses a series expansion about the branch point `z=-1/e`.
The loss of precision in `lambertw` is analogous to the loss of precision
in computing the `sqrt(1-x)` for `x` close to `1`.

### LambertW.finv(lambertw)

The functional inverse of the Lambert W function.
```
julia> finv(lambertw)(lambertw(1.0))
1.0

julia> finv(lambertw)(lambertw(1+im/2,3))
1.0 + 0.49999999999999944im
```

### omega constant

The [omega constant](http://en.wikipedia.org/wiki/Omega_constant)

```julia
julia> ω
ω = 0.5671432904097...

julia> omega
ω = 0.5671432904097...

julia> ω * exp(ω)
1.0

julia> big(ω)
5.67143290409783872999968662210355549753815787186512508135131079223045793086683e-01 with 256 bits of precision

julia> lambertw(1) == float(ω)
true
```
<!-- ### Notes -->

<!-- Both `lambertw` and `lambertwbp` throw `DomainErrors` rather than return `NaN`s. -->
<!-- This behavior is reversed by setting `LAMBERTW_USE_NAN=true` at the top of -->
<!-- the source file `lambertw.jl`. -->
 
<!--  LocalWords:  lambertw jacobisymbol julia ulia im eval LambertW
 -->
<!--  LocalWords:  lambertwbp lambertwm NaN bitstype Combinatorics
 -->
<!--  LocalWords:  BigInt imag sqrt
 -->
