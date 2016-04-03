# LambertW

This package implements the Lambert W function and associated omega constant

* lambertw
* lambertwbp
* omega constant

### lambertw

The [Lambert W function](http://en.wikipedia.org/wiki/Lambert_W_function)

```julia
lambertw(z,k)   # Lambert W function for argument z and branch index k
lambertw(z)     # the same as lambertw(z,0)
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

ulia> lambertw(-pi/2 + 0im)  / pi
4.6681174759251105e-18 + 0.5im
```

`lambertw` is vectorized, that is, it automatically maps over arrays.

### lambertwbp(x,k)

Returns `1 + W(-1/e + z)`, for `z` satisfying `0 <= abs(z) < 1/e`,
on the branch of index `k`, where `k` must be either `0` or `-1`. This
function is designed to minimize loss of precision near the branch point `z=-1/e`.
`lambertwbp(z,k)` converges to `Float64` precision for `abs(z) < 0.32`.

If `k=-1` and `imag(z) < 0`, the value on the branch `k=1` is returned.

`lambertwbp(z)` is equivalent to `lambertwbp(z,0)`.

```julia
julia> lambertwbp(1e-3,-1)
-0.07560894118662498

julia> lambertwbp(0)
-0.0
```

`lambertwbp` uses a series expansion about the branch point `z=-1/e`.
The loss of precision in `lambertw` is closely analogous to the loss of precision
in computing the floating point function `sqrt(1-x)` for `x` close to `1`.

`lambertwbp` is vectorized, that is automatically maps over arrays.

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
```

### Notes

Both `lambertw` and `lambertwbp` throw `DomainErrors` rather than return `NaN`s.
This behavior is reversed by setting `LAMBERTW_USE_NAN=true` at the top of
the source file `lambertw.jl`.
 
<!--  LocalWords:  lambertw jacobisymbol julia ulia im eval LambertW
 -->
<!--  LocalWords:  lambertwbp lambertwm NaN bitstype Combinatorics
 -->
<!--  LocalWords:  BigInt imag sqrt
 -->
