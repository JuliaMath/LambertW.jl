# SpecFun

These are a couple of mathematical functions that could be eventually
included in some other package.

* lambertw
* omega constant
* jacobisymbol

### lambertw

The [Lambert W function](http://en.wikipedia.org/wiki/Lambert_W_function)

```julia
lambertw(z,k)   # Lambert W function for argument z and branch index k
lambertw(z)     # the same as lambertw(z,0)
```

z may be Complex or Real.

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

### jacobisymbol

`jacobisymbol(a,n)` returns the Jacobi symbol. This is limited to bitstype integers.
This is faster than Combinatorics.jacobisymbol for bitstype inputs, but slower for
BigInt inputs. Thus, these methods are complementary.
