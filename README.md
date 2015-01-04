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

Examples:

```julia
julia> lambertw(10)
1.7455280027406994

julia> lambertw(e)
1

julia> lambertw(1)
ω = 0.5671432904097...

julia> lambertw(1.0)
0.5671432904097838
```

### lambertwm1(x)

Lambert W function of "-1/e + x" on the branch of index -1.
For Float64 x < 1e-7 (approximately) this is more accurate than
`lambertw(-1/e+x,-1)`.

```julia
julia> lambertwm1(1e-18)
-1.000000002331644

julia> lambertwm1(0)
-1.0
```

`lambertwm1` uses a series expansion about the pole `z=-/e`.

### Symbolic things

It's fun to see exact expressions:

```julia
julia> lambertw(:(-1/e))
-1

julia> lambertw(:(-pi/2))   # but, we are confusing real and complex domains.
:(complex(0,pi / 2))

julia> eval(ans)
0.0 + 1.5707963267948966im

julia> lambertw( :( x * exp(x) ) )
:x

julia> lambertw( :( (a+b*z) * exp((a+b*z)) ) )
:(a + b * z)

julia> lambertw( :( (a+b*z) * log((a+b*z)) ) )
:(log(a + b * z))
```

But, since we have no general mechanism to reduce expressions to normal form,
this symbolic knowledge is of limited use. For instance, we don't have a rule for this:

```julia
julia> lambertw(:(pi/-2))
:(lambertw(pi / -2,0))
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
```

### jacobisymbol

`jacobisymbol(a,n)` returns the Jacobi symbol. This is limited to bitstype integers.
This is faster than Combinatorics.jacobisymbol for bitstype inputs, but slower for
BigInt inputs. Thus, these methods are complementary.
