# SpecFun

### lambertw

The [Lambert W function](http://en.wikipedia.org/wiki/Lambert_W_function)

```julia
lambert(x,0)   # real domain and range,  branch index 0
lambert(x)     # the same as lambert(x,0)
lambert(x,-1)  # real domain and range,  branch index -1
lambert(z,k)   # complex domain and range,  branch index k
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
