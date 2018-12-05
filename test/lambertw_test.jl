
## integer arguments return floating point types
@test lambertw(0) isa AbstractFloat
@test lambertw(0) == 0

### math constant, MathConstants.e e

# could return math const e, but this would break type stability
@test lambertw(1) isa AbstractFloat
@test lambertw(MathConstants.e, 0) == 1


## convert irrationals to float

@test isapprox(lambertw(pi), 1.0736581947961492)
@test isapprox(lambertw(pi, 0), 1.0736581947961492)

## default branch is  k = 0
@test lambertw(1.0) == lambertw(1.0, 0)

## BigInt args return BigFloats
@test typeof(lambertw(BigInt(0))) == BigFloat
@test typeof(lambertw(BigInt(3))) == BigFloat

## Any Integer type allowed for second argument
@test lambertw(-0.2, -1) == lambertw(-0.2, BigInt(-1))

## BigInt for second arg does not promote the type
@test typeof(lambertw(-0.2, -1)) == typeof(lambertw(-0.2, BigInt(-1)))

for (z, k, res) in [ (0, 0 , 0), (complex(0, 0), 0 , 0),
                   (complex(0.0, 0), 0 , 0), (complex(1.0, 0), 0, 0.567143290409783873) ]
    if Int != Int32
        @test isapprox(lambertw(z, k), res)
        @test isapprox(lambertw(z), res)
    else
        @test isapprox(lambertw(z, k), res; rtol = 1e-14)
        @test isapprox(lambertw(z), res; rtol = 1e-14)
    end
end

for (z, k) in ((complex(1, 1), 2), (complex(1, 1), 0), (complex(.6, .6), 0),
              (complex(.6, -.6), 0))
    let w
        @test (w = lambertw(z, k) ; true)
        @test abs(w*exp(w) - z) < 1e-15
    end
end

@test abs(lambertw(complex(-3.0, -4.0), 0) - Complex(1.075073066569255, -1.3251023817343588)) < 1e-14
@test abs(lambertw(complex(-3.0, -4.0), 1) - Complex(0.5887666813694675, 2.7118802109452247)) < 1e-14
@test (lambertw(complex(.3, .3), 0); true)

# bug fix
# The routine will start at -1/e + eps * im, rather than -1/e + 0im,
# otherwise root finding will fail
if Int != Int32
    @test abs(lambertw(-1.0/MathConstants.e  + 0im, -1)) == 1
else
    @test abs(lambertw(-1.0/MathConstants.e  + 0im, -1) + 1) < 1e-7
end
# lambertw for BigFloat is more precise than Float64. Note
# that 70 digits in test is about 35 digits in W
let W
    for z in [ BigFloat(1), BigFloat(2), complex(BigFloat(1), BigFloat(1))]
        @test (W = lambertw(z); true)
        @test abs(z - W * exp(W)) < BigFloat(1)^(-70)
    end
end

###  ω constant

## get ω from recursion and compare to value from lambertw
let sp = precision(BigFloat)
    setprecision(512)
    @test lambertw(big(1)) == big(LambertW.omega)
    setprecision(sp)
end

@test lambertw(1) == float(LambertW.omega)
@test convert(Float16, LambertW.omega) == convert(Float16, 0.5674)
@test convert(Float32, LambertW.omega) == 0.56714326f0
@test lambertw(BigInt(1)) == big(LambertW.omega)
