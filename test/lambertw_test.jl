macro test_baddomain(expr)
    if LambertW.LAMBERTW_USE_NAN
        :(@test $(esc(expr)) === NaN)
    else
        :(@test_throws DomainError $(esc(expr)))
    end
end

### domain errors

@test_baddomain lambertw(-2.0,0)
@test_baddomain lambertw(-2.0,-1)
@test_baddomain lambertw(-2.0,1)
@test_baddomain lambertw(NaN)

## math constant e
@test_baddomain lambertw(e,1)
@test_baddomain lambertw(e,-1)

## integer arguments return floating point types
@test typeof(lambertw(0)) <: AbstractFloat
@test lambertw(0) == 0

### math constant, euler e

# could return math const e, but this would break type stability
@test typeof(lambertw(1)) <: AbstractFloat
@test lambertw(e,0) == 1

## value at branch point where real branches meet
@test lambertw(-1/e,0) == lambertw(-1/e,-1) == -1
@test typeof(lambertw(-1/e,0)) == typeof(lambertw(-1/e,-1)) <: AbstractFloat

### infinite args or return values

@test lambertw(0,-1) == lambertw(0.0,-1) == -Inf
@test lambertw(Inf,0) == Inf
@test lambertw(complex(Inf,1),0) == complex(Inf,1)
@test lambertw(complex(Inf,0),1) == complex(Inf,2pi)
@test lambertw(complex(-Inf,0),1) == complex(Inf,3pi)
@test lambertw(complex(0.0,0.0),-1) == complex(-Inf,0.0)

## default branch is  k = 0
@test lambertw(1.0) == lambertw(1.0,0)

## BigInt args return BigFloats
@test typeof(lambertw(BigInt(0))) == BigFloat
@test typeof(lambertw(BigInt(3))) == BigFloat

## Any Integer type allowed for second argument
@test lambertw(-0.2,-1) == lambertw(-0.2,BigInt(-1))

## BigInt for second arg does not promote the type
@test typeof(lambertw(-0.2,-1)) == typeof(lambertw(-0.2,BigInt(-1)))

for (z,k,res) in [ (0,0 ,0), (complex(0,0),0 ,0),
              (complex(0.0,0),0 ,0), (complex(1.0,0),0, 0.567143290409783873) ]
    @test_approx_eq  lambertw(z,k) res
    @test_approx_eq  lambertw(z) res
end

for (z,k) in ((complex(1,1),2), (complex(1,1),0),(complex(.6,.6),0),
              (complex(.6,-.6),0))
    let w
        @test (w = lambertw(z,k) ; true)
        @test abs(w*exp(w) - z) < 1e-15
    end
end

# bug fix
# The routine will start at -1/e + eps * im, rather than -1/e + 0im,
# otherwise root finding will fail
@test abs(lambertw(-1.0/e  + 0im,-1)) == 1

# lambertw for BigFloat is more precise than Float64. Note
# that 70 digits in test is about 35 digits in W
let W
    for z in [ BigFloat(1),  BigFloat(2), complex(BigFloat(1), BigFloat(1))]
        @test (W = lambertw(z); true)
        @test abs(z - W * exp(W)) < BigFloat(1)^(-70)
    end
end

###  ω constant

## get ω from recursion and compare to value from lambertw
let sp = precision(BigFloat)
    setprecision(512)
    @test lambertw(big(1)) == big(ω)
    setprecision(sp)
end

#@test lambertw(1) == ω
@test lambertw(1) == float(ω)
@test Float16(ω) == Float16(0.5674)
@test Float32(ω) == 0.56714326f0
@test lambertw(BigInt(1)) == big(ω)

###  expansion about branch point

# not a domain error, but not implemented
@test_throws ErrorException lambertwbp(1,1)

# Expansions about branch point converges almost to machine precision
# except near the radius of convergence.
# Complex args are not tested here.
let sp = precision(BigFloat), z = BigFloat(1)/10^12, wo, diff
    setprecision(2048)
    for i in 1:300
        # k = 0
        @test (wo = lambertwbp(Float64(z)); diff = abs(-1 + wo - lambertw(z-1/big(e))); true)
        if diff > 5e-16
            println(Float64(z), " ", Float64(diff))
        end
        @test diff < 5e-16
        #  k = -1
        @test (wo = lambertwbp(Float64(z),-1); diff = abs(-1 + wo - lambertw(z-1/big(e),-1)); true)
        if diff > 5e-16
            println(Float64(z), " ", Float64(diff))
        end
        @test diff < 5e-16
        z  *= 1.1
        if z > 0.23 break end
    end
    setprecision(sp)
end

# test the expansion about branch point for k=-1,
# by comparing to exact BigFloat calculation.
@test lambertwbp(1e-20,-1) - 1 - lambertw(-BigFloat(1)/big(e)+ BigFloat(1)/BigFloat(10)^BigFloat(20),-1) < 1e-16

## vectorization

@test lambertw([0.1,0.2]) == [lambertw(0.1),lambertw(0.2)]
@test lambertw([0.1+im ,0.2-im]) == [lambertw(0.1+im),lambertw(0.2-im)]
@test lambertw([0.1,-0.2],[0,-1]) == [lambertw(0.1,0),lambertw(-0.2,-1)]
