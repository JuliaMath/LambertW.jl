@test_throws DomainError lambertw(-2.0,0)
@test_throws DomainError lambertw(-2.0,-1)
@test_throws DomainError lambertw(-2.0,1)
@test lambertw(0,-1) == lambertw(0.0,-1) == -Inf
@test lambertw(-1/e,0) == lambertw(-1/e,-1) == -1
@test_throws DomainError lambertw(NaN) 

@test lambertw(complex(Inf,1),0) == complex(Inf,1)
@test lambertw(complex(Inf,0),1) == complex(Inf,2pi)
@test lambertw(complex(-Inf,0),1) == complex(Inf,3pi)
@test lambertw(1.0) == lambertw(1.0,0)
@test lambertw(complex(0.0,0.0),-1) == complex(-Inf,0.0)

@test typeof(lambertw(0)) <: FloatingPoint
@test lambertw(0) == 0

@test typeof(lambertw(1)) <: FloatingPoint
@test lambertw(1) == float(ω)
@test lambertw(BigInt(1)) == big(ω)
@test typeof(lambertw(BigInt(0))) == BigFloat
@test typeof(lambertw(BigInt(3))) == BigFloat

for tvals in [ (0,0,0), (complex(0,0),0,0),
              (0.0 + 0 * im ,0,0), (1.0 + 0 * im,0,0.567143290409783873) ]
    (z,k,res) = tvals
    @test_approx_eq  lambertw(z,k) res
end

for tvals in [ (0,0), (complex(0,0),0), (0.0,0), (complex(0.0,0),0) ]
    (z,res) = tvals
    @test_approx_eq  lambertw(z) res
end

for (z,k) in ( (complex(1,1),2), (complex(1,1),0),(complex(.6,.6),0),
     (complex(.6,-.6),0))
    let w
        @test (w = lambertw(z,k) ; true)  # so code-coverage see this
        @test abs(w*exp(w) - z) < 1e-15
    end
end    

@test lambertw(e,0) == 1
@test_throws DomainError lambertw(e,1)
@test_throws DomainError lambertw(e,-1)

let sp = get_bigfloat_precision()
    set_bigfloat_precision(512)
    @test lambertw(big(1)) == big(ω)
    set_bigfloat_precision(sp)
end

@test lambertw(1) == ω
@test_throws ErrorException lambertwbp(1,1)

let z = 1e-3, wo, diff
    @test (wo = lambertwbp(z); diff = abs(-1 + wo - lambertw(z-1/e)); true)
    @test 1e-9 < diff < 1e-6
end
    
let w
    for z in [ BigFloat(1),  BigFloat(2), complex(BigFloat(1), BigFloat(1))]
        @test (w = lambertw(z); true)
        @test abs(z - w * exp(w)) < BigFloat(1)^(-70)
    end
end

# test the expansion about branch point for k=-1,
# by comparing to exact BigFloat calculation.
@test lambertwbp(1e-20,-1) - 1 - lambertw(-BigFloat(1)/big(e)+ BigFloat(1)/BigFloat(10)^BigFloat(20),-1) < 1e-16

# Fails unless we offset the starting point slightly before root finding.
@test abs(lambertw(-1.0/e  + 0im,-1)) == 1
