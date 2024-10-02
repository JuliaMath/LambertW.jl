@testset "LambertW" begin

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
        @test isapprox(abs(lambertw(-1.0/MathConstants.e  + 0im, -1)), 1; atol=1e-15)
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

end # @testset "LambertW"

@testset "exceptions" begin
    ### domain errors

    @test_throws DomainError lambertw(-2.0, 0)
    @test_throws DomainError lambertw(-2.0, -1)
    @test_throws DomainError lambertw(-2.0, 1)
    @test isnan(lambertw(NaN))

    ## math constant e
    @test_throws DomainError lambertw(MathConstants.e, 1)
    @test_throws DomainError lambertw(MathConstants.e, -1)

    ###  expansion about branch point

    # not a domain error, but not implemented
    @test_throws ArgumentError lambertwbp(1, 1)
    @test_throws DomainError lambertw(.3, 2)

    @test_throws DomainError lambertwbp(1.1)
    @test_throws DomainError lambertwbp(complex(1.1))
end

@testset "branch point" begin
    # Expansions about branch point converges almost to machine precision
    # except near the radius of convergence.
    # Complex args are not tested here.

    if Int != Int32
        # Test double-precision expansion near branch point using BigFloats
        sp = precision(BigFloat)
        zinit = BigFloat(1)/10^12
        for z in (zinit, complex(zinit))
            setprecision(2048)
            for _ in 1:300 # We break from loop long before 300
                #            innerarg = z - 1 / big(MathConstants.e)
                innerarg = z + lambertwbranchpoint
                # branch k = 0
                wo = lambertwbp(Float64(z))
                xdiff = abs(-1 + wo - lambertw(innerarg))
                if xdiff > 5e-16
                    @warn(Float64(z), " ", Float64(xdiff))
                end
                @test xdiff < 5e-16

                #  branch k = -1
                wo = lambertwbp(Float64(z), -1)
                xdiff = abs(-1 + wo - lambertw(innerarg, -1))
                if xdiff > 5e-16
                    @warn(Float64(z), " ", Float64(xdiff))
                end
                @test xdiff < 5e-16
                z  *= 1.1
                if abs(z) > 0.23 break end
            end
            setprecision(sp)
        end
        # test the expansion about branch point for k=-1,
        # by comparing to exact BigFloat calculation.
        @test lambertwbp(1e-20, -1) - 1 - lambertw(-BigFloat(1)/big(MathConstants.e)+ BigFloat(1)/BigFloat(10)^BigFloat(20), -1) < 1e-16
        @test abs(lambertwbp(Complex(.01, .01), -1) - Complex(-0.2755038208041206, -0.1277888928494641)) < 1e-14

    else  # if Int != Int32
        @info "Skipping branch point tests"
    end
end

@testset "infinities" begin
    ### infinite args or return values
    @test lambertw(0, -1) == lambertw(0.0, -1) == -Inf
    @test lambertw(Inf, 0) == Inf
    @test lambertw(complex(Inf, 1), 0) == complex(Inf, 1)
    @test lambertw(complex(Inf, 0), 1) == complex(Inf, 2pi)
    @test lambertw(complex(-Inf, 0), 1) == complex(Inf, 3pi)
    @test lambertw(complex(0.0, 0.0), -1) == complex(-Inf, 0.0)
end

## value at branch point where real branches meet
@testset "at branch point" begin
    k0val = lambertw(lambertwbranchpoint, 0)
    km1val = lambertw(lambertwbranchpoint, -1)
    @test k0val == -1.0
    @test km1val == -1.0
    @test k0val isa Float64
    @test km1val isa Float64
end

@testset "complex code paths" begin
    z = complex(1/MathConstants.e - .01)
    @test isapprox(lambertw(z), 0.27251232622985155)
    @test isapprox(lambertw(z, -1), -2.6181060466381134 - 4.1495292474932475im)
end

@testset "finv" begin
    lambertw_inv = LambertW.finv(lambertw)
    z = 42.0
    w = lambertw(z)
    @test isapprox(z, lambertw_inv(w))
end

@testset "show" begin
    @test string(LambertW.Omega()) == "ω"
end

@testset "lambertw info" begin
    result = lambertw(1.0; info=true)
    @test result[1] == lambertw(1.0)
    @test result[2]
    @test result[3] > 1 && result[3] < 10

    for z in (10., complex(10), lambertwbranchpoint)
        res = lambertw(1.0; info=true)
        @test res[2]
        @test length(res) == 3
    end
end
