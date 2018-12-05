using LambertW
using LambertW: lambertwbranchpoint
using Test

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
end

@testset "LambertW.jl" begin
  include("lambertw_test.jl")
end

@testset "branch point" begin
    # Expansions about branch point converges almost to machine precision
    # except near the radius of convergence.
    # Complex args are not tested here.

    if Int != Int32
        # Test double-precision expansion near branch point using BigFloats
        sp = precision(BigFloat)
        z = BigFloat(1)/10^12
        setprecision(2048)
        for i in 1:300
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
            if z > 0.23 break end
        end
        setprecision(sp)

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
