using LambertW

@testset "does not error for small negative real complex" begin
  try lambertw(-0.4087512737796557 + 0.038455913083347504im, 0)
    @test true
  catch err
    @test false
  end
end
