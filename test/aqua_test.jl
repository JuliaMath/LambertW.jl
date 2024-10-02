using LambertW
using Aqua: Aqua

@testset "aqua deps compat" begin
    Aqua.test_deps_compat(LambertW)
end

# This often gives false positive
# @testset "aqua project toml formatting" begin
#     Aqua.test_project_toml_formatting(LambertW)
# end

@testset "aqua unbound_args" begin
    Aqua.test_unbound_args(LambertW)
end

@testset "aqua undefined exports" begin
    Aqua.test_undefined_exports(LambertW)
end

# Perhaps some of these should be fixed. Some are for combinations of types
# that make no sense.
@testset "aqua test ambiguities" begin
    Aqua.test_ambiguities([LambertW, Core, Base])
end

@testset "aqua piracies" begin
    Aqua.test_piracies(LambertW)
end

@testset "aqua project extras" begin
    Aqua.test_project_extras(LambertW)
end

@testset "aqua state deps" begin
    Aqua.test_stale_deps(LambertW)
end
