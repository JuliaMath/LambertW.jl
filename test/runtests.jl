using LambertW
using LambertW: lambertwbranchpoint
using Test

include("lambertw_test.jl")
include("aqua_test.jl")

if VERSION >= v"1.7"
    include("jet_test.jl")
end
