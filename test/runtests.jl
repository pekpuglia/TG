using Test
include("../src/TG.jl")
using .TG
using JuMP
using Ipopt
using LinearAlgebra
using SatelliteToolboxBase
using Setfield

@testset "TG tests" begin
    include("propagation_tests.jl")
end