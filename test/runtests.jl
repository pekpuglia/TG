using Test
include("../src/TG.jl")
using .TG
using JuMP
using Ipopt
using LinearAlgebra
using SatelliteToolboxBase
using Setfield

@testset "TG tests" begin
    
    ##########################################################################
    #
    #       Single Maneuver Tests
    #
    ##########################################################################
    include("single_maneuver_tests.jl")
    
    include("propagation_tests.jl")
end