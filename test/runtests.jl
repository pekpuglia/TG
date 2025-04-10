using Test
include("../src/TG.jl")
using .TG
using JuMP
using LinearAlgebra
using SatelliteToolboxBase

##########################################################################
#
#       Single Maneuver Tests
#
##########################################################################
@testset "Single Maneuver Models" begin
    LEOorb = KeplerianElements(
                  date_to_jd(2023, 1, 1, 0, 0, 0),
                  8000e3,
                  0.101111,
                  20 |> deg2rad,
                  50    |> deg2rad,
                  30     |> deg2rad,
                  70     |> deg2rad
    )

    LEOdeltaV = [1000, 0, 0]

    LEO_r_final, LEO_total_time, _, _ = final_position(LEOorb, LEOdeltaV, 0.25, 0.3)

    model, r_maneuver, v_post_maneuver, rf, vf = single_maneuver_model(LEOorb, LEO_r_final, LEO_total_time);

    optimize!(model)

    @test is_solved_and_feasible(model)

    @test norm(value.(model[:ΔV])) <= norm(LEOdeltaV)
end