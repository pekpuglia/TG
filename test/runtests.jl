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
@testset "Lambert test Curtis 5.2" begin
    r1 = [
        5000
        10000
        2100
    ]*1e3

    r2 = [
        -14600
        2500
        7000
    ]*1e3

    delta_t = 3600

    v1, v2 = lambert(r1, r2, delta_t)

    @test all(isapprox.(v1, [-5.9925, 1.9254, 3.2456]*1e3, rtol=1e-2))
    @test all(isapprox.(v2, [-3.3125, -4.1966, -0.38529]*1e3, rtol=1e-2))
end

@testset "Glandorf tests" begin
    r1 = [
        5000
        10000
        2100
    ]*1e3

    v1 = [
        -5992.495019799777
        1925.3667119503089
        3245.638049455475
    ]

    r2 = [
        -14600
        2500
        7000
    ]*1e3

    v2 = [
        -3312.458504221757
        -4196.619006657186
        -385.2890588564801
    ]

    delta_t = 3600

    P = P_glandorf(r1, v1, 0)
    Pinv = Pinv_glandorf(r1, v1, 0)
    @test all(isapprox.(eigen(P*Pinv).values, 1.0, atol=1e-4))

    Phi = P_glandorf(r2, v2, delta_t)*Pinv_glandorf(r1, v1, 0)
    Phi_other = P_glandorf(r2, v2, delta_t+10)*Pinv_glandorf(r1, v1, 10)

    @test all(isapprox.(Phi, Phi_other, rtol=1e-3))
end
# @testset "Single Maneuver Models" begin
#     @testset "LEO case" begin
#         LEOorb = KeplerianElements(
#                       date_to_jd(2023, 1, 1, 0, 0, 0),
#                       8000e3,
#                       0.101111,
#                       20 |> deg2rad,
#                       50    |> deg2rad,
#                       30     |> deg2rad,
#                       70     |> deg2rad
#         )
    
#         LEOdeltaV = [1000, 0, 0]
    
#         LEO_r_final, LEO_total_time, _, _, dir_v_final = final_position(LEOorb, LEOdeltaV, 0.25, 0.3)
    
#         model, r_maneuver, v_post_maneuver, rf, vf = single_maneuver_model(LEOorb, LEO_r_final, dir_v_final, LEO_total_time);
    
#         optimize!(model)
    
#         @test is_solved_and_feasible(model)
    
#         @test norm(value.(model[:ΔVmag])) <= norm(LEOdeltaV)
#     end

#     @testset "LEO high delta V" begin
#         LEOorb = KeplerianElements(
#             date_to_jd(2023, 1, 1, 0, 0, 0),
#             8000e3,
#             0.01,
#             1.5 |> deg2rad,
#             1.5    |> deg2rad,
#             1.5     |> deg2rad,
#             1.5     |> deg2rad
#         )
    
#         deltaV = [0, -2000, 0]
    
#         LEO_r_final, LEO_total_time, _, _, dir_v_final = final_position(LEOorb, LEOdeltaV, dir, 0.5, 0.7)
    
#         model, r_maneuver, v_post_maneuver, rf, vf = single_maneuver_model(LEOorb, LEO_r_final, dir_v_final, LEO_total_time);
    
#         optimize!(model)
    
#         @test is_solved_and_feasible(model)
    
#         @test norm(value.(model[:ΔVmag])) <= norm(LEOdeltaV)
#     end

#     @testset "GEO case" begin
#         orb = KeplerianElements(
#                       date_to_jd(2023, 1, 1, 0, 0, 0),
#                       30000e3,
#                       0.101111,
#                       20 |> deg2rad,
#                       50    |> deg2rad,
#                       30     |> deg2rad,
#                       70     |> deg2rad
#         )
    
#         deltaV = [-1000, 0, 0]
    
#         r_final, total_time, _, _, dir_v_final = final_position(
#             orb, 
#             deltaV,
#             0.25,
#             0.3)
    
#         model, r_maneuver, v_post_maneuver, rf, vf = single_maneuver_model(
#             orb, 
#             r_final, 
#             dir_v_final,
#             total_time
#         );
    
#         optimize!(model)
    
#         @test is_solved_and_feasible(model)
    
#         @test norm(value.(model[:ΔVmag])) <= norm(deltaV)
#     end
    
#     #remove? out of scope
#     @testset "High excentricity high deltaV" begin
#         orb0 = KeplerianElements(
#             date_to_jd(2023, 1, 1, 0, 0, 0),
#             30000e3,
#             0.6,
#             1.5 |> deg2rad,
#             1.5    |> deg2rad,
#             1.5     |> deg2rad,
#             1.5     |> deg2rad
#         )
#         deltaV = [0, -1500, 0]

#         r_final, total_time, orb_postman, orb_final, dir_v_final = final_position(
#             orb0, 
#             deltaV, 
#             0.5, 
#             0.7
#         )
            
#         ## solving maneuver
#         model, r_maneuver, v_post_maneuver, rf, vf = single_maneuver_model(
#             orb0, 
#             r_final, 
#             dir_v_final,
#             total_time);
#         ##
#         optimize!(model)

#         @test is_solved_and_feasible(model)
    
#         @test norm(value.(model[:ΔVmag])) <= norm(deltaV)
#     end

# end