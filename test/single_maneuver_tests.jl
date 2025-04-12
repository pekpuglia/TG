@testset "Single Maneuver Models" begin
    @testset "LEO case" begin
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
    
        LEO_r_final, LEO_total_time, _, _, dir_v_final = final_position(LEOorb, LEOdeltaV, 0.25, 0.3)
    
        model, r_maneuver, v_post_maneuver, rf, vf = single_maneuver_model(LEOorb, LEO_r_final, dir_v_final, LEO_total_time);
    
        optimize!(model)
    
        @test is_solved_and_feasible(model)
    
        @test norm(value.(model[:ΔVmag])) <= norm(LEOdeltaV)
    end

    @testset "LEO high delta V" begin
        LEOorb = KeplerianElements(
            date_to_jd(2023, 1, 1, 0, 0, 0),
            8000e3,
            0.01,
            1.5 |> deg2rad,
            1.5    |> deg2rad,
            1.5     |> deg2rad,
            1.5     |> deg2rad
        )
    
        deltaV = [0, -2000, 0]
    
        LEO_r_final, LEO_total_time, _, _, dir_v_final = final_position(LEOorb, LEOdeltaV, dir, 0.5, 0.7)
    
        model, r_maneuver, v_post_maneuver, rf, vf = single_maneuver_model(LEOorb, LEO_r_final, dir_v_final, LEO_total_time);
    
        optimize!(model)
    
        @test is_solved_and_feasible(model)
    
        @test norm(value.(model[:ΔVmag])) <= norm(LEOdeltaV)
    end

    @testset "GEO case" begin
        orb = KeplerianElements(
                      date_to_jd(2023, 1, 1, 0, 0, 0),
                      30000e3,
                      0.101111,
                      20 |> deg2rad,
                      50    |> deg2rad,
                      30     |> deg2rad,
                      70     |> deg2rad
        )
    
        deltaV = [-1000, 0, 0]
    
        r_final, total_time, _, _, dir_v_final = final_position(
            orb, 
            deltaV,
            0.25,
            0.3)
    
        model, r_maneuver, v_post_maneuver, rf, vf = single_maneuver_model(
            orb, 
            r_final, 
            dir_v_final,
            total_time
        );
    
        optimize!(model)
    
        @test is_solved_and_feasible(model)
    
        @test norm(value.(model[:ΔVmag])) <= norm(deltaV)
    end
    
    #remove? out of scope
    @testset "High excentricity high deltaV" begin
        orb0 = KeplerianElements(
            date_to_jd(2023, 1, 1, 0, 0, 0),
            30000e3,
            0.6,
            1.5 |> deg2rad,
            1.5    |> deg2rad,
            1.5     |> deg2rad,
            1.5     |> deg2rad
        )
        deltaV = [0, -1500, 0]

        r_final, total_time, orb_postman, orb_final, dir_v_final = final_position(
            orb0, 
            deltaV, 
            0.5, 
            0.7
        )
            
        ## solving maneuver
        model, r_maneuver, v_post_maneuver, rf, vf = single_maneuver_model(
            orb0, 
            r_final, 
            dir_v_final,
            total_time);
        ##
        optimize!(model)

        @test is_solved_and_feasible(model)
    
        @test norm(value.(model[:ΔVmag])) <= norm(deltaV)
    end

end