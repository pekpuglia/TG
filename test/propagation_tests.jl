@testset "Add orbital elements" begin
    rp = (6378+400)*1000.0
    ra = (6378+4000)*1000.0
    agiven = (rp + ra) / 2
    egiven = (ra - rp) / (ra + rp)
    orb0 = KeplerianElements(
        date_to_jd(2023, 1, 1, 0, 0, 0),
        agiven,
        egiven,
        30 |> deg2rad,
        15    |> deg2rad,
        60     |> deg2rad,
        180     |> deg2rad
    )
    r0, v0 = kepler_to_rv(orb0)

    model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
    "max_wall_time" => 30.0)
    )
    
    orbparams_i = add_orbital_elements!(model, false)
    r, v, a, e, i, Ω, ω, nu, M, E = getfield.(Ref(orbparams_i), fieldnames(FullOrbitalParameters))

    @constraint(model, a*(1-e) == rp)
    @constraint(model, a*(1+e) == ra)

    @constraint(model, i == orb0.i)

    @constraint(model, nu == orb0.f)

    @constraint(model, Ω == orb0.Ω)

    @constraint(model, ω == orb0.ω)

    optimize!(model)

    @test norm(value.(r) - r0) < 1
    @test norm(value.(v) - v0) < 1
    @test value(a) ≈ agiven
    @test value(e) ≈ egiven
    @test value(i) ≈ orb0.i atol = 1e-2
    @test mod(value(Ω), 2π) ≈ orb0.Ω atol = 1e-2
    @test mod(value(ω), 2π) ≈ orb0.ω atol = 1e-2

    model = Model(
        optimizer_with_attributes(Ipopt.Optimizer,
        "max_wall_time" => 30.0)
    )
    
    orbparams_i = add_orbital_elements!(model)
    r, v, a, e, i, Ω, ω, nu, M, E = getfield.(Ref(orbparams_i), fieldnames(FullOrbitalParameters))

    @constraint(model, r .== r0)
    @constraint(model, v .== v0)

    optimize!(model)

    @test norm(value.(r) - r0) < 1
    @test norm(value.(v) - v0) < 1
    @test value(a) ≈ agiven
    @test value(e) ≈ egiven
    @test value(i) ≈ orb0.i atol = 1e-2
    @test mod(value(Ω), 2π) ≈ orb0.Ω atol = 1e-2
    @test mod(value(ω), 2π) ≈ orb0.ω atol = 1e-2
end

@testset "Curtis Example 3.1" begin
    rp = 9600e3
    ra = 21000e3
    a = (rp + ra) / 2
    e = (ra - rp) / (ra + rp)
    orb0 = KeplerianElements(
        date_to_jd(2023, 1, 1, 0, 0, 0),
        a,
        e,
        30 |> deg2rad,
        0    |> deg2rad,
        0     |> deg2rad,
        0     |> deg2rad
    )
    r0, v0 = kepler_to_rv(orb0)
    orbf = @set orb0.f = 120 |> deg2rad
    r, v = kepler_to_rv(orbf)
    # T = orbital_period(orb0, GM_EARTH)

    model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
        "max_wall_time" => 30.0)
    )

    orbparams_i = add_orbital_elements!(model)
    ri, vi, ai, ei, ii, Ωi, ωi, nui, Mi, Ei = getfield.(Ref(orbparams_i), fieldnames(FullOrbitalParameters))
    orbparams_f = add_orbital_elements!(model)
    rf, vf, af, ef, i_f, Ωf, ωf, nuf, Mf, Ef = getfield.(Ref(orbparams_f), fieldnames(FullOrbitalParameters))

    @variable(model, Δt)

    add_coast_set_boundaries!(model, orbparams_i, orbparams_f, r0, v0, r, v, Δt)

    optimize!(model)

    @test is_solved_and_feasible(model)

    @test all(isapprox(value.(ri), r0, atol = 1e-3))
    @test all(isapprox(value.(vi), v0, atol = 1e-3))

    @test all(isapprox(value.(rf), r, atol = 1e-3))
    @test all(isapprox(value.(vf), v, atol = 1e-3))
    
    @test value(ai) ≈ a atol = 1e-3
    @test value(af) ≈ a atol = 1e-3

    @test value(ei) ≈ e atol = 1e-3
    @test value(ef) ≈ e atol = 1e-3

    @test value(Ei) ≈ 0.0 atol = 1e-6
    @test value(Ef) ≈ 1.7281 atol = 1e-3

    @test value(Mi) ≈ 0.0 atol = 1e-6
    @test value(Mf) ≈ 1.3601 atol = 1e-3

    @test value(model[:Δt]) ≈ 4077.0 atol = 1e-1
end

@testset "Curtis Example 3.2" begin
    rp = 9600e3
    ra = 21000e3
    a = (rp + ra) / 2
    e = (ra - rp) / (ra + rp)
    orbi = KeplerianElements(
        date_to_jd(2023, 1, 1, 0, 0, 0),
        a,
        e,
        30 |> deg2rad,
        0    |> deg2rad,
        0     |> deg2rad,
        0     |> deg2rad
    )
    ri, vi = kepler_to_rv(orbi)
    T = orbital_period(orbi, GM_EARTH)
    Δt = 10800.0

    model = Model(
        optimizer_with_attributes(Ipopt.Optimizer,
        "max_wall_time" => 30.0)
    )
    orbparams_i = add_orbital_elements!(model, true)
    r0, v0, a0, e0, i0, Ω0, ω0, nu0, M0, E0 = getfield.(Ref(orbparams_i), fieldnames(FullOrbitalParameters))

    orbparams_f = add_orbital_elements!(model, false)
    rf, vf, af, ef, i_f, Ωf, ωf, nuf, Mf, Ef = getfield.(Ref(orbparams_f), fieldnames(FullOrbitalParameters))

    @constraint(model, r0 .== ri)
    @constraint(model, v0 .== vi)

    @constraint(model, af == a0)
    @constraint(model, ef == e0)
    @constraint(model, i_f == i0)
    @constraint(model, Ωf == Ω0)
    @constraint(model, ωf == ω0)


    @constraint(model, Δt == (Mf - M0) / (2π) * T)

    optimize!(model)

    @test is_solved_and_feasible(model)

    @test all(isapprox(value.(r0), ri, atol = 1e-3))
    @test all(isapprox(value.(v0), vi, atol = 1e-3))
    
    @test value(a0) ≈ a atol = 1e-3
    @test value(af) ≈ a atol = 1e-3

    @test value(e0) ≈ e atol = 1e-3
    @test value(ef) ≈ e atol = 1e-3

    @test value(E0) ≈ 0.0 atol = 1e-6
    @test value(Ef) ≈ 3.4794 atol = 1e-3

    @test value(M0) ≈ 0.0 atol = 1e-6
    @test value(Mf) ≈ 3.6029 atol = 1e-3

    @test value(M0) ≈ 0.0 atol = 1e-6
    @test value(Mf) ≈ 3.6029 atol = 1e-3    

    @test value(nu0) ≈ 0.0 atol = 1e-6
    @test rad2deg(value(nuf)) ≈ 193.2 atol = 1e-1
end