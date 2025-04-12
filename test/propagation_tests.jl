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
    T = orbital_period(orb0, GM_EARTH)

    model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
        "max_wall_time" => 30.0)
    )

    ri, vi, ai, ei, ii, Ωi, ωi, nui, Mi, Ei = add_orbital_elements!(model)
    rf, vf, af, ef, i_f, Ωf, ωf, nuf, Mf, Ef = add_orbital_elements!(model)

    @constraint(model, ri .== r0)
    @constraint(model, vi .== v0)
    @constraint(model, rf .== r)
    @constraint(model, vf .== v)

    @variable(model, Δt)

    @constraint(model, Δt == (Mf - Mi) / (2π) * T)

    optimize!(model)

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