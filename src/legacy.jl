
include("./lambert.jl")
include("./coast_transfer.jl")

noNaNs(x::Real) = true
noNaNs(x::ForwardDiff.Dual) = !any(isnan, ForwardDiff.partials(x))

export propagate_coast
function propagate_coast(xi, yi, zi, vxi, vyi, vzi, ti, deltat)
    r = [xi; yi; zi]
    v = [vxi; vyi; vzi]
    
    #always shows 0 at the beginning of solve - why?
    (norm(r) < EARTH_EQUATORIAL_RADIUS) ? (@info "Passing through Earth r = $r") : nothing
    Vesc = √(2tbc_m0 / EARTH_EQUATORIAL_RADIUS)
    (norm(v) <= 100.0) ? (@info "Velocity very small v = $v") : nothing

    # @assert norm(v) > 1e3 "Velocity too small"
    # @assert norm(cross(r, v)) > 1e9 "Angular momentum too small"
    energy = -tbc_m0/norm(r) + (v'*v)/2
    @assert energy < 0 "Energy non negative: $energy, t = $ti, r = $r, v = $v"


    orbi = rv_to_kepler([xi, yi, zi], [vxi, vyi, vzi], ti)

    #low excentricity + f ≈ 0 causes NaN derivatives
    @assert noNaNs(orbi.a) "NaN found"
    @assert noNaNs(orbi.e) "NaN found"
    @assert noNaNs(orbi.f) "NaN found"
    @assert noNaNs(orbi.i) "NaN found"
    @assert noNaNs(orbi.t) "NaN found"
    @assert noNaNs(orbi.Ω) "NaN found"
    @assert noNaNs(orbi.ω) "NaN found"

    propagator = Propagators.init(Val(:TwoBody), orbi, propagation_type=Real)
    r, v = Propagators.propagate!(propagator, deltat)

    perigee = orbi.a * (1 - orbi.e)

    [r; v; propagator.tbd.orbk.t; perigee]
end

export memoize
# https://jump.dev/JuMP.jl/stable/tutorials/nonlinear/tips_and_tricks/#User-defined-operators-with-vector-outputs
"""
    memoize(foo::Function, n_outputs::Int)

Take a function `foo` and return a vector of length `n_outputs`, where element
`i` is a function that returns the equivalent of `foo(x...)[i]`.

To avoid duplication of work, cache the most-recent evaluations of `foo`.
Because `foo_i` is auto-differentiated with ForwardDiff, our cache needs to
work when `x` is a `Float64` and a `ForwardDiff.Dual`.
"""
function memoize(foo::Function, n_outputs::Int)
    last_x, last_f = nothing, nothing
    last_dx, last_dfdx = nothing, nothing
    function foo_i(i, x::T...) where {T<:Real}
        if T == Float64
            if x !== last_x
                last_x, last_f = x, foo(x...)
            end
            return last_f[i]::T
        else
            if x !== last_dx
                last_dx, last_dfdx = x, foo(x...)
            end
            return last_dfdx[i]::T
        end
    end
    return [(x...) -> foo_i(i, x...) for i in 1:n_outputs]
end

export add_coast_operators!
function add_coast_operators!(model)
    memoized_propagate_coast = memoize(propagate_coast, 8)
    
    #for some reason this is required to avoid error no method matching getindex(::Nothing, ::Int64)
    #no f clue why
    #maybe something to do with the memoization initialization?
    orb0 = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    8000e3,
    0.001111,
    2 |> deg2rad,
    50    |> deg2rad,
    30     |> deg2rad,
    1     |> deg2rad
)
    r0, v0 = kepler_to_rv(orb0)
    ForwardDiff.gradient(x -> memoized_propagate_coast[7](x...), [r0..., v0..., 0.0, 10.0])

    coast_position_x = @operator(model, coast_position_x, 8, memoized_propagate_coast[1])
    coast_position_y = @operator(model, coast_position_y, 8, memoized_propagate_coast[2])
    coast_position_z = @operator(model, coast_position_z, 8, memoized_propagate_coast[3])
    coast_velocity_x = @operator(model, coast_velocity_x, 8, memoized_propagate_coast[4])
    coast_velocity_y = @operator(model, coast_velocity_y, 8, memoized_propagate_coast[5])
    coast_velocity_z = @operator(model, coast_velocity_z, 8, memoized_propagate_coast[6])
    coast_time       = @operator(model, coast_time      , 8, memoized_propagate_coast[7])
    coast_perigee    = @operator(model, coast_perigee   , 8, memoized_propagate_coast[8])

    return (
        (r, v, t, dt) -> [
            coast_position_x(r..., v..., t, dt)
            coast_position_y(r..., v..., t, dt)
            coast_position_z(r..., v..., t, dt)
        ],
        (r, v, t, dt) -> [
            coast_velocity_x(r..., v..., t, dt)
            coast_velocity_y(r..., v..., t, dt)
            coast_velocity_z(r..., v..., t, dt)
        ],
        (r, v, t, dt) -> coast_time(r..., v..., t, dt),
        (r, v, t, dt) -> coast_perigee(r..., v..., t, dt)
    )
end

export single_maneuver_model
function single_maneuver_model(orb0, r_final, dir_v_final, total_time)
    ## auxiliary parameters
    Vesc = √(2tbc_m0 / EARTH_EQUATORIAL_RADIUS)
    Vmin = 100.0

    moon_distance = 384400.e2 

    model = Model(
        optimizer_with_attributes(Ipopt.Optimizer,
            "max_iter" => 10000,
            "max_wall_time" => 30.0)
    )

    #control variables
    Δt_maneuver = @variable(model, Δt_maneuver, start=0.5*total_time)
    @constraint(model, 0 <= Δt_maneuver <= total_time)  #set start value of constraints?

    #colocar no referencial local da v_pre_maneuver
    ΔVmag = @variable(model, 0 <= ΔVmag <= Vesc, start=1.0)
    
    ΔVdir = @variable(model, -1 <= ΔVdir[1:3] <= 1, start = 0.0)
    set_start_value(ΔVdir[1], 1.0)

    @constraint(model, ΔVdir' * ΔVdir == 1)

    # ΔV = @variable(model, -Vesc <= ΔV[1:3] <= Vesc, start = 1.0)

    coast_r, coast_v, coast_t, coast_perigee = add_coast_operators!(model)
    
    r0, v0 = kepler_to_rv(orb0)

    #first coast
    r_maneuver = coast_r(r0, v0, orb0.t, Δt_maneuver)

    @constraint(model, EARTH_EQUATORIAL_RADIUS^2 <= r_maneuver' * r_maneuver <= moon_distance^2)

    v_pre_maneuver = coast_v(r0, v0, orb0.t, Δt_maneuver)

    @constraints(model, begin
        Vesc >= v_pre_maneuver[1] >= -Vesc
        Vesc >= v_pre_maneuver[2] >= -Vesc 
        Vesc >= v_pre_maneuver[3] >= -Vesc 
    end)
    @constraint(model, Vesc^2 >= v_pre_maneuver' * v_pre_maneuver >= Vmin^2)
    @constraint(model, (v_pre_maneuver' * v_pre_maneuver) * √(r_maneuver' * r_maneuver) / (2tbc_m0) <= 1-1e-6)

    time_maneuver = coast_t(r0, v0, orb0.t, Δt_maneuver)
    
    #delta V is in tangential reference frame
    #x ∥ v
    #z ∥ h = r × v
    #y = z × x

    x_tang = v_pre_maneuver ./ √(v_pre_maneuver' * v_pre_maneuver)

    h = cross(r_maneuver, v_pre_maneuver)

    z_tang = h ./ √(h' * h)

    y_tang = cross(z_tang, x_tang)

    ΔV = ΔVmag * ΔVdir

    v_post_maneuver = v_pre_maneuver + [x_tang y_tang z_tang] * ΔV

    @constraints(model, begin
        Vesc >= v_post_maneuver[1] >= -Vesc
        Vesc >= v_post_maneuver[2] >= -Vesc 
        Vesc >= v_post_maneuver[3] >= -Vesc 
    end)
    @constraint(model, Vesc^2 >= v_post_maneuver' * v_post_maneuver >= Vmin^2)
    @constraint(model, (v_post_maneuver' * v_post_maneuver) * √(r_maneuver' * r_maneuver) / (2tbc_m0) <= 1-1e-6)
    
    #perigee outside of earth
    @constraint(model, coast_perigee(r_maneuver, v_post_maneuver, time_maneuver, total_time - Δt_maneuver) >= EARTH_EQUATORIAL_RADIUS)

    #second coast
    rf = coast_r(r_maneuver, v_post_maneuver, time_maneuver, total_time - Δt_maneuver)

    @constraint(model, EARTH_EQUATORIAL_RADIUS^2 <= rf' * rf <= moon_distance^2)

    vf = coast_v(r_maneuver, v_post_maneuver, time_maneuver, total_time - Δt_maneuver)
    
    @constraints(model, begin
        Vesc >= vf[1] >= -Vesc
        Vesc >= vf[2] >= -Vesc 
        Vesc >= vf[3] >= -Vesc 
    end)
    @constraint(model, Vesc^2 >= vf' * vf >= Vmin^2)
    @constraint(model, (vf' * vf) * √(rf' * rf) / (2tbc_m0) <= 1-1e-6)

    dir_vf = vf ./ √(vf' * vf)
    
    @constraints(model, begin
        -100.0 <= rf[1] - r_final[1] <= 100.0
        -100.0 <= rf[2] - r_final[2] <= 100.0
        -100.0 <= rf[3] - r_final[3] <= 100.0
        dot(dir_vf, dir_v_final) >= 0.97
    end)
    
    # @objective(model, MIN_SENSE, √(ΔV' * ΔV))
    @objective(model, MIN_SENSE, ΔVmag)

    model, r_maneuver, v_post_maneuver, rf, vf
end

#testing purposes only
export final_position
function final_position(orb0, deltaV, period_fraction1, period_fraction2)
    T0 = orbital_period(orb0, tbc_m0)

    r_preman, v_preman, orbp_preman = Propagators.propagate(Val(:TwoBody), period_fraction1*T0, orb0)
    orb_preman = orbp_preman.tbd.orbk

    v_postman = v_preman + deltaV
    orb_postman = rv_to_kepler(r_preman, v_postman, orb_preman.t)
    T_postman = orbital_period(orb_postman, tbc_m0)

    r_final, v_final, orbp_postman = Propagators.propagate(Val(:TwoBody), period_fraction2*T_postman, orb_postman)
    total_time = orbp_preman.tbd.Δt + orbp_postman.tbd.Δt
    r_final, total_time, orb_postman, orbp_postman.tbd.orbk, v_final
end

export FullOrbitalParameters
struct FullOrbitalParameters
    r
    v
    a
    e
    i
    Ω
    ω
    nu
    M
    E # excluir?
end

export add_orbital_elements!
function add_orbital_elements!(model)
    Vorb_sup = √(GM_EARTH/EARTH_EQUATORIAL_RADIUS)
    rscaled = @variable(model, [1:3])
    @constraint(model, rscaled' * rscaled >= 1)
    r = EARTH_EQUATORIAL_RADIUS * rscaled
    vscaled = @variable(model, [1:3])
    v = Vorb_sup*vscaled

    ascaled = @variable(model, lower_bound = 1.0)
    a = EARTH_EQUATORIAL_RADIUS * ascaled
    e = @variable(model, lower_bound = 0, upper_bound = 1) 
    i = @variable(model, lower_bound = 0, upper_bound = π, base_name = "i")
    Ω = @variable(model, base_name = "Ω")
    ω = @variable(model, base_name = "ω")
    nu = @variable(model, lower_bound = -2π, upper_bound = 2π, base_name = "nu")

    #rad!!!
    M = @variable(model, lower_bound = 0.0, base_name = "M")
    E = @variable(model, lower_bound = 0.0, base_name= "E")
    
    QxbarX = [
        -sin(Ω)*cos(i)*sin(ω)+cos(Ω)*cos(ω) -sin(Ω)*cos(i)*cos(ω)-cos(Ω)*sin(ω)  sin(Ω)*sin(i)
         cos(Ω)*cos(i)*sin(ω)+sin(Ω)*cos(ω)  cos(Ω)*cos(i)*cos(ω)-sin(Ω)*sin(ω) -cos(Ω)*sin(i)
         sin(i)*sin(ω)                               sin(i)*cos(ω)                              cos(i)
    ]

    #curtis chap 4
    #h^2/mu = p = a (1-e^2)
    r_perifocal = a*(1-e^2) * 1/(1+e*cos(nu)) * [cos(nu); sin(nu); 0]
    
    h = √(GM_EARTH*a*(1-e^2))
    v_perifocal = GM_EARTH / h * [-sin(nu); e + cos(nu); 0]


    @constraint(model, r .== QxbarX * r_perifocal)
    @constraint(model, v .== QxbarX * v_perifocal)
    
    @constraint(model, E - e*sin(E) == M)

    #curtis page 144 & 145
    @constraint(model, nu == 2 * atan(√(1+e)*sin(E/2), √(1-e)*cos(E/2)))

    FullOrbitalParameters(r, v, a, e, i, Ω, ω, nu, M, E)
end