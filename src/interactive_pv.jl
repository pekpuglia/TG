using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
include("TG.jl")
using .TG
using GLMakie
using LinearAlgebra
using Setfield
include("sample_orbits.jl")
##
function initial_guess_initorb!(orb1, tf, r, v)
    prop = Propagators.init(Val(:TwoBody), orb1)

    N = size(r)[2]

    for i = 1:N
        rguess, vguess = Propagators.propagate!(prop, (i-1)*tf/(N-1))
        set_start_value.(r[:, i], rguess)
        set_start_value.(v[:, i], vguess)
    end
end

struct Transfer
    X1::Vector
    X2::Vector
    mu::Float64
    transfer_time
    rcoast::Array
    vcoast::Array
    deltaVmags::Vector
    deltaVdirs::Matrix
    dts::Vector
end

struct ScaledTransfer
    transfer::Transfer
    L #always include?
    T
end

#scale options?
function lambert_transfer_model(X1, X2, deltat, MU, N)
    r1 = X1[1:3]
    v1 = X1[4:6]
    r2 = X2[1:3]
    v2 = X2[4:6]

    model = Model(Ipopt.Optimizer)

    r, v = add_coast_segment(model, deltat, N, "")

    dV = @variable(model, [1:2], lower_bound=0)
    dVdirs = @variable(model, [1:3, 1:2])
    @constraint(model, [i=1:2], dVdirs[:, i]' * dVdirs[:, i] == 1)

    @constraint(model, r[:, 1] .== r1)
    @constraint(model, r[:, end] .== r2)

    @constraint(model, dV[1] * dVdirs[:, 1] == v[:, 1] - v1)
    @constraint(model, dV[2] * dVdirs[:, 2] == (v[:, end] - v2))

    @objective(model, MIN_SENSE, sum(dV))

    model, Transfer(X1, X2, MU, deltat, r, v, dV, dVdirs, [deltat])
end
##
case_ind = 1
orb1, orb2 = ORBIT_STARTS[case_ind], ORBIT_ENDS[case_ind]
r1, v1 = kepler_to_rv(orb1)
r2, v2 = kepler_to_rv(orb2)

tf1 = orbital_period((orb1.a+orb2.a)/2, GM_EARTH) / 2

model = Model(Ipopt.Optimizer)

model, model_transfer = lambert_transfer_model([r1; v1], [r2; v2], tf1, GM_EARTH, 100)

initial_guess_initorb!(orb1, tf1, model_transfer.rcoast, model_transfer.vcoast)
##
optimize!(model)
##
solved_r = value.(model_transfer.rcoast)
solved_v = value.(model_transfer.vcoast)
solved_orb = rv_to_kepler(solved_r[:, 1], solved_v[:, 1])
solved_prop = Propagators.init(Val(:TwoBody), solved_orb)
##
f, ax3d = plot_orbit(orb1, orb2)
add_discretized_trajectory!(ax3d, solved_r)
f
##
save_with_views!(ax3d, f, "results/$(PREFIXES[case_ind])")
##
deltav1 = solved_v[:, 1] - v1
deltav2 = v2 - solved_v[:, end]
##
tspan, ppdot = ppdot_deltavs(solved_prop, deltav1, deltav2, tf1, 100)
normp = norm.(getindex.(ppdot, Ref(1:3)))
normpdot = [dot(ppdoti[1:3], ppdoti[4:6]) / norm(ppdoti[1:3]) for ppdoti in ppdot]
##
diagnose_ppdot(normp, normpdot)
##
f = Figure()
ax1 = Axis(f[1, 1], xlabel = "t (s)", ylabel = "|p|", title="Diagnostic: "*string(diagnose_ppdot(normp, normpdot)))
lines!(ax1, tspan, normp)
vlines!(ax1, tf1, linestyle=:dash, color=:gray)
ax2 = Axis(f[2, 1], xlabel = "t (s)", ylabel = L"d |p| / dt")
lines!(ax2, tspan, normpdot)
vlines!(ax2, tf1, linestyle=:dash, color=:gray)
f
##
save("results/"*PREFIXES[case_ind]*"_primer_vector.png", f)
##
