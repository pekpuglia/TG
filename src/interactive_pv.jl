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
##
case_ind = 1
orb1, orb2 = ORBIT_STARTS[case_ind], ORBIT_ENDS[case_ind]
r1, v1 = kepler_to_rv(orb1)
r2, v2 = kepler_to_rv(orb2)

tf1 = orbital_period((orb1.a+orb2.a)/2, GM_EARTH)

model = Model(Ipopt.Optimizer)

r, v = add_coast_segment(model, tf1, 200, "")

# dV = @variable(model, [1:2], lower_bound=0, start=0)

@constraint(model, r[:, 1] .== r1)
@constraint(model, r[:, end] .== r2)

# @constraint(model, dV[1]^2 == (v[:, 1] - v1)' * (v[:, 1] - v1))
# @constraint(model, dV[2]^2 == (v[:, end] - v2)' * (v[:, end] - v2))

# @objective(model, MIN_SENSE, sum(dV))

initial_guess_initorb!(orb1, tf1, r, v)
##
optimize!(model)
##
solved_r = value.(r)
solved_v = value.(v)
solved_orb = rv_to_kepler(solved_r[:, 1], solved_v[:, 1])
solved_prop = Propagators.init(Val(:TwoBody), solved_orb)
##
function add_discretized_trajectory!(ax3d, solved_r)
    scatter!(ax3d, solved_r[1, :], solved_r[2, :], solved_r[3, :], color="green")
end
f, ax3d = plot_orbit(orb1, orb2)
add_discretized_trajectory!(ax3d, solved_r)
f
##
function save_with_views!(ax3d, f, prefix)
    save(prefix*"_3d.png", f)
    
    az = [pi/2, 0, 0]
    el = [0, 0, pi/2]
    name = ["y+", "x+", "z+"]
    for (a, e, n) in zip(az, el, name)
        ax3d.azimuth = a
        ax3d.elevation = e
        save(prefix*"_"*n*".png", f)
    end
end
save_with_views!(ax3d, f, "results/$(PREFIXES[case_ind])")
##
function ppdot_dvs(deltav1, deltav2, delta_t, N)
    p0 = deltav1 / norm(deltav1)
    pf = deltav2 / norm(deltav2)
    p0dot = p0dot_tpbvp(p0, pf, delta_t, solved_prop)
    tspan = range(0, delta_t, N)
    ppdot = [TG.Phi_time(solved_prop, t) * [p0; p0dot] for t in tspan]
    tspan, ppdot
end
deltav1 = solved_v[:, 1] - v1
deltav2 = v2 - solved_v[:, end]
##
tspan, ppdot = ppdot_dvs(deltav1, deltav2, tf1, 100)
normp = norm.(getindex.(ppdot, Ref(1:3)))
normpdot = [dot(ppdoti[1:3], ppdoti[4:6]) / norm(ppdoti[1:3]) for ppdoti in ppdot]
##
#doesn't account for continuity yet
@enum PRIMER_DIAGNOSTIC IC_FC IC_LA ED_FC ED_LA MID OPT
function diagnose_ppdot(normp, normp_dot)
    normp_dot0 = normp_dot[1]
    normp_dotf = normp_dot[end]
    
    max_normp_dot = maximum(abs.(normp_dot))
    tol = 1e-4*max_normp_dot

    if normp_dot0 > tol && normp_dotf < -tol
        IC_FC
    elseif normp_dot0 > tol && normp_dotf > tol
        IC_LA
    elseif normp_dot0 < -tol && normp_dotf < -tol
        ED_FC
    elseif normp_dot0 < -tol && normp_dotf > tol
        ED_LA
    elseif any(normp .> 1)
        MID
    else
        OPT
    end
end
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