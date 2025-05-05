using SatelliteToolboxBase
using SatelliteToolboxPropagators
include("TG.jl")
using .TG
using Setfield
using LinearAlgebra
using GLMakie
using JuMP
using Ipopt
## estimate of transfer_time
extra_phase = +60
hohmann_start_phase_frac = 1
a1 = 7000e3
a2 = 8000e3

hohmann_time = orbital_period((a1+a2)/2, GM_EARTH) / 2
initial_coast_time = hohmann_start_phase_frac * extra_phase / 360 * orbital_period(a1, GM_EARTH)
terminal_coast_time = (extra_phase*(1 - hohmann_start_phase_frac)) / 360 * orbital_period(a2, GM_EARTH)

transfer_time = initial_coast_time + hohmann_time + terminal_coast_time

##
orb1 = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    a1,
    0.0,
    51 |> deg2rad,
    0    |> deg2rad,
    0     |> deg2rad,
    0     |> deg2rad
)

orb2 = KeplerianElements(
    orb1.t + transfer_time / 86400,
    a2,
    0.0,
    orb1.i,
    orb1.Ω,
    orb1.ω,
    deg2rad(180+extra_phase) + orb1.f
)
##
r1, v1             = kepler_to_rv(orb1)
r2, v2             = kepler_to_rv(orb2)
lambsol = lambert(r1, r2, (orb2.t - orb1.t)*86400, setsilent=false, prograde=false, RAAN = orb1.Ω, i=orb1.i)
##
plot_orbit(
    rv_to_kepler(r1, lambsol.v1),
    rv_to_kepler(r1, v1),
    rv_to_kepler(r2, lambsol.v2),
    rv_to_kepler(r2, v2),
)
##
hohmann_start = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    (orb1.a + orb2.a) / 2,
    (orb2.a - orb1.a) / (orb1.a + orb2.a),
    orb1.i,
    orb1.Ω,
    orb1.ω + deg2rad(hohmann_start_phase_frac*extra_phase),
    0 |> deg2rad
)
hohmann_end = @set hohmann_start.f = π
##
transfer = CoastTransfer(orb1, orb2, prograde=true)
##
f = plot_orbit(orb1, orb2,transfer_start(transfer))
##
# save("./src/late_hohmann.png", f)
##
t_list = transfer_time * (0.0:0.01:1.0)
p_history = p.(Ref(transfer), t_list)
##
pnorm_history = norm.(p_history)
##
#analisar
#descobrir por que N não é inversível no caso Hohmann e o que fazer
f = lines(t_list, pnorm_history)
##
# save("./src/late_hohmann_p_history.png", f)
##
p_x = getindex.(p_history, 1)
p_y = getindex.(p_history, 2)
p_z = getindex.(p_history, 3)
##
lines(p_x, p_y, p_z)
## porkchop plot over 1 synodic period
T1 = orbital_period(a1, GM_EARTH)
T2 = orbital_period(a2, GM_EARTH)
Tsyn = T1*T2 / abs(T2 - T1)
##
_, _, prop = Propagators.propagate(Val(:TwoBody), -transfer_time, orb2)
orb2_init = prop.tbd.orbk
##
N = 100
time_offsets1 =  range(0, 2T2, N)
time_offsets2 = range(0, 2T2, N)
orb1_porkchop = [Propagators.propagate(Val(:TwoBody), t, orb1)[3].tbd.orbk for t in time_offsets1]
orb2_porkchop = [Propagators.propagate(Val(:TwoBody), t, orb2_init)[3].tbd.orbk for t in time_offsets2]
##
function porkchop_transfer(orb1, orb2)
    if orb1.t >= orb2.t
        return nothing
    else
        return LambertTransfer2(orb1, orb2)
    end
end
##
transfer_grid = porkchop_transfer.(orb1_porkchop, permutedims(orb2_porkchop))
##
f = Figure()
ax = Axis(f[1, 1], xlabel= "departure", ylabel = "arrival", title="Cost of transfer", aspect=1.0)
cont = contourf!(ax, 1:N, 1:N, cost.(transfer_grid))
hyperbolic = findall(x -> !isnothing(x) && !x.lambert.is_elliptic, transfer_grid)
scatter!(ax, getfield.(hyperbolic, :I) .|> first, getfield.(hyperbolic, :I) .|> last, color=:red)
Colorbar(f[1, 2], cont)
f
## TODO
# add retrograde maneuver?/exclude weird retrograde cases
# optimize initial coast
# better time frames
# check sigma > sigma_par
#adicionar ifs nas funções de stumpff
##
# save("./src/better_porkchop.png", f)