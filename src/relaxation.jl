using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using TG
using GLMakie
using LinearAlgebra
using Setfield
##
function plot_orbit_fix(orbs::KeplerianElements...)
    N = 100
    θ = LinRange(0, 2π, N)

    fig = Figure()
    ax3d = Axis3(fig[1, 1])
    
    for (orb, c) in zip(orbs, Makie.wong_colors())
        orbit = ((@set orb.f = theta) for theta in θ) .|> kepler_to_rv .|> first |> stack
        current_pos, current_velocity = kepler_to_rv(orb)
        velocity_arrow_base = current_pos
        velocity_arrow_tip = velocity_arrow_base + current_velocity / √sum(current_velocity .^ 2) * 0.4 * √sum(current_pos .^ 2)
        velocity_arrow_data = [velocity_arrow_base velocity_arrow_tip] |> Matrix
        lines!(ax3d, orbit[1, :], orbit[2, :], orbit[3, :], color=c)
        scatter!(ax3d, current_pos, markersize=20, color=c)
        lines!(ax3d, velocity_arrow_data[1, :], velocity_arrow_data[2, :], velocity_arrow_data[3, :], color=c)
    end
    
    wireframe!(ax3d, Sphere(Point3(0.0), EARTH_EQUATORIAL_RADIUS), color=:cyan, alpha=0.3)
    fig, ax3d
end
##
a1 = 7000e3
a2 = 9000e3

hohmann_time = orbital_period((a1+a2)/2, GM_EARTH) / 2

transfer_time = hohmann_time

orb1 = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    a1,
    0.0,
    51 |> deg2rad,
    0    |> deg2rad,
    0     |> deg2rad,
    0     |> deg2rad
)
prop1 = Propagators.init(Val(:TwoBody), orb1)

orb2 = KeplerianElements(
    orb1.t + transfer_time / 86400,
    a2,
    0.0,
    orb1.i,
    orb1.Ω,
    orb1.ω,
    deg2rad(180) + orb1.f
)

hohmann_orb = KeplerianElements(
    orb1.t,
    (a1+a2)/2,
    (a2-a1) / (a1+a2),
    orb1.i,
    orb1.Ω,
    orb1.ω,
    deg2rad(60)
)

prop2 = Propagators.init(Val(:TwoBody), orb2)
f, ax3d = plot_orbit_fix(orb1, orb2, hohmann_orb)
f
##
save("./report/img/hohmann_condition.png", f)
##
ax3d.azimuth = -π/2
ax3d.elevation = orb1.i
save("./report/img/hohmann_condition_in_plane.png", f)
##
r1, v1 = kepler_to_rv(orb1)
r2, v2 = kepler_to_rv(orb2)
##
function two_body_dyn(X, mu)
    r = X[1:3]
    v = X[4:6]
    [
        v
        -mu/(√(r'*r)^3)*r
    ]
end

function euler(f, X, dt)
    X + f(X)*dt
end

function RK4(f, X, dt)
    k1 = dt*f(X)
    k2 = dt*f(X+k1/2)
    k3 = dt*f(X+k2/2)
    k4 = dt*f(X+k3)
    X + (k1+2k2+2k3+k4)/6
end

function RK8(f, X, dt)
    c = zeros(Float64, (16,))
    c[1] =  0
    c[2] =  .4e-1
    c[3] =  .9648736013787361245235039379666356743708e-1
    c[4] =  .1447310402068104186785255906949953511556
    c[5] =  .576
    c[6] =  .2272326564618766017153738192188229509142
    c[7] =  .5407673435381233982846261807811770490858
    c[8] =  .64
    c[9] =  .48
    c[10] =  .6754e-1
    c[11] =  .25
    c[12] =  .6770920153543242682384311058159603931192
    c[13] =  .8115
    c[14] =  .906
    c[15] =  1
    c[16] =  1

    a = zeros(Float64, (16,16))
    a[2,1] =  .4e-1

    a[3,1] = -.198852731918229097650241511466089129345e-1
    a[3,2] =  .1163726333296965222173745449432724803716

    a[4,1] =  .3618276005170260466963139767374883778890e-1
    a[4,2] =  0
    a[4,3] =  .1085482801551078140088941930212465133667

    a[5,1] =  2.272114264290177409193144938921415409241
    a[5,2] =  0
    a[5,3] = -8.526886447976398578316416192982602292786
    a[5,4] =  6.830772183686221169123271254061186883545

    a[6,1] =  .5094385535389374394512668566783434123978e-1
    a[6,2] =  0
    a[6,3] =  0
    a[6,4] =  .1755865049809071110203693328749561646990
    a[6,5] =  .70229612707574674987780067603244497535e-3

    a[7,1] =  .1424783668683284782770955365543878809824
    a[7,2] =  0
    a[7,3] =  0
    a[7,4] = -.3541799434668684104094753917518523845155
    a[7,5] =  .7595315450295100889001534202778550159932e-1
    a[7,6] =  .6765157656337123215269906939508560510196

    a[8,1] =  .7111111111111111111111111111111111111111e-1
    a[8,2] =  0
    a[8,3] =  0
    a[8,4] =  0
    a[8,5] =  0
    a[8,6] =  .3279909287605898328568406057725491803016
    a[8,7] =  .2408979601282990560320482831163397085872

    a[9,1] =  .7125e-1
    a[9,2] =  0
    a[9,3] =  0
    a[9,4] =  0
    a[9,5] =  0
    a[9,6] =  .3268842451575245554847578757216915662785
    a[9,7] =  .1156157548424754445152421242783084337215
    a[9,8] = -.3375e-1

    a[10,1] =  .4822677322465810178387112087673611111111e-1
    a[10,2] =  0
    a[10,3] =  0
    a[10,4] =  0
    a[10,5] =  0
    a[10,6] =  .3948559980495400110769549704186108167677e-1
    a[10,7] =  .1058851161934658144373823566907778072121
    a[10,8] = -.2152006320474309346664428710937500000000e-1
    a[10,9] = -.1045374260183348238623046875000000000000

    a[11,1] = -.2609113435754923412210928689962011065179e-1
    a[11,2] =  0
    a[11,3] =  0
    a[11,4] =  0
    a[11,5] =  0
    a[11,6] =  .3333333333333333333333333333333333333333e-1
    a[11,7] = -.1652504006638105086724681598195267241410
    a[11,8] =  .3434664118368616658319419895678838776647e-1
    a[11,9] =  .1595758283215209043195814910843067811951
    a[11,10] =  .2140857321828193385584684233447183324979

    a[12,1] = -.362842339625565859076509979091267105528e-1
    a[12,2] =  0
    a[12,3] =  0
    a[12,4] =  0
    a[12,5] =  0
    a[12,6] = -1.096167597427208807028761474420297770752
    a[12,7] =  .1826035504321331052308236240517254331348
    a[12,8] =  .708225444417068325613028685455625123741e-1
    a[12,9] = -.231364701848243126999929738482630407146e-1
    a[12,10] =  .2711204726320932916455631550463654973432
    a[12,11] =  1.308133749422980744437146904349994472286

    a[13,1] = -.5074635056416974879347823927726392374259
    a[13,2] =  0
    a[13,3] =  0
    a[13,4] =  0
    a[13,5] =  0
    a[13,6] = -6.631342198657237090355284142048733580937
    a[13,7] = -.252748010090880105270020973014860316405
    a[13,8] = -.4952612380036095562991116175550167835424
    a[13,9] =  .293252554525388690285739720360003594753
    a[13,10] =  1.440108693768280908474851998204423941413
    a[13,11] =  6.237934498647055877243623886838802127716
    a[13,12] =  .7270192054526987638549835199880202544289

    a[14,1] =  .6130118256955931701496387847232542148725
    a[14,2] =  0
    a[14,3] =  0
    a[14,4] =  0
    a[14,5] =  0
    a[14,6] =  9.088803891640463313341034206647776279557
    a[14,7] = -.407378815629344868103315381138325162923
    a[14,8] =  1.790733389490374687043894756399015035977
    a[14,9] =  .714927166761755073724875250629602731782
    a[14,10] = -1.438580857841722850237810322456327208949
    a[14,11] = -8.263329312064740580595954649844133476994
    a[14,12] = -1.537570570808865115231450725068826856201
    a[14,13] =  .3453832827564871699090880801079644428793

    a[15,1] = -1.211697910343873872490625222495537087293
    a[15,2] =  0
    a[15,3] =  0
    a[15,4] =  0
    a[15,5] =  0
    a[15,6] = -19.05581871559595277753334676575234493500
    a[15,7] =  1.26306067538987510135943101851905310045
    a[15,8] = -6.913916969178458046793476128409110926069
    a[15,9] = -.676462266509498065300115641383621209887
    a[15,10] =  3.367860445026607887090352785684064242560
    a[15,11] =  18.00675164312590810020103216906571965203
    a[15,12] =  6.838828926794279896350389904990814350968
    a[15,13] = -1.031516451921950498420447675652291096155
    a[15,14] =  .4129106232130622755368055554332539084021

    a[16,1] =  2.157389007494053627033175177985666660692
    a[16,2] =  0
    a[16,3] =  0
    a[16,4] =  0
    a[16,5] =  0
    a[16,6] =  23.80712219809580523172312179815279712750
    a[16,7] =  .88627792492165554903036801415266308369
    a[16,8] =  13.13913039759876381480201677314222971522
    a[16,9] = -2.604415709287714883747369630937415176632
    a[16,10] = -5.193859949783872300189266203049579105962
    a[16,11] = -20.41234071154150778768154893536134356354
    a[16,12] = -12.30085625250572261314889445241581039623
    a[16,13] =  1.521553095008539362178397458330791655267
    a[16,14] =  0
    a[16,15] =  0

    b = zeros(Float64, (16,))
    b[1] =  .1458885278405539719101539582255752917034e-1
    b[2] =  0
    b[3] =  0
    b[4] =  0
    b[5] =  0
    b[6] =  0
    b[7] =  0
    b[8] =  .2024197887889332650566666683195656097825e-2
    b[9] =  .2178047084569716646796256135839225745895
    b[10] =  .1274895340854389692868677968654808668201
    b[11] =  .2244617745463131861258531547137348031621
    b[12] =  .1787254491259903095100090833796054447157
    b[13] =  .7594344758096557172908303416513173076283e-1
    b[14] =  .1294845879197561516869001434704642286297
    b[15] =  .2947744761261941714007911131590716605202e-1
    b[16] =  0

    F = zeros(Float64, (length(X),16))

    for i = 1:16
        Y = X + dt * sum(a[i, 1:(i-1)]' .* F[:, 1:(i-1)], dims=2)
        F[:, i] = f(Y)
    end

    X + dt*sum(b' .* F, dims=2)
end
#test with xdot(x) = x
# RK8(xdot, 1, 1) = e

function Xs(f, X0, dt, N, integrator=euler)
    sol = [X0]
    for i = 1:N
        push!(sol, integrator(f, sol[i], dt))
    end
    sol
end
##
function add_coast_segment(model, deltat, N, ind)
    #[(x, y, z), time instants]
    r = @variable(model, [1:3, 1:N+1], base_name="r_coast_$ind")
    v = @variable(model, [1:3, 1:N+1], base_name="v_coast_$ind")


    for i = 1:(N+1)
        # @constraint(model, r[:, i]' * r[:, i] >= EARTH_EQUATORIAL_RADIUS^2)
        @constraint(model, cross(r[:, i], v[:, i])[3] >=0)
    end

    for i = 1:N
        @constraint(model, RK4(X -> two_body_dyn(X, GM_EARTH), [r[:, i]; v[:, i]], deltat/N) .== [r[:, i+1]; v[:, i+1]])
    end

    r, v
end
##
#solve lambert problem to use as initial condition
N = 20

model = Model(Ipopt.Optimizer)

r, v = add_coast_segment(model, transfer_time, N, "lamb")

dt10 = 0/3*transfer_time
dt20 = 3/3*transfer_time

r1_lamb, v1_lamb = Propagators.propagate!(prop1, dt10)
r2_lamb = Propagators.propagate!(prop2, dt10+dt20-transfer_time)[1]

@constraint(model, r[:, 1] .== r1_lamb)
@constraint(model, r[:, end] .== r2_lamb)

#follow natural orbit
for i = 1:size(r)[2]
    rlamb_guess, vlamb_guess = Propagators.propagate!(prop1, dt10 + dt20 * (i/(N+1)))
    set_start_value.(r[:, i], rlamb_guess)
    set_start_value.(v[:, i], vlamb_guess)
end
##
optimize!(model)
##
lamb1 = rv_to_kepler(value.(r[:, 1]), value.(v[:, 1]))
lamb_prop = Propagators.init(Val(:TwoBody), lamb1)
f, ax3d = plot_orbit_fix(orb1, orb2, lamb1)
f
# save("./report/img/hohmann_lambert_guess.png", f)
##
ax3d.azimuth = -π/2
ax3d.elevation = orb1.i
save("./report/img/hohmann_lambert_guess_in_plane.png", f)
##
N=50
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
    "max_iter" => 10_000
    )
)

dt1 = @variable(model, base_name="dt1", lower_bound=0, upper_bound=transfer_time)
dt2 = @variable(model, base_name="dt2", lower_bound=0, upper_bound=transfer_time)

@constraint(model, dt1+dt2 <= transfer_time)

deltaV1mag = @variable(model, lower_bound = 0, base_name="dV1mag")
deltaV2mag = @variable(model, lower_bound = 0, base_name="dV2mag")

deltaV1dir = @variable(model, [1:3], base_name="dV1dir")
@constraint(model, deltaV1dir' * deltaV1dir == 1)
deltaV2dir = @variable(model, [1:3], base_name="dV2dir")
@constraint(model, deltaV2dir' * deltaV2dir == 1)


rcoast1, vcoast1 = add_coast_segment(model, dt1, N, 1)
rcoast2, vcoast2 = add_coast_segment(model, dt2, N, 2)
rcoast3, vcoast3 = add_coast_segment(model, transfer_time-dt1-dt2, N, 3)

@constraint(model, rcoast1[:, 1] .== r1)
@constraint(model, vcoast1[:, 1] .== v1)

@constraint(model, rcoast2[:, 1] .== rcoast1[:, end])
deltaV1 = deltaV1mag * deltaV1dir
@constraint(model, vcoast2[:, 1] .== vcoast1[:, end] + deltaV1)


@constraint(model, rcoast3[:, 1] .== rcoast2[:, end])
deltaV2 = deltaV2mag * deltaV2dir
@constraint(model, vcoast3[:, 1] .== vcoast2[:, end] + deltaV2)


@constraint(model, rcoast3[:, end] .== r2)
@constraint(model, vcoast3[:, end] .== v2)

@objective(model, MIN_SENSE, deltaV1mag + deltaV2mag)

#start condition - set to maneuvers along lambert between start and end
#lambert maneuver
set_start_value(dt1, dt10)
set_start_value(dt2, dt20)

dv1 = value.(v[:, 1]) - v1
set_start_value(deltaV1mag, norm(dv1))
set_start_value.(deltaV1dir, dv1/norm(dv1))

dv2 = value.(v[:, end]) - v2
set_start_value(deltaV2mag, norm(dv2))
set_start_value.(deltaV2dir, dv2/norm(dv2))

for i = 1:(N+1)
    r1i, v1i = Propagators.propagate!(prop1, dt10*(i-1)/N)
    set_start_value.(rcoast1[:, i], r1i)
    set_start_value.(vcoast1[:, i], v1i)
    
    r2i, v2i = Propagators.propagate!(lamb_prop, dt20*(i-1)/N)
    set_start_value.(rcoast2[:, i], r2i)
    set_start_value.(vcoast2[:, i], v2i)
    
    r3i, v3i = Propagators.propagate!(prop2, (transfer_time-dt10-dt20)*((i-1)/N - 1))
    set_start_value.(rcoast3[:, i], r3i)
    set_start_value.(vcoast3[:, i], v3i)
end
##
optimize!(model)
##
deltaV1_solved = value.(deltaV1)
deltaV2_solved = value.(deltaV2)
##
deltaVmodel = objective_value(model)
##
transf_orb1 = rv_to_kepler(value.(rcoast1[:, 1]), value.(vcoast1[:, 1]))
transf_orb2 = rv_to_kepler(value.(rcoast2[:, 1]), value.(vcoast2[:, 1]))
transf_orb3 = rv_to_kepler(value.(rcoast3[:, 1]), value.(vcoast3[:, 1]))
##
v1n = norm(v1)
v2n = norm(v2)
vtransf1 = √(2*GM_EARTH*(1/a1 - 1/(a1+a2)))
vtransf2 = √(2*GM_EARTH*(1/a2 - 1/(a1+a2)))
deltaV_hohmann = v2n - vtransf2 + vtransf1 - v1n
##
solved_r = [value.(rcoast1) value.(rcoast2) value.(rcoast3)]
##
# f = Figure()
# ax3d = Axis3(f[1, 1])

f, ax3d = plot_orbit_fix(orb1, orb2)

scatter!(ax3d, solved_r[1, :], solved_r[2, :], solved_r[3, :], color="green")
scatter!(ax3d, solved_r[1, (N+1):(N+1):end], solved_r[2, (N+1):(N+1):end], solved_r[3, (N+1):(N+1):end], markersize=20, color=:red)
wireframe!(ax3d, Sphere(Point3(0.0), EARTH_EQUATORIAL_RADIUS), color=:cyan, alpha=0.3)

f
##
save("./report/img/hohmann_solved.png", f)
ax3d.azimuth = -π/2
ax3d.elevation = orb1.i
save("./report/img/hohmann_solved_in_plane.png", f)
##
plot_orbit(orb1, orb2, transf_orb1, transf_orb2, transf_orb3)
## primer vector 

p0 = deltaV1_solved / norm(deltaV1_solved)
pf = deltaV2_solved / norm(deltaV2_solved)

inter_impulse_prop = Propagators.init(Val(:TwoBody), transf_orb2)
Phi = TG.Phi_time(inter_impulse_prop, value(dt2))

M = Phi[1:3, 1:3]
N = Phi[1:3, 4:6]

if abs(det(N)) <= 1e-10
    #CHECK THIS IS TRUE
    #only works for coplanar transfers?
    A = N * [r1 v1]
    b = pf - M * p0
    p0dot = [r1 v1] * (A \ b)
else
    p0dot = N \ (pf - M * p0)
end
## plot p and pdot
tspan = range(0, value(dt2), 100)
##
ppdot = [TG.Phi_time(inter_impulse_prop, t) * [p0; p0dot] for t in tspan]
##
normpdot = [dot(ppdoti[1:3], ppdoti[4:6]) / norm(ppdoti[1:3]) for ppdoti in ppdot]
##
f = Figure()
ax1 = Axis(f[1, 1], xlabel = "t (s)", ylabel = "|p|")
lines!(ax1, tspan, norm.(getindex.(ppdot, Ref(1:3))))

ax2 = Axis(f[2, 1], xlabel = "t (s)", ylabel = L"d |p| / dt")
lines!(ax2, tspan, normpdot)

# save("./report/img/hohmann_primer_vector_history.png", f)