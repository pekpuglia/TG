using CasADi
##
x = SX("x")
y = SX("y")
α = 1
b = 100
f = (α - x)^2 + b*(y - x^2)^2

nlp = Dict(
    "x" => vcat([x ; y]), 
    "f" => f,
    "g" => vcat([
        x - y
    ]),

    );
S = casadi.nlpsol("S", "ipopt", nlp);

sol = S(x0 = [0, 0], lbg = [1.0]);

println("Optimal solution: x = ", sol["x"].toarray()[1], ", y = ", sol["x"].toarray()[2])