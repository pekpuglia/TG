using JuMP
using Ipopt
##
model = Model(Ipopt.Optimizer)
@variable(model, x)
@variable(model, y)
@constraint(model, x+y == 1)
@objective(model, Min, (1 - x)^2 + 100 * (y - x^2)^2)
##
optimize!(model)
##
objective_function_type(model)
##
objective_value(model)
##
value(x), value(y)