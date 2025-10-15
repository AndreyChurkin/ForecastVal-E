"""
This testing script launches 'Nash_Bargaining.jl' for several allocation problems.

Andrey Churkin https://andreychurkin.ru/

"""

"""
Note: 
This script can take the full characteristic function (all 2^N coalitions) to derive the feasible set as the Core of the game.
The coalitional structure and the order of coalitions in the characteristic function `V` cannot be arbitrary.
The order must match the structure defined in the `Shapley_matrix` function.
"""

"""
Note: 
This Nash bargaining formulation considers value games and expects positive input values to be allocated.
That is, it is expected that the players achieve some positive profits or savings together.
"""



cd(dirname(@__FILE__))
pwd()

using JuMP, Ipopt
using Suppressor



# # Functions for defining the coalitional structure and calculating the Nucleolus:
include("../functions/valuation/Shapley_matrix.jl")
include("../functions/valuation/Nash_Bargaining.jl")





# # Example #1: Symmetric bargaining problem with 2 players

n = 2 # (number of players)
d = [0,0] # (disagreement points)
V = [0,10,0] # (characteristic function)
X_correct = [5.0,5.0] # (the correct solution)

println()
println("Example #1: Symmetric bargaining problem with 2 players")

X_Nash_bargaining, Nash_bargaining_model = Nash_bargaining(n, d; total_value = V[n], objective_formulation = "log")
# X_Nash_bargaining, Nash_bargaining_model = Nash_bargaining(n, d; total_value = V[n], objective_formulation = "prod")

println("The correct Nash bargaining solution:")
println(X_correct')
println("The computed solution:")
println(X_Nash_bargaining)
if isapprox(vec(X_Nash_bargaining), vec(X_correct); atol=1e-4)
    printstyled("✅ The Nash bargaining result is correct for this case."; color = :green)
else
    printstyled("❌ WARNING: The computed Nash bargaining solution is not correct for this case!"; color = :magenta)
end
println()




# # Example #2: Asymmetric bargaining problem with 2 players

n = 2 # (number of players)
d = [2,1] # (disagreement points)
V = [2,10,1] # (characteristic function)
X_correct = [5.5,4.5] # (the correct solution)

println()
println("Example #2: Asymmetric bargaining problem with 2 players")

X_Nash_bargaining, Nash_bargaining_model = Nash_bargaining(n, d; total_value = V[n], objective_formulation = "log")
# X_Nash_bargaining, Nash_bargaining_model = Nash_bargaining(n, d; total_value = V[n], objective_formulation = "prod")

println("The correct Nash bargaining solution:")
println(X_correct')
println("The computed solution:")
println(X_Nash_bargaining)
if isapprox(vec(X_Nash_bargaining), vec(X_correct); atol=1e-4)
    printstyled("✅ The Nash bargaining result is correct for this case."; color = :green)
else
    printstyled("❌ WARNING: The computed Nash bargaining solution is not correct for this case!"; color = :magenta)
end
println()
 



# # Example #2: Asymmetric bargaining problem with 3 players

n = 3 # (number of players)
d = [2,1,0] # (disagreement points)
total_value = 30
X_correct = [11.0,10.0,9.0] # (the correct solution)

println()
println("Example #3: Asymmetric bargaining problem with 3 players")

X_Nash_bargaining, Nash_bargaining_model = Nash_bargaining(n, d; total_value = total_value, objective_formulation = "log")
# X_Nash_bargaining, Nash_bargaining_model = Nash_bargaining(n, d; total_value = total_value, objective_formulation = "prod")

println("The correct Nash bargaining solution:")
println(X_correct')
println("The computed solution:")
println(X_Nash_bargaining)
if isapprox(vec(X_Nash_bargaining), vec(X_correct); atol=1e-4)
    printstyled("✅ The Nash bargaining result is correct for this case."; color = :green)
else
    printstyled("❌ WARNING: The computed Nash bargaining solution is not correct for this case!"; color = :magenta)
end
println()




# # Example #3: Asymmetric cooperative with 3 players and the Core constraints

n = 3 # (number of players)
d = [2.0,1.0,0.5] # (disagreement points)
# V = [2.0,5.0,7.0,4.0,1.0,3.5,0.5] # <-- The core constraints are not binding
V = [2.0,5.0,7.0,4.0,1.0,4.0,0.5] # <-- constraint for coalition v({2,3}) is made binding
X_correct = [3.0,2.25,1.75] # (the correct solution)

println()
println("Example #3: Asymmetric cooperative with 3 players and the Core constraints")

X_Nash_bargaining, Nash_bargaining_model = Nash_bargaining(n, d; V = V, objective_formulation = "log")
# X_Nash_bargaining, Nash_bargaining_model = Nash_bargaining(n, d; total_value = total_value, objective_formulation = "prod")

println("The correct Nash bargaining solution:")
println(X_correct')
println("The computed solution:")
println(X_Nash_bargaining)
if isapprox(vec(X_Nash_bargaining), vec(X_correct); atol=1e-3)
    printstyled("✅ The Nash bargaining result is correct for this case."; color = :green)
else
    printstyled("❌ WARNING: The computed Nash bargaining solution is not correct for this case!"; color = :magenta)
end
println()