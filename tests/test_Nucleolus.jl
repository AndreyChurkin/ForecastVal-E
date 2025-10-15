"""
This test script launches 'Nucleolus.jl' for several well-known cooperative games.
It verifies the correctness of the computed solutions against known results.
The purpose is to check that the linear programs within 'Nucleolus.jl' work correctly and produce valid Nucleolus allocations.

The test cases are taken from: Guajardo, Mario, and Kurt Jörnsten, "Common mistakes in computing the nucleolus," European Journal of Operational Research 241, no. 3 (2015).

Andrey Churkin https://andreychurkin.ru/

"""


cd(dirname(@__FILE__))
pwd()

using JuMP, Gurobi
using Suppressor



# # Functions for defining the coalitional structure and calculating the Nucleolus:
include("../functions/valuation/Shapley_matrix.jl")
include("../functions/valuation/Nucleolus.jl")

"""
Note: 
The coalitional structure and the order of coalitions in the characteristic function `V` cannot be arbitrary.
The order must match the structure defined in the `Shapley_matrix` function.
"""


"""
Note: 
Games with positive characteristic functions are cost games.
In such games, players achieve synergy in cost reduction (lower values due to cooperation) and allocate the total cost.

Games with negative characteristic functions are value games.
In such games, players achieve higher values (e.g., profits) when cooperating.

It is important to select the right sign of the characteristic function.
Otherwise, Nucleolus calculations will be incorrect.

Nucleolus.jl uses the cost formulation for cooperative games.
In cost games, the excess is a measure of how satisfied a coalition S is with the cost allocation X. 
The larger the excess of S, the more satisfied the coalition S is.
"""


# # Example #1: Insurances (Lemaire, 1991):

V = -[46125,69187.5,90000,53812.5,17437.5,30750,5812.5] # (characteristic function)
n = 3 # (number of players)
X_correct = -[52687.5,24468.75,12843.75] # (the correct solution)

println()
println("Example #1: Insurances (Lemaire, 1991), with 3 players")

X_Nucleolus = Nucleolus(n,V)

println("The correct Nucleolus solution:")
println(X_correct')
println("The computed Nucleolus:")
println(X_Nucleolus)
if isapprox(vec(X_Nucleolus), vec(X_correct); atol=1e-8)
    printstyled("✅ The Nucleolus result is correct for this case."; color = :green)
else
    printstyled("❌ WARNING: The computed Nucleolus is not correct for this case!"; color = :magenta)
end
println()




# # Example #2: Mobiles in broadcast transmission (Hasan et al., 2011):

V = [8,8,10,20,19,10,20,19,1,10,13,12,10,12,11] # (characteristic function)
n = 4 # (number of players)
X_correct = [7.5,0.5,1.5,10.5] # (the correct solution)

println()
println("Example #2: Mobiles in broadcast transmission (Hasan et al., 2011), with 4 players")

X_Nucleolus = Nucleolus(n,V)

println("The correct Nucleolus solution:")
println(X_correct')
println("The computed Nucleolus:")
println(X_Nucleolus)
if isapprox(vec(X_Nucleolus), vec(X_correct); atol=1e-8)
    printstyled("✅ The Nucleolus result is correct for this case."; color = :green)
else
    printstyled("❌ WARNING: The computed Nucleolus is not correct for this case!"; color = :magenta)
end
println()




# # Example #3: Manufacturing (Oh and Shin, 2012):

V = [375144,568232,772500,551055,245280,452411,211239] # (characteristic function)
n = 3 # (number of players)
X_correct = [331695.25,233051.25,207753.5] # (the correct solution)

println()
println("Example #3: Manufacturing (Oh and Shin, 2012), with 3 players")

X_Nucleolus = Nucleolus(n,V)

println("The correct Nucleolus solution:")
println(X_correct')
println("The computed Nucleolus:")
println(X_Nucleolus)
if isapprox(vec(X_Nucleolus), vec(X_correct); atol=1e-8)
    printstyled("✅ The Nucleolus result is correct for this case."; color = :green)
else
    printstyled("❌ WARNING: The computed Nucleolus is not correct for this case!"; color = :magenta)
end
println()




# # Example #4: Production and transportation planning by Sakawa et al. (2001):

V = -[0.060,0.378,0.536,0.906,1,0.573,0.747,0.874,0.485,0.144,0.546,0.634,0.255,0.408,0.561,0.182,0.168,0.337,0.706,0.801,0.379,0.538,0.668,0.279,0.03,0.383,0.473,0.083,0.249,0.386,0,0] # (characteristic function)
n = 5 # (number of players)
X_correct = -[0.165000,0.320500,0.084500,0.374500,0.055500] # (the correct solution)

println()
println("Example #4: Production and transportation planning by Sakawa et al. (2001), with 5 players")

X_Nucleolus = Nucleolus(n,V)

println("The correct Nucleolus solution:")
println(X_correct')
println("The computed Nucleolus:")
println(X_Nucleolus)
if isapprox(vec(X_Nucleolus), vec(X_correct); atol=1e-8)
    printstyled("✅ The Nucleolus result is correct for this case."; color = :green)
else
    printstyled("❌ WARNING: The computed Nucleolus is not correct for this case!"; color = :magenta)
end
println()