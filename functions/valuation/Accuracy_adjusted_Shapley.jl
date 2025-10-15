"""
This function calculates the accuracy-adjusted Shapley value for a given forecast-sharing game.

Input: 
    1) value of the grand coalition to allocate
    2) the original Shapley value
    3) forecast errors for all players
    4) parameter α – determines the share of the grand coalition value that is allocated based on forecast accuracy
    5) parameter ρ – controls the shape of the error penalisation function, 5) parameter τ – sets the error threshold
Output: 
    1) accuracy-adjusted Shapley value
    2) accuracy scores
    3) accuracy-based payments

Andrey Churkin https://andreychurkin.ru/

"""



function AA_Shapley(; V_GC, original_Shapley, errors, α, ρ = 2.0, τ = 0.05)

    error_scores = 1 .- (abs.(errors)./τ).^ρ # Computing absolute-error scores:
    accuracy_pot = V_GC * α
    Shapley_pot = V_GC * (1 - α)
    accuracy_bonuses = accuracy_pot * error_scores/sum(error_scores)
    Shapley_residual = original_Shapley * (1 - α)
    accuracy_adjusted_Shapley = Shapley_residual + accuracy_bonuses

    only_accuracy_based_payments = V_GC * error_scores/sum(error_scores)

    if !isapprox(sum(accuracy_adjusted_Shapley), V_GC; atol=1e-5)
    println()
    printstyled("Warning: Accuracy-adjusted Shapley does not allocate the grand coalition value exactly!"; color = :magenta)
    println()
    println("sum(accuracy_adjusted_Shapley) = ",sum(accuracy_adjusted_Shapley),", V_GC = ",V_GC)
    end

    return accuracy_adjusted_Shapley, error_scores, only_accuracy_based_payments

end