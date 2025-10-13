"""
This script contains functions for rescaling and normalising vectors of allocations.
It can be useful to avoid negative payments in allocation solutions.

Andrey Churkin https://andreychurkin.ru/

"""



function rescale_allocation(X_input)

    X_rescaled_softmax = zeros(length(X_input))
    X_rescaled_additive_shift = zeros(length(X_input))
    X_rescaled_discard_negatives = zeros(length(X_input))
    
    X_input_softmax = X_input / maximum(abs.(X_input)) # Normalising input for softmax
    exp_sum_softmax = sum(exp(xx) for xx in X_input_softmax) # Sum of all exponentials
    for (i,x) in enumerate(X_input_softmax)
        X_rescaled_softmax[i] = exp(x)/exp_sum_softmax .* sum(X_input)
    end

    if any(x -> x < 0, X_input)
        X_rescaled_additive_shift = X_input .- minimum(X_input)
        X_rescaled_additive_shift = X_rescaled_additive_shift ./ sum(X_rescaled_additive_shift) .* sum(X_input)

        X_rescaled_discard_negatives = max.(X_input, 0.0)
        X_rescaled_discard_negatives = X_rescaled_discard_negatives ./ sum(X_rescaled_discard_negatives) .* sum(X_input)
    else # Do nothing:
        X_rescaled_additive_shift = X_input
        X_rescaled_discard_negatives = X_input
    end

    return X_rescaled_softmax, X_rescaled_additive_shift, X_rescaled_discard_negatives

end

# # # Quick test:
# X = [-3.333 14.774 16.883]
# X2_softmax, X2_additive_shift, X2_discard_negatives = rescale_allocation(X)