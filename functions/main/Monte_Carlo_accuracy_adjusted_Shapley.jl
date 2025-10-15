"""
This script calculates accuracy-adjusted Shapley values for the pre-calculated Monte Carlo simulations.
The accuracy filter uses the ex post forecast errors to adjust the original Shapley value allocations.

Andrey Churkin https://andreychurkin.ru/

"""

using JLD2
using Plots
using Plots.PlotMeasures
using StatsPlots

cd(dirname(@__FILE__))
pwd()


# # read Monte Carlo simulations data, define the players:

data = load("../../results/Monte_Carlo_coalitional_analysis_results_1000_positive_runs.jld2")
Monte_Carlo_coalitional_results = data["Monte_Carlo_results_for_export"]

Players = [22 32 14]


N_iterations = length(Monte_Carlo_coalitional_results)
N_players = length(Players)


# # Set up parameters for the accuracy filter:
τ = 0.05 # Positive scale factor for forecast errors (MW)
ρ = 2 # Smoothness used to calculate the forecast error weights
α = 1.0 # Accuracy-protection share (fraction of the grand coalition value reserved for the accuracy-based payments). 30% of payments will be based on the forecast accuracy, and 70% based on the Shapley value



plot_results_Shapley = zeros(Float64, 0, N_players)
plot_EVA_results_errors = zeros(Float64, 0, N_players)
plot_results_accuracy_filtered_Shapley = zeros(Float64, 0, N_players)

for iteration = 1:N_iterations
    Shapley_value = Monte_Carlo_coalitional_results[iteration][:coalitional_analysis][:Shapley]'
    global plot_results_Shapley = vcat(plot_results_Shapley, Shapley_value)
   
    forecasts = Monte_Carlo_coalitional_results[iteration][:EVA_forecasts]'
    observations = Monte_Carlo_coalitional_results[iteration][:EVA_observations]'
    errors = forecasts .- observations

    global plot_EVA_results_errors = vcat(plot_EVA_results_errors, errors)

    GC_avoided_cost = Monte_Carlo_coalitional_results[iteration][:GC_avoided_cost]

    # # Now we compute the accuracy-filtered Shapley value:

    error_scores = 1 .- (abs.(errors)./τ).^ρ # Computing absolute-error scores:
    accuracy_pot = GC_avoided_cost * α
    Shapley_pot = GC_avoided_cost * (1 - α)
    accuracy_bonuses = accuracy_pot * error_scores/sum(error_scores)
    Shapley_residual = Shapley_value * (1 - α)

    accuracy_filtered_Shapley = Shapley_residual + accuracy_bonuses
    global plot_results_accuracy_filtered_Shapley = vcat(plot_results_accuracy_filtered_Shapley, accuracy_filtered_Shapley)

    if !isapprox(sum(accuracy_filtered_Shapley), GC_avoided_cost; atol=1e-5)
        println()
        printstyled("Warning: Accuracy-filtered Shapley does not allocate the grand coalition value exactly!"; color = :magenta)
        println()
        println("iteration = ",iteration)
        println("sum(accuracy_filtered_Shapley) = ",sum(accuracy_filtered_Shapley),", GC_avoided_cost = ",GC_avoided_cost)

    end
    
end




""" Plotting starts here """


# fz = 18
fz = 22

# allocation_plot_range_max = maximum([maximum(plot_results_Shapley),maximum(plot_results_accuracy_filtered_Shapley)])
# allocation_plot_range_min = minimum([minimum(plot_results_Shapley),minimum(plot_results_accuracy_filtered_Shapley)])

allocation_plot_range_max = 50
allocation_plot_range_min = -10

plt_Monte_Carlo_accuracy_filtered_Shapley = plot(
    title = "Monte Carlo coalitional analysis \n(accuracy-filtered Shapley value)",
    xlabel = "Forecast error, MW",
    ylabel = "Allocated avoided cost, EUR/h",

    size = (1200,1200), # width and height of the whole plot (in px)

    # xlim = (1,33),
    # xlim = (minimum(plot_EVA_results_errors) - 0.05, maximum(plot_EVA_results_errors) + 0.05),
    # xlim = (minimum(plot_EVA_results_errors) - 0.01, maximum(plot_EVA_results_errors) + 0.01),
    xlim = (-0.06, 0.06),

    # ylim = (400, 750),
    # ylim = (minimum(plot_results_Shapley) - 10, maximum(plot_results_Shapley) + 10),
    # ylim = (minimum(plot_results_Shapley) - 5, maximum(plot_results_Shapley) + 5),
    ylim = (allocation_plot_range_min - 5, allocation_plot_range_max + 5),


    xtickfontsize = fz, ytickfontsize = fz,
    fontfamily = "Courier", 
    titlefontsize = fz,
    xguidefontsize = fz,
    yguidefontsize = fz,
    legendfontsize = fz,

    # legend = false,
    
    framestyle = :box,
    # margin = 10mm,

    left_margin=19mm, top_margin=10mm, right_margin=10mm, bottom_margin=10mm,

    # minorgrid = :true,
    # minorgrid = :false,

    # xticks = -0.05:0.01:0.05

)

plt_density_accuracy_filtered_Shapley = plot(
    title = "Density of the allocation solutions \n(accuracy-filtered Shapley value)",
    size = (1200,1200),
    bandwidth = 0.7, # kernel bandwidth in kernel density estimation (KDE)
    xtickfontsize = fz, ytickfontsize = fz,
    fontfamily = "Courier", 
    titlefontsize = fz,
    xguidefontsize = fz,
    yguidefontsize = fz,
    legendfontsize = fz,

    framestyle = :box,
    # margin = 10mm,

    left_margin=19mm, top_margin=10mm, right_margin=10mm, bottom_margin=10mm,

    xlim = (allocation_plot_range_min - 5, allocation_plot_range_max + 5),
    ylim = (0,0.15)
)


iterative_player_colors = [palette(:tab10)[2], palette(:tab10)[8], palette(:tab10)[5]]


for player_i = 2
# for player_i = 1:length(Players)
    global player_plot = player_i

    scatter!(plt_Monte_Carlo_accuracy_filtered_Shapley,
            plot_EVA_results_errors[:,player_i],
            plot_results_accuracy_filtered_Shapley[:,player_i],
            label = "EVA "*string(Players[player_i]),
            # color = :black,
            # color = :gray,
            # color = palette(:tab10)[5],
            color = iterative_player_colors[player_i],
            markersize = 8,
            markerstrokewidth = 0,
            alpha = 0.5
    )


    density!(plt_density_accuracy_filtered_Shapley,
        bandwidth = 0.7, # kernel bandwidth in kernel density estimation (KDE)

        [plot_results_accuracy_filtered_Shapley[:,player_i]], 
        label = "EVA "*string(Players[player_i]),
        xlabel = "Allocated avoided cost (EUR/h)",
        ylabel = "",
        color = iterative_player_colors[player_i],
        w = 6,

        fillrange = 0,
        fill = (iterative_player_colors[player_i], 0.4),
    )

end

display(plt_Monte_Carlo_accuracy_filtered_Shapley)

savefig(plt_Monte_Carlo_accuracy_filtered_Shapley,"../../results/plt_Monte_Carlo_accuracy_filtered_Shapley_alpha_"*string(α)*"palyer"*string(player_plot)*".svg")
savefig(plt_Monte_Carlo_accuracy_filtered_Shapley,"../../results/plt_Monte_Carlo_accuracy_filtered_Shapley_alpha_"*string(α)*"palyer"*string(player_plot)*".png")
savefig(plt_Monte_Carlo_accuracy_filtered_Shapley,"../../results/plt_Monte_Carlo_accuracy_filtered_Shapley_alpha_"*string(α)*"palyer"*string(player_plot)*".pdf")

display(plt_density_accuracy_filtered_Shapley)

savefig(plt_density_accuracy_filtered_Shapley,"../../results/plt_density_accuracy_filtered_Shapley_alpha_"*string(α)*"palyer"*string(player_plot)*".svg")
savefig(plt_density_accuracy_filtered_Shapley,"../../results/plt_density_accuracy_filtered_Shapley_alpha_"*string(α)*"palyer"*string(player_plot)*".png")
savefig(plt_density_accuracy_filtered_Shapley,"../../results/plt_density_accuracy_filtered_Shapley_alpha_"*string(α)*"palyer"*string(player_plot)*".pdf")
