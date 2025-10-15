"""
This script visualises the results of the Monte Carlo coalitional analysis.

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



""" Preparing data for plotting """

# N_iterations = length(Monte_Carlo_coalitional_results)
N_iterations = 1000
N_players = length(Players)
N_buses = length(Monte_Carlo_coalitional_results[1][:GC_BM_voltages])

plot_results_EV_load_alpha = []
plot_results_Shapley = zeros(Float64, 0, N_players)
plot_results_Nucleolus = zeros(Float64, 0, N_players)
plot_results_leave_one_out = zeros(Float64, 0, N_players)

plot_results_BM_voltages = zeros(Float64, 0, N_buses)
plot_results_DA_p = []
plot_results_DA_p_flex = []
plot_results_BM_p_gen = []
plot_results_BM_p_flex = []

plot_results_EVA_forecasts = zeros(Float64, 0, N_players)
plot_results_EVA_observations = zeros(Float64, 0, N_players)

plot_results_GC_avoided_costs = []


for iteration = 1:N_iterations
    global plot_results_Shapley = vcat(plot_results_Shapley, Monte_Carlo_coalitional_results[iteration][:coalitional_analysis][:Shapley]')
    global plot_results_Nucleolus = vcat(plot_results_Nucleolus, Monte_Carlo_coalitional_results[iteration][:coalitional_analysis][:Nucleolus])
    global plot_results_leave_one_out = vcat(plot_results_leave_one_out, Monte_Carlo_coalitional_results[iteration][:coalitional_analysis][:X_leave_one_out]')


    global plot_results_BM_voltages = vcat(plot_results_BM_voltages, Monte_Carlo_coalitional_results[iteration][:GC_BM_voltages]')
    
    p_gen_DA = Monte_Carlo_coalitional_results[iteration][:GC_pg_DA][1]
    global plot_results_DA_p = vcat(plot_results_DA_p, p_gen_DA)
    p_flex_DA = Monte_Carlo_coalitional_results[iteration][:GC_pg_flex_DA_up][7] - Monte_Carlo_coalitional_results[iteration][:GC_pg_flex_DA_down][7]
    global plot_results_DA_p_flex = vcat(plot_results_DA_p_flex, p_flex_DA)

    p_gen_BM_regultion = Monte_Carlo_coalitional_results[iteration][:GC_pg_BM_up][1] - Monte_Carlo_coalitional_results[iteration][:GC_pg_BM_down][1]
    global plot_results_BM_p_gen = vcat(plot_results_BM_p_gen, p_gen_BM_regultion)
    p_flex_BM_regultion = Monte_Carlo_coalitional_results[iteration][:GC_pg_flex_BM_up][7] - Monte_Carlo_coalitional_results[iteration][:GC_pg_flex_BM_down][7]
    global plot_results_BM_p_flex = vcat(plot_results_BM_p_flex, p_flex_BM_regultion)

    global plot_results_EVA_forecasts = vcat(plot_results_EVA_forecasts, Monte_Carlo_coalitional_results[iteration][:EVA_forecasts]')
    global plot_results_EVA_observations = vcat(plot_results_EVA_observations, Monte_Carlo_coalitional_results[iteration][:EVA_observations]')

    global plot_results_GC_avoided_costs = vcat(plot_results_GC_avoided_costs, Monte_Carlo_coalitional_results[iteration][:GC_avoided_cost])

end

plot_EVA_results_errors = plot_results_EVA_forecasts .- plot_results_EVA_observations


count_negative_avoided_costs = count(i->(i<0), plot_results_GC_avoided_costs)
println()
if count_negative_avoided_costs == 0
    printstyled("Total avoided costs (savings) are positive in all ",N_iterations," simulations."; color = :green)
    println()
else
    printstyled("WARNING: Total avoided costs (savings) are negative in ",count_negative_avoided_costs," out of ",N_iterations," simulations!"; color = :magenta)
    println()
end

neg_avoided_costs_indices = findall(x -> x < 0, plot_results_GC_avoided_costs)

neg_plot_results_Shapley = plot_results_Shapley[neg_avoided_costs_indices,:]
neg_plot_EVA_results_errors = plot_EVA_results_errors[neg_avoided_costs_indices,:]




""" Plotting starts here """


# fz = 18
fz = 22

allocation_plot_range_max = maximum([maximum(plot_results_Shapley),maximum(plot_results_Nucleolus),maximum(plot_results_leave_one_out)])
allocation_plot_range_min = minimum([minimum(plot_results_Shapley),minimum(plot_results_Nucleolus),minimum(plot_results_leave_one_out)])


plt_Monte_Carlo_coalitional_analysis_Shapley = plot(
    title = "Monte Carlo coalitional analysis \n(Shapley value)",
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

plt_density_Shapley = plot(
    title = "Density of the allocation solutions \n(Shapley value)",
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

plt_Monte_Carlo_coalitional_analysis_Nucleolus = plot(
    title = "Monte Carlo coalitional analysis \n(Nucleolus)",
    xlabel = "Forecast error, MW",
    ylabel = "Allocated avoided cost, EUR/h",

    size = (1200,1200), # width and height of the whole plot (in px)
#     aspect_ratio = :equal,

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

plt_density_Nucleolus = plot(
    title = "Density of the allocation solutions \n(Nucleolus)",
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

plt_Monte_Carlo_coalitional_analysis_leave_one_out = plot(
    title = "Monte Carlo coalitional analysis \n(Leave-one-out)",
    xlabel = "Forecast error, MW",
    ylabel = "Allocated avoided cost, EUR/h",

    size = (1200,1200), # width and height of the whole plot (in px)
#     aspect_ratio = :equal,

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

plt_density_leave_one_out = plot(
    title = "Density of the allocation solutions \n(Leave-one-out)",
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


for player_i = 3
# for player_i = 1:length(Players)
    global player_plot = player_i

    scatter!(plt_Monte_Carlo_coalitional_analysis_Shapley,
            plot_EVA_results_errors[:,player_i],
            plot_results_Shapley[:,player_i],
            label = "EVA "*string(Players[player_i]),
            # color = :black,
            # color = :gray,
            # color = palette(:tab10)[5],
            color = iterative_player_colors[player_i],
            markersize = 8,
            markerstrokewidth = 0,
            alpha = 0.5
    )

    if !isempty(neg_plot_results_Shapley)
        scatter!(plt_Monte_Carlo_coalitional_analysis,
            neg_plot_EVA_results_errors[:,player_i],
            neg_plot_results_Shapley[:,player_i],
            label = "Markets with negative savings",
            # color = :black,
            # color = :gray,
            color = palette(:tab10)[4],
            markersize = 7,
            markerstrokewidth = 1,
            # alpha = 0.5,
            m = :x
        )
    end

    density!(plt_density_Shapley,
        bandwidth = 0.7, # kernel bandwidth in kernel density estimation (KDE)

        [plot_results_Shapley[:,player_i]], 
        label = "EVA "*string(Players[player_i]),
        xlabel = "Allocated avoided cost (EUR/h)",
        ylabel = "",
        color = iterative_player_colors[player_i],
        w = 6,

        fillrange = 0,
        fill = (iterative_player_colors[player_i], 0.4),
    )
    
    scatter!(plt_Monte_Carlo_coalitional_analysis_Nucleolus,
        plot_EVA_results_errors[:,player_i],
        plot_results_Nucleolus[:,player_i],
        label = "EVA "*string(Players[player_i]),
        # color = :black,
        # color = :gray,
        # color = palette(:tab10)[5],
        color = iterative_player_colors[player_i],
        markersize = 8,
        markerstrokewidth = 0,
        alpha = 0.5
    )

    density!(plt_density_Nucleolus,
        bandwidth = 0.7, # kernel bandwidth in kernel density estimation (KDE)

        [plot_results_Nucleolus[:,player_i]], 
        label = "EVA "*string(Players[player_i]),
        xlabel = "Allocated avoided cost (EUR/h)",
        ylabel = "",
        color = iterative_player_colors[player_i],
        w = 6,

        fillrange = 0,
        fill = (iterative_player_colors[player_i], 0.4),
    )

        scatter!(plt_Monte_Carlo_coalitional_analysis_leave_one_out,
        plot_EVA_results_errors[:,player_i],
        plot_results_leave_one_out[:,player_i],
        label = "EVA "*string(Players[player_i]),
        # color = :black,
        # color = :gray,
        # color = palette(:tab10)[5],
        color = iterative_player_colors[player_i],
        markersize = 8,
        markerstrokewidth = 0,
        alpha = 0.5
    )

    density!(plt_density_leave_one_out,
        bandwidth = 0.7, # kernel bandwidth in kernel density estimation (KDE)

        [plot_results_leave_one_out[:,player_i]], 
        label = "EVA "*string(Players[player_i]),
        xlabel = "Allocated avoided cost (EUR/h)",
        ylabel = "",
        color = iterative_player_colors[player_i],
        w = 6,

        fillrange = 0,
        fill = (iterative_player_colors[player_i], 0.4),
    )

end


display(plt_Monte_Carlo_coalitional_analysis_Shapley)

savefig(plt_Monte_Carlo_coalitional_analysis_Shapley,"../../results/plt_Monte_Carlo_coalitional_analysis_Shapley_v1_palyer"*string(player_plot)*".svg")
savefig(plt_Monte_Carlo_coalitional_analysis_Shapley,"../../results/plt_Monte_Carlo_coalitional_analysis_Shapley_v1_palyer"*string(player_plot)*".png")
savefig(plt_Monte_Carlo_coalitional_analysis_Shapley,"../../results/plt_Monte_Carlo_coalitional_analysis_Shapley_v1_palyer"*string(player_plot)*".pdf")


display(plt_density_Shapley)

savefig(plt_density_Shapley,"../../results/plt_Monte_Carlo_density_Shapley_v1_palyer"*string(player_plot)*".svg")
savefig(plt_density_Shapley,"../../results/plt_Monte_Carlo_density_Shapley_v1_palyer"*string(player_plot)*".png")
savefig(plt_density_Shapley,"../../results/plt_Monte_Carlo_density_Shapley_v1_palyer"*string(player_plot)*".pdf")



display(plt_Monte_Carlo_coalitional_analysis_Nucleolus)

savefig(plt_Monte_Carlo_coalitional_analysis_Nucleolus,"../../results/plt_Monte_Carlo_coalitional_analysis_Nucleolus_v1_palyer"*string(player_plot)*".svg")
savefig(plt_Monte_Carlo_coalitional_analysis_Nucleolus,"../../results/plt_Monte_Carlo_coalitional_analysis_Nucleolus_v1_palyer"*string(player_plot)*".png")
savefig(plt_Monte_Carlo_coalitional_analysis_Nucleolus,"../../results/plt_Monte_Carlo_coalitional_analysis_Nucleolus_v1_palyer"*string(player_plot)*".pdf")

display(plt_density_Nucleolus)

savefig(plt_density_Nucleolus,"../../results/plt_Monte_Carlo_density_Nucleolus_v1_palyer"*string(player_plot)*".svg")
savefig(plt_density_Nucleolus,"../../results/plt_Monte_Carlo_density_Nucleolus_v1_palyer"*string(player_plot)*".png")
savefig(plt_density_Nucleolus,"../../results/plt_Monte_Carlo_density_Nucleolus_v1_palyer"*string(player_plot)*".pdf")



display(plt_Monte_Carlo_coalitional_analysis_leave_one_out)

savefig(plt_Monte_Carlo_coalitional_analysis_leave_one_out,"../../results/plt_Monte_Carlo_coalitional_analysis_leave_one_out_v1_palyer"*string(player_plot)*".svg")
savefig(plt_Monte_Carlo_coalitional_analysis_leave_one_out,"../../results/plt_Monte_Carlo_coalitional_analysis_leave_one_out_v1_palyer"*string(player_plot)*".png")
savefig(plt_Monte_Carlo_coalitional_analysis_leave_one_out,"../../results/plt_Monte_Carlo_coalitional_analysis_leave_one_out_v1_palyer"*string(player_plot)*".pdf")

display(plt_density_leave_one_out)

savefig(plt_density_leave_one_out,"../../results/plt_Monte_Carlo_density_leave_one_out_v1_palyer"*string(player_plot)*".svg")
savefig(plt_density_leave_one_out,"../../results/plt_Monte_Carlo_density_leave_one_out_v1_palyer"*string(player_plot)*".png")
savefig(plt_density_leave_one_out,"../../results/plt_Monte_Carlo_density_leave_one_out_v1_palyer"*string(player_plot)*".pdf")






plt_avoided_costs = plot(
    title = "Density of the total avoided costs",

    size = (1200,1200),

    xtickfontsize = fz, ytickfontsize = fz,
    fontfamily = "Courier", 
    titlefontsize = fz,
    xguidefontsize = fz,
    yguidefontsize = fz,
    legendfontsize = fz,

    legend = false,
    
    framestyle = :box,
    margin = 5mm,

    # ylim = (0, 0.01)

)

density!(plt_avoided_costs,

    bandwidth = 0.7, # kernel bandwidth in kernel density estimation (KDE)

    [plot_results_GC_avoided_costs], 
    # title = "Distribution of avoided costs over all simulations",
    xlabel = "Total avoided cost (EUR/h)",
    ylabel = "Density",
    color = :grey,
    w = 6,
    # fillalpha = 0.3,
    # fillcolor = :blue,
    fillrange = 0,
    fill = (:grey, 0.4),


)

display(plt_avoided_costs)




# # Let's analyse the worst market instances with negative savings:

results_to_sort = hcat(plot_results_GC_avoided_costs, plot_results_EVA_forecasts, plot_results_EVA_observations, plot_results_Shapley)
println("Sorted results: avoided cost, [forecasts x3], [observations x3], [Shapley x3]")
sorted_results = sortslices(results_to_sort, dims = 1, by = row -> row[1])
