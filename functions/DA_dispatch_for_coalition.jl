"""
This code formulates the OPF optimisation model corresponding to the day-ahead dispatch for a coalition of EVAs.
EVAs within a certain coalition are assumed to provide their forecasts.
For other EVAs, DSO uses the available “status quo” forecast, which is typically less accurate.

Andrey Churkin https://andreychurkin.ru/

"""



if optimize_with == "Gurobi"
    optimizer = Gurobi.Optimizer
#     optimizer = optimizer_with_attributes(() -> Gurobi.Optimizer(), "NonConvex" => 2, "OutputFlag" => 0)
elseif optimize_with == "Ipopt"
    optimizer = Ipopt.Optimizer
else
    println()
    printstyled("Warning: optimizer was not specified correctly!"; color = :magenta)
    println()
end


global OPFmodel_DA = Model(optimizer)

set_silent(OPFmodel_DA) # Suppressing all solver output (including parameter changes messages). Use "unset_silent(model)" to restore output

# set_optimizer_attributes(OPFmodel_DA,"tol" => 1e-05) # Ipopt
set_optimizer_attributes(OPFmodel_DA,"OutputFlag" => 0) # Gurobi
set_optimizer_attributes(OPFmodel_DA,"NonConvex" => 2) # Gurobi
# set_optimizer_attributes(OPFmodel_DA,"MIPGap" => 10^-6) # Gurobi

@variable(OPFmodel_DA, network_data["gen"][string(i)]["pmin"] <= pg[i = 1:length(network_data["gen"])] <= network_data["gen"][string(i)]["pmax"]) # generation power (real) at node i
@variable(OPFmodel_DA, network_data["gen"][string(i)]["qmin"] <= qg[i = 1:length(network_data["gen"])] <= network_data["gen"][string(i)]["qmax"]) # generation power (imaginary) at node i
@variable(OPFmodel_DA, p[i = 1:length(network_data["bus"])]) # active load at node i
@variable(OPFmodel_DA, q[i = 1:length(network_data["bus"])]) # reactive load at node i
@variable(OPFmodel_DA, network_data["bus"][string(i)]["vmin"]^2 <= w[i = 1:length(network_data["bus"])] <= network_data["bus"][string(i)]["vmax"]^2) # squared voltage magnitude of node i
@variable(OPFmodel_DA, IL[l = 1:length(network_data["branch"])]) # current magnitute squared
@variable(OPFmodel_DA, pl[l = 1:length(network_data["branch"])]) # active power flow from i to j
@variable(OPFmodel_DA, ql[l = 1:length(network_data["branch"])]) # reactive power flow from i to j

@variable(OPFmodel_DA, 0 <= pg_flex_down[i = 1:length(flex_data["bus"])] <= flex_data["bus"][string(i)]["flex_p_down"])
@variable(OPFmodel_DA, 0 <= pg_flex_up[i = 1:length(flex_data["bus"])] <= flex_data["bus"][string(i)]["flex_p_up"])
@variable(OPFmodel_DA, 0 <= qg_flex_down[i = 1:length(flex_data["bus"])] <= flex_data["bus"][string(i)]["flex_q_down"])
@variable(OPFmodel_DA, 0 <= qg_flex_up[i = 1:length(flex_data["bus"])] <= flex_data["bus"][string(i)]["flex_q_up"])

bus_loads_p = zeros(length(network_data["bus"])) # a vector of nodal actve loads
bus_loads_q = zeros(length(network_data["bus"])) # a vector of nodal reactve loads

for i = 1:length(network_data["load"])
        bus_loads_p[network_data["load"][string(i)]["load_bus"]] = network_data["load"][string(i)]["pd"]
        bus_loads_q[network_data["load"][string(i)]["load_bus"]] = network_data["load"][string(i)]["qd"]
end

@constraint(OPFmodel_DA, [i = 1:length(network_data["bus"])], p[i] == bus_loads_p[i]) # active loads
@constraint(OPFmodel_DA, [i = 1:length(network_data["bus"])], q[i] == bus_loads_q[i]) # reactive loads

# Branch flow equations:
@constraint(OPFmodel_DA, [l = 1:length(network_data["branch"]); network_data["branch"][string(l)]["br_status"] == 1],
        pl[l]
        ==
        p[network_data["branch"][string(l)]["t_bus"]]
        -
        pg_flex_up[network_data["branch"][string(l)]["t_bus"]]
        +
        pg_flex_down[network_data["branch"][string(l)]["t_bus"]]
        +
        network_data["branch"][string(l)]["br_r"]*IL[l]
        +
        sum(pl[k] for k = 1:1:length(network_data["branch"]) if network_data["branch"][string(k)]["f_bus"] == network_data["branch"][string(l)]["t_bus"])
        +
        # # If in the coalition, then use "load_forecast":
        EVA_data["bus"][string(network_data["branch"][string(l)]["t_bus"])]["load_forecast"] * u_buses[coalition, network_data["branch"][string(l)]["t_bus"]]
        +
        # # If not, use “status_quo_forecast”
        EVA_data["bus"][string(network_data["branch"][string(l)]["t_bus"])]["status_quo_forecast"] * (u_buses[coalition, network_data["branch"][string(l)]["t_bus"]] == 0)
        )

@constraint(OPFmodel_DA, [l = 1:length(network_data["branch"]); network_data["branch"][string(l)]["br_status"] == 1],
        ql[l]
        ==
        q[network_data["branch"][string(l)]["t_bus"]]
        -
        qg_flex_up[network_data["branch"][string(l)]["t_bus"]]
        +
        qg_flex_down[network_data["branch"][string(l)]["t_bus"]]
        +
        network_data["branch"][string(l)]["br_x"]*IL[l]
        +
        sum(ql[k] for k = 1:1:length(network_data["branch"]) if network_data["branch"][string(k)]["f_bus"] == network_data["branch"][string(l)]["t_bus"]) 
)

# Voltage relation:
@constraint(OPFmodel_DA, [l = 1:length(network_data["branch"]); network_data["branch"][string(l)]["br_status"] == 1],
        w[network_data["branch"][string(l)]["t_bus"]]
        ==
        w[network_data["branch"][string(l)]["f_bus"]]
        +
        (network_data["branch"][string(l)]["br_r"]^2 + network_data["branch"][string(l)]["br_x"]^2)*IL[l]
        -
        2*(network_data["branch"][string(l)]["br_r"]*pl[l] + network_data["branch"][string(l)]["br_x"]*ql[l]) 
)

# Power flows and current relation (the exact formulation):
@constraint(OPFmodel_DA, [l = 1:length(network_data["branch"]); network_data["branch"][string(l)]["br_status"] == 1],
        pl[l]^2 + ql[l]^2
        ==
        IL[l]*w[network_data["branch"][string(l)]["f_bus"]] 
)

# Special conditions for the interface node
if size(gen_line) == ()
        @constraint(OPFmodel_DA, pg[1] == pl[gen_line])
        @constraint(OPFmodel_DA, qg[1] == ql[gen_line])
else
        @constraint(OPFmodel_DA, pg[1] == sum(pl[line] for line in gen_line)) # !? check if this is correct for such cases !?
        @constraint(OPFmodel_DA, qg[1] == sum(ql[line] for line in gen_line))
end

if primary_substation_voltage_fixed == true
    @constraint(OPFmodel_DA, w[gen_bus] == 1) # <--- fix the voltage (this may change network's flexibility significantly)
end

if include_line_limits == true # imposing MVA limits for branches
    @constraint(OPFmodel_DA, [l = 1:length(network_data["branch"]); network_data["branch"][string(l)]["rate_a"] >= 10^-4],
                    pl[l]^2 + ql[l]^2 <= network_data["branch"][string(l)]["rate_a"]^2)
end

@objective(OPFmodel_DA, Min, pg[1]*TSO_data["cost_DA"]
            + sum(pg_flex_up[i]*flex_data["bus"][string(i)]["cost_DA_up"]
                    + pg_flex_down[i]*flex_data["bus"][string(i)]["cost_DA_down"]
                    for i = 1:length(flex_data["bus"]))
)
