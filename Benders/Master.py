from load_data import *
import gurobipy as gp
from gurobipy import GRB
import time
import math
from sub_results import *

def master(data,res):

# MODEL
    model = gp.Model("Master problem")

# VARIABLES

    # Capacities
    p_cap = model.addMVar(10, lb=0, name="p_cap")
    # 0:wind, 1:solar, 2:elec, 3:reac, 4: dest, 5: batt, 6:h2_sto, 7:raw_MeOH_sto, 8:CO2_sto, 9:MeOH_sto

    # q auxiliary variable
    q = model.addMVar(1, lb=0, name="q")


# CONSTRAINTS

    # Max capacities
    max_cap = model.addConstr( p_cap <= data.max_capacities, name='max_cap')

    min_cap = model.addConstr( p_cap >= data.max_capacities*0.7, name='min_cap')

    # Objective function of the sub problem
    q_constr = model.addConstr(
        q >= 

        res.p_flow @ (data.Wind * p_cap[0] + data.Solar * p_cap[1])  # Power
        - sum(res.max_grid * data.grid_max)                         # Grid
        + res.init_bat * p_cap[5] * data.Battery[2] + res.fin_bat * p_cap[5] * data.Battery[3]
        + gp.quicksum(
            res.bat_min * p_cap[5] * data.Battery[4]
            - res.bat_max * p_cap[5] * data.Battery[5]
            - res.C_bat_1 * p_cap[5] * data.Battery[6]
            - res.C_bat_2 * p_cap[5] * data.Battery[6]
        )                                                           # Battery
        + gp.quicksum(
            -res.elec_max * p_cap[2]
            + res.elec_min * p_cap[2] * data.op_elec[2]
            - res.elec_ramp_up * p_cap[2] * data.op_elec[0]
            - res.elec_ramp_down * p_cap[2] * data.op_elec[1]
        )                                                           # Electrolyzer
        + res.balance @ data.Demand                                 # Demand
        + res.init_h2 * p_cap[6] * data.H2storage[2]
        + res.fin_h2 * p_cap[6] * data.H2storage[3]
        + gp.quicksum(
            res.h2_min * p_cap[6] * data.H2storage[4]
            - res.h2_max * p_cap[6] * data.H2storage[5]
            - res.C_rate_h2_1 * p_cap[6] * data.H2storage[6]
            - res.C_rate_h2_2 * p_cap[6] * data.H2storage[6]
        )                                                           # H2 Storage
        + res.init_CO2 * p_cap[8] * data.CO2storage[2]
        + res.fin_CO2 * p_cap[8] * data.CO2storage[3]
        + gp.quicksum(
            res.CO2_min * p_cap[8] * data.CO2storage[4]
            - res.CO2_max * p_cap[8] * data.CO2storage[5]
            + res.CO2_str_min * data.CO2storage[6]
            - res.CO2_str_max * data.CO2storage[7]
        )                                                           # CO2 Storage
        + gp.quicksum(
            -res.reac_max * p_cap[3]
            + res.reac_min * p_cap[3] * data.op_reac[2]
            - res.Reac_ramp_up * p_cap[3] * data.op_reac[0]
            - res.Reac_ramp_down * p_cap[3] * data.op_reac[1]
        )                                                           # Reactor
        + res.init_MeOH_sto * p_cap[7] * data.raw_MeOH_sto[2]
        + res.fin_MeOH_sto * p_cap[7] * data.raw_MeOH_sto[3]
        + gp.quicksum(
            res.raw_MeOH_min * p_cap[7] * data.raw_MeOH_sto[4]
            - res.raw_MeOH_max * p_cap[7] * data.raw_MeOH_sto[5]
            - res.C_raw_MeOH_1 * p_cap[7] * data.raw_MeOH_sto[6]
            - res.C_raw_MeOH_2 * p_cap[7] * data.raw_MeOH_sto[6]
        )                                                           # Raw Storage
        + gp.quicksum(
            - res.dest_max * p_cap[4]
            + res.dest_min * p_cap[4] * data.op_dest[2]
        )
        + gp.quicksum(
            - res.dest_ramp_up * p_cap[4] * data.op_dest[0]
            - res.dest_ramp_down * p_cap[4] * data.op_dest[1]
        )                                                           # Distillator
        + res.init_pure * p_cap[9] * data.Pure_MeOH_sto[2]
        + res.fin_pure * p_cap[9] * data.Pure_MeOH_sto[3]
        + gp.quicksum(
            res.MeOH_min * p_cap[9] * data.Pure_MeOH_sto[4]
            - res.MeOH_max * p_cap[9] * data.Pure_MeOH_sto[5]
            - res.C_rate_MeOH_1 * p_cap[9] * data.Pure_MeOH_sto[6]
            - res.C_rate_MeOH_2 * p_cap[9] * data.Pure_MeOH_sto[6]
        )
        )

# OBJECTIVE FUNCTION

    objective = p_cap @ data.cost_cap + q
     
    model.setObjective(objective, GRB.MINIMIZE)

    # Start timer
    start_time = time.time()

    # Write the equations
    # model.update()

    # model.write('master.lp')

    # Optimize the Gurobi model
    #model.setParam('BarConvTol', 1e-4) #(to reduce Crossover time)
    model.setParam('Method', 2)  # Barrier method
    model.setParam('Crossover', 0)
    model.setParam('MIPGap', 0.1)  # Allow a 0.1% gap
    model.optimize()

    # End timer
    end_time = time.time()

    # Calculate elapsed time
    elapsed_time = end_time - start_time
    
    # Print optimal objective value
    optimal_objective = model.objVal

# Extract P_cap

    P_cap = p_cap.x

    return optimal_objective, elapsed_time, P_cap 
