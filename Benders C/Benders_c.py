from load_data import *
from load_data_one_year import *

import gurobipy as gp
from gurobipy import GRB
import math

from Sub_c import *
from Ray_c import *

import psutil
import os

def get_memory_usage():
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / (1024 ** 2)  # Memory in MB

# Before running the optimization
memory_before = get_memory_usage()
print(f"Memory before optimization: {memory_before:.2f} MB")


data = load_data_one()

# MODEL
model = gp.Model("Master problem")

# VARIABLES

# Capacities
p_cap = model.addMVar(10, lb=0, name="p_cap")  #  vtype=GRB.INTEGER,
# 0:wind, 1:solar, 2:elec, 3:reac, 4: dest, 5: batt, 6:h2_sto, 7:raw_MeOH_sto, 8:CO2_sto, 9:MeOH_sto
#p_aux = model.addMVar(10, lb=0, vtype=GRB.INTEGER, name="p_aux")  
intervals = np.array([10, 10, 5, 0.1, 0.1, 1, 1, 20, 20, 20])

# q auxiliary variable
q = model.addMVar(1, lb=0, name="q")

# CONSTRAINTS

# Max capacities
max_cap = model.addConstr( p_cap <= data.max_capacities, name='max_cap')

#aux_cap = model.addConstr( p_cap == p_aux*intervals, name='aux_cap')

# Accelerating convergence value
# min_cap_1 = model.addConstr( p_cap[2] >= 60, name='min_cap_1') # Check these values manually
# min_cap_2 = model.addConstr( p_cap[3] >= 8, name='min_cap_2')
# min_cap_3 = model.addConstr( p_cap[4] >= 5, name='min_cap_3')
 
# OBJECTIVE FUNCTION

objective = p_cap @ data.cost_cap + q
    
model.setObjective(objective, GRB.MINIMIZE)


def solve_Mas(data,res,opt_cut):

    if opt_cut:  
        # Add extreme POINTS
        q_constr = model.addConstr(
            q >= 
            res.p_flow @ (data.Wind * p_cap[0] + data.Solar * p_cap[1])  # Power
            - gp.quicksum(res.max_grid * data.grid_max)                         # Grid
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
            )                                                           # Battery
            + gp.quicksum(    
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
            )                                                           # Battery
            + gp.quicksum(
                - res.Reac_ramp_up * p_cap[3] * data.op_reac[0]
                - res.Reac_ramp_down * p_cap[3] * data.op_reac[1]
            )                                                           # Reactor
            + res.init_MeOH_sto * p_cap[7] * data.raw_MeOH_sto[2]
            + res.fin_MeOH_sto * p_cap[7] * data.raw_MeOH_sto[3]
            + gp.quicksum(
                res.raw_MeOH_min * p_cap[7] * data.raw_MeOH_sto[4]
                - res.raw_MeOH_max * p_cap[7] * data.raw_MeOH_sto[5]
                # - res.C_raw_MeOH_1 * p_cap[7] * data.raw_MeOH_sto[6]
                # - res.C_raw_MeOH_2 * p_cap[7] * data.raw_MeOH_sto[6]
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
                # - res.C_rate_MeOH_1 * p_cap[9] * data.Pure_MeOH_sto[6]
                # - res.C_rate_MeOH_2 * p_cap[9] * data.Pure_MeOH_sto[6]
            )                                                           # Pure Storage
            )
        
    else:
        # Add extreme RAYS
        ray_constr = model.addConstr(
            0 >= 

            res.p_flow @ (data.Wind * p_cap[0] + data.Solar * p_cap[1])  # Power
            - gp.quicksum(res.max_grid * data.grid_max)                         # Grid
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
            )
            + gp.quicksum(    
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
            )
            + gp.quicksum(   
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
            )                                                           # Pure Storage
            )                                                                                                      

    # Write the equations
    model.update()

    model.write('master.lp')
    # model.setParam("OutputFlag", 0)

    model.optimize()

    optimal_objective = model.objVal

    return optimal_objective, p_cap.x


def execute():
    UB = math.inf
    LB = -math.inf
    Delta = 10000
    P_cap = data.max_capacities*0.8
    it = 1

    UBs = []
    LBs = []
    PCs = []

    while UB - LB > Delta:  # it <= 20 :

        print(f"------------------- Run Sub {it}-----------------------")
        opt_cut, sub_obj, duals = sub(data,P_cap)

        if opt_cut:
            UB = min(UB, sub_obj + P_cap @ data.cost_cap)
        else:
            print(f"------------------- Run Ray {it}-----------------------")
            ray_obj, duals = ray(data,P_cap)

        print(f"------------------- Run Master {it}--------------------")
        mas_obj, P_cap = solve_Mas(data, duals, opt_cut)
        print("P_cap: ",P_cap)

        LB = mas_obj

        if opt_cut:
            print(f" SUB It: {it} UB: {UB} LB: {LB} Sub: {sub_obj}")
        else:
            print(f" RAY It: {it} UB: {UB} LB: {LB} Sub: {ray_obj}")

        it += 1

        UBs.append(UB)
        LBs.append(LB)
        PCs.append(P_cap)

    print("Correct Ending")

    val = [UBs,LBs,PCs]

    return val,duals
import time
start_time = time.time()

values, duals = execute()
end_time = time.time()
print(f"Array saved to Excel in {end_time - start_time:.2f} seconds!")

# After running the optimization
memory_after = get_memory_usage()
print(f"Memory after optimization: {memory_after:.2f} MB")

# Memory difference
memory_difference = memory_after - memory_before
print(f"Memory used by optimization: {memory_difference:.2f} MB")


process = psutil.Process(os.getpid())  # Get current process
rss_memory = process.memory_info().rss / (1024 ** 2)  # Convert to MB
print(f"Resident Memory Usage: {rss_memory:.2f} MB")

psutil.virtual_memory()