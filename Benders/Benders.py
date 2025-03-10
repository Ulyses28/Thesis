import numpy as np
from load_data import *
from load_data_one_year import *

import gurobipy as gp
from gurobipy import GRB
import math

from Sub import *
from Ray import *


data = load_data_one()

def Benders(data):

    # MODEL
    model = gp.Model("Master problem")

    # VARIABLES

    # Capacities
    p_cap = model.addMVar(10, lb=0, name="p_cap")  #  vtype=GRB.INTEGER,
    # 0:wind, 1:solar, 2:elec, 3:reac, 4: dest, 5: batt, 6:h2_sto, 7:raw_MeOH_sto, 8:CO2_sto, 9:MeOH_sto
    # p_aux = model.addMVar(10, lb=0, vtype=GRB.INTEGER, name="p_aux")  
    intervals = np.array([10, 10, 5, 1, 1, 10, 1, 120, 120, 120])

    # q auxiliary variable
    q = model.addMVar(1, lb=0, name="q")

    # CONSTRAINTS

    # Max capacities
    max_cap = model.addConstr( p_cap <= data.max_capacities, name='max_cap')

    # aux_cap = model.addConstr( p_cap == p_aux*intervals, name='aux_cap') # Force the capacities to be integer values
    
    # OBJECTIVE FUNCTION

    objective = p_cap @ data.cost_cap+ q
        
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
                    - res.reac_max * p_cap[3]
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
        P_cap = data.max_capacities
        it = 1

        UBs = []
        LBs = []
        PCs = []

        while it <= 50: #UB - LB > Delta:  # it <= 20 :

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

        bounds = [UBs,LBs,PCs]

        # Create DataFrame for UBs, LBs, and PCs
        df_values = pd.DataFrame({'UBs': UBs, 'LBs': LBs})
        df_PCs = pd.DataFrame(PCs)  # Convert list of lists into DataFrame

        # Merge everything into one table
        df_final = pd.concat([df_values, df_PCs], axis=1)

        # Save to Excel
        df_final.to_excel("output.xlsx", sheet_name="Values", index=False)

        # Fixing duals issue
        if isinstance(duals, list):  
            df_duals = pd.DataFrame({'Duals': duals})
            print("HI")
        else:
            df_duals = pd.DataFrame()  # Empty DataFrame in case of unexpected format

        # Save duals only if it's not empty
        if not df_duals.empty:
            with pd.ExcelWriter("output.xlsx", mode='a', engine='openpyxl') as writer:
                df_duals.to_excel(writer, sheet_name="Duals", index=False)

        print("Results saved to output.xlsx")


        return bounds,duals
    import time
    start_time = time.time()

    bounds, duals = execute()
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Array saved to Excel in {end_time - start_time:.2f} seconds!")


# Extract results

    results = {
        "Capacities": p_cap.x,   
        "CAPEX":p_cap.x * data.cost_cap,
        }

    return elapsed_time, results, bounds, duals