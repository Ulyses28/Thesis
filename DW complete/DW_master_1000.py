import gurobipy as gp
from gurobipy import GRB
import math
from DW_X1_results import *
from DW_duals_1000 import *

def DW_master(data,P_cap,Xs):

# Master model

    # Master model NOTE (setupMASTER)
    master_DW = gp.Model("DW_master")

    # VARIABLES
    # lambda 
    lambdas = []
    for k in range(len(Xs)):
        lambdas.append(master_DW.addMVar(1, lb=0, name=f"lambda_{k}"))

# MASTER CONSTRAINTS

    lambdas_total = master_DW.addConstr( gp.quicksum(lambdas) == 1)

    # X1 = Xs[k]
    e_bat_SOC = master_DW.addConstr(gp.quicksum((Xs[k].bat_sto[4379] - Xs[k].bat_sto[4380] + Xs[k].bat_min[4380] - Xs[k].bat_max[4380]) * lambdas[k] for k in range(len(Xs))) <= 0)

    e_elec = master_DW.addConstr(gp.quicksum((Xs[k].p_flow[4380] - Xs[k].elec_max[4380] + Xs[k].elec_min[4380] - Xs[k].elec_op[4380] - Xs[k].elec_ramp_up[4379] + Xs[k].elec_ramp_down[4379] + Xs[k].elec_ramp_up[4380] - Xs[k].elec_ramp_down[4380]) * lambdas[k] for k in range(len(Xs))) <= 0, name='e_elec')

    m_h2_sto_SOC = master_DW.addConstr(gp.quicksum((Xs[k].h2_storage[4379] - Xs[k].h2_storage[4380] + Xs[k].h2_min[4380] - Xs[k].h2_max[4380]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_h2_sto_SOC')

    m_CO2_sto_SOC = master_DW.addConstr(gp.quicksum((Xs[k].CO2_sto[4379] - Xs[k].CO2_sto[4380] + Xs[k].CO2_min[4380] - Xs[k].CO2_max[4380]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_CO2_sto_SOC')

    m_rawMeOH_reac = master_DW.addConstr(gp.quicksum((-Xs[k].reac_max[4380] + Xs[k].reac_min[4380] + Xs[k].reac_H2[4380]*data.Reactor[0] + Xs[k].reac_CO2[4380]*data.Reactor[1] + Xs[k].reac_el[4380]*data.Reactor[2] + Xs[k].reac_heat[4380]*data.Reactor[3] - Xs[k].Reac_ramp_up[4379] + Xs[k].Reac_ramp_down[4379] + Xs[k].Reac_ramp_up[4380] - Xs[k].Reac_ramp_down[4380] + Xs[k].compr_MeOH[4380]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_rawMeOH_reac')

    m_rawMeOH_sto_SOC = master_DW.addConstr(gp.quicksum((-Xs[k].rawMeOH_sto[4380] + Xs[k].rawMeOH_sto[4379] + Xs[k].raw_MeOH_min[4380] - Xs[k].raw_MeOH_max[4380]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_rawMeOH_sto_SOC')

    m_MeOH_dest = master_DW.addConstr(gp.quicksum((- Xs[k].dest_max[4380] + Xs[k].dest_min[4380] + Xs[k].dest_rawMeOH[4380]*data.Destilator[0] + Xs[k].dest_water[4380]*data.Destilator[1] + Xs[k].dest_heat[4380]*data.Destilator[3] + Xs[k].dest_sto[4380] - Xs[k].dest_ramp_up[4379] + Xs[k].dest_ramp_down[4379] + Xs[k].dest_ramp_up[4380] - Xs[k].dest_ramp_down[4380]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_MeOH_dest')

    m_MeOH_SOC = master_DW.addConstr(gp.quicksum((Xs[k].pure_sto[4379] - Xs[k].pure_sto[4380]  + Xs[k].MeOH_min[4380] - Xs[k].MeOH_max[4380]) * lambdas[k] for k in range(len(Xs))) <= 0, name="m_MeOH_SOC")

    # OBJECTIVE FUNCTION MASTER

    objective = []
    for k in range(len(Xs)): 
        dual_objective = (
            Xs[k].p_flow @ (data.Wind * P_cap[0] + data.Solar * P_cap[1])
            - sum(Xs[k].max_grid * data.grid_max)                         # Grid
            + Xs[k].init_bat * P_cap[5] * data.Battery[2] + Xs[k].fin_bat * P_cap[5] * data.Battery[3]
            + sum(
                Xs[k].bat_min * P_cap[5] * data.Battery[4]
                - Xs[k].bat_max * P_cap[5] * data.Battery[5]
                - Xs[k].C_bat_1 * P_cap[5] * data.Battery[6]
                - Xs[k].C_bat_2 * P_cap[5] * data.Battery[6]
            )                                                           # Battery
            + sum(
                -Xs[k].elec_max * P_cap[2]
                + Xs[k].elec_min * P_cap[2] * data.op_elec[2]
            )
            + sum(
                - Xs[k].elec_ramp_up * P_cap[2] * data.op_elec[0]
                - Xs[k].elec_ramp_down * P_cap[2] * data.op_elec[1]
            )                                                           # Electrolyzer
            + Xs[k].balance @ data.Demand                                     # Demand
            + Xs[k].init_h2 * P_cap[6] * data.H2storage[2]
            + Xs[k].fin_h2 * P_cap[6] * data.H2storage[3]
            + sum(
                Xs[k].h2_min * P_cap[6] * data.H2storage[4]
                - Xs[k].h2_max * P_cap[6] * data.H2storage[5]
                - Xs[k].C_rate_h2_1 * P_cap[6] * data.H2storage[6]
                - Xs[k].C_rate_h2_2 * P_cap[6] * data.H2storage[6]
            )                                                           # H2 Storage
            + Xs[k].init_CO2 * P_cap[8] * data.CO2storage[2]
            + Xs[k].fin_CO2 * P_cap[8] * data.CO2storage[3]
            + sum(
                Xs[k].CO2_min * P_cap[8] * data.CO2storage[4]
                - Xs[k].CO2_max * P_cap[8] * data.CO2storage[5]
                + Xs[k].CO2_str_min * data.CO2storage[6]
                - Xs[k].CO2_str_max * data.CO2storage[7]
            )                                                           # CO2 Storage
            + sum(
                -Xs[k].reac_max * P_cap[3]
                + Xs[k].reac_min * P_cap[3] * data.op_reac[2]
            )
            + sum(
                - Xs[k].Reac_ramp_up * P_cap[3] * data.op_reac[0]
                - Xs[k].Reac_ramp_down * P_cap[3] * data.op_reac[1]
            )                                                           # Reactor
            + Xs[k].init_MeOH_sto * P_cap[7] * data.raw_MeOH_sto[2]
            + Xs[k].fin_MeOH_sto * P_cap[7] * data.raw_MeOH_sto[3]
            + sum(
                Xs[k].raw_MeOH_min * P_cap[7] * data.raw_MeOH_sto[4]
                - Xs[k].raw_MeOH_max * P_cap[7] * data.raw_MeOH_sto[5]
                - Xs[k].C_raw_MeOH_1 * P_cap[7] * data.raw_MeOH_sto[6]
                - Xs[k].C_raw_MeOH_2 * P_cap[7] * data.raw_MeOH_sto[6]
            )                                                           # Raw Storage
            + sum(
                - Xs[k].dest_max * P_cap[4]
                + Xs[k].dest_min * P_cap[4] * data.op_dest[2]
            )
            + sum(
                - Xs[k].dest_ramp_up * P_cap[4] * data.op_dest[0]
                - Xs[k].dest_ramp_down * P_cap[4] * data.op_dest[1]
            )                                                           # Distillator
            + Xs[k].init_pure * P_cap[9] * data.Pure_MeOH_sto[2]
            + Xs[k].fin_pure * P_cap[9] * data.Pure_MeOH_sto[3]
            + sum(
                Xs[k].MeOH_min * P_cap[9] * data.Pure_MeOH_sto[4]
                - Xs[k].MeOH_max * P_cap[9] * data.Pure_MeOH_sto[5]
                - Xs[k].C_rate_MeOH_1 * P_cap[9] * data.Pure_MeOH_sto[6]
                - Xs[k].C_rate_MeOH_2 * P_cap[9] * data.Pure_MeOH_sto[6]
            )                                                           # Pure Storage
        )

        objective.append(dual_objective*lambdas[k])

    master_DW.setObjective(gp.quicksum(objective), GRB.MAXIMIZE)

    master_DW.optimize()
    
    master_DW.update()
    master_DW.write('DW_master_1000.lp')

    my_pi = DW_duals(
                e_bat_SOC=e_bat_SOC.Pi,
                e_elec=e_elec.Pi,
                m_h2_sto_SOC=m_h2_sto_SOC.Pi,
                m_CO2_sto_SOC=m_CO2_sto_SOC.Pi,
                m_rawMeOH_reac=m_rawMeOH_reac.Pi,
                m_rawMeOH_sto_SOC=m_rawMeOH_sto_SOC.Pi,
                m_MeOH_dest=m_MeOH_dest.Pi,
                m_MeOH_SOC=m_MeOH_SOC.Pi,
            )
    my_kappa = lambdas_total.Pi
    for k in range(len(Xs)):
        print("---------------------------------------lambda",lambdas[k].x)
    return master_DW.objVal, my_pi, my_kappa