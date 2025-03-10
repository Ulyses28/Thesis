import gurobipy as gp
from gurobipy import GRB
import math
from DW_X1_results import *
from DW_duals import *

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
    e_bat_SOC_0 = master_DW.addConstr(gp.quicksum((Xs[k].init_bat - Xs[k].bat_sto[0] + Xs[k].bat_min[0] - Xs[k].bat_max[0]) * lambdas[k] for k in range(len(Xs))) <= 0)
    e_bat_SOC = master_DW.addConstr(gp.quicksum((Xs[k].bat_sto[:-1] - Xs[k].bat_sto[1:] + Xs[k].bat_min[1:-1] - Xs[k].bat_max[1:-1]) * lambdas[k] for k in range(len(Xs))) <= 0)
    e_bat_SOC_h = master_DW.addConstr(gp.quicksum((Xs[k].fin_bat + Xs[k].bat_sto[data.hours-2] + Xs[k].bat_min[data.hours-1] - Xs[k].bat_max[data.hours-1]) * lambdas[k] for k in range(len(Xs))) <= 0)

    e_elec_0 = master_DW.addConstr(gp.quicksum((Xs[k].p_flow[0] - Xs[k].elec_max[0] + Xs[k].elec_min[0] - Xs[k].elec_op[0] + Xs[k].elec_ramp_up[0] - Xs[k].elec_ramp_down[0]) * lambdas[k] for k in range(len(Xs))) <= 0, name='e_elec_0')
    e_elec = master_DW.addConstr(gp.quicksum((Xs[k].p_flow[1:-1] - Xs[k].elec_max[1:-1] + Xs[k].elec_min[1:-1] - Xs[k].elec_op[1:-1] - Xs[k].elec_ramp_up[:-1] + Xs[k].elec_ramp_down[:-1] + Xs[k].elec_ramp_up[1:] - Xs[k].elec_ramp_down[1:]) * lambdas[k] for k in range(len(Xs))) <= 0, name='e_elec')
    e_elec_h = master_DW.addConstr(gp.quicksum((Xs[k].p_flow[data.hours-1] - Xs[k].elec_max[data.hours-1] + Xs[k].elec_min[data.hours-1] - Xs[k].elec_op[data.hours-1] - Xs[k].elec_ramp_up[data.hours-2] + Xs[k].elec_ramp_down[data.hours-2]) * lambdas[k] for k in range(len(Xs))) <= 0, name='e_elec_h')

    m_h2_sto_SOC_0 = master_DW.addConstr(gp.quicksum((Xs[k].init_h2 - Xs[k].h2_storage[0] + Xs[k].h2_min[0] - Xs[k].h2_max[0]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_h2_sto_SOC_0')
    m_h2_sto_SOC = master_DW.addConstr(gp.quicksum((Xs[k].h2_storage[:-1] - Xs[k].h2_storage[1:] + Xs[k].h2_min[1:-1] - Xs[k].h2_max[1:-1]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_h2_sto_SOC')
    m_h2_sto_SOC_h = master_DW.addConstr(gp.quicksum((Xs[k].fin_h2 + Xs[k].h2_storage[data.hours-2] + Xs[k].h2_min[data.hours-1] - Xs[k].h2_max[data.hours-1]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_h2_sto_SOC_h')

    m_CO2_sto_SOC_0 = master_DW.addConstr(gp.quicksum((Xs[k].init_CO2 - Xs[k].CO2_sto[0] + Xs[k].CO2_min[0] - Xs[k].CO2_max[0]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_CO2_sto_SOC_0')
    m_CO2_sto_SOC = master_DW.addConstr(gp.quicksum((Xs[k].CO2_sto[:-1] - Xs[k].CO2_sto[1:] + Xs[k].CO2_min[1:-1] - Xs[k].CO2_max[1:-1]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_CO2_sto_SOC')
    m_CO2_sto_SOC_h = master_DW.addConstr(gp.quicksum((Xs[k].fin_CO2 + Xs[k].CO2_sto[data.hours-2] + Xs[k].CO2_min[data.hours-1] - Xs[k].CO2_max[data.hours-1]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_CO2_sto_SOC_h') 

    m_rawMeOH_reac_0 = master_DW.addConstr(gp.quicksum((-Xs[k].reac_max[0] + Xs[k].reac_min[0] + Xs[k].reac_H2[0]*data.Reactor[0] + Xs[k].reac_CO2[0]*data.Reactor[1] + Xs[k].reac_el[0]*data.Reactor[2] + Xs[k].reac_heat[0]*data.Reactor[3] + Xs[k].Reac_ramp_up[0] - Xs[k].Reac_ramp_down[0] + Xs[k].compr_MeOH[0]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_rawMeOH_reac_0')
    m_rawMeOH_reac = master_DW.addConstr(gp.quicksum((-Xs[k].reac_max[1:-1] + Xs[k].reac_min[1:-1] + Xs[k].reac_H2[1:-1]*data.Reactor[0] + Xs[k].reac_CO2[1:-1]*data.Reactor[1] + Xs[k].reac_el[1:-1]*data.Reactor[2] + Xs[k].reac_heat[1:-1]*data.Reactor[3] - Xs[k].Reac_ramp_up[:-1] + Xs[k].Reac_ramp_down[:-1] + Xs[k].Reac_ramp_up[1:] - Xs[k].Reac_ramp_down[1:] + Xs[k].compr_MeOH[1:-1]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_rawMeOH_reac')
    m_rawMeOH_reac_h = master_DW.addConstr(gp.quicksum((-Xs[k].reac_max[data.hours-1] + Xs[k].reac_min[data.hours-1] + Xs[k].reac_H2[data.hours-1]*data.Reactor[0] + Xs[k].reac_CO2[data.hours-1]*data.Reactor[1] + Xs[k].reac_el[data.hours-1]*data.Reactor[2] + Xs[k].reac_heat[data.hours-1]*data.Reactor[3] - Xs[k].Reac_ramp_up[data.hours-2] + Xs[k].Reac_ramp_down[data.hours-2] + Xs[k].compr_MeOH[data.hours-1]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_rawMeOH_reac_h')

    m_rawMeOH_sto_SOC_0 = master_DW.addConstr(gp.quicksum((Xs[k].init_MeOH_sto - Xs[k].rawMeOH_sto[0] + Xs[k].raw_MeOH_min[0] - Xs[k].raw_MeOH_max[0]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_rawMeOH_sto_SOC_0')
    m_rawMeOH_sto_SOC = master_DW.addConstr(gp.quicksum((-Xs[k].rawMeOH_sto[1:] + Xs[k].rawMeOH_sto[:-1] + Xs[k].raw_MeOH_min[1:-1] - Xs[k].raw_MeOH_max[1:-1]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_rawMeOH_sto_SOC')
    m_rawMeOH_sto_SOC_h = master_DW.addConstr(gp.quicksum((Xs[k].fin_MeOH_sto + Xs[k].rawMeOH_sto[data.hours-2] + Xs[k].raw_MeOH_min[data.hours-1] - Xs[k].raw_MeOH_max[data.hours-1]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_rawMeOH_sto_SOC_h')

    m_MeOH_dest_0 = master_DW.addConstr(gp.quicksum((- Xs[k].dest_max[0] + Xs[k].dest_min[0] + Xs[k].dest_rawMeOH[0]*data.Destilator[0] + Xs[k].dest_water[0]*data.Destilator[1] + Xs[k].dest_heat[0]*data.Destilator[3] + Xs[k].dest_sto[0] + Xs[k].dest_ramp_up[0] - Xs[k].dest_ramp_down[0]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_MeOH_dest_0')
    m_MeOH_dest = master_DW.addConstr(gp.quicksum((- Xs[k].dest_max[1:-1] + Xs[k].dest_min[1:-1] + Xs[k].dest_rawMeOH[1:-1]*data.Destilator[0] + Xs[k].dest_water[1:-1]*data.Destilator[1] + Xs[k].dest_heat[1:-1]*data.Destilator[3] + Xs[k].dest_sto[1:-1] - Xs[k].dest_ramp_up[:-1] + Xs[k].dest_ramp_down[:-1] + Xs[k].dest_ramp_up[1:] - Xs[k].dest_ramp_down[1:]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_MeOH_dest')
    m_MeOH_dest_h = master_DW.addConstr(gp.quicksum((- Xs[k].dest_max[data.hours-1] + Xs[k].dest_min[data.hours-1] + Xs[k].dest_rawMeOH[data.hours-1]*data.Destilator[0] + Xs[k].dest_water[data.hours-1]*data.Destilator[1] + Xs[k].dest_heat[data.hours-1]*data.Destilator[3] + Xs[k].dest_sto[data.hours-1] - Xs[k].dest_ramp_up[data.hours-2] + Xs[k].dest_ramp_down[data.hours-2]) * lambdas[k] for k in range(len(Xs))) <= 0, name='m_MeOH_dest_h')

    m_MeOH_SOC_0 = master_DW.addConstr(gp.quicksum((Xs[k].init_pure - Xs[k].pure_sto[0] + Xs[k].MeOH_min[0] - Xs[k].MeOH_max[0]) * lambdas[k] for k in range(len(Xs))) <= 0, name="m_MeOH_SOC_0")
    m_MeOH_SOC = master_DW.addConstr(gp.quicksum((Xs[k].pure_sto[:-1] - Xs[k].pure_sto[1:]  + Xs[k].MeOH_min[1:-1] - Xs[k].MeOH_max[1:-1]) * lambdas[k] for k in range(len(Xs))) <= 0, name="m_MeOH_SOC")
    m_MeOH_SOC_h = master_DW.addConstr(gp.quicksum((Xs[k].fin_pure + Xs[k].pure_sto[data.hours-2] + Xs[k].MeOH_min[data.hours-1] - Xs[k].MeOH_max[data.hours-1]) * lambdas[k] for k in range(len(Xs))) <= 0, name="m_MeOH_SOC_h")

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

    my_pi = DW_duals(
                e_bat_SOC_0=e_bat_SOC_0.Pi,
                e_bat_SOC=e_bat_SOC.Pi,
                e_bat_SOC_h=e_bat_SOC_h.Pi,
                e_elec_0=e_elec_0.Pi,
                e_elec=e_elec.Pi,
                e_elec_h=e_elec_h.Pi,
                m_h2_sto_SOC_0=m_h2_sto_SOC_0.Pi,
                m_h2_sto_SOC=m_h2_sto_SOC.Pi,
                m_h2_sto_SOC_h=m_h2_sto_SOC_h.Pi,
                m_CO2_sto_SOC_0=m_CO2_sto_SOC_0.Pi,
                m_CO2_sto_SOC=m_CO2_sto_SOC.Pi,
                m_CO2_sto_SOC_h=m_CO2_sto_SOC_h.Pi,
                m_rawMeOH_reac_0=m_rawMeOH_reac_0.Pi,
                m_rawMeOH_reac=m_rawMeOH_reac.Pi,
                m_rawMeOH_reac_h=m_rawMeOH_reac_h.Pi,
                m_rawMeOH_sto_SOC_0=m_rawMeOH_sto_SOC_0.Pi,
                m_rawMeOH_sto_SOC=m_rawMeOH_sto_SOC.Pi,
                m_rawMeOH_sto_SOC_h=m_rawMeOH_sto_SOC_h.Pi,
                m_MeOH_dest_0=m_MeOH_dest_0.Pi,
                m_MeOH_dest=m_MeOH_dest.Pi,
                m_MeOH_dest_h=m_MeOH_dest_h.Pi,
                m_MeOH_SOC_0=m_MeOH_SOC_0.Pi,
                m_MeOH_SOC=m_MeOH_SOC.Pi,
                m_MeOH_SOC_h=m_MeOH_SOC_h.Pi,
            )
    my_kappa = lambdas_total.Pi
    for k in range(len(Xs)):
        print("---------------------------------------lambda",lambdas[k].x)
    return master_DW.objVal, my_pi, my_kappa