import gurobipy as gp
from gurobipy import GRB
import math
from DW_X1_results import X1_Results
from sub_results import *

def DW_sub(data,P_cap,my_pi,my_kappa):

#   Sub model NOTE (setupSUB)
    sub_DW = gp.Model("DW_sub")

# VARIABLES
    # Power flow
    p_flow = sub_DW.addMVar(data.hours, lb=-10000, ub=10000, name="p_flow")

    # Max grid
    max_grid = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="max_grid")

    # Electrolyzer operation
    elec_max = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="Elec_max")
    elec_min = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="Elec_min")
    elec_op = sub_DW.addMVar(data.hours, lb=-10000, ub=10000, name="Elec_op")
    elec_H2O = sub_DW.addMVar(data.hours, lb=-10000, ub=10000, name="elec_H2O")
    elec_ramp_up = sub_DW.addMVar(data.hours-1, lb=0, ub=10000, name="elec_ramp_up")
    elec_ramp_down = sub_DW.addMVar(data.hours-1, lb=0, ub=10000, name="elec_ramp_down")

    # Battery
    init_bat = sub_DW.addMVar(1, lb=-10000, ub=10000, name="init_bat")
    fin_bat = sub_DW.addMVar(1, lb=0, ub=10000, name="fin_bat")
    bat_sto = sub_DW.addMVar(data.hours-1, lb=-10000, ub=10000, name="bat_sto")
    bat_min = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="bat_min")
    bat_max = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="bat_max")
    C_bat_1 = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="C_bat_1")
    C_bat_2 = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="C_bat_2")

    # Desalinator
    des_H2O = sub_DW.addMVar(data.hours, lb=-10000, ub=10000, name="des_H2O")

    # Hydrogen storage
    init_h2 = sub_DW.addMVar(1, lb=-10000, ub=10000, name="initialize_h2")
    fin_h2 = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_h2")
    h2_storage = sub_DW.addMVar(data.hours-1, lb=-10000, ub=10000, name="h2_storage")
    h2_min = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="h2_min")
    h2_max = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="h2_max")
    C_rate_h2_1 = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="C_rate_h2_1")
    C_rate_h2_2 = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="C_rate_h2_2")

    # H2 Compressor
    compressor_h2 = sub_DW.addMVar(data.hours, lb=-10000, ub=10000, name="compressor_h2")

    # CO2 Storage
    init_CO2 = sub_DW.addMVar(1, lb=-10000, ub=10000, name="initialize_CO2")
    fin_CO2 = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_CO2")
    CO2_sto = sub_DW.addMVar(data.hours-1, lb=-10000, ub=10000, name="CO2_storage")
    CO2_min = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="CO2_min")
    CO2_max = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="CO2_max")
    CO2_str_min = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="CO2_stream_min")
    CO2_str_max = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="CO2_stream_max")
    CO2_comp = sub_DW.addMVar(data.hours, lb=-10000, ub=10000, name="CO2_comp")

    # Feed syngas Compressor
    compressor_e = sub_DW.addMVar(data.hours, lb=-10000, ub=10000, name="compressor_e")

    # Methanol reactor operation
    reac_max = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="Reac_max")
    reac_min = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="Reac_min")
    reac_H2 =  sub_DW.addMVar(data.hours, lb=-10000, ub=10000, name="Reac_H2")
    reac_CO2 = sub_DW.addMVar(data.hours, lb=-10000, ub=10000, name="reac_CO2")
    reac_el = sub_DW.addMVar(data.hours, lb=-10000, ub=10000, name="Reac_el")
    reac_heat = sub_DW.addMVar(data.hours, lb=-10000, ub=10000, name="Reac_heat") 
    Reac_ramp_up = sub_DW.addMVar(data.hours-1, lb=0, ub=10000, name="Reac_ramp_up")
    Reac_ramp_down = sub_DW.addMVar(data.hours-1, lb=0, ub=10000, name="Reac_ramp_down")

    # Methanol compressor
    compr_MeOH = sub_DW.addMVar(data.hours, lb=-10000, ub=10000, name="compr_MeOH")

    # Raw Methanol storage
    init_MeOH_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="initialize_MeOH_sto")
    fin_MeOH_sto = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_MeOH_sto")
    rawMeOH_sto = sub_DW.addMVar(data.hours-1, lb=-10000, ub=10000, name="rawMeOH_storage")
    raw_MeOH_min = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="raw_MeOH_min")
    raw_MeOH_max = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="raw_MeOH_max")
    C_raw_MeOH_1 = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="C_rate_raw_MeOH_1")
    C_raw_MeOH_2 = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="C_rate_raw_MeOH_2")

    # Methanol distillator
    dest_max = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="Dest_max")
    dest_min = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="Dest_min")
    dest_rawMeOH = sub_DW.addMVar(data.hours, lb=-10000, ub=10000, name="Dest_rawMeOH")
    dest_water = sub_DW.addMVar(data.hours, lb=-10000, ub=10000, name="Dest_water")
    dest_heat = sub_DW.addMVar(data.hours, lb=-10000, ub=10000, name="Dest_heat")
    dest_sto = sub_DW.addMVar(data.hours, lb=-10000, ub=10000, name="Dest_sto")
    dest_ramp_up = sub_DW.addMVar(data.hours-1, lb=0, ub=10000, name="Dest_ramp_up")
    dest_ramp_down = sub_DW.addMVar(data.hours-1, lb=0, ub=10000, name="Dest_ramp_down")

    # Pure Methanol storage
    init_pure = sub_DW.addMVar(1, lb=-10000, ub=10000, name="initialize_pure")
    fin_pure = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_pure")
    pure_sto = sub_DW.addMVar(data.hours-1, lb=-10000, ub=10000, name="pure_storage")
    MeOH_min = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="MeOH_min")
    MeOH_max = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="MeOH_max")
    C_rate_MeOH_1 = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="C_rate_MeOH_1")
    C_rate_MeOH_2 = sub_DW.addMVar(data.hours, lb=0, ub=10000, name="C_rate_MeOH_2")

    # Meet demand
    balance = sub_DW.addMVar(data.hours, lb=0, name="balance")

#   SUB CONSTRAINTS
    # Energy
    e_buy = sub_DW.addConstr(-p_flow - max_grid <= data.Prices, name='e_buy')
    e_sell = sub_DW.addConstr(p_flow <= 0, name='e_sell')

    # Battery
    e_bat_cha = sub_DW.addConstr(- bat_sto*data.Battery[0] + p_flow[:-1] - C_bat_1[:-1] <= 0)
    e_bat_dis = sub_DW.addConstr(+ bat_sto/data.Battery[1] - p_flow[:-1] - C_bat_2[:-1] <= 0)
    e_bat_cha_h = sub_DW.addConstr(+ p_flow[data.hours-1] - C_bat_1[data.hours-1] <= 0)
    e_bat_dis_h = sub_DW.addConstr(- p_flow[data.hours-1] - C_bat_2[data.hours-1] <= 0)

    # Desalinator
    e_des = sub_DW.addConstr(des_H2O + p_flow <= 0, name='e_des')

    # Electrolyzer
    m_h2_elec = sub_DW.addConstr(elec_op*data.Electrolyzer[0] + elec_H2O*data.Electrolyzer[1] + compressor_h2 <= 0, name='m_h2_elec')
    m_H2O_elec = sub_DW.addConstr(- elec_H2O - des_H2O*data.Desalinator[0] <= 0, name='m_H2O_elec')

    # H2 Storage
    m_h2_sto_cha = sub_DW.addConstr(- h2_storage*data.H2storage[0] - C_rate_h2_1[:-1] - compressor_h2[:-1] <= 0, name='m_h2_sto_cha')
    m_h2_sto_dis = sub_DW.addConstr(h2_storage/data.H2storage[1] - C_rate_h2_2[:-1] - reac_H2[:-1] - compressor_e[:-1]*data.Compressor[0] <= 0, name='m_h2_sto_dis')
    m_h2_sto_cha_h = sub_DW.addConstr(- C_rate_h2_1[data.hours-1] - compressor_h2[data.hours-1] <= 0, name='m_h2_sto_cha_h')
    m_h2_sto_dis_h = sub_DW.addConstr(- C_rate_h2_2[data.hours-1] - reac_H2[data.hours-1] - compressor_e[data.hours-1]*data.Compressor[0] <= 0, name='m_h2_sto_dis_h')

    # CO2 Storage
    m_CO2_buy = sub_DW.addConstr(CO2_str_min - CO2_str_max + CO2_comp <= data.CO2_prices, name='m_CO2_buy')
    m_CO2_sto_cha = sub_DW.addConstr(- CO2_sto*data.CO2storage[0] - CO2_comp[:-1] <= 0, name='m_CO2_sto_cha')
    m_CO2_sto_dis = sub_DW.addConstr(CO2_sto/data.CO2storage[1] - compressor_e[:-1]*data.Compressor[0] - reac_CO2[:-1] <= 0, name='m_CO2_sto_dis')
    m_CO2_sto_cha_h = sub_DW.addConstr(- CO2_comp[data.hours-1] <= 0, name='m_CO2_sto_cha_h')
    m_CO2_sto_dis_h = sub_DW.addConstr(- compressor_e[data.hours-1]*data.Compressor[0] - reac_CO2[data.hours-1] <= 0, name='m_CO2_sto_dis_h')

    # Syngas compressor
    e_compr = sub_DW.addConstr(p_flow + compressor_e <= 0, name='e_compr')

    # Methanol reactor
    m_h2_reac = sub_DW.addConstr(- compressor_h2 - reac_H2 - compressor_e*data.Compressor[0] <= 0, name='m_h2_reac')
    m_CO2_reac = sub_DW.addConstr(- CO2_comp - reac_CO2 - compressor_e*data.Compressor[0] <= 0, name='m_CO2_reac') 
    e_reac = sub_DW.addConstr(p_flow - reac_el <= 0, name='e_reac')
    q_reac = sub_DW.addConstr(- reac_heat + p_flow*data.Heater[0] <= 0, name='q_reac')

    # MeOH Storage
    m_rawMeOH_sto_cha = sub_DW.addConstr(-rawMeOH_sto*data.raw_MeOH_sto[0] - C_raw_MeOH_1[:-1] - compr_MeOH[:-1] <= 0, name='m_rawMeOH_sto_cha')
    m_rawMeOH_sto_dis = sub_DW.addConstr(+rawMeOH_sto/data.raw_MeOH_sto[1] - C_raw_MeOH_2[:-1] - dest_rawMeOH[:-1] <= 0, name='m_rawMeOH_sto_dis')
    m_rawMeOH_sto_cha_h = sub_DW.addConstr( - C_raw_MeOH_1[data.hours-1] - compr_MeOH[data.hours-1] <= 0, name='m_rawMeOH_sto_cha')
    m_rawMeOH_sto_dis_h = sub_DW.addConstr( - C_raw_MeOH_2[data.hours-1] - dest_rawMeOH[data.hours-1] <= 0, name='m_rawMeOH_sto_dis')

    # Methanol distilator
    m_rawMeOH_dest = sub_DW.addConstr(- compr_MeOH - dest_rawMeOH <= 0, name='m_rawMeOH_dest')
    q_dest = sub_DW.addConstr(- dest_heat + p_flow*data.Heater[0] <= 0, name='q_dest')
    m_h2O_dest = sub_DW.addConstr(- dest_water <= 0, name='m_h2O_dest')

    # Pure MeOH Storagess
    m_MeOH_cha = sub_DW.addConstr(- dest_sto[:-1] - C_rate_MeOH_1[:-1] - pure_sto*data.Pure_MeOH_sto[0] <= 0, name='m_MeOH_cha')
    m_MeOH_dis = sub_DW.addConstr(pure_sto/data.Pure_MeOH_sto[1] - C_rate_MeOH_2[:-1] + balance[:-1] <= 0, name="m_MeOH_dis") 
    m_MeOH_cha_h = sub_DW.addConstr(- dest_sto[data.hours-1] - C_rate_MeOH_1[data.hours-1] <= 0, name='m_MeOH_cha_h')
    m_MeOH_dis_h = sub_DW.addConstr(- C_rate_MeOH_2[data.hours-1] + balance[data.hours-1] <= 0, name="m_MeOH_dis_h")


    # Sub objective

    c_x = (
    p_flow @ (data.Wind * P_cap[0] + data.Solar * P_cap[1])
    - gp.quicksum(max_grid * data.grid_max)                         # Grid
    + init_bat * P_cap[5] * data.Battery[2] + fin_bat * P_cap[5] * data.Battery[3]
    + gp.quicksum(
        bat_min * P_cap[5] * data.Battery[4]
        - bat_max * P_cap[5] * data.Battery[5]
        - C_bat_1 * P_cap[5] * data.Battery[6]
        - C_bat_2 * P_cap[5] * data.Battery[6]
    )                                                           # Battery
    + gp.quicksum(
        -elec_max * P_cap[2]
        + elec_min * P_cap[2] * data.op_elec[2]
    )
    + gp.quicksum(
        - elec_ramp_up * P_cap[2] * data.op_elec[0]
        - elec_ramp_down * P_cap[2] * data.op_elec[1]
    )                                                           # Electrolyzer
    + balance @ data.Demand                                     # Demand
    + init_h2 * P_cap[6] * data.H2storage[2]
    + fin_h2 * P_cap[6] * data.H2storage[3]
    + gp.quicksum(
        h2_min * P_cap[6] * data.H2storage[4]
        - h2_max * P_cap[6] * data.H2storage[5]
        - C_rate_h2_1 * P_cap[6] * data.H2storage[6]
        - C_rate_h2_2 * P_cap[6] * data.H2storage[6]
    )                                                           # H2 Storage
    + init_CO2 * P_cap[8] * data.CO2storage[2]
    + fin_CO2 * P_cap[8] * data.CO2storage[3]
    + gp.quicksum(
        CO2_min * P_cap[8] * data.CO2storage[4]
        - CO2_max * P_cap[8] * data.CO2storage[5]
        + CO2_str_min * data.CO2storage[6]
        - CO2_str_max * data.CO2storage[7]
    )                                                           # CO2 Storage
    + gp.quicksum(
        -reac_max * P_cap[3]
        + reac_min * P_cap[3] * data.op_reac[2]
    )
    + gp.quicksum(
        - Reac_ramp_up * P_cap[3] * data.op_reac[0]
        - Reac_ramp_down * P_cap[3] * data.op_reac[1]
    )                                                           # Reactor
    + init_MeOH_sto * P_cap[7] * data.raw_MeOH_sto[2]
    + fin_MeOH_sto * P_cap[7] * data.raw_MeOH_sto[3]
    + gp.quicksum(
        raw_MeOH_min * P_cap[7] * data.raw_MeOH_sto[4]
        - raw_MeOH_max * P_cap[7] * data.raw_MeOH_sto[5]
        - C_raw_MeOH_1 * P_cap[7] * data.raw_MeOH_sto[6]
        - C_raw_MeOH_2 * P_cap[7] * data.raw_MeOH_sto[6]
    )                                                           # Raw Storage
    + gp.quicksum(
        - dest_max * P_cap[4]
        + dest_min * P_cap[4] * data.op_dest[2]
    )
    + gp.quicksum(
        - dest_ramp_up * P_cap[4] * data.op_dest[0]
        - dest_ramp_down * P_cap[4] * data.op_dest[1]
    )                                                           # Distillator
    + init_pure * P_cap[9] * data.Pure_MeOH_sto[2]
    + fin_pure * P_cap[9] * data.Pure_MeOH_sto[3]
    + gp.quicksum(
        MeOH_min * P_cap[9] * data.Pure_MeOH_sto[4]
        - MeOH_max * P_cap[9] * data.Pure_MeOH_sto[5]
        - C_rate_MeOH_1 * P_cap[9] * data.Pure_MeOH_sto[6]
        - C_rate_MeOH_2 * P_cap[9] * data.Pure_MeOH_sto[6]
    )                                                           # Pure Storage
    )
        
    pi_A0_x = (
    # Battery SOC terms

    my_pi.e_bat_SOC_0 * (init_bat - bat_sto[0] + bat_min[0] - bat_max[0]) +
    gp.quicksum(my_pi.e_bat_SOC[t] * (bat_sto[t] - bat_sto[t + 1] + bat_min[t + 1] - bat_max[t + 1])
            for t in range(data.hours - 2)) +
    my_pi.e_bat_SOC_h * (fin_bat + bat_sto[data.hours - 2] + bat_min[data.hours - 1] - bat_max[data.hours - 1]) +

    # Electricity flow terms
    my_pi.e_elec_0 * (p_flow[0] - elec_max[0] + elec_min[0] - elec_op[0] + elec_ramp_up[0] - elec_ramp_down[0]) +
    gp.quicksum(my_pi.e_elec[t] * (p_flow[t + 1] - elec_max[t + 1] + elec_min[t + 1] - elec_op[t + 1] -
                            elec_ramp_up[t] + elec_ramp_down[t] + elec_ramp_up[t + 1] - elec_ramp_down[t + 1])
            for t in range(data.hours - 2)) +
    my_pi.e_elec_h * (p_flow[data.hours - 1] - elec_max[data.hours - 1] + elec_min[data.hours - 1] -
                elec_op[data.hours - 1] - elec_ramp_up[data.hours - 2] + elec_ramp_down[data.hours - 2]) +

    # Hydrogen storage SOC terms
    my_pi.m_h2_sto_SOC_0 * (init_h2 - h2_storage[0] + h2_min[0] - h2_max[0]) +
    gp.quicksum(my_pi.m_h2_sto_SOC[t] * (h2_storage[t] - h2_storage[t + 1] + h2_min[t + 1] - h2_max[t + 1])
            for t in range(data.hours - 2)) +
    my_pi.m_h2_sto_SOC_h * (fin_h2 + h2_storage[data.hours - 2] + h2_min[data.hours - 1] - h2_max[data.hours - 1]) +

    # CO2 storage SOC terms
    my_pi.m_CO2_sto_SOC_0 * (init_CO2 - CO2_sto[0] + CO2_min[0] - CO2_max[0]) +
    gp.quicksum(my_pi.m_CO2_sto_SOC[t] * (CO2_sto[t] - CO2_sto[t + 1] + CO2_min[t + 1] - CO2_max[t + 1])
            for t in range(data.hours - 2)) +
    my_pi.m_CO2_sto_SOC_h * (fin_CO2 + CO2_sto[data.hours - 2] + CO2_min[data.hours - 1] - CO2_max[data.hours - 1]) +

    # Methanol reactor terms
    my_pi.m_rawMeOH_reac_0 * (-reac_max[0] + reac_min[0] + reac_H2[0] * data.Reactor[0] +
                        reac_CO2[0] * data.Reactor[1] + reac_el[0] * data.Reactor[2] +
                        reac_heat[0] * data.Reactor[3] + Reac_ramp_up[0] - Reac_ramp_down[0] +
                        compr_MeOH[0]) +
    gp.quicksum(my_pi.m_rawMeOH_reac[t] * (-reac_max[t + 1] + reac_min[t + 1] + reac_H2[t + 1] * data.Reactor[0] +
                                    reac_CO2[t + 1] * data.Reactor[1] + reac_el[t + 1] * data.Reactor[2] +
                                    reac_heat[t + 1] * data.Reactor[3] - Reac_ramp_up[t] + Reac_ramp_down[t] +
                                    Reac_ramp_up[t + 1] - Reac_ramp_down[t + 1] + compr_MeOH[t + 1])
            for t in range(data.hours - 2)) +
    my_pi.m_rawMeOH_reac_h * (-reac_max[data.hours - 1] + reac_min[data.hours - 1] +
                        reac_H2[data.hours - 1] * data.Reactor[0] +
                        reac_CO2[data.hours - 1] * data.Reactor[1] +
                        reac_el[data.hours - 1] * data.Reactor[2] +
                        reac_heat[data.hours - 1] * data.Reactor[3] - Reac_ramp_up[data.hours - 2] +
                        Reac_ramp_down[data.hours - 2] + compr_MeOH[data.hours - 1]) +

    # Raw methanol storage SOC terms
    my_pi.m_rawMeOH_sto_SOC_0 * (init_MeOH_sto - rawMeOH_sto[0] + raw_MeOH_min[0] - raw_MeOH_max[0]) +
    gp.quicksum(my_pi.m_rawMeOH_sto_SOC[t] * (rawMeOH_sto[t + 1] - rawMeOH_sto[t] + raw_MeOH_min[t + 1] - raw_MeOH_max[t + 1])
            for t in range(data.hours - 2)) +
    my_pi.m_rawMeOH_sto_SOC_h * (fin_MeOH_sto + rawMeOH_sto[data.hours - 2] + raw_MeOH_min[data.hours - 1] -
                            raw_MeOH_max[data.hours - 1]) +

    # Methanol distillation terms
    my_pi.m_MeOH_dest_0 * (-dest_max[0] + dest_min[0] + dest_rawMeOH[0] * data.Destilator[0] +
                    dest_water[0] * data.Destilator[1] + dest_heat[0] * data.Destilator[3] +
                    dest_sto[0] + dest_ramp_up[0] - dest_ramp_down[0]) +
    gp.quicksum(my_pi.m_MeOH_dest[t] * (-dest_max[t + 1] + dest_min[t + 1] +
                                dest_rawMeOH[t + 1] * data.Destilator[0] +
                                dest_water[t + 1] * data.Destilator[1] +
                                dest_heat[t + 1] * data.Destilator[3] + dest_sto[t + 1] -
                                dest_ramp_up[t] + dest_ramp_down[t] + dest_ramp_up[t + 1] - dest_ramp_down[t + 1])
            for t in range(data.hours - 2)) +
    my_pi.m_MeOH_dest_h * (-dest_max[data.hours - 1] + dest_min[data.hours - 1] +
                    dest_rawMeOH[data.hours - 1] * data.Destilator[0] +
                    dest_water[data.hours - 1] * data.Destilator[1] +
                    dest_heat[data.hours - 1] * data.Destilator[3] + dest_sto[data.hours - 1] -
                    dest_ramp_up[data.hours - 2] + dest_ramp_down[data.hours - 2]) +

    # Pure methanol storage SOC terms
    my_pi.m_MeOH_SOC_0 * (init_pure - pure_sto[0] + MeOH_min[0] - MeOH_max[0]) +
    gp.quicksum(my_pi.m_MeOH_SOC[t] * (pure_sto[t] - pure_sto[t + 1] + MeOH_min[t + 1] - MeOH_max[t + 1])
            for t in range(data.hours - 2)) +
    my_pi.m_MeOH_SOC_h * (fin_pure + pure_sto[data.hours - 2] + MeOH_min[data.hours - 1] - MeOH_max[data.hours - 1])
    )

# Set the objective in the model
    sub_DW.setObjective(c_x - pi_A0_x - my_kappa, GRB.MAXIMIZE)

    sub_DW.update()

    sub_DW.write('wolfe.lp')

    sub_DW.optimize()


    if sub_DW.Status == GRB.OPTIMAL:

        # Print optimal objective value
        sub_val = sub_DW.objVal

        X1 = X1_Results(
            p_flow=p_flow.X,
            max_grid=max_grid.X,

            elec_max=elec_max.X,
            elec_min=elec_min.X,
            elec_op=elec_op.X,
            elec_H2O=elec_H2O.X,
            elec_ramp_up=elec_ramp_up.X,
            elec_ramp_down=elec_ramp_down.X,

            init_bat=init_bat.X,
            fin_bat=fin_bat.X,
            bat_sto=bat_sto.X,
            bat_min=bat_min.X,
            bat_max=bat_max.X,
            C_bat_1=C_bat_1.X,
            C_bat_2=C_bat_2.X,
            
            des_H2O=des_H2O.X,

            init_h2=init_h2.X,
            fin_h2=fin_h2.X,
            h2_storage=h2_storage.X,
            h2_min=h2_min.X,
            h2_max=h2_max.X,
            C_rate_h2_1=C_rate_h2_1.X,
            C_rate_h2_2=C_rate_h2_2.X,

            compressor_h2=compressor_h2.X,

            init_CO2=init_CO2.X,
            fin_CO2=fin_CO2.X,
            CO2_sto=CO2_sto.X,
            CO2_min=CO2_min.X,
            CO2_max=CO2_max.X,
            CO2_str_min=CO2_str_min.X,
            CO2_str_max=CO2_str_max.X,
            CO2_comp=CO2_comp.X,

            compressor_e=compressor_e.X,

            reac_max=reac_max.X,
            reac_min=reac_min.X,
            reac_H2=reac_H2.X,
            reac_CO2=reac_CO2.X,
            reac_el=reac_el.X,
            reac_heat=reac_heat.X,
            Reac_ramp_up=Reac_ramp_up.X,
            Reac_ramp_down=Reac_ramp_down.X,

            compr_MeOH=compr_MeOH.X,

            init_MeOH_sto=init_MeOH_sto.X,
            fin_MeOH_sto=fin_MeOH_sto.X,
            rawMeOH_sto=rawMeOH_sto.X,
            raw_MeOH_min=raw_MeOH_min.X,
            raw_MeOH_max=raw_MeOH_max.X,
            C_raw_MeOH_1=C_raw_MeOH_1.X,
            C_raw_MeOH_2=C_raw_MeOH_2.X,

            dest_max=dest_max.X,
            dest_min=dest_min.X,
            dest_rawMeOH=dest_rawMeOH.X,
            dest_water=dest_water.X,
            dest_heat=dest_heat.X,
            dest_sto=dest_sto.X,
            dest_ramp_up=dest_ramp_up.X,
            dest_ramp_down=dest_ramp_down.X,

            init_pure=init_pure.X,
            fin_pure=fin_pure.X,
            pure_sto=pure_sto.X,
            MeOH_min=MeOH_min.X,
            MeOH_max=MeOH_max.X,
            C_rate_MeOH_1=C_rate_MeOH_1.X,
            C_rate_MeOH_2=C_rate_MeOH_2.X,

            balance=balance.X,
        )

    return sub_val, X1