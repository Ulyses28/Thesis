import gurobipy as gp
from gurobipy import GRB
import math
from sub_results import *
import scipy.sparse as sp

def sub(data,P_cap):

    # Create model
    dual_model = gp.Model("Dual_Methanol_Plant")

    # Power flow 
    p_flow = dual_model.addMVar(data.hours, lb=-math.inf, name="p_flow")

    max_grid = dual_model.addMVar(data.hours, lb=0, name="max_grid")

    # Electrolyzer operation
    elec_max = dual_model.addMVar(data.hours, lb=0, name="Elec_max")
    elec_min = dual_model.addMVar(data.hours, lb=0, name="Elec_min")
    elec_op = dual_model.addMVar(data.hours, lb=-math.inf, name="Elec_op")
    elec_H2O = dual_model.addMVar(data.hours, lb=-math.inf, name="elec_H2O")
    elec_ramp_up = dual_model.addMVar(data.hours-1, lb=0, name="elec_ramp_up")
    elec_ramp_down = dual_model.addMVar(data.hours-1, lb=0, name="elec_ramp_down")

    # Battery
    init_bat = dual_model.addMVar(1, lb=-math.inf, name="init_bat")
    fin_bat = dual_model.addMVar(1, lb=0, name="fin_bat")
    bat_sto = dual_model.addMVar(data.hours-1, lb=-math.inf, name="bat_sto")
    bat_min = dual_model.addMVar(data.hours, lb=0, name="bat_min")
    bat_max = dual_model.addMVar(data.hours, lb=0, name="bat_max")
    C_bat_1 = dual_model.addMVar(data.hours, lb=0, name="C_bat_1")
    C_bat_2 = dual_model.addMVar(data.hours, lb=0, name="C_bat_2")

    # Desalinator
    des_H2O = dual_model.addMVar(data.hours, lb=-math.inf, name="des_H2O")

    # Hydrogen storage
    init_h2 = dual_model.addMVar(1, lb=-math.inf, name="initialize_h2")
    fin_h2 = dual_model.addMVar(1, lb=0, name="finalize_h2")
    h2_storage = dual_model.addMVar(data.hours-1, lb=-math.inf, name="h2_storage")
    h2_min = dual_model.addMVar(data.hours, lb=0, name="h2_min")
    h2_max = dual_model.addMVar(data.hours, lb=0, name="h2_max")
    C_rate_h2_1 = dual_model.addMVar(data.hours, lb=0, name="C_rate_h2_1")
    C_rate_h2_2 = dual_model.addMVar(data.hours, lb=0, name="C_rate_h2_2")

    # H2 Compressor
    compressor_h2 = dual_model.addMVar(data.hours, lb=-math.inf, name="compressor_h2")

    # CO2 Storage
    init_CO2 = dual_model.addMVar(1, lb=-math.inf, name="initialize_CO2")
    fin_CO2 = dual_model.addMVar(1, lb=0, name="finalize_CO2")
    CO2_sto = dual_model.addMVar(data.hours-1, lb=-math.inf, name="CO2_storage")
    CO2_min = dual_model.addMVar(data.hours, lb=0, name="CO2_min")
    CO2_max = dual_model.addMVar(data.hours, lb=0, name="CO2_max")
    CO2_str_min = dual_model.addMVar(data.hours, lb=0, name="CO2_stream_min")
    CO2_str_max = dual_model.addMVar(data.hours, lb=0, name="CO2_stream_max")
    CO2_comp = dual_model.addMVar(data.hours, lb=-math.inf, name="CO2_comp")

    # Feed syngas Compressor
    compressor_e = dual_model.addMVar(data.hours, lb=-math.inf, name="compressor_e")

    # Methanol reactor operation
    reac_max = dual_model.addMVar(data.hours, lb=0, name="Reac_max")
    reac_min = dual_model.addMVar(data.hours, lb=0, name="Reac_min")
    reac_H2 =  dual_model.addMVar(data.hours, lb=-math.inf, name="Reac_H2")
    reac_CO2 = dual_model.addMVar(data.hours, lb=-math.inf, name="reac_CO2")
    reac_el = dual_model.addMVar(data.hours, lb=-math.inf, name="Reac_el")
    reac_heat = dual_model.addMVar(data.hours, lb=-math.inf, name="Reac_heat") 
    Reac_ramp_up = dual_model.addMVar(data.hours-1, lb=0, name="Reac_ramp_up")
    Reac_ramp_down = dual_model.addMVar(data.hours-1, lb=0, name="Reac_ramp_down")

    # Methanol compressor
    compr_MeOH = dual_model.addMVar(data.hours, lb=-math.inf, name="compr_MeOH")

    # Raw Methanol storage
    init_MeOH_sto = dual_model.addMVar(1, lb=-math.inf, name="initialize_MeOH_sto")
    fin_MeOH_sto = dual_model.addMVar(1, lb=0, name="finalize_MeOH_sto")
    rawMeOH_sto = dual_model.addMVar(data.hours-1, lb=-math.inf, name="rawMeOH_storage")
    raw_MeOH_min = dual_model.addMVar(data.hours, lb=0, name="raw_MeOH_min")
    raw_MeOH_max = dual_model.addMVar(data.hours, lb=0, name="raw_MeOH_max")
    C_raw_MeOH_1 = dual_model.addMVar(data.hours, lb=0, name="C_rate_raw_MeOH_1")
    C_raw_MeOH_2 = dual_model.addMVar(data.hours, lb=0, name="C_rate_raw_MeOH_2")

    # Methanol distillator
    dest_max = dual_model.addMVar(data.hours, lb=0, name="Dest_max")
    dest_min = dual_model.addMVar(data.hours, lb=0, name="Dest_min")
    dest_rawMeOH = dual_model.addMVar(data.hours, lb=-math.inf, name="Dest_rawMeOH")
    dest_water = dual_model.addMVar(data.hours, lb=-math.inf, name="Dest_water")
    dest_heat = dual_model.addMVar(data.hours, lb=-math.inf, name="Dest_heat")
    dest_sto = dual_model.addMVar(data.hours, lb=-math.inf, name="Dest_sto")
    dest_ramp_up = dual_model.addMVar(data.hours-1, lb=0, name="Dest_ramp_up")
    dest_ramp_down = dual_model.addMVar(data.hours-1, lb=0, name="Dest_ramp_down")

    # Pure Methanol storage
    init_pure = dual_model.addMVar(1, lb=-math.inf, name="initialize_pure")
    fin_pure = dual_model.addMVar(1, lb=0, name="finalize_pure")
    pure_sto = dual_model.addMVar(data.hours-1, lb=-math.inf, name="pure_storage")
    MeOH_min = dual_model.addMVar(data.hours, lb=0, name="MeOH_min")
    MeOH_max = dual_model.addMVar(data.hours, lb=0, name="MeOH_max")
    C_rate_MeOH_1 = dual_model.addMVar(data.hours, lb=0, name="C_rate_MeOH_1")
    C_rate_MeOH_2 = dual_model.addMVar(data.hours, lb=0, name="C_rate_MeOH_2")   

    # Meet demand
    balance = dual_model.addMVar(data.hours, lb=0, name="balance")

# CONSTRAINTS

    # Energy
    e_buy = dual_model.addConstr(-p_flow - max_grid <= data.Prices*15, name='e_buy')
    e_sell = dual_model.addConstr(p_flow <= 0, name='e_sell')

    # Battery
    e_bat_SOC_0 = dual_model.addConstr(init_bat - bat_sto[0] + bat_min[0] - bat_max[0] <= 0)
    e_bat_SOC = dual_model.addConstr(bat_sto[:-1] - bat_sto[1:] + bat_min[1:-1] - bat_max[1:-1] <= 0)
    e_bat_SOC_h = dual_model.addConstr(fin_bat + bat_sto[data.hours-2] + bat_min[data.hours-1] - bat_max[data.hours-1] <= 0)
    e_bat_cha = dual_model.addConstr(- bat_sto*data.Battery[0] + p_flow[:-1] - C_bat_1[:-1] <= 0)
    e_bat_dis = dual_model.addConstr(+ bat_sto/data.Battery[1] - p_flow[:-1] - C_bat_2[:-1] <= 0)
    e_bat_cha_h = dual_model.addConstr(+ p_flow[data.hours-1] - C_bat_1[data.hours-1] <= 0)
    e_bat_dis_h = dual_model.addConstr(- p_flow[data.hours-1] - C_bat_2[data.hours-1] <= 0)

    # Desalinator
    e_des = dual_model.addConstr(des_H2O + p_flow <= 0, name='e_des')

    # Electrolyzer
    e_elec_0 = dual_model.addConstr(p_flow[0] - elec_max[0] + elec_min[0] - elec_op[0] + elec_ramp_up[0] - elec_ramp_down[0] <= 0, name='e_elec_0')
    e_elec = dual_model.addConstr(p_flow[1:-1] - elec_max[1:-1] + elec_min[1:-1] - elec_op[1:-1] - elec_ramp_up[:-1] + elec_ramp_down[:-1] + elec_ramp_up[1:] - elec_ramp_down[1:] <= 0, name='e_elec')
    e_elec_h = dual_model.addConstr(p_flow[data.hours-1] - elec_max[data.hours-1] + elec_min[data.hours-1] - elec_op[data.hours-1] - elec_ramp_up[data.hours-2] + elec_ramp_down[data.hours-2] <= 0, name='e_elec_h')
    m_h2_elec = dual_model.addConstr(elec_op*data.Electrolyzer[0] + elec_H2O*data.Electrolyzer[1] + compressor_h2 <= 0, name='m_h2_elec')
    m_H2O_elec = dual_model.addConstr(- elec_H2O - des_H2O*data.Desalinator[0] <= 0, name='m_H2O_elec')
    
    # H2 Storage
    m_h2_sto_SOC_0 = dual_model.addConstr(init_h2 - h2_storage[0] + h2_min[0] - h2_max[0] <= 0, name='m_h2_sto_SOC_0')
    m_h2_sto_SOC = dual_model.addConstr(h2_storage[:-1] - h2_storage[1:] + h2_min[1:-1] - h2_max[1:-1] <= 0, name='m_h2_sto_SOC')
    m_h2_sto_SOC_h = dual_model.addConstr(fin_h2 + h2_storage[data.hours-2] + h2_min[data.hours-1] - h2_max[data.hours-1] <= 0, name='m_h2_sto_SOC_h')
    m_h2_sto_cha = dual_model.addConstr(- h2_storage*data.H2storage[0] - C_rate_h2_1[:-1] - compressor_h2[:-1] <= 0, name='m_h2_sto_cha')
    m_h2_sto_dis = dual_model.addConstr(h2_storage/data.H2storage[1] - C_rate_h2_2[:-1] - reac_H2[:-1] - compressor_e[:-1]*data.Compressor[0] <= 0, name='m_h2_sto_dis')
    m_h2_sto_cha_h = dual_model.addConstr(- C_rate_h2_1[data.hours-1] - compressor_h2[data.hours-1] <= 0, name='m_h2_sto_cha_h')
    m_h2_sto_dis_h = dual_model.addConstr(- C_rate_h2_2[data.hours-1] - reac_H2[data.hours-1] - compressor_e[data.hours-1]*data.Compressor[0] <= 0, name='m_h2_sto_dis_h')

    # CO2 Storage
    m_CO2_sto_SOC_0 = dual_model.addConstr(init_CO2 - CO2_sto[0] + CO2_min[0] - CO2_max[0] <= 0, name='m_CO2_sto_SOC_0')
    m_CO2_sto_SOC = dual_model.addConstr(CO2_sto[:-1] - CO2_sto[1:] + CO2_min[1:-1] - CO2_max[1:-1] <= 0, name='m_CO2_sto_SOC')
    m_CO2_sto_SOC_h = dual_model.addConstr(fin_CO2 + CO2_sto[data.hours-2] + CO2_min[data.hours-1] - CO2_max[data.hours-1] <= 0, name='m_CO2_sto_SOC_h') 
    m_CO2_buy = dual_model.addConstr(CO2_str_min - CO2_str_max + CO2_comp <= data.CO2_prices*30, name='m_CO2_buy')
    m_CO2_sto_cha = dual_model.addConstr(- CO2_sto*data.CO2storage[0] - CO2_comp[:-1] <= 0, name='m_CO2_sto_cha')
    m_CO2_sto_dis = dual_model.addConstr(CO2_sto/data.CO2storage[1] - compressor_e[:-1]*data.Compressor[0] - reac_CO2[:-1] <= 0, name='m_CO2_sto_dis')
    m_CO2_sto_cha_h = dual_model.addConstr(- CO2_comp[data.hours-1] <= 0, name='m_CO2_sto_cha_h')
    m_CO2_sto_dis_h = dual_model.addConstr(- compressor_e[data.hours-1]*data.Compressor[0] - reac_CO2[data.hours-1] <= 0, name='m_CO2_sto_dis_h')

    # Syngas compressor
    e_compr = dual_model.addConstr(p_flow + compressor_e <= 0, name='e_compr')

    # Methanol reactor
    m_h2_reac = dual_model.addConstr(- compressor_h2 - reac_H2 - compressor_e*data.Compressor[0] <= 0, name='m_h2_reac')
    m_CO2_reac = dual_model.addConstr(- CO2_comp - reac_CO2 - compressor_e*data.Compressor[0] <= 0, name='m_CO2_reac') 
    e_reac = dual_model.addConstr(p_flow - reac_el <= 0, name='e_reac')
    q_reac = dual_model.addConstr(- reac_heat + p_flow*data.Heater[0] <= 0, name='q_reac')
    m_rawMeOH_reac_0 = dual_model.addConstr(-reac_max[0] + reac_min[0] + reac_H2[0]*data.Reactor[0] + reac_CO2[0]*data.Reactor[1] + reac_el[0]*data.Reactor[2] + reac_heat[0]*data.Reactor[3] + Reac_ramp_up[0] - Reac_ramp_down[0] + compr_MeOH[0] <= 0, name='m_rawMeOH_reac_0')
    m_rawMeOH_reac = dual_model.addConstr(-reac_max[1:-1] + reac_min[1:-1] + reac_H2[1:-1]*data.Reactor[0] + reac_CO2[1:-1]*data.Reactor[1] + reac_el[1:-1]*data.Reactor[2] + reac_heat[1:-1]*data.Reactor[3] - Reac_ramp_up[:-1] + Reac_ramp_down[:-1] + Reac_ramp_up[1:] - Reac_ramp_down[1:] + compr_MeOH[1:-1] <= 0, name='m_rawMeOH_reac')
    m_rawMeOH_reac_h = dual_model.addConstr(-reac_max[data.hours-1] + reac_min[data.hours-1] + reac_H2[data.hours-1]*data.Reactor[0] + reac_CO2[data.hours-1]*data.Reactor[1] + reac_el[data.hours-1]*data.Reactor[2] + reac_heat[data.hours-1]*data.Reactor[3] - Reac_ramp_up[data.hours-2] + Reac_ramp_down[data.hours-2] + compr_MeOH[data.hours-1] <= 0, name='m_rawMeOH_reac_h')

    # MeOH Storage
    m_rawMeOH_sto_SOC_0 = dual_model.addConstr(init_MeOH_sto - rawMeOH_sto[0] + raw_MeOH_min[0] - raw_MeOH_max[0] <= 0, name='m_rawMeOH_sto_SOC_0')
    m_rawMeOH_sto_SOC = dual_model.addConstr(-rawMeOH_sto[1:] + rawMeOH_sto[:-1] + raw_MeOH_min[1:-1] - raw_MeOH_max[1:-1] <= 0, name='m_rawMeOH_sto_SOC')
    m_rawMeOH_sto_SOC_h = dual_model.addConstr(fin_MeOH_sto + rawMeOH_sto[data.hours-2] + raw_MeOH_min[data.hours-1] - raw_MeOH_max[data.hours-1] <= 0, name='m_rawMeOH_sto_SOC_h')
    m_rawMeOH_sto_cha = dual_model.addConstr(-rawMeOH_sto*data.raw_MeOH_sto[0] - C_raw_MeOH_1[:-1] - compr_MeOH[:-1] <= 0, name='m_rawMeOH_sto_cha')
    m_rawMeOH_sto_dis = dual_model.addConstr(+rawMeOH_sto/data.raw_MeOH_sto[1] - C_raw_MeOH_2[:-1] - dest_rawMeOH[:-1] <= 0, name='m_rawMeOH_sto_dis')
    m_rawMeOH_sto_cha_h = dual_model.addConstr( - C_raw_MeOH_1[data.hours-1] - compr_MeOH[data.hours-1] <= 0, name='m_rawMeOH_sto_cha')
    m_rawMeOH_sto_dis_h = dual_model.addConstr( - C_raw_MeOH_2[data.hours-1] - dest_rawMeOH[data.hours-1] <= 0, name='m_rawMeOH_sto_dis')

    # Methanol distilator
    m_rawMeOH_dest = dual_model.addConstr(- compr_MeOH - dest_rawMeOH <= 0, name='m_rawMeOH_dest')
    q_dest = dual_model.addConstr(- dest_heat + p_flow*data.Heater[0] <= 0, name='q_dest')
    m_h2O_dest = dual_model.addConstr(- dest_water <= 0, name='m_h2O_dest')
    m_MeOH_dest_0 = dual_model.addConstr(- dest_max[0] + dest_min[0] + dest_rawMeOH[0]*data.Destilator[0] + dest_water[0]*data.Destilator[1] + dest_heat[0]*data.Destilator[3] + dest_sto[0] + dest_ramp_up[0] - dest_ramp_down[0] <= 0, name='m_MeOH_dest_0')
    m_MeOH_dest = dual_model.addConstr(- dest_max[1:-1] + dest_min[1:-1] + dest_rawMeOH[1:-1]*data.Destilator[0] + dest_water[1:-1]*data.Destilator[1] + dest_heat[1:-1]*data.Destilator[3] + dest_sto[1:-1] - dest_ramp_up[:-1] + dest_ramp_down[:-1] + dest_ramp_up[1:] - dest_ramp_down[1:] <= 0, name='m_MeOH_dest')
    m_MeOH_dest_h = dual_model.addConstr(- dest_max[data.hours-1] + dest_min[data.hours-1] + dest_rawMeOH[data.hours-1]*data.Destilator[0] + dest_water[data.hours-1]*data.Destilator[1] + dest_heat[data.hours-1]*data.Destilator[3] + dest_sto[data.hours-1] - dest_ramp_up[data.hours-2] + dest_ramp_down[data.hours-2] <= 0, name='m_MeOH_dest_h')

    # Pure MeOH Storage
    m_MeOH_SOC_0 = dual_model.addConstr(init_pure - pure_sto[0] + MeOH_min[0] - MeOH_max[0] <= 0, name="m_MeOH_SOC_0")
    m_MeOH_SOC = dual_model.addConstr(- pure_sto[1:] + pure_sto[:-1] + MeOH_min[1:-1] - MeOH_max[1:-1] <= 0, name="m_MeOH_SOC")
    m_MeOH_SOC_h = dual_model.addConstr(fin_pure + pure_sto[data.hours-2] + MeOH_min[data.hours-1] - MeOH_max[data.hours-1] <= 0, name="m_MeOH_SOC_h")
    m_MeOH_cha = dual_model.addConstr(- dest_sto[:-1] - C_rate_MeOH_1[:-1] - pure_sto*data.Pure_MeOH_sto[0] <= 0, name='m_MeOH_cha')
    m_MeOH_dis = dual_model.addConstr(pure_sto/data.Pure_MeOH_sto[1] - C_rate_MeOH_2[:-1] + balance[:-1] <= 0, name="m_MeOH_dis") 
    m_MeOH_cha_h = dual_model.addConstr(- dest_sto[data.hours-1] - C_rate_MeOH_1[data.hours-1] <= 0, name='m_MeOH_cha_h')
    m_MeOH_dis_h = dual_model.addConstr(- C_rate_MeOH_2[data.hours-1] + balance[data.hours-1] <= 0, name="m_MeOH_dis_h") 
    
    
# OBJECTIVE FUNCTION

    # Dual objective

    dual_objective = (
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
        
    dual_model.setObjective(dual_objective, GRB.MAXIMIZE)
    # Write lp model
    # dual_model.update()
    # dual_model.write('sub.lp')

    # dual_model.setParam("OutputFlag", 0)
    # Solve dual model
    # dual_model.setParam('Method', 2)  # Barrier method
    # dual_model.setParam('Crossover', 0)
    # dual_model.setParam('MIPGap', 0.1)  # Allow a 0.1% gap
    # dual_model.Params.InfUnbdInfo = 1
    # dual_model.computeIIS()
    # dual_model.write("dual_model.ilp")
    dual_model.optimize()
    A = dual_model.getA()  # Returns a sparse matrix

    if dual_model.Status == GRB.OPTIMAL:

        # Print optimal objective value
        objective_val = dual_model.objVal

        duals = Sub_Results(
            p_flow=p_flow.X,
            max_grid=max_grid.X,
            init_bat=init_bat.X,
            fin_bat=fin_bat.X,
            bat_min=bat_min.X,
            bat_max=bat_max.X,
            C_bat_1=C_bat_1.X,
            C_bat_2=C_bat_2.X,
            elec_min=elec_min.X,
            elec_max=elec_max.X,
            elec_ramp_up=elec_ramp_up.X,
            elec_ramp_down=elec_ramp_down.X,
            init_h2=init_h2.X,
            fin_h2=fin_h2.X,
            h2_min=h2_min.X,
            h2_max=h2_max.X,
            C_rate_h2_1=C_rate_h2_1.X,
            C_rate_h2_2=C_rate_h2_2.X,
            init_CO2=init_CO2.X,
            fin_CO2=fin_CO2.X,
            CO2_min=CO2_min.X,
            CO2_max=CO2_max.X,
            CO2_str_min=CO2_str_min.X,
            CO2_str_max=CO2_str_max.X,
            reac_min=reac_min.X,
            reac_max=reac_max.X,
            Reac_ramp_up=Reac_ramp_up.X,
            Reac_ramp_down=Reac_ramp_down.X,
            init_MeOH_sto=init_MeOH_sto.X,
            fin_MeOH_sto=fin_MeOH_sto.X,
            raw_MeOH_min=raw_MeOH_min.X,
            raw_MeOH_max=raw_MeOH_max.X,
            C_raw_MeOH_1=C_raw_MeOH_1.X,
            C_raw_MeOH_2=C_raw_MeOH_2.X,
            dest_min=dest_min.X,
            dest_max=dest_max.X,
            dest_ramp_up=dest_ramp_up.X,
            dest_ramp_down=dest_ramp_down.X,
            init_pure=init_pure.X,
            fin_pure=fin_pure.X,
            MeOH_min=MeOH_min.X,
            MeOH_max=MeOH_max.X,
            C_rate_MeOH_1=C_rate_MeOH_1.X,
            C_rate_MeOH_2=C_rate_MeOH_2.X,
            balance=balance.X,
        )

        return True, objective_val, duals
    
    elif dual_model.Status == GRB.INF_OR_UNBD:
        return False, None, None
    else:
        # Unknown or unhandled status
        raise ValueError(f"Subproblem not solved. Solver status: {dual_model.Status}")