import gurobipy as gp
from gurobipy import GRB
import math
from DW_X1_results import X1_Results
from sub_results import *

def DW_sub_1(data,P_cap,my_pi,my_kappa,res):

#   Sub model NOTE (setupSUB)
    sub_DW = gp.Model("DW_sub")

# VARIABLES
    # Power flow 
    p_flow = sub_DW.addMVar(4380, lb=-math.inf, ub=10000, name="p_flow")

    max_grid = sub_DW.addMVar(4380, lb=0, name="max_grid")

    # Electrolyzer operation
    elec_max = sub_DW.addMVar(4380, lb=0, name="Elec_max")
    elec_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="Elec_min")
    elec_op = sub_DW.addMVar(4380, lb=-math.inf, name="Elec_op")
    elec_H2O = sub_DW.addMVar(4380, lb=-math.inf, name="elec_H2O")
    elec_ramp_up = sub_DW.addMVar(4380, lb=0, ub=10000, name="elec_ramp_up")
    elec_ramp_down = sub_DW.addMVar(4380, lb=0, name="elec_ramp_down")

    # Battery
    init_bat = sub_DW.addMVar(1, lb=-math.inf, name="init_bat")
    bat_sto = sub_DW.addMVar(4380, lb=-10000, ub=10000, name="bat_sto")
    bat_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="bat_min")
    bat_max = sub_DW.addMVar(4380, lb=0, name="bat_max")
    C_bat_1 = sub_DW.addMVar(4380, lb=0, name="C_bat_1")
    C_bat_2 = sub_DW.addMVar(4380, lb=0, name="C_bat_2")

    # Desalinator
    des_H2O = sub_DW.addMVar(4380, lb=-math.inf, name="des_H2O")

    # Hydrogen storage
    init_h2 = sub_DW.addMVar(1, lb=-math.inf, name="initialize_h2")
    h2_storage = sub_DW.addMVar(4380, lb=-10000, name="h2_storage")
    h2_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="h2_min")
    h2_max = sub_DW.addMVar(4380, lb=0, name="h2_max")
    C_rate_h2_1 = sub_DW.addMVar(4380, lb=0, name="C_rate_h2_1")
    C_rate_h2_2 = sub_DW.addMVar(4380, lb=0, name="C_rate_h2_2")

    # H2 Compressor
    compressor_h2 = sub_DW.addMVar(4380, lb=-math.inf, name="compressor_h2")

    # CO2 Storage
    init_CO2 = sub_DW.addMVar(1, lb=-math.inf, name="initialize_CO2")
    CO2_sto = sub_DW.addMVar(4380, lb=-10000, name="CO2_storage")
    CO2_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="CO2_min")
    CO2_max = sub_DW.addMVar(4380, lb=0, name="CO2_max")
    CO2_str_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="CO2_stream_min")
    CO2_str_max = sub_DW.addMVar(4380, lb=0, name="CO2_stream_max")
    CO2_comp = sub_DW.addMVar(4380, lb=-math.inf, name="CO2_comp")

    # Feed syngas Compressor
    compressor_e = sub_DW.addMVar(4380, lb=-math.inf, name="compressor_e")

    # Methanol reactor operation
    reac_max = sub_DW.addMVar(4380, lb=0, name="Reac_max")
    reac_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="Reac_min")
    reac_H2 =  sub_DW.addMVar(4380, lb=-math.inf, name="Reac_H2")
    reac_CO2 = sub_DW.addMVar(4380, lb=-math.inf, name="reac_CO2")
    reac_el = sub_DW.addMVar(4380, lb=-math.inf, name="Reac_el")
    reac_heat = sub_DW.addMVar(4380, lb=-math.inf, name="Reac_heat") 
    Reac_ramp_up = sub_DW.addMVar(4380, lb=0, ub=10000, name="Reac_ramp_up")
    Reac_ramp_down = sub_DW.addMVar(4380, lb=0, name="Reac_ramp_down")

    # Methanol compressor
    compr_MeOH = sub_DW.addMVar(4380, lb=-math.inf, name="compr_MeOH")

    # Raw Methanol storage
    init_MeOH_sto = sub_DW.addMVar(1, lb=-math.inf, name="initialize_MeOH_sto")
    rawMeOH_sto = sub_DW.addMVar(4380, lb=-10000, ub=10000, name="rawMeOH_storage")
    raw_MeOH_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="raw_MeOH_min")
    raw_MeOH_max = sub_DW.addMVar(4380, lb=0, name="raw_MeOH_max")
    C_raw_MeOH_1 = sub_DW.addMVar(4380, lb=0, name="C_rate_raw_MeOH_1")
    C_raw_MeOH_2 = sub_DW.addMVar(4380, lb=0, name="C_rate_raw_MeOH_2")

    # Methanol distillator
    dest_max = sub_DW.addMVar(4380, lb=0, name="Dest_max")
    dest_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="Dest_min")
    dest_rawMeOH = sub_DW.addMVar(4380, lb=-math.inf, name="Dest_rawMeOH")
    dest_water = sub_DW.addMVar(4380, lb=-math.inf, name="Dest_water")
    dest_heat = sub_DW.addMVar(4380, lb=-math.inf, name="Dest_heat")
    dest_sto = sub_DW.addMVar(4380, lb=-math.inf, name="Dest_sto")
    dest_ramp_up = sub_DW.addMVar(4380, lb=0, ub=10000, name="Dest_ramp_up")
    dest_ramp_down = sub_DW.addMVar(4380, lb=0, name="Dest_ramp_down")

    # Pure Methanol storage
    init_pure = sub_DW.addMVar(1, lb=-math.inf, name="initialize_pure")
    pure_sto = sub_DW.addMVar(4380, lb=-10000, name="pure_storage")
    MeOH_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="MeOH_min")
    MeOH_max = sub_DW.addMVar(4380, lb=0, name="MeOH_max")
    C_rate_MeOH_1 = sub_DW.addMVar(4380, lb=0, name="C_rate_MeOH_1")
    C_rate_MeOH_2 = sub_DW.addMVar(4380, lb=0, name="C_rate_MeOH_2")   

    # Meet demand
    balance = sub_DW.addMVar(4380, lb=0, ub=10000, name="balance")

#   SUB CONSTRAINTS
    # Energy
    e_buy = sub_DW.addConstr(-p_flow - max_grid <= data.Prices[:4380], name='e_buy')
    e_sell = sub_DW.addConstr(p_flow <= 0, name='e_sell')

    # Battery
    e_bat_SOC_0 = sub_DW.addConstr(init_bat - bat_sto[0] + bat_min[0] - bat_max[0] <= 0)
    e_bat_SOC = sub_DW.addConstr(bat_sto[:-1] - bat_sto[1:] + bat_min[1:] - bat_max[1:] <= 0)
    e_bat_cha = sub_DW.addConstr(- bat_sto*data.Battery[0] + p_flow - C_bat_1 <= 0)
    e_bat_dis = sub_DW.addConstr(+ bat_sto/data.Battery[1] - p_flow - C_bat_2 <= 0)

    # Desalinator
    e_des = sub_DW.addConstr(des_H2O + p_flow <= 0, name='e_des')

    # Electrolyzer
    e_elec_0 = sub_DW.addConstr(p_flow[0] - elec_max[0] + elec_min[0] - elec_op[0] + elec_ramp_up[0] - elec_ramp_down[0] <= 0, name='e_elec_0')
    e_elec = sub_DW.addConstr(p_flow[1:] - elec_max[1:] + elec_min[1:] - elec_op[1:] - elec_ramp_up[:-1] + elec_ramp_down[:-1] + elec_ramp_up[1:] - elec_ramp_down[1:] <= 0, name='e_elec')
    m_h2_elec = sub_DW.addConstr(elec_op*data.Electrolyzer[0] + elec_H2O*data.Electrolyzer[1] + compressor_h2 <= 0, name='m_h2_elec')
    m_H2O_elec = sub_DW.addConstr(- elec_H2O - des_H2O*data.Desalinator[0] <= 0, name='m_H2O_elec')
    
    # H2 Storage
    m_h2_sto_SOC_0 = sub_DW.addConstr(init_h2 - h2_storage[0] + h2_min[0] - h2_max[0] <= 0, name='m_h2_sto_SOC_0')
    m_h2_sto_SOC = sub_DW.addConstr(h2_storage[:-1] - h2_storage[1:] + h2_min[1:] - h2_max[1:] <= 0, name='m_h2_sto_SOC')
    m_h2_sto_cha = sub_DW.addConstr(- h2_storage*data.H2storage[0] - C_rate_h2_1 - compressor_h2 <= 0, name='m_h2_sto_cha')
    m_h2_sto_dis = sub_DW.addConstr(h2_storage/data.H2storage[1] - C_rate_h2_2 - reac_H2 - compressor_e*data.Compressor[0] <= 0, name='m_h2_sto_dis')

    # CO2 Storage
    m_CO2_sto_SOC_0 = sub_DW.addConstr(init_CO2 - CO2_sto[0] + CO2_min[0] - CO2_max[0] <= 0, name='m_CO2_sto_SOC_0')
    m_CO2_sto_SOC = sub_DW.addConstr(CO2_sto[:-1] - CO2_sto[1:] + CO2_min[1:] - CO2_max[1:] <= 0, name='m_CO2_sto_SOC')
    m_CO2_buy = sub_DW.addConstr(CO2_str_min - CO2_str_max + CO2_comp <= data.CO2_prices[:4380], name='m_CO2_buy')
    m_CO2_sto_cha = sub_DW.addConstr(- CO2_sto*data.CO2storage[0] - CO2_comp <= 0, name='m_CO2_sto_cha')
    m_CO2_sto_dis = sub_DW.addConstr(CO2_sto/data.CO2storage[1] - compressor_e*data.Compressor[0] - reac_CO2 <= 0, name='m_CO2_sto_dis')

    # Syngas compressor
    e_compr = sub_DW.addConstr(p_flow + compressor_e <= 0, name='e_compr')

    # Methanol reactor
    m_h2_reac = sub_DW.addConstr(- compressor_h2 - reac_H2 - compressor_e*data.Compressor[0] <= 0, name='m_h2_reac')
    m_CO2_reac = sub_DW.addConstr(- CO2_comp - reac_CO2 - compressor_e*data.Compressor[0] <= 0, name='m_CO2_reac') 
    e_reac = sub_DW.addConstr(p_flow - reac_el <= 0, name='e_reac')
    q_reac = sub_DW.addConstr(- reac_heat + p_flow*data.Heater[0] <= 0, name='q_reac')
    m_rawMeOH_reac_0 = sub_DW.addConstr(-reac_max[0] + reac_min[0] + reac_H2[0]*data.Reactor[0] + reac_CO2[0]*data.Reactor[1] + reac_el[0]*data.Reactor[2] + reac_heat[0]*data.Reactor[3] + Reac_ramp_up[0] - Reac_ramp_down[0] + compr_MeOH[0] <= 0, name='m_rawMeOH_reac_0')
    m_rawMeOH_reac = sub_DW.addConstr(-reac_max[1:] + reac_min[1:] + reac_H2[1:]*data.Reactor[0] + reac_CO2[1:]*data.Reactor[1] + reac_el[1:]*data.Reactor[2] + reac_heat[1:]*data.Reactor[3] - Reac_ramp_up[:-1] + Reac_ramp_down[:-1] + Reac_ramp_up[1:] - Reac_ramp_down[1:] + compr_MeOH[1:] <= 0, name='m_rawMeOH_reac')

    # MeOH Storage
    m_rawMeOH_sto_SOC_0 = sub_DW.addConstr(init_MeOH_sto - rawMeOH_sto[0] + raw_MeOH_min[0] - raw_MeOH_max[0] <= 0, name='m_rawMeOH_sto_SOC_0')
    m_rawMeOH_sto_SOC = sub_DW.addConstr(-rawMeOH_sto[1:] + rawMeOH_sto[:-1] + raw_MeOH_min[1:] - raw_MeOH_max[1:] <= 0, name='m_rawMeOH_sto_SOC')
    m_rawMeOH_sto_cha = sub_DW.addConstr(-rawMeOH_sto*data.raw_MeOH_sto[0] - C_raw_MeOH_1 - compr_MeOH <= 0, name='m_rawMeOH_sto_cha')
    m_rawMeOH_sto_dis = sub_DW.addConstr(+rawMeOH_sto/data.raw_MeOH_sto[1] - C_raw_MeOH_2 - dest_rawMeOH <= 0, name='m_rawMeOH_sto_dis')

    # Methanol distilator
    m_rawMeOH_dest = sub_DW.addConstr(- compr_MeOH - dest_rawMeOH <= 0, name='m_rawMeOH_dest')
    q_dest = sub_DW.addConstr(- dest_heat + p_flow*data.Heater[0] <= 0, name='q_dest')
    m_h2O_dest = sub_DW.addConstr(- dest_water <= 0, name='m_h2O_dest')
    m_MeOH_dest_0 = sub_DW.addConstr(- dest_max[0] + dest_min[0] + dest_rawMeOH[0]*data.Destilator[0] + dest_water[0]*data.Destilator[1] + dest_heat[0]*data.Destilator[3] + dest_sto[0] + dest_ramp_up[0] - dest_ramp_down[0] <= 0, name='m_MeOH_dest_0')
    m_MeOH_dest = sub_DW.addConstr(- dest_max[1:] + dest_min[1:] + dest_rawMeOH[1:]*data.Destilator[0] + dest_water[1:]*data.Destilator[1] + dest_heat[1:]*data.Destilator[3] + dest_sto[1:] - dest_ramp_up[:-1] + dest_ramp_down[:-1] + dest_ramp_up[1:] - dest_ramp_down[1:] <= 0, name='m_MeOH_dest')

    # Pure MeOH Storage
    m_MeOH_SOC_0 = sub_DW.addConstr(init_pure - pure_sto[0] + MeOH_min[0] - MeOH_max[0] <= 0, name="m_MeOH_SOC_0")
    m_MeOH_SOC = sub_DW.addConstr(- pure_sto[1:] + pure_sto[:-1] + MeOH_min[1:] - MeOH_max[1:] <= 0, name="m_MeOH_SOC")
    m_MeOH_cha = sub_DW.addConstr(- dest_sto - C_rate_MeOH_1 - pure_sto*data.Pure_MeOH_sto[0] <= 0, name='m_MeOH_cha')
    m_MeOH_dis = sub_DW.addConstr(pure_sto/data.Pure_MeOH_sto[1] - C_rate_MeOH_2 + balance <= 0, name="m_MeOH_dis")




    # Sub objective

    c_x = (
    p_flow @ (data.Wind[:4380] * P_cap[0] + data.Solar[:4380] * P_cap[1])
    - gp.quicksum(max_grid * data.grid_max)                     # Grid
    + init_bat * P_cap[5] * data.Battery[2] 
    + gp.quicksum(
        bat_min * P_cap[5] * data.Battery[4]
        - bat_max * P_cap[5] * data.Battery[5]
        - C_bat_1 * P_cap[5] * data.Battery[6]
        - C_bat_2 * P_cap[5] * data.Battery[6]
    )                                                           # Battery
    + gp.quicksum(
        - elec_max * P_cap[2]
        + elec_min * P_cap[2] * data.op_elec[2]
    )
    + gp.quicksum(
        - elec_ramp_up * P_cap[2] * data.op_elec[0]
        - elec_ramp_down * P_cap[2] * data.op_elec[1]
    )                                                           # Electrolyzer
    + balance @ data.Demand[:4380]                                    # Demand
    + init_h2 * P_cap[6] * data.H2storage[2]
    + gp.quicksum(
        h2_min * P_cap[6] * data.H2storage[4]
        - h2_max * P_cap[6] * data.H2storage[5]
        - C_rate_h2_1 * P_cap[6] * data.H2storage[6]
        - C_rate_h2_2 * P_cap[6] * data.H2storage[6]
    )                                                           # H2 Storage
    + init_CO2 * P_cap[8] * data.CO2storage[2]
    + gp.quicksum(
        CO2_min * P_cap[8] * data.CO2storage[4]
        - CO2_max * P_cap[8] * data.CO2storage[5]
        + CO2_str_min * data.CO2storage[6]
        - CO2_str_max * data.CO2storage[7]
    )                                                           # CO2 Storage
    + gp.quicksum(
        - reac_max * P_cap[3]
        + reac_min * P_cap[3] * data.op_reac[2]
    )
    + gp.quicksum(
        - Reac_ramp_up * P_cap[3] * data.op_reac[0]
        - Reac_ramp_down * P_cap[3] * data.op_reac[1]
    )                                                           # Reactor
    + init_MeOH_sto * P_cap[7] * data.raw_MeOH_sto[2]
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
    + gp.quicksum(
        MeOH_min * P_cap[9] * data.Pure_MeOH_sto[4]
        - MeOH_max * P_cap[9] * data.Pure_MeOH_sto[5]
        - C_rate_MeOH_1 * P_cap[9] * data.Pure_MeOH_sto[6]
        - C_rate_MeOH_2 * P_cap[9] * data.Pure_MeOH_sto[6]
    )                                                           # Pure Storage
    )
        
    pi_A0_x = (
    # Battery SOC terms

    my_pi.e_bat_SOC * bat_sto[4379] +

    # Electricity flow terms
    my_pi.e_elec * ( - elec_ramp_up[4379] + elec_ramp_down[4379] )+

    # Hydrogen storage SOC terms
    my_pi.m_h2_sto_SOC * h2_storage[4379] +

    # CO2 storage SOC terms
    my_pi.m_CO2_sto_SOC * CO2_sto[4379] +

    # Methanol reactor terms
    my_pi.m_rawMeOH_reac *  (- Reac_ramp_up[4379] + Reac_ramp_down[4379]) +

    # Raw methanol storage SOC terms
    my_pi.m_rawMeOH_sto_SOC * (- rawMeOH_sto[4379]) +

    # Methanol distillation terms
    my_pi.m_MeOH_dest * (- dest_ramp_up[4379] + dest_ramp_down[4379] ) +

    # Pure methanol storage SOC terms
    my_pi.m_MeOH_SOC * (pure_sto[4379]) )

# Set the objective in the model
    sub_DW.setObjective(c_x - pi_A0_x - my_kappa, GRB.MAXIMIZE)

    sub_DW.update()
    sub_DW.write('sub1.lp')

    sub_DW.optimize()


    if sub_DW.Status == GRB.OPTIMAL:

        # Print optimal objective value
        sub_val = sub_DW.objVal

        res.p_flow[:4380]=p_flow.X
        res.max_grid[:4380]=max_grid.X

        res.elec_max[:4380]=elec_max.X
        res.elec_min[:4380]=elec_min.X
        res.elec_op[:4380]=elec_op.X
        res.elec_H2O[:4380]=elec_H2O.X
        res.elec_ramp_up[:4380]=elec_ramp_up.X
        res.elec_ramp_down[:4380]=elec_ramp_down.X

        res.init_bat=init_bat.X
        res.bat_sto[:4380]=bat_sto.X
        res.bat_min[:4380]=bat_min.X
        res.bat_max[:4380]=bat_max.X
        res.C_bat_1[:4380]=C_bat_1.X
        res.C_bat_2[:4380]=C_bat_2.X
        
        res.des_H2O[:4380]=des_H2O.X

        res.init_h2=init_h2.X
        res.h2_storage[:4380]=h2_storage.X
        res.h2_min[:4380]=h2_min.X
        res.h2_max[:4380]=h2_max.X
        res.C_rate_h2_1[:4380]=C_rate_h2_1.X
        res.C_rate_h2_2[:4380]=C_rate_h2_2.X

        res.compressor_h2[:4380]=compressor_h2.X

        res.init_CO2=init_CO2.X
        res.CO2_sto[:4380]=CO2_sto.X
        res.CO2_min[:4380]=CO2_min.X
        res.CO2_max[:4380]=CO2_max.X
        res.CO2_str_min[:4380]=CO2_str_min.X
        res.CO2_str_max[:4380]=CO2_str_max.X
        res.CO2_comp[:4380]=CO2_comp.X

        res.compressor_e[:4380]=compressor_e.X

        res.reac_max[:4380]=reac_max.X
        res.reac_min[:4380]=reac_min.X
        res.reac_H2[:4380]=reac_H2.X
        res.reac_CO2[:4380]=reac_CO2.X
        res.reac_el[:4380]=reac_el.X
        res.reac_heat[:4380]=reac_heat.X
        res.Reac_ramp_up[:4380]=Reac_ramp_up.X
        res.Reac_ramp_down[:4380]=Reac_ramp_down.X

        res.compr_MeOH[:4380]=compr_MeOH.X

        res.init_MeOH_sto=init_MeOH_sto.X
        res.rawMeOH_sto[:4380]=rawMeOH_sto.X
        res.raw_MeOH_min[:4380]=raw_MeOH_min.X
        res.raw_MeOH_max[:4380]=raw_MeOH_max.X
        res.C_raw_MeOH_1[:4380]=C_raw_MeOH_1.X
        res.C_raw_MeOH_2[:4380]=C_raw_MeOH_2.X

        res.dest_max[:4380]=dest_max.X
        res.dest_min[:4380]=dest_min.X
        res.dest_rawMeOH[:4380]=dest_rawMeOH.X
        res.dest_water[:4380]=dest_water.X
        res.dest_heat[:4380]=dest_heat.X
        res.dest_sto[:4380]=dest_sto.X
        res.dest_ramp_up[:4380]=dest_ramp_up.X
        res.dest_ramp_down[:4380]=dest_ramp_down.X

        res.init_pure=init_pure.X
        res.pure_sto[:4380]=pure_sto.X
        res.MeOH_min[:4380]=MeOH_min.X
        res.MeOH_max[:4380]=MeOH_max.X
        res.C_rate_MeOH_1[:4380]=C_rate_MeOH_1.X
        res.C_rate_MeOH_2[:4380]=C_rate_MeOH_2.X

        res.balance[:4380]=balance.X

    return sub_val, res



def DW_sub_2(data,P_cap,my_pi,my_kappa,res):

#   Sub model NOTE (setupSUB)
    sub_DW = gp.Model("DW_sub")

# VARIABLES
    # Power flow 
    p_flow = sub_DW.addMVar(4380, lb=-10000, ub=10000, name="p_flow")

    max_grid = sub_DW.addMVar(4380, lb=0, name="max_grid")

    # Electrolyzer operation
    elec_max = sub_DW.addMVar(4380, lb=0, ub=10000, name="Elec_max")
    elec_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="Elec_min")
    elec_op = sub_DW.addMVar(4380, lb=-10000, ub=10000, name="Elec_op")
    elec_H2O = sub_DW.addMVar(4380, lb=-10000, name="elec_H2O")
    elec_ramp_up = sub_DW.addMVar(4380-1, lb=0, name="elec_ramp_up")
    elec_ramp_down = sub_DW.addMVar(4380-1, lb=0, ub=10000, name="elec_ramp_down")

    # Battery
    fin_bat = sub_DW.addMVar(1, lb=0, name="fin_bat")
    bat_sto = sub_DW.addMVar(4380-1, lb=-10000, ub=10000, name="bat_sto")
    bat_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="bat_min")
    bat_max = sub_DW.addMVar(4380, lb=0, ub=10000, name="bat_max")
    C_bat_1 = sub_DW.addMVar(4380, lb=0, name="C_bat_1")
    C_bat_2 = sub_DW.addMVar(4380, lb=0, name="C_bat_2")

    # Desalinator
    des_H2O = sub_DW.addMVar(4380, lb=-math.inf, name="des_H2O")

    # Hydrogen storage
    fin_h2 = sub_DW.addMVar(1, lb=0, name="finalize_h2")
    h2_storage = sub_DW.addMVar(4380-1, lb=-10000, ub=10000, name="h2_storage")
    h2_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="h2_min")
    h2_max = sub_DW.addMVar(4380, lb=0, ub=10000, name="h2_max")
    C_rate_h2_1 = sub_DW.addMVar(4380, lb=0, name="C_rate_h2_1")
    C_rate_h2_2 = sub_DW.addMVar(4380, lb=0, name="C_rate_h2_2")

    # H2 Compressor
    compressor_h2 = sub_DW.addMVar(4380, lb=-math.inf, name="compressor_h2")

    # CO2 Storage
    fin_CO2 = sub_DW.addMVar(1, lb=0, name="finalize_CO2")
    CO2_sto = sub_DW.addMVar(4380-1, lb=-10000, name="CO2_storage")
    CO2_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="CO2_min")
    CO2_max = sub_DW.addMVar(4380, lb=0, ub=10000, name="CO2_max")
    CO2_str_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="CO2_stream_min")
    CO2_str_max = sub_DW.addMVar(4380, lb=0, ub=10000, name="CO2_stream_max")
    CO2_comp = sub_DW.addMVar(4380, lb=-10000, name="CO2_comp")

    # Feed syngas Compressor
    compressor_e = sub_DW.addMVar(4380, lb=-math.inf, name="compressor_e")

    # Methanol reactor operation
    reac_max = sub_DW.addMVar(4380, lb=0, ub=10000, name="Reac_max")
    reac_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="Reac_min")
    reac_H2 =  sub_DW.addMVar(4380, lb=-10000, ub=10000, name="Reac_H2")
    reac_CO2 = sub_DW.addMVar(4380, lb=-10000, ub=10000, name="reac_CO2")
    reac_el = sub_DW.addMVar(4380, lb=-10000, ub=10000, name="Reac_el")
    reac_heat = sub_DW.addMVar(4380, lb=-10000, ub=10000, name="Reac_heat") 
    Reac_ramp_up = sub_DW.addMVar(4380-1, lb=0, ub=10000, name="Reac_ramp_up")
    Reac_ramp_down = sub_DW.addMVar(4380-1, lb=0, ub=10000, name="Reac_ramp_down")

    # Methanol compressor
    compr_MeOH = sub_DW.addMVar(4380, lb=-math.inf, name="compr_MeOH")

    # Raw Methanol storage
    fin_MeOH_sto = sub_DW.addMVar(1, lb=0, name="finalize_MeOH_sto")
    rawMeOH_sto = sub_DW.addMVar(4380-1, lb=-10000, ub=10000, name="rawMeOH_storage")
    raw_MeOH_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="raw_MeOH_min")
    raw_MeOH_max = sub_DW.addMVar(4380, lb=0, ub=10000, name="raw_MeOH_max")
    C_raw_MeOH_1 = sub_DW.addMVar(4380, lb=0, name="C_rate_raw_MeOH_1")
    C_raw_MeOH_2 = sub_DW.addMVar(4380, lb=0, name="C_rate_raw_MeOH_2")

    # Methanol distillator
    dest_max = sub_DW.addMVar(4380, lb=0, ub=10000, name="Dest_max")
    dest_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="Dest_min")
    dest_rawMeOH = sub_DW.addMVar(4380, lb=-math.inf, ub=10000, name="Dest_rawMeOH")
    dest_water = sub_DW.addMVar(4380, lb=-math.inf, ub=10000, name="Dest_water")
    dest_heat = sub_DW.addMVar(4380, lb=-math.inf, ub=10000, name="Dest_heat")
    dest_sto = sub_DW.addMVar(4380, lb=-math.inf, ub=10000, name="Dest_sto")
    dest_ramp_up = sub_DW.addMVar(4380-1, lb=0, ub=10000, name="Dest_ramp_up")
    dest_ramp_down = sub_DW.addMVar(4380-1, lb=0, ub=10000, name="Dest_ramp_down")

    # Pure Methanol storage
    fin_pure = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_pure")
    pure_sto = sub_DW.addMVar(4380-1, lb=-10000, ub=10000, name="pure_storage")
    MeOH_min = sub_DW.addMVar(4380, lb=0, ub=10000, name="MeOH_min")
    MeOH_max = sub_DW.addMVar(4380, lb=0, ub=10000, name="MeOH_max")
    C_rate_MeOH_1 = sub_DW.addMVar(4380, lb=0, ub=10000, name="C_rate_MeOH_1")
    C_rate_MeOH_2 = sub_DW.addMVar(4380, lb=0, ub=10000, name="C_rate_MeOH_2")   

    # Meet demand
    balance = sub_DW.addMVar(4380, lb=0, ub=10000, name="balance")

#   SUB CONSTRAINTS
    # Energy
    e_buy = sub_DW.addConstr(-p_flow - max_grid <= data.Prices[4380:], name='e_buy')
    e_sell = sub_DW.addConstr(p_flow <= 0, name='e_sell')

    # Battery
    e_bat_SOC = sub_DW.addConstr(bat_sto[:-1] - bat_sto[1:] + bat_min[1:-1] - bat_max[1:-1] <= 0)
    e_bat_SOC_h = sub_DW.addConstr(fin_bat + bat_sto[4378] + bat_min[4379] - bat_max[4379] <= 0)
    e_bat_cha = sub_DW.addConstr(- bat_sto*data.Battery[0] + p_flow[:-1] - C_bat_1[:-1] <= 0)
    e_bat_dis = sub_DW.addConstr(+ bat_sto/data.Battery[1] - p_flow[:-1] - C_bat_2[:-1] <= 0)
    e_bat_cha_h = sub_DW.addConstr(+ p_flow[4379] - C_bat_1[4379] <= 0)
    e_bat_dis_h = sub_DW.addConstr(- p_flow[4379] - C_bat_2[4379] <= 0)

    # Desalinator
    e_des = sub_DW.addConstr(des_H2O + p_flow <= 0, name='e_des')

    # Electrolyzer
    e_elec = sub_DW.addConstr(p_flow[1:-1] - elec_max[1:-1] + elec_min[1:-1] - elec_op[1:-1] - elec_ramp_up[:-1] + elec_ramp_down[:-1] + elec_ramp_up[1:] - elec_ramp_down[1:] <= 0, name='e_elec')
    e_elec_h = sub_DW.addConstr(p_flow[4379] - elec_max[4379] + elec_min[4379] - elec_op[4379] - elec_ramp_up[4378] + elec_ramp_down[4378] <= 0, name='e_elec_h')
    m_h2_elec = sub_DW.addConstr(elec_op*data.Electrolyzer[0] + elec_H2O*data.Electrolyzer[1] + compressor_h2 <= 0, name='m_h2_elec')
    m_H2O_elec = sub_DW.addConstr(- elec_H2O - des_H2O*data.Desalinator[0] <= 0, name='m_H2O_elec')
    
    # H2 Storage
    m_h2_sto_SOC = sub_DW.addConstr(h2_storage[:-1] - h2_storage[1:] + h2_min[1:-1] - h2_max[1:-1] <= 0, name='m_h2_sto_SOC')
    m_h2_sto_SOC_h = sub_DW.addConstr(fin_h2 + h2_storage[4378] + h2_min[4379] - h2_max[4379] <= 0, name='m_h2_sto_SOC_h')
    m_h2_sto_cha = sub_DW.addConstr(- h2_storage*data.H2storage[0] - C_rate_h2_1[:-1] - compressor_h2[:-1] <= 0, name='m_h2_sto_cha')
    m_h2_sto_dis = sub_DW.addConstr(h2_storage/data.H2storage[1] - C_rate_h2_2[:-1] - reac_H2[:-1] - compressor_e[:-1]*data.Compressor[0] <= 0, name='m_h2_sto_dis')
    m_h2_sto_cha_h = sub_DW.addConstr(- C_rate_h2_1[4379] - compressor_h2[4379] <= 0, name='m_h2_sto_cha_h')
    m_h2_sto_dis_h = sub_DW.addConstr(- C_rate_h2_2[4379] - reac_H2[4379] - compressor_e[4379]*data.Compressor[0] <= 0, name='m_h2_sto_dis_h')

    # CO2 Storage
    m_CO2_sto_SOC = sub_DW.addConstr(CO2_sto[:-1] - CO2_sto[1:] + CO2_min[1:-1] - CO2_max[1:-1] <= 0, name='m_CO2_sto_SOC')
    m_CO2_sto_SOC_h = sub_DW.addConstr(fin_CO2 + CO2_sto[4378] + CO2_min[4379] - CO2_max[4379] <= 0, name='m_CO2_sto_SOC_h') 
    m_CO2_buy = sub_DW.addConstr(CO2_str_min - CO2_str_max + CO2_comp <= data.CO2_prices[4380:], name='m_CO2_buy')
    m_CO2_sto_cha = sub_DW.addConstr(- CO2_sto*data.CO2storage[0] - CO2_comp[:-1] <= 0, name='m_CO2_sto_cha')
    m_CO2_sto_dis = sub_DW.addConstr(CO2_sto/data.CO2storage[1] - compressor_e[:-1]*data.Compressor[0] - reac_CO2[:-1] <= 0, name='m_CO2_sto_dis')
    m_CO2_sto_cha_h = sub_DW.addConstr(- CO2_comp[4379] <= 0, name='m_CO2_sto_cha_h')
    m_CO2_sto_dis_h = sub_DW.addConstr(- compressor_e[4379]*data.Compressor[0] - reac_CO2[4379] <= 0, name='m_CO2_sto_dis_h')

    # Syngas compressor
    e_compr = sub_DW.addConstr(p_flow + compressor_e <= 0, name='e_compr')

    # Methanol reactor
    m_h2_reac = sub_DW.addConstr(- compressor_h2 - reac_H2 - compressor_e*data.Compressor[0] <= 0, name='m_h2_reac')
    m_CO2_reac = sub_DW.addConstr(- CO2_comp - reac_CO2 - compressor_e*data.Compressor[0] <= 0, name='m_CO2_reac') 
    e_reac = sub_DW.addConstr(p_flow - reac_el <= 0, name='e_reac')
    q_reac = sub_DW.addConstr(- reac_heat + p_flow*data.Heater[0] <= 0, name='q_reac')
    m_rawMeOH_reac = sub_DW.addConstr(-reac_max[1:-1] + reac_min[1:-1] + reac_H2[1:-1]*data.Reactor[0] + reac_CO2[1:-1]*data.Reactor[1] + reac_el[1:-1]*data.Reactor[2] + reac_heat[1:-1]*data.Reactor[3] - Reac_ramp_up[:-1] + Reac_ramp_down[:-1] + Reac_ramp_up[1:] - Reac_ramp_down[1:] + compr_MeOH[1:-1] <= 0, name='m_rawMeOH_reac')
    m_rawMeOH_reac_h = sub_DW.addConstr(-reac_max[4379] + reac_min[4379] + reac_H2[4379]*data.Reactor[0] + reac_CO2[4379]*data.Reactor[1] + reac_el[4379]*data.Reactor[2] + reac_heat[4379]*data.Reactor[3] - Reac_ramp_up[4378] + Reac_ramp_down[4378] + compr_MeOH[4379] <= 0, name='m_rawMeOH_reac_h')

    # MeOH Storage
    m_rawMeOH_sto_SOC = sub_DW.addConstr(-rawMeOH_sto[1:] + rawMeOH_sto[:-1] + raw_MeOH_min[1:-1] - raw_MeOH_max[1:-1] <= 0, name='m_rawMeOH_sto_SOC')
    m_rawMeOH_sto_SOC_h = sub_DW.addConstr(fin_MeOH_sto + rawMeOH_sto[4378] + raw_MeOH_min[4379] - raw_MeOH_max[4379] <= 0, name='m_rawMeOH_sto_SOC_h')
    m_rawMeOH_sto_cha = sub_DW.addConstr(-rawMeOH_sto*data.raw_MeOH_sto[0] - C_raw_MeOH_1[:-1] - compr_MeOH[:-1] <= 0, name='m_rawMeOH_sto_cha')
    m_rawMeOH_sto_dis = sub_DW.addConstr(+rawMeOH_sto/data.raw_MeOH_sto[1] - C_raw_MeOH_2[:-1] - dest_rawMeOH[:-1] <= 0, name='m_rawMeOH_sto_dis')
    m_rawMeOH_sto_cha_h = sub_DW.addConstr( - C_raw_MeOH_1[4379] - compr_MeOH[4379] <= 0, name='m_rawMeOH_sto_cha')
    m_rawMeOH_sto_dis_h = sub_DW.addConstr( - C_raw_MeOH_2[4379] - dest_rawMeOH[4379] <= 0, name='m_rawMeOH_sto_dis')

    # Methanol distilator
    m_rawMeOH_dest = sub_DW.addConstr(- compr_MeOH - dest_rawMeOH <= 0, name='m_rawMeOH_dest')
    q_dest = sub_DW.addConstr(- dest_heat + p_flow*data.Heater[0] <= 0, name='q_dest')
    m_h2O_dest = sub_DW.addConstr(- dest_water <= 0, name='m_h2O_dest')
    m_MeOH_dest = sub_DW.addConstr(- dest_max[1:-1] + dest_min[1:-1] + dest_rawMeOH[1:-1]*data.Destilator[0] + dest_water[1:-1]*data.Destilator[1] + dest_heat[1:-1]*data.Destilator[3] + dest_sto[1:-1] - dest_ramp_up[:-1] + dest_ramp_down[:-1] + dest_ramp_up[1:] - dest_ramp_down[1:] <= 0, name='m_MeOH_dest')
    m_MeOH_dest_h = sub_DW.addConstr(- dest_max[4379] + dest_min[4379] + dest_rawMeOH[4379]*data.Destilator[0] + dest_water[4379]*data.Destilator[1] + dest_heat[4379]*data.Destilator[3] + dest_sto[4379] - dest_ramp_up[4378] + dest_ramp_down[4378] <= 0, name='m_MeOH_dest_h')

    # Pure MeOH Storage
    m_MeOH_SOC = sub_DW.addConstr(- pure_sto[1:] + pure_sto[:-1] + MeOH_min[1:-1] - MeOH_max[1:-1] <= 0, name="m_MeOH_SOC")
    m_MeOH_SOC_h = sub_DW.addConstr(fin_pure + pure_sto[4378] + MeOH_min[4379] - MeOH_max[4379] <= 0, name="m_MeOH_SOC_h")
    m_MeOH_cha = sub_DW.addConstr(- dest_sto[:-1] - C_rate_MeOH_1[:-1] - pure_sto*data.Pure_MeOH_sto[0] <= 0, name='m_MeOH_cha')
    m_MeOH_dis = sub_DW.addConstr(pure_sto/data.Pure_MeOH_sto[1] - C_rate_MeOH_2[:-1] + balance[:-1] <= 0, name="m_MeOH_dis") 
    m_MeOH_cha_h = sub_DW.addConstr(- dest_sto[4379] - C_rate_MeOH_1[4379] <= 0, name='m_MeOH_cha_h')
    m_MeOH_dis_h = sub_DW.addConstr(- C_rate_MeOH_2[4379] + balance[4379] <= 0, name="m_MeOH_dis_h") 


    # Sub objective

    c_x = (
    p_flow @ (data.Wind[4380:] * P_cap[0] + data.Solar[4380:] * P_cap[1])
    - gp.quicksum(max_grid * data.grid_max)                         # Grid
    + fin_bat * P_cap[5] * data.Battery[3]
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
    + balance @ data.Demand[4380:]                                     # Demand
    + fin_h2 * P_cap[6] * data.H2storage[3]
    + gp.quicksum(
        h2_min * P_cap[6] * data.H2storage[4]
        - h2_max * P_cap[6] * data.H2storage[5]
        - C_rate_h2_1 * P_cap[6] * data.H2storage[6]
        - C_rate_h2_2 * P_cap[6] * data.H2storage[6]
    )                                                           # H2 Storage
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

    my_pi.e_bat_SOC * (- bat_sto[0] + bat_min[0] - bat_max[0]) +

    # Electricity flow terms
    my_pi.e_elec * (p_flow[0] - elec_max[0] + elec_min[0] - elec_op[0] + elec_ramp_up[0] - elec_ramp_down[0]) +

    # Hydrogen storage SOC terms
    my_pi.m_h2_sto_SOC * ( - h2_storage[0] + h2_min[0] - h2_max[0]) +

    # CO2 storage SOC terms
    my_pi.m_CO2_sto_SOC * ( - CO2_sto[0] + CO2_min[0] - CO2_max[0]) +

    # Methanol reactor terms
    my_pi.m_rawMeOH_reac * (-reac_max[0] + reac_min[0] + reac_H2[0] * data.Reactor[0] +
                                    reac_CO2[0] * data.Reactor[1] + reac_el[0] * data.Reactor[2] +
                                    reac_heat[0] * data.Reactor[3] + Reac_ramp_up[0] - Reac_ramp_down[0] + compr_MeOH[0]) +

    # Raw methanol storage SOC terms
    my_pi.m_rawMeOH_sto_SOC * (rawMeOH_sto[0] + raw_MeOH_min[0] - raw_MeOH_max[0]) +

    # Methanol distillation terms
    my_pi.m_MeOH_dest * (-dest_max[0] + dest_min[0] +
                                dest_rawMeOH[0] * data.Destilator[0] + dest_water[0] * data.Destilator[1] +
                                dest_heat[0] * data.Destilator[3] + dest_sto[0] + dest_ramp_up[0] - dest_ramp_down[0]) +

    # Pure methanol storage SOC terms
    my_pi.m_MeOH_SOC * ( - pure_sto[0] + MeOH_min[0] - MeOH_max[0])
    )

# Set the objective in the model
    sub_DW.setObjective(c_x - pi_A0_x - my_kappa, GRB.MAXIMIZE)


    sub_DW.update()
    sub_DW.write('sub2.lp')

    sub_DW.optimize()


    if sub_DW.Status == GRB.OPTIMAL:

        # Print optimal objective value
        sub_val = sub_DW.objVal

        res.p_flow[4380:]=p_flow.X
        res.max_grid[4380:]=max_grid.X

        res.elec_max[4380:]=elec_max.X
        res.elec_min[4380:]=elec_min.X
        res.elec_op[4380:]=elec_op.X
        res.elec_H2O[4380:]=elec_H2O.X
        res.elec_ramp_up[4380:]=elec_ramp_up.X
        res.elec_ramp_down[4380:]=elec_ramp_down.X

        res.fin_bat=fin_bat.X
        res.bat_sto[4380:]=bat_sto.X
        res.bat_min[4380:]=bat_min.X
        res.bat_max[4380:]=bat_max.X
        res.C_bat_1[4380:]=C_bat_1.X
        res.C_bat_2[4380:]=C_bat_2.X
        
        res.des_H2O[4380:]=des_H2O.X

        res.fin_h2=fin_h2.X
        res.h2_storage[4380:]=h2_storage.X
        res.h2_min[4380:]=h2_min.X
        res.h2_max[4380:]=h2_max.X
        res.C_rate_h2_1[4380:]=C_rate_h2_1.X
        res.C_rate_h2_2[4380:]=C_rate_h2_2.X

        res.compressor_h2[4380:]=compressor_h2.X

        res.fin_CO2=fin_CO2.X
        res.CO2_sto[4380:]=CO2_sto.X
        res.CO2_min[4380:]=CO2_min.X
        res.CO2_max[4380:]=CO2_max.X
        res.CO2_str_min[4380:]=CO2_str_min.X
        res.CO2_str_max[4380:]=CO2_str_max.X
        res.CO2_comp[4380:]=CO2_comp.X

        res.compressor_e[4380:]=compressor_e.X

        res.reac_max[4380:]=reac_max.X
        res.reac_min[4380:]=reac_min.X
        res.reac_H2[4380:]=reac_H2.X
        res.reac_CO2[4380:]=reac_CO2.X
        res.reac_el[4380:]=reac_el.X
        res.reac_heat[4380:]=reac_heat.X
        res.Reac_ramp_up[4380:]=Reac_ramp_up.X
        res.Reac_ramp_down[4380:]=Reac_ramp_down.X

        res.compr_MeOH[4380:]=compr_MeOH.X

        res.fin_MeOH_sto=fin_MeOH_sto.X
        res.rawMeOH_sto[4380:]=rawMeOH_sto.X
        res.raw_MeOH_min[4380:]=raw_MeOH_min.X
        res.raw_MeOH_max[4380:]=raw_MeOH_max.X
        res.C_raw_MeOH_1[4380:]=C_raw_MeOH_1.X
        res.C_raw_MeOH_2[4380:]=C_raw_MeOH_2.X

        res.dest_max[4380:]=dest_max.X
        res.dest_min[4380:]=dest_min.X
        res.dest_rawMeOH[4380:]=dest_rawMeOH.X
        res.dest_water[4380:]=dest_water.X
        res.dest_heat[4380:]=dest_heat.X
        res.dest_sto[4380:]=dest_sto.X
        res.dest_ramp_up[4380:]=dest_ramp_up.X
        res.dest_ramp_down[4380:]=dest_ramp_down.X

        res.fin_pure=fin_pure.X
        res.pure_sto[4380:]=pure_sto.X
        res.MeOH_min[4380:]=MeOH_min.X
        res.MeOH_max[4380:]=MeOH_max.X
        res.C_rate_MeOH_1[4380:]=C_rate_MeOH_1.X
        res.C_rate_MeOH_2[4380:]=C_rate_MeOH_2.X

        res.balance[4380:]=balance.X


    return sub_val, res