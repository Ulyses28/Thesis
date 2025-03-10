import gurobipy as gp
from gurobipy import GRB

def DW_sub_hour(data,P_cap,my_pi,my_kappa,results,h):

#   Sub model NOTE (setupSUB)
    sub_DW = gp.Model("DW_sub")

# VARIABLES
    # Power flow
    p_flow = sub_DW.addMVar(1, lb=-10000, ub=10000, name="p_flow")

    # Max grid
    max_grid = sub_DW.addMVar(1, lb=0, ub=10000, name="max_grid")

    # Electrolyzer operation
    elec_max = sub_DW.addMVar(1, lb=0, ub=10000, name="Elec_max")
    elec_min = sub_DW.addMVar(1, lb=0, ub=10000, name="Elec_min")
    elec_op = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Elec_op")
    elec_H2O = sub_DW.addMVar(1, lb=-10000, ub=10000, name="elec_H2O")
    elec_ramp_up = sub_DW.addMVar(1, lb=0, ub=10000, name="elec_ramp_up")
    elec_ramp_down = sub_DW.addMVar(1, lb=0, ub=10000, name="elec_ramp_down")

    # Battery
    init_bat = sub_DW.addMVar(1, lb=-1000, ub=100, name="init_bat")
    # fin_bat = sub_DW.addMVar(1, lb=0, ub=1000, name="fin_bat")
    bat_sto = sub_DW.addMVar(1, lb=-1000, ub=100, name="bat_sto")
    bat_min = sub_DW.addMVar(1, lb=0, ub=1000, name="bat_min")
    bat_max = sub_DW.addMVar(1, lb=0, ub=1000, name="bat_max")
    C_bat_1 = sub_DW.addMVar(1, lb=0, ub=1000, name="C_bat_1")
    C_bat_2 = sub_DW.addMVar(1, lb=0, ub=1000, name="C_bat_2")

    # Desalinator
    des_H2O = sub_DW.addMVar(1, lb=-1000, ub=1000, name="des_H2O")

    # Hydrogen storage
    init_h2 = sub_DW.addMVar(1, lb=-10000, ub=10000, name="initialize_h2")
    # fin_h2 = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_h2")
    h2_storage = sub_DW.addMVar(1, lb=-10000, ub=10000, name="h2_storage")
    h2_min = sub_DW.addMVar(1, lb=0, ub=10000, name="h2_min")
    h2_max = sub_DW.addMVar(1, lb=0, ub=10000, name="h2_max")
    C_rate_h2_1 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_h2_1")
    C_rate_h2_2 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_h2_2")

    # H2 Compressor
    compressor_h2 = sub_DW.addMVar(1, lb=-10000, ub=10000, name="compressor_h2")

    # CO2 Storage
    init_CO2 = sub_DW.addMVar(1, lb=-10000, ub=10000, name="initialize_CO2")
    # fin_CO2 = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_CO2")
    CO2_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="CO2_storage")
    CO2_min = sub_DW.addMVar(1, lb=0, ub=10000, name="CO2_min")
    CO2_max = sub_DW.addMVar(1, lb=0, ub=10000, name="CO2_max")
    CO2_str_min = sub_DW.addMVar(1, lb=0, ub=10000, name="CO2_stream_min")
    CO2_str_max = sub_DW.addMVar(1, lb=0, ub=10000, name="CO2_stream_max")
    CO2_comp = sub_DW.addMVar(1, lb=-10000, ub=10000, name="CO2_comp")

    # Feed syngas Compressor
    compressor_e = sub_DW.addMVar(1, lb=-10000, ub=10000, name="compressor_e")

    # Methanol reactor operation
    reac_max = sub_DW.addMVar(1, lb=0, ub=10000, name="Reac_max")
    reac_min = sub_DW.addMVar(1, lb=0, ub=10000, name="Reac_min")
    reac_H2 =  sub_DW.addMVar(1, lb=-10000, ub=10000, name="Reac_H2")
    reac_CO2 = sub_DW.addMVar(1, lb=-10000, ub=10000, name="reac_CO2")
    reac_el = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Reac_el")
    reac_heat = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Reac_heat") 
    Reac_ramp_up = sub_DW.addMVar(1, lb=0, ub=10000, name="Reac_ramp_up")
    Reac_ramp_down = sub_DW.addMVar(1, lb=0, ub=10000, name="Reac_ramp_down")

    # Methanol compressor
    compr_MeOH = sub_DW.addMVar(1, lb=-10000, ub=10000, name="compr_MeOH")

    # Raw Methanol storage
    init_MeOH_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="initialize_MeOH_sto")
    # fin_MeOH_sto = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_MeOH_sto")
    rawMeOH_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="rawMeOH_storage")
    raw_MeOH_min = sub_DW.addMVar(1, lb=0, ub=10000, name="raw_MeOH_min")
    raw_MeOH_max = sub_DW.addMVar(1, lb=0, ub=10000, name="raw_MeOH_max")
    C_raw_MeOH_1 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_raw_MeOH_1")
    C_raw_MeOH_2 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_raw_MeOH_2")

    # Methanol distillator
    dest_max = sub_DW.addMVar(1, lb=0, ub=10000, name="Dest_max")
    dest_min = sub_DW.addMVar(1, lb=0, ub=10000, name="Dest_min")
    dest_rawMeOH = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Dest_rawMeOH")
    dest_water = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Dest_water")
    dest_heat = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Dest_heat")
    dest_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Dest_sto")
    dest_ramp_up = sub_DW.addMVar(1, lb=0, ub=10000, name="Dest_ramp_up")
    dest_ramp_down = sub_DW.addMVar(1, lb=0, ub=10000, name="Dest_ramp_down")

    # Pure Methanol storage
    init_pure = sub_DW.addMVar(1, lb=-10000, ub=10000, name="initialize_pure")
    # fin_pure = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_pure")
    pure_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="pure_storage")
    MeOH_min = sub_DW.addMVar(1, lb=0, ub=10000, name="MeOH_min")
    MeOH_max = sub_DW.addMVar(1, lb=0, ub=10000, name="MeOH_max")
    C_rate_MeOH_1 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_MeOH_1")
    C_rate_MeOH_2 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_MeOH_2")

    # Meet demand
    balance = sub_DW.addMVar(1, lb=0, name="balance")

#   SUB CONSTRAINTS
    # Energy
    e_buy = sub_DW.addConstr(-p_flow - max_grid <= data.Prices[h], name='e_buy')
    e_sell = sub_DW.addConstr(p_flow <= 0, name='e_sell')

    # Battery
    e_bat_cha = sub_DW.addConstr(- bat_sto*data.Battery[0] + p_flow - C_bat_1 <= 0)
    e_bat_dis = sub_DW.addConstr(+ bat_sto/data.Battery[1] - p_flow - C_bat_2 <= 0)

    # Desalinator
    e_des = sub_DW.addConstr(des_H2O + p_flow <= 0, name='e_des')

    # Electrolyzer
    m_h2_elec = sub_DW.addConstr(elec_op*data.Electrolyzer[0] + elec_H2O*data.Electrolyzer[1] + compressor_h2 <= 0, name='m_h2_elec')
    m_H2O_elec = sub_DW.addConstr(- elec_H2O - des_H2O*data.Desalinator[0] <= 0, name='m_H2O_elec')

    # H2 Storage
    m_h2_sto_cha = sub_DW.addConstr(- h2_storage*data.H2storage[0] - C_rate_h2_1 - compressor_h2 <= 0, name='m_h2_sto_cha')
    m_h2_sto_dis = sub_DW.addConstr(h2_storage/data.H2storage[1] - C_rate_h2_2 - reac_H2 - compressor_e*data.Compressor[0] <= 0, name='m_h2_sto_dis')

    # CO2 Storage
    m_CO2_buy = sub_DW.addConstr(CO2_str_min - CO2_str_max + CO2_comp <= data.CO2_prices[h], name='m_CO2_buy')
    m_CO2_sto_cha = sub_DW.addConstr(- CO2_sto*data.CO2storage[0] - CO2_comp <= 0, name='m_CO2_sto_cha')
    m_CO2_sto_dis = sub_DW.addConstr(CO2_sto/data.CO2storage[1] - compressor_e*data.Compressor[0] - reac_CO2 <= 0, name='m_CO2_sto_dis')
    
    # Syngas compressor
    e_compr = sub_DW.addConstr(p_flow + compressor_e <= 0, name='e_compr')

    # Methanol reactor
    m_h2_reac = sub_DW.addConstr(- compressor_h2 - reac_H2 - compressor_e*data.Compressor[0] <= 0, name='m_h2_reac')
    m_CO2_reac = sub_DW.addConstr(- CO2_comp - reac_CO2 - compressor_e*data.Compressor[0] <= 0, name='m_CO2_reac') 
    e_reac = sub_DW.addConstr(p_flow - reac_el <= 0, name='e_reac')
    q_reac = sub_DW.addConstr(- reac_heat + p_flow*data.Heater[0] <= 0, name='q_reac')

    # MeOH Storage
    m_rawMeOH_sto_cha = sub_DW.addConstr(-rawMeOH_sto*data.raw_MeOH_sto[0] - C_raw_MeOH_1 - compr_MeOH <= 0, name='m_rawMeOH_sto_cha')
    m_rawMeOH_sto_dis = sub_DW.addConstr(+rawMeOH_sto/data.raw_MeOH_sto[1] - C_raw_MeOH_2 - dest_rawMeOH <= 0, name='m_rawMeOH_sto_dis')

    # Methanol distilator
    m_rawMeOH_dest = sub_DW.addConstr(- compr_MeOH - dest_rawMeOH <= 0, name='m_rawMeOH_dest')
    q_dest = sub_DW.addConstr(- dest_heat + p_flow*data.Heater[0] <= 0, name='q_dest')
    m_h2O_dest = sub_DW.addConstr(- dest_water <= 0, name='m_h2O_dest')

    # Pure MeOH Storagess
    m_MeOH_cha = sub_DW.addConstr(- dest_sto - C_rate_MeOH_1 - pure_sto*data.Pure_MeOH_sto[0] <= 0, name='m_MeOH_cha')
    m_MeOH_dis = sub_DW.addConstr(pure_sto/data.Pure_MeOH_sto[1] - C_rate_MeOH_2 + balance <= 0, name="m_MeOH_dis")


    # Sub objective

    c_x = (
    p_flow * (data.Wind[h] * P_cap[0] + data.Solar[h] * P_cap[1])
    - max_grid * data.grid_max                        # Grid
    + bat_min * P_cap[5] * data.Battery[4]
    - bat_max * P_cap[5] * data.Battery[5]
    - C_bat_1 * P_cap[5] * data.Battery[6]
    - C_bat_2 * P_cap[5] * data.Battery[6]            # Battery
    - elec_max * P_cap[2]
    + elec_min * P_cap[2] * data.op_elec[2]
    - elec_ramp_up * P_cap[2] * data.op_elec[0]
    - elec_ramp_down * P_cap[2] * data.op_elec[1]           # Electrolyzer
    + balance * data.Demand[h]                                     # Demand
    + h2_min * P_cap[6] * data.H2storage[4]
    - h2_max * P_cap[6] * data.H2storage[5]
    - C_rate_h2_1 * P_cap[6] * data.H2storage[6]
    - C_rate_h2_2 * P_cap[6] * data.H2storage[6]                     # H2 Storage
    + CO2_min * P_cap[8] * data.CO2storage[4]
    - CO2_max * P_cap[8] * data.CO2storage[5]
    + CO2_str_min * data.CO2storage[6]
    - CO2_str_max * data.CO2storage[7]
    - reac_max * P_cap[3]
    + reac_min * P_cap[3] * data.op_reac[2]
    - Reac_ramp_up * P_cap[3] * data.op_reac[0]
    - Reac_ramp_down * P_cap[3] * data.op_reac[1]                   # Reactor
    + raw_MeOH_min * P_cap[7] * data.raw_MeOH_sto[4]
    - raw_MeOH_max * P_cap[7] * data.raw_MeOH_sto[5]
    - C_raw_MeOH_1 * P_cap[7] * data.raw_MeOH_sto[6]
    - C_raw_MeOH_2 * P_cap[7] * data.raw_MeOH_sto[6]
    - dest_max * P_cap[4]
    + dest_min * P_cap[4] * data.op_dest[2]
    - dest_ramp_up * P_cap[4] * data.op_dest[0]
    - dest_ramp_down * P_cap[4] * data.op_dest[1]         # Distillator
    + MeOH_min * P_cap[9] * data.Pure_MeOH_sto[4]
    - MeOH_max * P_cap[9] * data.Pure_MeOH_sto[5]
    - C_rate_MeOH_1 * P_cap[9] * data.Pure_MeOH_sto[6]
    - C_rate_MeOH_2 * P_cap[9] * data.Pure_MeOH_sto[6]
    )                                                   # Pure Storage
        
    pi_A0_x = (
    # Battery SOC terms
    my_pi.e_bat_SOC[h-1] * (- bat_sto + bat_min - bat_max) + my_pi.e_bat_SOC[h] * bat_sto +

    # Electricity flow terms
    my_pi.e_elec[h-1] * (p_flow - elec_max + elec_min - elec_op + elec_ramp_up - elec_ramp_down) + 
    my_pi.e_elec[h] * (- elec_ramp_up + elec_ramp_down ) +

    # Hydrogen storage SOC terms
    my_pi.m_h2_sto_SOC[h-1] * (- h2_storage + h2_min - h2_max) + my_pi.m_h2_sto_SOC[h] * h2_storage + 

    # CO2 storage SOC terms
    my_pi.m_CO2_sto_SOC[h-1] * (- CO2_sto + CO2_min - CO2_max) + my_pi.m_CO2_sto_SOC[h] * CO2_sto +

    # Methanol reactor terms
    my_pi.m_rawMeOH_reac[h-1] * (-reac_max + reac_min + reac_H2 * data.Reactor[0] + reac_CO2 * data.Reactor[1] + 
    reac_el * data.Reactor[2] + reac_heat * data.Reactor[3] + Reac_ramp_up - Reac_ramp_down + compr_MeOH) +
    my_pi.m_rawMeOH_reac[h] * ( - Reac_ramp_up + Reac_ramp_down) +

    # Raw methanol storage SOC terms
    my_pi.m_rawMeOH_sto_SOC[h-1] * (rawMeOH_sto + raw_MeOH_min - raw_MeOH_max) - my_pi.m_rawMeOH_sto_SOC[h] * rawMeOH_sto +

    # Methanol distillation terms
    my_pi.m_MeOH_dest[h-1] * (-dest_max + dest_min + dest_rawMeOH * data.Destilator[0] + dest_water * data.Destilator[1] +
    dest_heat * data.Destilator[3] + dest_sto + dest_ramp_up - dest_ramp_down) + my_pi.m_MeOH_dest[h] * (- dest_ramp_up + dest_ramp_down) +

    # Pure methanol storage SOC terms
    my_pi.m_MeOH_SOC[h-1] * (- pure_sto + MeOH_min - MeOH_max) +  my_pi.m_MeOH_SOC[h] * pure_sto )

# Set the objective in the model
    sub_DW.setObjective(c_x - pi_A0_x - my_kappa, GRB.MAXIMIZE)


    sub_DW.setParam("OutputFlag", 0)

    sub_DW.optimize()


    if sub_DW.Status == GRB.OPTIMAL:

        sub_val = sub_DW.objVal

        results.p_flow[h]=p_flow.X
        results.max_grid[h]=max_grid.X

        results.elec_max[h]=elec_max.X
        results.elec_min[h]=elec_min.X
        results.elec_op[h]=elec_op.X
        results.elec_H2O[h]=elec_H2O.X
        results.elec_ramp_up[h]=elec_ramp_up.X
        results.elec_ramp_down[h]=elec_ramp_down.X

        results.bat_sto[h]=bat_sto.X
        results.bat_min[h]=bat_min.X
        results.bat_max[h]=bat_max.X
        results.C_bat_1[h]=C_bat_1.X
        results.C_bat_2[h]=C_bat_2.X
    
        results.des_H2O[h]=des_H2O.X

        results.h2_storage[h]=h2_storage.X
        results.h2_min[h]=h2_min.X
        results.h2_max[h]=h2_max.X
        results.C_rate_h2_1[h]=C_rate_h2_1.X
        results.C_rate_h2_2[h]=C_rate_h2_2.X

        results.compressor_h2[h]=compressor_h2.X

        results.CO2_sto[h]=CO2_sto.X
        results.CO2_min[h]=CO2_min.X
        results. CO2_max[h]=CO2_max.X
        results.CO2_str_min[h]=CO2_str_min.X
        results.CO2_str_max[h]=CO2_str_max.X
        results.CO2_comp[h]=CO2_comp.X

        results.compressor_e[h]=compressor_e.X

        results.reac_max[h]=reac_max.X
        results.reac_min[h]=reac_min.X
        results.reac_H2[h]=reac_H2.X
        results.reac_CO2[h]=reac_CO2.X
        results.reac_el[h]=reac_el.X
        results.reac_heat[h]=reac_heat.X
        results.Reac_ramp_up[h]=Reac_ramp_up.X
        results.Reac_ramp_down[h]=Reac_ramp_down.X

        results.compr_MeOH[h]=compr_MeOH.X

        results.rawMeOH_sto[h]=rawMeOH_sto.X
        results. raw_MeOH_min[h]=raw_MeOH_min.X
        results.raw_MeOH_max[h]=raw_MeOH_max.X
        results.C_raw_MeOH_1[h]=C_raw_MeOH_1.X
        results.C_raw_MeOH_2[h]=C_raw_MeOH_2.X

        results.dest_max[h]=dest_max.X
        results.dest_min[h]=dest_min.X
        results.dest_rawMeOH[h]=dest_rawMeOH.X
        results.dest_water[h]=dest_water.X
        results.dest_heat[h]=dest_heat.X
        results.dest_sto[h]=dest_sto.X
        results.dest_ramp_up[h]=dest_ramp_up.X
        results.dest_ramp_down[h]=dest_ramp_down.X

        results.pure_sto[h]=pure_sto.X
        results.MeOH_min[h]=MeOH_min.X
        results.MeOH_max[h]=MeOH_max.X
        results.C_rate_MeOH_1[h]=C_rate_MeOH_1.X
        results.C_rate_MeOH_2[h]=C_rate_MeOH_2.X

        results.balance[h]=balance.X

    print("IT: ", h)

    return sub_val, results



def DW_sub_hour_0(data,P_cap,my_pi,my_kappa,results):

#   Sub model NOTE (setupSUB)
    sub_DW = gp.Model("DW_sub")

# VARIABLES
    # Power flow
    p_flow = sub_DW.addMVar(1, lb=-10000, ub=10000, name="p_flow")

    # Max grid
    max_grid = sub_DW.addMVar(1, lb=0, ub=10000, name="max_grid")

    # Electrolyzer operation
    elec_max = sub_DW.addMVar(1, lb=0, ub=10000, name="Elec_max")
    elec_min = sub_DW.addMVar(1, lb=0, ub=10000, name="Elec_min")
    elec_op = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Elec_op")
    elec_H2O = sub_DW.addMVar(1, lb=-10000, ub=10000, name="elec_H2O")
    elec_ramp_up = sub_DW.addMVar(1, lb=0, ub=10000, name="elec_ramp_up")
    elec_ramp_down = sub_DW.addMVar(1, lb=0, ub=10000, name="elec_ramp_down")

    # Battery
    init_bat = sub_DW.addMVar(1, lb=-1000, ub=100, name="init_bat")
    # fin_bat = sub_DW.addMVar(1, lb=0, ub=1000, name="fin_bat")
    bat_sto = sub_DW.addMVar(1, lb=-1000, ub=100, name="bat_sto")
    bat_min = sub_DW.addMVar(1, lb=0, ub=1000, name="bat_min")
    bat_max = sub_DW.addMVar(1, lb=0, ub=1000, name="bat_max")
    C_bat_1 = sub_DW.addMVar(1, lb=0, ub=1000, name="C_bat_1")
    C_bat_2 = sub_DW.addMVar(1, lb=0, ub=1000, name="C_bat_2")

    # Desalinator
    des_H2O = sub_DW.addMVar(1, lb=-1000, ub=1000, name="des_H2O")

    # Hydrogen storage
    init_h2 = sub_DW.addMVar(1, lb=-10000, ub=10000, name="initialize_h2")
    # fin_h2 = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_h2")
    h2_storage = sub_DW.addMVar(1, lb=-10000, ub=10000, name="h2_storage")
    h2_min = sub_DW.addMVar(1, lb=0, ub=10000, name="h2_min")
    h2_max = sub_DW.addMVar(1, lb=0, ub=10000, name="h2_max")
    C_rate_h2_1 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_h2_1")
    C_rate_h2_2 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_h2_2")

    # H2 Compressor
    compressor_h2 = sub_DW.addMVar(1, lb=-10000, ub=10000, name="compressor_h2")

    # CO2 Storage
    init_CO2 = sub_DW.addMVar(1, lb=-10000, ub=10000, name="initialize_CO2")
    # fin_CO2 = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_CO2")
    CO2_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="CO2_storage")
    CO2_min = sub_DW.addMVar(1, lb=0, ub=10000, name="CO2_min")
    CO2_max = sub_DW.addMVar(1, lb=0, ub=10000, name="CO2_max")
    CO2_str_min = sub_DW.addMVar(1, lb=0, ub=10000, name="CO2_stream_min")
    CO2_str_max = sub_DW.addMVar(1, lb=0, ub=10000, name="CO2_stream_max")
    CO2_comp = sub_DW.addMVar(1, lb=-10000, ub=10000, name="CO2_comp")

    # Feed syngas Compressor
    compressor_e = sub_DW.addMVar(1, lb=-10000, ub=10000, name="compressor_e")

    # Methanol reactor operation
    reac_max = sub_DW.addMVar(1, lb=0, ub=10000, name="Reac_max")
    reac_min = sub_DW.addMVar(1, lb=0, ub=10000, name="Reac_min")
    reac_H2 =  sub_DW.addMVar(1, lb=-10000, ub=10000, name="Reac_H2")
    reac_CO2 = sub_DW.addMVar(1, lb=-10000, ub=10000, name="reac_CO2")
    reac_el = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Reac_el")
    reac_heat = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Reac_heat") 
    Reac_ramp_up = sub_DW.addMVar(1, lb=0, ub=10000, name="Reac_ramp_up")
    Reac_ramp_down = sub_DW.addMVar(1, lb=0, ub=10000, name="Reac_ramp_down")

    # Methanol compressor
    compr_MeOH = sub_DW.addMVar(1, lb=-10000, ub=10000, name="compr_MeOH")

    # Raw Methanol storage
    init_MeOH_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="initialize_MeOH_sto")
    # fin_MeOH_sto = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_MeOH_sto")
    rawMeOH_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="rawMeOH_storage")
    raw_MeOH_min = sub_DW.addMVar(1, lb=0, ub=10000, name="raw_MeOH_min")
    raw_MeOH_max = sub_DW.addMVar(1, lb=0, ub=10000, name="raw_MeOH_max")
    C_raw_MeOH_1 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_raw_MeOH_1")
    C_raw_MeOH_2 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_raw_MeOH_2")

    # Methanol distillator
    dest_max = sub_DW.addMVar(1, lb=0, ub=10000, name="Dest_max")
    dest_min = sub_DW.addMVar(1, lb=0, ub=10000, name="Dest_min")
    dest_rawMeOH = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Dest_rawMeOH")
    dest_water = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Dest_water")
    dest_heat = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Dest_heat")
    dest_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Dest_sto")
    dest_ramp_up = sub_DW.addMVar(1, lb=0, ub=10000, name="Dest_ramp_up")
    dest_ramp_down = sub_DW.addMVar(1, lb=0, ub=10000, name="Dest_ramp_down")

    # Pure Methanol storage
    init_pure = sub_DW.addMVar(1, lb=-10000, ub=10000, name="initialize_pure")
    # fin_pure = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_pure")
    pure_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="pure_storage")
    MeOH_min = sub_DW.addMVar(1, lb=0, ub=10000, name="MeOH_min")
    MeOH_max = sub_DW.addMVar(1, lb=0, ub=10000, name="MeOH_max")
    C_rate_MeOH_1 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_MeOH_1")
    C_rate_MeOH_2 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_MeOH_2")

    # Meet demand
    balance = sub_DW.addMVar(1, lb=0, name="balance")

#   SUB CONSTRAINTS
    # Energy
    e_buy = sub_DW.addConstr(-p_flow - max_grid <= data.Prices[0], name='e_buy')
    e_sell = sub_DW.addConstr(p_flow <= 0, name='e_sell')

    # Battery
    e_bat_cha = sub_DW.addConstr(- bat_sto*data.Battery[0] + p_flow - C_bat_1 <= 0)
    e_bat_dis = sub_DW.addConstr(+ bat_sto/data.Battery[1] - p_flow - C_bat_2 <= 0)

    # Desalinator
    e_des = sub_DW.addConstr(des_H2O + p_flow <= 0, name='e_des')

    # Electrolyzer
    m_h2_elec = sub_DW.addConstr(elec_op*data.Electrolyzer[0] + elec_H2O*data.Electrolyzer[1] + compressor_h2 <= 0, name='m_h2_elec')
    m_H2O_elec = sub_DW.addConstr(- elec_H2O - des_H2O*data.Desalinator[0] <= 0, name='m_H2O_elec')

    # H2 Storage
    m_h2_sto_cha = sub_DW.addConstr(- h2_storage*data.H2storage[0] - C_rate_h2_1 - compressor_h2 <= 0, name='m_h2_sto_cha')
    m_h2_sto_dis = sub_DW.addConstr(h2_storage/data.H2storage[1] - C_rate_h2_2 - reac_H2 - compressor_e*data.Compressor[0] <= 0, name='m_h2_sto_dis')

    # CO2 Storage
    m_CO2_buy = sub_DW.addConstr(CO2_str_min - CO2_str_max + CO2_comp <= data.CO2_prices[0], name='m_CO2_buy')
    m_CO2_sto_cha = sub_DW.addConstr(- CO2_sto*data.CO2storage[0] - CO2_comp <= 0, name='m_CO2_sto_cha')
    m_CO2_sto_dis = sub_DW.addConstr(CO2_sto/data.CO2storage[1] - compressor_e*data.Compressor[0] - reac_CO2 <= 0, name='m_CO2_sto_dis')

    # Syngas compressor
    e_compr = sub_DW.addConstr(p_flow + compressor_e <= 0, name='e_compr')

    # Methanol reactor
    m_h2_reac = sub_DW.addConstr(- compressor_h2 - reac_H2 - compressor_e*data.Compressor[0] <= 0, name='m_h2_reac')
    m_CO2_reac = sub_DW.addConstr(- CO2_comp - reac_CO2 - compressor_e*data.Compressor[0] <= 0, name='m_CO2_reac') 
    e_reac = sub_DW.addConstr(p_flow - reac_el <= 0, name='e_reac')
    q_reac = sub_DW.addConstr(- reac_heat + p_flow*data.Heater[0] <= 0, name='q_reac')

    # MeOH Storage
    m_rawMeOH_sto_cha = sub_DW.addConstr(-rawMeOH_sto*data.raw_MeOH_sto[0] - C_raw_MeOH_1 - compr_MeOH <= 0, name='m_rawMeOH_sto_cha')
    m_rawMeOH_sto_dis = sub_DW.addConstr(+rawMeOH_sto/data.raw_MeOH_sto[1] - C_raw_MeOH_2 - dest_rawMeOH <= 0, name='m_rawMeOH_sto_dis')

    # Methanol distilator
    m_rawMeOH_dest = sub_DW.addConstr(- compr_MeOH - dest_rawMeOH <= 0, name='m_rawMeOH_dest')
    q_dest = sub_DW.addConstr(- dest_heat + p_flow*data.Heater[0] <= 0, name='q_dest')
    m_h2O_dest = sub_DW.addConstr(- dest_water <= 0, name='m_h2O_dest')

    # Pure MeOH Storagess
    m_MeOH_cha = sub_DW.addConstr(- dest_sto - C_rate_MeOH_1 - pure_sto*data.Pure_MeOH_sto[0] <= 0, name='m_MeOH_cha')
    m_MeOH_dis = sub_DW.addConstr(pure_sto/data.Pure_MeOH_sto[1] - C_rate_MeOH_2 + balance <= 0, name="m_MeOH_dis")

    # Sub objective

    c_x = (
    p_flow * (data.Wind[0] * P_cap[0] + data.Solar[0] * P_cap[1])
    - max_grid * data.grid_max                                # Grid
    + init_bat * P_cap[5] * data.Battery[2]
    + bat_min * P_cap[5] * data.Battery[4]
    - bat_max * P_cap[5] * data.Battery[5]
    - C_bat_1 * P_cap[5] * data.Battery[6]
    - C_bat_2 * P_cap[5] * data.Battery[6]                   # Battery
    - elec_max * P_cap[2]
    + elec_min * P_cap[2] * data.op_elec[2]
    - elec_ramp_up * P_cap[2] * data.op_elec[0]
    - elec_ramp_down * P_cap[2] * data.op_elec[1]               # Electrolyzer
    + balance * data.Demand[0]                                  # Demand
    + init_h2 * P_cap[6] * data.H2storage[2]
    + h2_min * P_cap[6] * data.H2storage[4]
    - h2_max * P_cap[6] * data.H2storage[5]
    - C_rate_h2_1 * P_cap[6] * data.H2storage[6]
    - C_rate_h2_2 * P_cap[6] * data.H2storage[6]                 # H2 Storage
    + init_CO2 * P_cap[8] * data.CO2storage[2]
    + CO2_min * P_cap[8] * data.CO2storage[4]
    - CO2_max * P_cap[8] * data.CO2storage[5]
    + CO2_str_min * data.CO2storage[6]
    - CO2_str_max * data.CO2storage[7]                    # CO2 Storage
    - reac_max * P_cap[3]
    + reac_min * P_cap[3] * data.op_reac[2]
    - Reac_ramp_up * P_cap[3] * data.op_reac[0]
    - Reac_ramp_down * P_cap[3] * data.op_reac[1]                   # Reactor
    + init_MeOH_sto * P_cap[7] * data.raw_MeOH_sto[2]
    + raw_MeOH_min * P_cap[7] * data.raw_MeOH_sto[4]
    - raw_MeOH_max * P_cap[7] * data.raw_MeOH_sto[5]
    - C_raw_MeOH_1 * P_cap[7] * data.raw_MeOH_sto[6]
    - C_raw_MeOH_2 * P_cap[7] * data.raw_MeOH_sto[6]          # Raw Storage
    - dest_max * P_cap[4]
    + dest_min * P_cap[4] * data.op_dest[2]
    - dest_ramp_up * P_cap[4] * data.op_dest[0]
    - dest_ramp_down * P_cap[4] * data.op_dest[1]             # Distillator
    + init_pure * P_cap[9] * data.Pure_MeOH_sto[2]
    + MeOH_min * P_cap[9] * data.Pure_MeOH_sto[4]
    - MeOH_max * P_cap[9] * data.Pure_MeOH_sto[5]
    - C_rate_MeOH_1 * P_cap[9] * data.Pure_MeOH_sto[6]
    - C_rate_MeOH_2 * P_cap[9] * data.Pure_MeOH_sto[6]       # Pure Storage
    )                                                           
        
    pi_A0_x = (
    # Battery SOC terms
    my_pi.e_bat_SOC_0 * (init_bat - bat_sto[0] + bat_min[0] - bat_max[0]) + my_pi.e_bat_SOC[0] * bat_sto[0] +

    # Electricity flow terms
    my_pi.e_elec_0 * (p_flow[0] - elec_max[0] + elec_min[0] - elec_op[0] + elec_ramp_up[0] - elec_ramp_down[0]) +
    my_pi.e_elec[0] * ( - elec_ramp_up[0] + elec_ramp_down[0] ) +

    # Hydrogen storage SOC terms
    my_pi.m_h2_sto_SOC_0 * (init_h2 - h2_storage[0] + h2_min[0] - h2_max[0]) + my_pi.m_h2_sto_SOC[0] * h2_storage[0] +

    # CO2 storage SOC terms
    my_pi.m_CO2_sto_SOC_0 * (init_CO2 - CO2_sto[0] + CO2_min[0] - CO2_max[0]) + my_pi.m_CO2_sto_SOC[0] * (CO2_sto[0] ) +

    # Methanol reactor terms
    my_pi.m_rawMeOH_reac_0 * (-reac_max[0] + reac_min[0] + reac_H2[0] * data.Reactor[0] + reac_CO2[0] * data.Reactor[1] + 
    reac_el[0] * data.Reactor[2] + reac_heat[0] * data.Reactor[3] + Reac_ramp_up[0] - Reac_ramp_down[0] + compr_MeOH[0]) + 
    my_pi.m_rawMeOH_reac[0] * (- Reac_ramp_up[0] + Reac_ramp_down[0]) +

    # Raw methanol storage SOC terms
    my_pi.m_rawMeOH_sto_SOC_0 * (init_MeOH_sto - rawMeOH_sto[0] + raw_MeOH_min[0] - raw_MeOH_max[0]) +
    my_pi.m_rawMeOH_sto_SOC[0] * (- rawMeOH_sto[0]) +

    # Methanol distillation terms
    my_pi.m_MeOH_dest_0 * (-dest_max[0] + dest_min[0] + dest_rawMeOH[0] * data.Destilator[0] + dest_water[0] * data.Destilator[1] + 
    dest_heat[0] * data.Destilator[3] + dest_sto[0] + dest_ramp_up[0] - dest_ramp_down[0]) +
    my_pi.m_MeOH_dest[0] * (- dest_ramp_up[0] + dest_ramp_down[0] ) +

    # Pure methanol storage SOC terms
    my_pi.m_MeOH_SOC_0 * (init_pure - pure_sto[0] + MeOH_min[0] - MeOH_max[0]) +
    my_pi.m_MeOH_SOC[0] * (pure_sto[0])
    )

# Set the objective in the model
    sub_DW.setObjective(c_x - pi_A0_x - my_kappa, GRB.MAXIMIZE)

    sub_DW.update()

    sub_DW.write('wolfe.lp')

    sub_DW.optimize()


    if sub_DW.Status == GRB.OPTIMAL:

        results.p_flow[0]=p_flow.X
        results.max_grid[0]=max_grid.X

        results.elec_max[0]=elec_max.X
        results.elec_min[0]=elec_min.X
        results.elec_op[0]=elec_op.X
        results.elec_H2O[0]=elec_H2O.X
        results.elec_ramp_up[0]=elec_ramp_up.X
        results.elec_ramp_down[0]=elec_ramp_down.X

        results.init_bat=init_bat.X,
        results.bat_sto[0]=bat_sto.X
        results.bat_min[0]=bat_min.X
        results.bat_max[0]=bat_max.X
        results.C_bat_1[0]=C_bat_1.X
        results.C_bat_2[0]=C_bat_2.X
            
        results.des_H2O[0]=des_H2O.X

        results.init_h2=init_h2.X,
        results.h2_storage[0]=h2_storage.X
        results.h2_min[0]=h2_min.X
        results.h2_max[0]=h2_max.X
        results.C_rate_h2_1[0]=C_rate_h2_1.X
        results.C_rate_h2_2[0]=C_rate_h2_2.X

        results.compressor_h2[0]=compressor_h2.X

        results.init_CO2=init_CO2.X,
        results.CO2_sto[0]=CO2_sto.X
        results.CO2_min[0]=CO2_min.X
        results. CO2_max[0]=CO2_max.X
        results.CO2_str_min[0]=CO2_str_min.X
        results.CO2_str_max[0]=CO2_str_max.X
        results.CO2_comp[0]=CO2_comp.X

        results.compressor_e[0]=compressor_e.X

        results.reac_max[0]=reac_max.X
        results.reac_min[0]=reac_min.X
        results.reac_H2[0]=reac_H2.X
        results.reac_CO2[0]=reac_CO2.X
        results.reac_el[0]=reac_el.X
        results.reac_heat[0]=reac_heat.X
        results.Reac_ramp_up[0]=Reac_ramp_up.X
        results.Reac_ramp_down[0]=Reac_ramp_down.X

        results.compr_MeOH[0]=compr_MeOH.X

        results.init_MeOH_sto=init_MeOH_sto.X,
        results.rawMeOH_sto[0]=rawMeOH_sto.X
        results. raw_MeOH_min[0]=raw_MeOH_min.X
        results.raw_MeOH_max[0]=raw_MeOH_max.X
        results.C_raw_MeOH_1[0]=C_raw_MeOH_1.X
        results.C_raw_MeOH_2[0]=C_raw_MeOH_2.X

        results.dest_max[0]=dest_max.X
        results.dest_min[0]=dest_min.X
        results.dest_rawMeOH[0]=dest_rawMeOH.X
        results.dest_water[0]=dest_water.X
        results.dest_heat[0]=dest_heat.X
        results.dest_sto[0]=dest_sto.X
        results.dest_ramp_up[0]=dest_ramp_up.X
        results.dest_ramp_down[0]=dest_ramp_down.X

        results.init_pure=init_pure.X,
        results.pure_sto[0]=pure_sto.X
        results.MeOH_min[0]=MeOH_min.X
        results.MeOH_max[0]=MeOH_max.X
        results.C_rate_MeOH_1[0]=C_rate_MeOH_1.X
        results.C_rate_MeOH_2[0]=C_rate_MeOH_2.X

        results.balance[0]=balance.X


    return sub_DW.objVal, results


def DW_sub_hour_h(data,P_cap,my_pi,my_kappa,results):

#   Sub model NOTE (setupSUB)
    sub_DW = gp.Model("DW_sub")

# VARIABLES
    # Power flow
    p_flow = sub_DW.addMVar(1, lb=-10000, ub=10000, name="p_flow")

    # Max grid
    max_grid = sub_DW.addMVar(1, lb=0, ub=10000, name="max_grid")

    # Electrolyzer operation
    elec_max = sub_DW.addMVar(1, lb=0, ub=10000, name="Elec_max")
    elec_min = sub_DW.addMVar(1, lb=0, ub=10000, name="Elec_min")
    elec_op = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Elec_op")
    elec_H2O = sub_DW.addMVar(1, lb=-10000, ub=10000, name="elec_H2O")
    elec_ramp_up = sub_DW.addMVar(1, lb=0, ub=10000, name="elec_ramp_up")
    elec_ramp_down = sub_DW.addMVar(1, lb=0, ub=10000, name="elec_ramp_down")

    # Battery
    fin_bat = sub_DW.addMVar(1, lb=0, ub=1000, name="fin_bat")
    bat_sto = sub_DW.addMVar(1, lb=-1000, ub=100, name="bat_sto")
    bat_min = sub_DW.addMVar(1, lb=0, ub=1000, name="bat_min")
    bat_max = sub_DW.addMVar(1, lb=0, ub=1000, name="bat_max")
    C_bat_1 = sub_DW.addMVar(1, lb=0, ub=1000, name="C_bat_1")
    C_bat_2 = sub_DW.addMVar(1, lb=0, ub=1000, name="C_bat_2")

    # Desalinator
    des_H2O = sub_DW.addMVar(1, lb=-1000, ub=1000, name="des_H2O")

    # Hydrogen storage
    fin_h2 = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_h2")
    h2_storage = sub_DW.addMVar(1, lb=-10000, ub=10000, name="h2_storage")
    h2_min = sub_DW.addMVar(1, lb=0, ub=10000, name="h2_min")
    h2_max = sub_DW.addMVar(1, lb=0, ub=10000, name="h2_max")
    C_rate_h2_1 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_h2_1")
    C_rate_h2_2 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_h2_2")

    # H2 Compressor
    compressor_h2 = sub_DW.addMVar(1, lb=-10000, ub=10000, name="compressor_h2")

    # CO2 Storage
    fin_CO2 = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_CO2")
    CO2_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="CO2_storage")
    CO2_min = sub_DW.addMVar(1, lb=0, ub=10000, name="CO2_min")
    CO2_max = sub_DW.addMVar(1, lb=0, ub=10000, name="CO2_max")
    CO2_str_min = sub_DW.addMVar(1, lb=0, ub=10000, name="CO2_stream_min")
    CO2_str_max = sub_DW.addMVar(1, lb=0, ub=10000, name="CO2_stream_max")
    CO2_comp = sub_DW.addMVar(1, lb=-10000, ub=10000, name="CO2_comp")

    # Feed syngas Compressor
    compressor_e = sub_DW.addMVar(1, lb=-10000, ub=10000, name="compressor_e")

    # Methanol reactor operation
    reac_max = sub_DW.addMVar(1, lb=0, ub=10000, name="Reac_max")
    reac_min = sub_DW.addMVar(1, lb=0, ub=10000, name="Reac_min")
    reac_H2 =  sub_DW.addMVar(1, lb=-10000, ub=10000, name="Reac_H2")
    reac_CO2 = sub_DW.addMVar(1, lb=-10000, ub=10000, name="reac_CO2")
    reac_el = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Reac_el")
    reac_heat = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Reac_heat") 
    Reac_ramp_up = sub_DW.addMVar(1, lb=0, ub=10000, name="Reac_ramp_up")
    Reac_ramp_down = sub_DW.addMVar(1, lb=0, ub=10000, name="Reac_ramp_down")

    # Methanol compressor
    compr_MeOH = sub_DW.addMVar(1, lb=-10000, ub=10000, name="compr_MeOH")

    # Raw Methanol storage
    fin_MeOH_sto = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_MeOH_sto")
    rawMeOH_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="rawMeOH_storage")
    raw_MeOH_min = sub_DW.addMVar(1, lb=0, ub=10000, name="raw_MeOH_min")
    raw_MeOH_max = sub_DW.addMVar(1, lb=0, ub=10000, name="raw_MeOH_max")
    C_raw_MeOH_1 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_raw_MeOH_1")
    C_raw_MeOH_2 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_raw_MeOH_2")

    # Methanol distillator
    dest_max = sub_DW.addMVar(1, lb=0, ub=10000, name="Dest_max")
    dest_min = sub_DW.addMVar(1, lb=0, ub=10000, name="Dest_min")
    dest_rawMeOH = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Dest_rawMeOH")
    dest_water = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Dest_water")
    dest_heat = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Dest_heat")
    dest_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Dest_sto")
    dest_ramp_up = sub_DW.addMVar(1, lb=0, ub=10000, name="Dest_ramp_up")
    dest_ramp_down = sub_DW.addMVar(1, lb=0, ub=10000, name="Dest_ramp_down")

    # Pure Methanol storage
    fin_pure = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_pure")
    pure_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="pure_storage")
    MeOH_min = sub_DW.addMVar(1, lb=0, ub=10000, name="MeOH_min")
    MeOH_max = sub_DW.addMVar(1, lb=0, ub=10000, name="MeOH_max")
    C_rate_MeOH_1 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_MeOH_1")
    C_rate_MeOH_2 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_MeOH_2")

    # Meet demand
    balance = sub_DW.addMVar(1, lb=0, name="balance")

#   SUB CONSTRAINTS
    # Energy
    e_buy = sub_DW.addConstr(-p_flow - max_grid <= data.Prices[data.hours-1], name='e_buy')
    e_sell = sub_DW.addConstr(p_flow <= 0, name='e_sell')

    # Battery
    e_bat_cha_h = sub_DW.addConstr(+ p_flow - C_bat_1 <= 0)
    e_bat_dis_h = sub_DW.addConstr(- p_flow - C_bat_2 <= 0)

    # Desalinator
    e_des = sub_DW.addConstr(des_H2O + p_flow <= 0, name='e_des')

    # Electrolyzer
    m_h2_elec = sub_DW.addConstr(elec_op*data.Electrolyzer[0] + elec_H2O*data.Electrolyzer[1] + compressor_h2 <= 0, name='m_h2_elec')
    m_H2O_elec = sub_DW.addConstr(- elec_H2O - des_H2O*data.Desalinator[0] <= 0, name='m_H2O_elec')

    # H2 Storage
    m_h2_sto_cha_h = sub_DW.addConstr(- C_rate_h2_1 - compressor_h2 <= 0, name='m_h2_sto_cha_h')
    m_h2_sto_dis_h = sub_DW.addConstr(- C_rate_h2_2 - reac_H2 - compressor_e*data.Compressor[0] <= 0, name='m_h2_sto_dis_h')

    # CO2 Storage
    m_CO2_buy = sub_DW.addConstr(CO2_str_min - CO2_str_max + CO2_comp <= data.CO2_prices[data.hours-1], name='m_CO2_buy')
    m_CO2_sto_cha_h = sub_DW.addConstr(- CO2_comp <= 0, name='m_CO2_sto_cha_h')
    m_CO2_sto_dis_h = sub_DW.addConstr(- compressor_e*data.Compressor[0] - reac_CO2 <= 0, name='m_CO2_sto_dis_h')

    # Syngas compressor
    e_compr = sub_DW.addConstr(p_flow + compressor_e <= 0, name='e_compr')

    # Methanol reactor
    m_h2_reac = sub_DW.addConstr(- compressor_h2 - reac_H2 - compressor_e*data.Compressor[0] <= 0, name='m_h2_reac')
    m_CO2_reac = sub_DW.addConstr(- CO2_comp - reac_CO2 - compressor_e*data.Compressor[0] <= 0, name='m_CO2_reac') 
    e_reac = sub_DW.addConstr(p_flow - reac_el <= 0, name='e_reac')
    q_reac = sub_DW.addConstr(- reac_heat + p_flow*data.Heater[0] <= 0, name='q_reac')

    # MeOH Storage
    m_rawMeOH_sto_cha_h = sub_DW.addConstr( - C_raw_MeOH_1 - compr_MeOH <= 0, name='m_rawMeOH_sto_cha')
    m_rawMeOH_sto_dis_h = sub_DW.addConstr( - C_raw_MeOH_2 - dest_rawMeOH <= 0, name='m_rawMeOH_sto_dis')

    # Methanol distilator
    m_rawMeOH_dest = sub_DW.addConstr(- compr_MeOH - dest_rawMeOH <= 0, name='m_rawMeOH_dest')
    q_dest = sub_DW.addConstr(- dest_heat + p_flow*data.Heater[0] <= 0, name='q_dest')
    m_h2O_dest = sub_DW.addConstr(- dest_water <= 0, name='m_h2O_dest')

    # Pure MeOH Storagess
    m_MeOH_cha_h = sub_DW.addConstr(- dest_sto - C_rate_MeOH_1 <= 0, name='m_MeOH_cha_h')
    m_MeOH_dis_h = sub_DW.addConstr(- C_rate_MeOH_2 + balance <= 0, name="m_MeOH_dis_h")


    # Sub objective

    c_x = (
    p_flow * (data.Wind[data.hours-1] * P_cap[0] + data.Solar[data.hours-1] * P_cap[1])
    - max_grid * data.grid_max                       # Grid
    + fin_bat * P_cap[5] * data.Battery[3]
    + bat_min * P_cap[5] * data.Battery[4]
    - bat_max * P_cap[5] * data.Battery[5]
    - C_bat_1 * P_cap[5] * data.Battery[6]
    - C_bat_2 * P_cap[5] * data.Battery[6]           # Battery
    - elec_max * P_cap[2] + elec_min * P_cap[2] * data.op_elec[2]     # Electrolyzer
    + balance * data.Demand[data.hours-1]                                   # Demand
    + fin_h2 * P_cap[6] * data.H2storage[3]
    + h2_min * P_cap[6] * data.H2storage[4]
    - h2_max * P_cap[6] * data.H2storage[5]
    - C_rate_h2_1 * P_cap[6] * data.H2storage[6]
    - C_rate_h2_2 * P_cap[6] * data.H2storage[6]    # H2 Storage
    + fin_CO2 * P_cap[8] * data.CO2storage[3]
    + CO2_min * P_cap[8] * data.CO2storage[4]
    - CO2_max * P_cap[8] * data.CO2storage[5]
    + CO2_str_min * data.CO2storage[6]
    - CO2_str_max * data.CO2storage[7]          # CO2 Storage
    - reac_max * P_cap[3]
    + reac_min * P_cap[3] * data.op_reac[2]       # Reactor
    + fin_MeOH_sto * P_cap[7] * data.raw_MeOH_sto[3]
    + raw_MeOH_min * P_cap[7] * data.raw_MeOH_sto[4]
    - raw_MeOH_max * P_cap[7] * data.raw_MeOH_sto[5]
    - C_raw_MeOH_1 * P_cap[7] * data.raw_MeOH_sto[6]
    - C_raw_MeOH_2 * P_cap[7] * data.raw_MeOH_sto[6]       # Raw Storage
    - dest_max * P_cap[4]
    + dest_min * P_cap[4] * data.op_dest[2]       # Distillator
    + fin_pure * P_cap[9] * data.Pure_MeOH_sto[3]
    + MeOH_min * P_cap[9] * data.Pure_MeOH_sto[4]
    - MeOH_max * P_cap[9] * data.Pure_MeOH_sto[5]
    - C_rate_MeOH_1 * P_cap[9] * data.Pure_MeOH_sto[6]
    - C_rate_MeOH_2 * P_cap[9] * data.Pure_MeOH_sto[6]
    )                                                           # Pure Storage
        
    pi_A0_x = (

    # Battery SOC terms
    my_pi.e_bat_SOC_h * (fin_bat + bat_min - bat_max) +

    # Electricity flow terms
    my_pi.e_elec_h * (p_flow - elec_max + elec_min - elec_op) +

    # Hydrogen storage SOC terms
    my_pi.m_h2_sto_SOC_h * (fin_h2 + h2_min - h2_max) +

    # CO2 storage SOC terms
    my_pi.m_CO2_sto_SOC_h * (fin_CO2 + CO2_min - CO2_max) +

    # Methanol reactor terms
    my_pi.m_rawMeOH_reac_h * (-reac_max + reac_min + reac_H2 * data.Reactor[0] + reac_CO2 * data.Reactor[1] + reac_el * data.Reactor[2] + reac_heat * data.Reactor[3] + compr_MeOH) +

    # Raw methanol storage SOC terms
    my_pi.m_rawMeOH_sto_SOC_h * (fin_MeOH_sto + raw_MeOH_min - raw_MeOH_max) +

    # Methanol distillation terms
    my_pi.m_MeOH_dest_h * (-dest_max + dest_min + dest_rawMeOH * data.Destilator[0] + dest_water * data.Destilator[1] + dest_heat * data.Destilator[3] + dest_sto) +

    # Pure methanol storage SOC terms
    my_pi.m_MeOH_SOC_h * (fin_pure + MeOH_min - MeOH_max)
    )

# Set the objective in the model
    sub_DW.setObjective(c_x - pi_A0_x - my_kappa, GRB.MAXIMIZE)

    sub_DW.update()

    sub_DW.write('wolfe.lp')

    sub_DW.optimize()


    if sub_DW.Status == GRB.OPTIMAL:

        # Print optimal objective value
        sub_val = sub_DW.objVal

        results.p_flow[data.hours-1]=p_flow.X
        results.max_grid[data.hours-1]=max_grid.X

        results.elec_max[data.hours-1]=elec_max.X
        results.elec_min[data.hours-1]=elec_min.X
        results.elec_op[data.hours-1]=elec_op.X
        results.elec_H2O[data.hours-1]=elec_H2O.X
        results.elec_ramp_up[data.hours-1]=elec_ramp_up.X
        results.elec_ramp_down[data.hours-1]=elec_ramp_down.X

        results.fin_bat=fin_bat.X
        results.bat_sto[data.hours-1]=bat_sto.X
        results.bat_min[data.hours-1]=bat_min.X
        results.bat_max[data.hours-1]=bat_max.X
        results.C_bat_1[data.hours-1]=C_bat_1.X
        results.C_bat_2[data.hours-1]=C_bat_2.X
            
        results.des_H2O[data.hours-1]=des_H2O.X

        results.fin_h2 = fin_h2.X
        results.h2_storage[data.hours-1]=h2_storage.X
        results.h2_min[data.hours-1]=h2_min.X
        results.h2_max[data.hours-1]=h2_max.X
        results.C_rate_h2_1[data.hours-1]=C_rate_h2_1.X
        results.C_rate_h2_2[data.hours-1]=C_rate_h2_2.X

        results.compressor_h2[data.hours-1]=compressor_h2.X

        results.fin_CO2 = fin_CO2.X
        results.CO2_sto[data.hours-1]=CO2_sto.X
        results.CO2_min[data.hours-1]=CO2_min.X
        results. CO2_max[data.hours-1]=CO2_max.X
        results.CO2_str_min[data.hours-1]=CO2_str_min.X
        results.CO2_str_max[data.hours-1]=CO2_str_max.X
        results.CO2_comp[data.hours-1]=CO2_comp.X

        results.compressor_e[data.hours-1]=compressor_e.X

        results.reac_max[data.hours-1]=reac_max.X
        results.reac_min[data.hours-1]=reac_min.X
        results.reac_H2[data.hours-1]=reac_H2.X
        results.reac_CO2[data.hours-1]=reac_CO2.X
        results.reac_el[data.hours-1]=reac_el.X
        results.reac_heat[data.hours-1]=reac_heat.X
        results.Reac_ramp_up[data.hours-1]=Reac_ramp_up.X
        results.Reac_ramp_down[data.hours-1]=Reac_ramp_down.X

        results.compr_MeOH[data.hours-1]=compr_MeOH.X

        results.fin_MeOH_sto = fin_MeOH_sto.X
        results.rawMeOH_sto[data.hours-1]=rawMeOH_sto.X
        results. raw_MeOH_min[data.hours-1]=raw_MeOH_min.X
        results.raw_MeOH_max[data.hours-1]=raw_MeOH_max.X
        results.C_raw_MeOH_1[data.hours-1]=C_raw_MeOH_1.X
        results.C_raw_MeOH_2[data.hours-1]=C_raw_MeOH_2.X

        results.dest_max[data.hours-1]=dest_max.X
        results.dest_min[data.hours-1]=dest_min.X
        results.dest_rawMeOH[data.hours-1]=dest_rawMeOH.X
        results.dest_water[data.hours-1]=dest_water.X
        results.dest_heat[data.hours-1]=dest_heat.X
        results.dest_sto[data.hours-1]=dest_sto.X
        results.dest_ramp_up[data.hours-1]=dest_ramp_up.X
        results.dest_ramp_down[data.hours-1]=dest_ramp_down.X

        results.fin_pure = fin_pure.X
        results.pure_sto[data.hours-1]=pure_sto.X
        results.MeOH_min[data.hours-1]=MeOH_min.X
        results.MeOH_max[data.hours-1]=MeOH_max.X
        results.C_rate_MeOH_1[data.hours-1]=C_rate_MeOH_1.X
        results.C_rate_MeOH_2[data.hours-1]=C_rate_MeOH_2.X

        results.balance[data.hours-1]=balance.X

    return sub_val, results


def DW_sub_hour_h_1(data,P_cap,my_pi,my_kappa,results):

#   Sub model NOTE (setupSUB)
    sub_DW = gp.Model("DW_sub")

# VARIABLES
    # Power flow
    p_flow = sub_DW.addMVar(1, lb=-10000, ub=10000, name="p_flow")

    # Max grid
    max_grid = sub_DW.addMVar(1, lb=0, ub=10000, name="max_grid")

    # Electrolyzer operation
    elec_max = sub_DW.addMVar(1, lb=0, ub=10000, name="Elec_max")
    elec_min = sub_DW.addMVar(1, lb=0, ub=10000, name="Elec_min")
    elec_op = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Elec_op")
    elec_H2O = sub_DW.addMVar(1, lb=-10000, ub=10000, name="elec_H2O")
    elec_ramp_up = sub_DW.addMVar(1, lb=0, ub=10000, name="elec_ramp_up")
    elec_ramp_down = sub_DW.addMVar(1, lb=0, ub=10000, name="elec_ramp_down")

    # Battery
    init_bat = sub_DW.addMVar(1, lb=-1000, ub=100, name="init_bat")
    # fin_bat = sub_DW.addMVar(1, lb=0, ub=1000, name="fin_bat")
    bat_sto = sub_DW.addMVar(1, lb=-1000, ub=100, name="bat_sto")
    bat_min = sub_DW.addMVar(1, lb=0, ub=1000, name="bat_min")
    bat_max = sub_DW.addMVar(1, lb=0, ub=1000, name="bat_max")
    C_bat_1 = sub_DW.addMVar(1, lb=0, ub=1000, name="C_bat_1")
    C_bat_2 = sub_DW.addMVar(1, lb=0, ub=1000, name="C_bat_2")

    # Desalinator
    des_H2O = sub_DW.addMVar(1, lb=-1000, ub=1000, name="des_H2O")

    # Hydrogen storage
    init_h2 = sub_DW.addMVar(1, lb=-10000, ub=10000, name="initialize_h2")
    # fin_h2 = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_h2")
    h2_storage = sub_DW.addMVar(1, lb=-10000, ub=10000, name="h2_storage")
    h2_min = sub_DW.addMVar(1, lb=0, ub=10000, name="h2_min")
    h2_max = sub_DW.addMVar(1, lb=0, ub=10000, name="h2_max")
    C_rate_h2_1 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_h2_1")
    C_rate_h2_2 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_h2_2")

    # H2 Compressor
    compressor_h2 = sub_DW.addMVar(1, lb=-10000, ub=10000, name="compressor_h2")

    # CO2 Storage
    init_CO2 = sub_DW.addMVar(1, lb=-10000, ub=10000, name="initialize_CO2")
    # fin_CO2 = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_CO2")
    CO2_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="CO2_storage")
    CO2_min = sub_DW.addMVar(1, lb=0, ub=10000, name="CO2_min")
    CO2_max = sub_DW.addMVar(1, lb=0, ub=10000, name="CO2_max")
    CO2_str_min = sub_DW.addMVar(1, lb=0, ub=10000, name="CO2_stream_min")
    CO2_str_max = sub_DW.addMVar(1, lb=0, ub=10000, name="CO2_stream_max")
    CO2_comp = sub_DW.addMVar(1, lb=-10000, ub=10000, name="CO2_comp")

    # Feed syngas Compressor
    compressor_e = sub_DW.addMVar(1, lb=-10000, ub=10000, name="compressor_e")

    # Methanol reactor operation
    reac_max = sub_DW.addMVar(1, lb=0, ub=10000, name="Reac_max")
    reac_min = sub_DW.addMVar(1, lb=0, ub=10000, name="Reac_min")
    reac_H2 =  sub_DW.addMVar(1, lb=-10000, ub=10000, name="Reac_H2")
    reac_CO2 = sub_DW.addMVar(1, lb=-10000, ub=10000, name="reac_CO2")
    reac_el = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Reac_el")
    reac_heat = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Reac_heat") 
    Reac_ramp_up = sub_DW.addMVar(1, lb=0, ub=10000, name="Reac_ramp_up")
    Reac_ramp_down = sub_DW.addMVar(1, lb=0, ub=10000, name="Reac_ramp_down")

    # Methanol compressor
    compr_MeOH = sub_DW.addMVar(1, lb=-10000, ub=10000, name="compr_MeOH")

    # Raw Methanol storage
    init_MeOH_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="initialize_MeOH_sto")
    # fin_MeOH_sto = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_MeOH_sto")
    rawMeOH_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="rawMeOH_storage")
    raw_MeOH_min = sub_DW.addMVar(1, lb=0, ub=10000, name="raw_MeOH_min")
    raw_MeOH_max = sub_DW.addMVar(1, lb=0, ub=10000, name="raw_MeOH_max")
    C_raw_MeOH_1 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_raw_MeOH_1")
    C_raw_MeOH_2 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_raw_MeOH_2")

    # Methanol distillator
    dest_max = sub_DW.addMVar(1, lb=0, ub=10000, name="Dest_max")
    dest_min = sub_DW.addMVar(1, lb=0, ub=10000, name="Dest_min")
    dest_rawMeOH = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Dest_rawMeOH")
    dest_water = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Dest_water")
    dest_heat = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Dest_heat")
    dest_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="Dest_sto")
    dest_ramp_up = sub_DW.addMVar(1, lb=0, ub=10000, name="Dest_ramp_up")
    dest_ramp_down = sub_DW.addMVar(1, lb=0, ub=10000, name="Dest_ramp_down")

    # Pure Methanol storage
    init_pure = sub_DW.addMVar(1, lb=-10000, ub=10000, name="initialize_pure")
    # fin_pure = sub_DW.addMVar(1, lb=0, ub=10000, name="finalize_pure")
    pure_sto = sub_DW.addMVar(1, lb=-10000, ub=10000, name="pure_storage")
    MeOH_min = sub_DW.addMVar(1, lb=0, ub=10000, name="MeOH_min")
    MeOH_max = sub_DW.addMVar(1, lb=0, ub=10000, name="MeOH_max")
    C_rate_MeOH_1 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_MeOH_1")
    C_rate_MeOH_2 = sub_DW.addMVar(1, lb=0, ub=10000, name="C_rate_MeOH_2")

    # Meet demand
    balance = sub_DW.addMVar(1, lb=0, name="balance")

#   SUB CONSTRAINTS
    # Energy
    e_buy = sub_DW.addConstr(-p_flow - max_grid <= data.Prices[data.hours-2], name='e_buy')
    e_sell = sub_DW.addConstr(p_flow <= 0, name='e_sell')

    # Battery
    e_bat_cha = sub_DW.addConstr(- bat_sto*data.Battery[0] + p_flow - C_bat_1 <= 0)
    e_bat_dis = sub_DW.addConstr(+ bat_sto/data.Battery[1] - p_flow - C_bat_2 <= 0)

    # Desalinator
    e_des = sub_DW.addConstr(des_H2O + p_flow <= 0, name='e_des')

    # Electrolyzer
    m_h2_elec = sub_DW.addConstr(elec_op*data.Electrolyzer[0] + elec_H2O*data.Electrolyzer[1] + compressor_h2 <= 0, name='m_h2_elec')
    m_H2O_elec = sub_DW.addConstr(- elec_H2O - des_H2O*data.Desalinator[0] <= 0, name='m_H2O_elec')

    # H2 Storage
    m_h2_sto_cha = sub_DW.addConstr(- h2_storage*data.H2storage[0] - C_rate_h2_1 - compressor_h2 <= 0, name='m_h2_sto_cha')
    m_h2_sto_dis = sub_DW.addConstr(h2_storage/data.H2storage[1] - C_rate_h2_2 - reac_H2 - compressor_e*data.Compressor[0] <= 0, name='m_h2_sto_dis')

    # CO2 Storage
    m_CO2_buy = sub_DW.addConstr(CO2_str_min - CO2_str_max + CO2_comp <= data.CO2_prices[data.hours-2], name='m_CO2_buy')
    m_CO2_sto_cha = sub_DW.addConstr(- CO2_sto*data.CO2storage[0] - CO2_comp <= 0, name='m_CO2_sto_cha')
    m_CO2_sto_dis = sub_DW.addConstr(CO2_sto/data.CO2storage[1] - compressor_e*data.Compressor[0] - reac_CO2 <= 0, name='m_CO2_sto_dis')
    
    # Syngas compressor
    e_compr = sub_DW.addConstr(p_flow + compressor_e <= 0, name='e_compr')

    # Methanol reactor
    m_h2_reac = sub_DW.addConstr(- compressor_h2 - reac_H2 - compressor_e*data.Compressor[0] <= 0, name='m_h2_reac')
    m_CO2_reac = sub_DW.addConstr(- CO2_comp - reac_CO2 - compressor_e*data.Compressor[0] <= 0, name='m_CO2_reac') 
    e_reac = sub_DW.addConstr(p_flow - reac_el <= 0, name='e_reac')
    q_reac = sub_DW.addConstr(- reac_heat + p_flow*data.Heater[0] <= 0, name='q_reac')

    # MeOH Storage
    m_rawMeOH_sto_cha = sub_DW.addConstr(-rawMeOH_sto*data.raw_MeOH_sto[0] - C_raw_MeOH_1 - compr_MeOH <= 0, name='m_rawMeOH_sto_cha')
    m_rawMeOH_sto_dis = sub_DW.addConstr(+rawMeOH_sto/data.raw_MeOH_sto[1] - C_raw_MeOH_2 - dest_rawMeOH <= 0, name='m_rawMeOH_sto_dis')

    # Methanol distilator
    m_rawMeOH_dest = sub_DW.addConstr(- compr_MeOH - dest_rawMeOH <= 0, name='m_rawMeOH_dest')
    q_dest = sub_DW.addConstr(- dest_heat + p_flow*data.Heater[0] <= 0, name='q_dest')
    m_h2O_dest = sub_DW.addConstr(- dest_water <= 0, name='m_h2O_dest')

    # Pure MeOH Storagess
    m_MeOH_cha = sub_DW.addConstr(- dest_sto - C_rate_MeOH_1 - pure_sto*data.Pure_MeOH_sto[0] <= 0, name='m_MeOH_cha')
    m_MeOH_dis = sub_DW.addConstr(pure_sto/data.Pure_MeOH_sto[1] - C_rate_MeOH_2 + balance <= 0, name="m_MeOH_dis")


    # Sub objective

    c_x = (
    p_flow * (data.Wind[data.hours-2] * P_cap[0] + data.Solar[data.hours-2] * P_cap[1])
    - max_grid * data.grid_max                        # Grid
    + bat_min * P_cap[5] * data.Battery[4]
    - bat_max * P_cap[5] * data.Battery[5]
    - C_bat_1 * P_cap[5] * data.Battery[6]
    - C_bat_2 * P_cap[5] * data.Battery[6]            # Battery
    - elec_max * P_cap[2]
    + elec_min * P_cap[2] * data.op_elec[2]
    - elec_ramp_up * P_cap[2] * data.op_elec[0]
    - elec_ramp_down * P_cap[2] * data.op_elec[1]           # Electrolyzer
    + balance * data.Demand[data.hours-2]                                     # Demand
    + h2_min * P_cap[6] * data.H2storage[4]
    - h2_max * P_cap[6] * data.H2storage[5]
    - C_rate_h2_1 * P_cap[6] * data.H2storage[6]
    - C_rate_h2_2 * P_cap[6] * data.H2storage[6]                     # H2 Storage
    + CO2_min * P_cap[8] * data.CO2storage[4]
    - CO2_max * P_cap[8] * data.CO2storage[5]
    + CO2_str_min * data.CO2storage[6]
    - CO2_str_max * data.CO2storage[7]
    - reac_max * P_cap[3]
    + reac_min * P_cap[3] * data.op_reac[2]
    - Reac_ramp_up * P_cap[3] * data.op_reac[0]
    - Reac_ramp_down * P_cap[3] * data.op_reac[1]                   # Reactor
    + raw_MeOH_min * P_cap[7] * data.raw_MeOH_sto[4]
    - raw_MeOH_max * P_cap[7] * data.raw_MeOH_sto[5]
    - C_raw_MeOH_1 * P_cap[7] * data.raw_MeOH_sto[6]
    - C_raw_MeOH_2 * P_cap[7] * data.raw_MeOH_sto[6]
    - dest_max * P_cap[4]
    + dest_min * P_cap[4] * data.op_dest[2]
    - dest_ramp_up * P_cap[4] * data.op_dest[0]
    - dest_ramp_down * P_cap[4] * data.op_dest[1]         # Distillator
    + MeOH_min * P_cap[9] * data.Pure_MeOH_sto[4]
    - MeOH_max * P_cap[9] * data.Pure_MeOH_sto[5]
    - C_rate_MeOH_1 * P_cap[9] * data.Pure_MeOH_sto[6]
    - C_rate_MeOH_2 * P_cap[9] * data.Pure_MeOH_sto[6]
    )                                                   # Pure Storage
    
        
    pi_A0_x = (
    # Battery SOC terms
    my_pi.e_bat_SOC[data.hours-3] * (- bat_sto + bat_min - bat_max) + my_pi.e_bat_SOC_h * bat_sto +

    # Electricity flow terms
    my_pi.e_elec[data.hours-3] * (p_flow - elec_max + elec_min - elec_op + elec_ramp_up - elec_ramp_down) +
    my_pi.e_elec_h * (- elec_ramp_up + elec_ramp_down ) +

    # Hydrogen storage SOC terms
    my_pi.m_h2_sto_SOC[data.hours-3] * (- h2_storage + h2_min - h2_max) + my_pi.m_h2_sto_SOC_h * h2_storage +

    # CO2 storage SOC terms
    my_pi.m_CO2_sto_SOC[data.hours-3] * (CO2_sto - CO2_sto + CO2_min - CO2_max) + my_pi.m_CO2_sto_SOC_h * CO2_sto +

    # Methanol reactor terms

    my_pi.m_rawMeOH_reac[data.hours-3] * (-reac_max + reac_min + reac_H2 * data.Reactor[0] + reac_CO2 * data.Reactor[1] +
    reac_el * data.Reactor[2] + reac_heat * data.Reactor[3]  + Reac_ramp_up - Reac_ramp_down + compr_MeOH) 
    + my_pi.m_rawMeOH_reac_h * (-Reac_ramp_up + Reac_ramp_down) +

    # Raw methanol storage SOC terms
    my_pi.m_rawMeOH_sto_SOC[data.hours-3] * (rawMeOH_sto  + raw_MeOH_min - raw_MeOH_max) 
    - my_pi.m_rawMeOH_sto_SOC_h * rawMeOH_sto +

    # Methanol distillation terms
    my_pi.m_MeOH_dest[data.hours-3] * (-dest_max + dest_min + dest_rawMeOH * data.Destilator[0] + dest_water * data.Destilator[1] + 
    dest_heat * data.Destilator[3] + dest_sto + dest_ramp_up - dest_ramp_down) + my_pi.m_MeOH_dest_h * (- dest_ramp_up + dest_ramp_down) +

    # Pure methanol storage SOC terms
    my_pi.m_MeOH_SOC[data.hours-3] * (- pure_sto + MeOH_min - MeOH_max) + my_pi.m_MeOH_SOC_h * pure_sto )

# Set the objective in the model
    sub_DW.setObjective(c_x - pi_A0_x - my_kappa, GRB.MAXIMIZE)

    sub_DW.update()

    sub_DW.write('wolfe.lp')

    sub_DW.optimize()


    if sub_DW.Status == GRB.OPTIMAL:

        # Print optimal objective value
        sub_val = sub_DW.objVal

        results.p_flow[data.hours-2]=p_flow.X
        results.max_grid[data.hours-2]=max_grid.X

        results.elec_max[data.hours-2]=elec_max.X
        results.elec_min[data.hours-2]=elec_min.X
        results.elec_op[data.hours-2]=elec_op.X
        results.elec_H2O[data.hours-2]=elec_H2O.X
        results.elec_ramp_up[data.hours-2]=elec_ramp_up.X
        results.elec_ramp_down[data.hours-2]=elec_ramp_down.X

        results.bat_sto[data.hours-2]=bat_sto.X
        results.bat_min[data.hours-2]=bat_min.X
        results.bat_max[data.hours-2]=bat_max.X
        results.C_bat_1[data.hours-2]=C_bat_1.X
        results.C_bat_2[data.hours-2]=C_bat_2.X
            
        results.des_H2O[data.hours-2]=des_H2O.X

        results.h2_storage[data.hours-2]=h2_storage.X
        results.h2_min[data.hours-2]=h2_min.X
        results.h2_max[data.hours-2]=h2_max.X
        results.C_rate_h2_1[data.hours-2]=C_rate_h2_1.X
        results.C_rate_h2_2[data.hours-2]=C_rate_h2_2.X

        results.compressor_h2[data.hours-2]=compressor_h2.X

        results.CO2_sto[data.hours-2]=CO2_sto.X
        results.CO2_min[data.hours-2]=CO2_min.X
        results. CO2_max[data.hours-2]=CO2_max.X
        results.CO2_str_min[data.hours-2]=CO2_str_min.X
        results.CO2_str_max[data.hours-2]=CO2_str_max.X
        results.CO2_comp[data.hours-2]=CO2_comp.X

        results.compressor_e[data.hours-2]=compressor_e.X

        results.reac_max[data.hours-2]=reac_max.X
        results.reac_min[data.hours-2]=reac_min.X
        results.reac_H2[data.hours-2]=reac_H2.X
        results.reac_CO2[data.hours-2]=reac_CO2.X
        results.reac_el[data.hours-2]=reac_el.X
        results.reac_heat[data.hours-2]=reac_heat.X
        results.Reac_ramp_up[data.hours-2]=Reac_ramp_up.X
        results.Reac_ramp_down[data.hours-2]=Reac_ramp_down.X

        results.compr_MeOH[data.hours-2]=compr_MeOH.X

        results.rawMeOH_sto[data.hours-2]=rawMeOH_sto.X
        results. raw_MeOH_min[data.hours-2]=raw_MeOH_min.X
        results.raw_MeOH_max[data.hours-2]=raw_MeOH_max.X
        results.C_raw_MeOH_1[data.hours-2]=C_raw_MeOH_1.X
        results.C_raw_MeOH_2[data.hours-2]=C_raw_MeOH_2.X

        results.dest_max[data.hours-2]=dest_max.X
        results.dest_min[data.hours-2]=dest_min.X
        results.dest_rawMeOH[data.hours-2]=dest_rawMeOH.X
        results.dest_water[data.hours-2]=dest_water.X
        results.dest_heat[data.hours-2]=dest_heat.X
        results.dest_sto[data.hours-2]=dest_sto.X
        results.dest_ramp_up[data.hours-2]=dest_ramp_up.X
        results.dest_ramp_down[data.hours-2]=dest_ramp_down.X

        results.pure_sto[data.hours-2]=pure_sto.X
        results.MeOH_min[data.hours-2]=MeOH_min.X
        results.MeOH_max[data.hours-2]=MeOH_max.X
        results.C_rate_MeOH_1[data.hours-2]=C_rate_MeOH_1.X
        results.C_rate_MeOH_2[data.hours-2]=C_rate_MeOH_2.X

        results.balance[data.hours-2]=balance.X

    return sub_val, results