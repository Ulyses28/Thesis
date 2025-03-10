from load_data import *
import gurobipy as gp
from gurobipy import GRB
import time

def energy_model_maint(data,grid):

# MODEL
    model = gp.Model("Methanol Plant")

# VARIABLES

    # Capacities
    p_cap = model.addMVar(10, lb=0, name="p_cap") #
    # 0:wind, 1:solar, 2:elec, 3:reac, 4: dest, 5: batt, 6:h2_sto, 7:raw_MeOH_sto, 8:CO2_sto, 9:MeOH_sto
    
    # Energy
    e_buy = model.addMVar(data.hours, lb=0, ub=grid, name="e_buy") # Wind energy bought
    e_sell = model.addMVar(data.hours, lb=0, name="e_sell") # Wind energy sold to the grid, currently no revenue (acts like spilled energy)

    # Battery
    e_bat_SOC = model.addMVar(data.hours, lb=0, name="e_bat_SOC") # Storage SOC
    e_bat_cha = model.addMVar(data.hours, lb=0, name="e_bat_cha") # Storage inflow
    e_bat_dis = model.addMVar(data.hours, lb=0, name="e_bat_dis") # Storage outflow

    # Desalinator
    e_des = model.addMVar(data.hours, lb=0, name="e_des") # Energy input to desalinate water

    # Electrolyzer
    e_elec = model.addMVar(data.hours, lb=0, name="e_elec") # Energy input
    m_h2_elec = model.addMVar(data.hours, lb=0, name="m_h2_elec") # Hydrogen production in kg
    m_H2O_elec = model.addMVar(data.hours, lb=0, name="m_H2O") # H2O input
    

    # H2 Storage
    m_h2_sto_SOC = model.addMVar(data.hours, lb=0, name="m_h2_sto_SOC")  # Storage SOC
    m_h2_sto_cha = model.addMVar(data.hours, lb=0, name="m_h2_sto_cha")  # Storage inflow
    m_h2_sto_dis = model.addMVar(data.hours, lb=0, name="m_h2_sto_dis")  # Storage outflow

    # CO2 Storage
    m_CO2_sto_SOC = model.addMVar(data.hours, lb=0, name="m_CO2_sto_SOC")  # Storage SOC
    m_CO2_buy = model.addMVar(data.hours, lb=0, name="m_CO2_buy")  # Storage inflow
    m_CO2_sto_cha = model.addMVar(data.hours, lb=0, name="m_CO2_sto_cha")  # Storage inflow
    m_CO2_sto_dis = model.addMVar(data.hours, lb=0, name="m_CO2_sto_dis")  # Storage outflow

    # Syngas compressor
    e_compr = model.addMVar(data.hours, lb=0, name="e_compr") # Energy input

    # Methanol reactor
    m_h2_reac = model.addMVar(data.hours, lb=0, name="m_h2_reac") # Direct h2 input
    m_CO2_reac = model.addMVar(data.hours, lb=0, name="m_CO2_reac") # Direct CO2 inflow / CO2 storage outflow
    e_reac = model.addMVar(data.hours, lb=0, name="e_reac")  # Energy input
    q_reac = model.addMVar(data.hours, lb=0, name="heat_reac")  # Heat input
    m_rawMeOH_reac = model.addMVar(data.hours, lb=0, name="m_MeOH_reac")  # Methanol production in kg

    # Meth Storage
    m_rawMeOH_sto_SOC = model.addMVar(data.hours, lb=0, name="m_MeOH_sto_SOC") # Storage SOC
    m_rawMeOH_sto_cha = model.addMVar(data.hours, lb=0, name="m_MeOH_sto_cha")# Storage inflow
    m_rawMeOH_sto_dis = model.addMVar(data.hours, lb=0, name="m_MeOH_sto_dis")# Storage outflow

    # Methanol destilator
    m_rawMeOH_dest = model.addMVar(data.hours, lb=0, name="m_MeOH_dest") # Direct meth input
    m_MeOH_dest = model.addMVar(data.hours, lb=0, name="m_pure_dest") # Pure Methanol production in kg
    q_dest = model.addMVar(data.hours, lb=0, name="heat_dest")  # Heat input
    m_h2O_dest = model.addMVar(data.hours, lb=0, name="water_dest")  # Heat input

    # Pure Meth Storage
    m_MeOH_SOC = model.addMVar(data.hours, lb=0, name="m_pure_SOC") # Storage volume
    m_MeOH_cha = model.addMVar(data.hours, lb=0, name="m_pure_cha") # Storage inflow / destilator output
    m_MeOH_dis = model.addMVar(data.hours, lb=0, name="m_pure_dis") # Storage outflow

    # Maintenance 
    m = model.addMVar(data.hours, lb=0) # maintenance variable
    z_elec = model.addMVar(data.hours, lb=0) # Auxiliary variable
    z_reac = model.addMVar(data.hours, lb=0) # Auxiliary variable
    z_dest = model.addMVar(data.hours, lb=0) # Auxiliary variable


# CONSTRAINTS

    # Power flow 
    Power_flow = model.addConstr( data.Wind*p_cap[0] + data.Solar*p_cap[1]  + e_buy + e_bat_dis == e_bat_cha + e_elec + e_reac + e_compr + e_des + (q_reac+q_dest)*data.Heater[0] + e_sell, name='Power_flow')
   
    # Max capacities
    max_cap = model.addConstr( p_cap <= data.max_capacities, name='max_cap')

    # Battery
    initialize_bat = model.addConstr( e_bat_SOC[0] == p_cap[5]*data.Battery[2], name='initialize_bat') # Initial SOC
    finalize_bat = model.addConstr( e_bat_SOC[data.hours-1] == p_cap[5]*data.Battery[3], name='finalize_bat') # Final SOC
    bat_storage = model.addConstr( e_bat_SOC[1:] == e_bat_SOC[:-1] + e_bat_cha[:-1]*data.Battery[0] - e_bat_dis[:-1]/data.Battery[1], name='bat_storage')
    # Battery storage capacity
    bat_min = model.addConstr( e_bat_SOC >= p_cap[5]*data.Battery[4], name='bat_min')
    bat_max = model.addConstr( e_bat_SOC <= p_cap[5]*data.Battery[5], name='bat_max')
    # C rate bat
    C_rate_1 = model.addConstr( e_bat_cha <= p_cap[5]*data.Battery[6], name='C_rate_1')
    C_rate_2 = model.addConstr( e_bat_dis <= p_cap[5]*data.Battery[6], name='C_rate_2')

    # Desalinator
    Des_H2O = model.addConstr(e_des == data.Desalinator[0]*m_H2O_elec, name='Des_H2O')

    # Electrolyzer operation
    Elec_max = model.addConstr(e_elec <= z_elec, name='Elec_max')
    Elec_min = model.addConstr(e_elec >= z_elec*data.op_elec[2], name='Elec_min')
    Elec_op =  model.addConstr(m_h2_elec*data.Electrolyzer[0] == e_elec, name='Elec_op')
    Elec_H2O = model.addConstr(m_h2_elec*data.Electrolyzer[1] == m_H2O_elec, name='Elec_H2O')
    Elec_ramp_up = model.addConstr(e_elec[1:]-e_elec[:-1] <= p_cap[2]*data.op_elec[0], name='Elec_ramp_up')
    Elec_ramp_down = model.addConstr(e_elec[:-1]-e_elec[1:] <= p_cap[2]*data.op_elec[1], name='Elec_ramp_down')

    # Hydrogen storage
    initialize_h2 = model.addConstr( m_h2_sto_SOC[0] == p_cap[6]*data.H2storage[2], name='initialize_h2') # Initial SOC
    finalize_h2 = model.addConstr( m_h2_sto_SOC[data.hours-1] == p_cap[6]*data.H2storage[3], name='finalize_h2') # Final SOC
    h2_storage = model.addConstr( m_h2_sto_SOC[1:] == m_h2_sto_SOC[:-1] + m_h2_sto_cha[:-1]*data.H2storage[0] - m_h2_sto_dis[:-1]/data.H2storage[1], name='h2_storage')
    # Hydrogen storage capacity
    h2_min = model.addConstr( m_h2_sto_SOC >= p_cap[6]*data.H2storage[4], name='h2_min')
    h2_max = model.addConstr( m_h2_sto_SOC <= p_cap[6]*data.H2storage[5], name='h2_max')
    # C rate h2
    C_rate_h2_1 = model.addConstr( m_h2_sto_cha <= p_cap[6]*data.H2storage[6], name='C_rate_h2_1')
    C_rate_h2_2 = model.addConstr( m_h2_sto_dis <= p_cap[6]*data.H2storage[6], name='C_rate_h2_2')

    # H2 Compressor
    compressor_h2 = model.addConstr( m_h2_elec == m_h2_reac + m_h2_sto_cha , name='compressor_h2')

    # CO2 storage
    initialize_CO2 = model.addConstr( m_CO2_sto_SOC[0] == p_cap[8]*data.CO2storage[2], name='initialize_CO2') # Initial SOC
    finalize_CO2 = model.addConstr( m_CO2_sto_SOC[data.hours-1] == p_cap[8]*data.CO2storage[3], name='finalize_CO2') # Final SOC
    CO2_storage = model.addConstr( m_CO2_sto_SOC[1:] == m_CO2_sto_SOC[:-1] + m_CO2_sto_cha[:-1]*data.CO2storage[0] - m_CO2_sto_dis[:-1]/data.CO2storage[1], name='CO2_storage')
    # CO2 storage capacity
    CO2_min = model.addConstr( m_CO2_sto_SOC >= p_cap[8]*data.CO2storage[4], name='CO2_min')
    CO2_max = model.addConstr( m_CO2_sto_SOC <= p_cap[8]*data.CO2storage[5], name='CO2_max')
    # CO2 as a stream
    CO2_stream_min = model.addConstr( m_CO2_buy >= data.CO2storage[6]* m, name='CO2_stream_min')
    CO2_stream_max = model.addConstr( m_CO2_buy <= data.CO2storage[7]* m, name='CO2_stream_max')
    # CO2 Compressor
    CO2_comp = model.addConstr( m_CO2_buy == m_CO2_reac + m_CO2_sto_cha, name='compressor_CO2')


    # Feed syngas Compressor
    compressor_e = model.addConstr(e_compr == (m_CO2_reac + m_CO2_sto_dis + m_h2_reac + m_h2_sto_dis)*data.Compressor[0], name='compressor_e')

    # Methanol reactor operation
    Reac_max = model.addConstr(m_rawMeOH_reac <= z_reac, name='Reac_max')
    Reac_min = model.addConstr(m_rawMeOH_reac >= z_reac*data.op_reac[2], name='Reac_min')
    Reac_H2 =  model.addConstr(m_rawMeOH_reac*data.Reactor[0] == m_h2_reac + m_h2_sto_dis, name='Raw_MeOH_H2')
    Reac_CO2 =  model.addConstr(m_rawMeOH_reac*data.Reactor[1] == m_CO2_reac + m_CO2_sto_dis, name='Raw_MeOH_CO2')
    Reac_el = model.addConstr(m_rawMeOH_reac*data.Reactor[2] == e_reac , name='Reac_el')
    Reac_heat =  model.addConstr(m_rawMeOH_reac*data.Reactor[3] == q_reac, name='Raw_MeOH_heat')
    Reac_ramp_up = model.addConstr(m_rawMeOH_reac[1:]-m_rawMeOH_reac[:-1] <= p_cap[3]*data.op_reac[0], name='Reac_ramp_up')
    Reac_ramp_down = model.addConstr(m_rawMeOH_reac[:-1]-m_rawMeOH_reac[1:] <= p_cap[3]*data.op_reac[1], name='Reac_ramp_down')

    # Methanol compressor
    compressor_MeOH = model.addConstr( m_rawMeOH_reac == m_rawMeOH_dest + m_rawMeOH_sto_cha , name='compressor_MeOH')
        
    # Raw Methanol storage
    initialize_MeOH_sto = model.addConstr( m_rawMeOH_sto_SOC[0] == p_cap[7]*data.raw_MeOH_sto[2], name='initialize_MeOH_sto') # Initial SOC
    finalize_MeOH_sto = model.addConstr( m_rawMeOH_sto_SOC[data.hours-1] == p_cap[7]*data.raw_MeOH_sto[3], name='finalize_MeOH_sto') # Final SOC
    rawMeOH_storage = model.addConstr( m_rawMeOH_sto_SOC[1:] == m_rawMeOH_sto_SOC[:-1] + m_rawMeOH_sto_cha[:-1]*data.raw_MeOH_sto[0] - m_rawMeOH_sto_dis[:-1]/data.raw_MeOH_sto[1], name='rawMeOH_storage')
    # Raw Methanol storage capacity
    raw_MeOH_min = model.addConstr( m_rawMeOH_sto_SOC >= p_cap[7]*data.raw_MeOH_sto[4], name='raw_MeOH_min')    
    raw_MeOH_max = model.addConstr( m_rawMeOH_sto_SOC <= p_cap[7]*data.raw_MeOH_sto[5], name='raw_MeOH_max')
    # C rate raw MeOH
    C_rate_raw_MeOH_1 = model.addConstr( m_rawMeOH_sto_cha <= p_cap[7]*data.raw_MeOH_sto[6], name='C_rate_raw_MeOH_1')
    C_rate_raw_MeOH_2 = model.addConstr( m_rawMeOH_sto_dis <= p_cap[7]*data.raw_MeOH_sto[6], name='C_rate_raw_MeOH_2')
    
    # Methanol distillator
    Dest_max = model.addConstr(m_MeOH_dest <= z_dest, name='Dest_max')
    Dest_min = model.addConstr(m_MeOH_dest >= z_dest*data.op_dest[2], name='Dest_min')
    Dest_rawMeOH =  model.addConstr(m_MeOH_dest*data.Destilator[0] == m_rawMeOH_dest + m_rawMeOH_sto_dis, name='Dest_op')
    Dest_water = model.addConstr( m_MeOH_dest*data.Destilator[1] == m_h2O_dest , name='Dest_water')
    Dest_heat = model.addConstr( m_MeOH_dest*data.Destilator[3] == q_dest , name='Dest_en')
    Dest_sto = model.addConstr(m_MeOH_dest == m_MeOH_cha , name='Dest_sto')
    Dest_ramp_up = model.addConstr(m_MeOH_dest[1:]-m_MeOH_dest[:-1] <= p_cap[4]*data.op_dest[0], name='Dest_ramp_up')
    Dest_ramp_down = model.addConstr(m_MeOH_dest[:-1]-m_MeOH_dest[1:] <= p_cap[4]*data.op_dest[1], name='Dest_ramp_down')
        
    # Pure Methanol storage
    initialize_pure = model.addConstr( m_MeOH_SOC[0] == p_cap[9]*data.Pure_MeOH_sto[2], name='initialize_pure') # Initial SOC
    finalize_pure = model.addConstr( m_MeOH_SOC[data.hours-1] == p_cap[9]*data.Pure_MeOH_sto[3], name='finalize_pure') # Final SOC
    pure_storage = model.addConstr( m_MeOH_SOC[1:] == m_MeOH_SOC[:-1] + m_MeOH_cha[:-1]*data.Pure_MeOH_sto[0] - m_MeOH_dis[:-1]/data.Pure_MeOH_sto[1], name='pure_storage')
    # Pure Methanol storage capacity
    MeOH_min = model.addConstr( m_MeOH_SOC >= p_cap[9]*data.Pure_MeOH_sto[4], name='MeOH_min')    
    MeOH_max = model.addConstr( m_MeOH_SOC <= p_cap[9]*data.Pure_MeOH_sto[5], name='MeOH_max')
     # C rate raw MeOH
    C_rate_MeOH_1 = model.addConstr( m_MeOH_cha <= p_cap[9]*data.Pure_MeOH_sto[6], name='C_rate_MeOH_1')
    C_rate_MeOH_2 = model.addConstr( m_MeOH_dis <= p_cap[9]*data.Pure_MeOH_sto[6], name='C_rate_MeOH_2')

    # Maintenance
    model.addConstr(m == data.Maintenance)
    # Electrolyzer
    model.addConstr(z_elec <= p_cap[2])         
    model.addConstr(z_elec <= data.max_capacities[2] * m)     
    model.addConstr(z_elec >= p_cap[2] - (1 - m) * data.max_capacities[2])
    # Reactor
    model.addConstr(z_reac <= p_cap[3])         
    model.addConstr(z_reac <= data.max_capacities[3] * m)     
    model.addConstr(z_reac >= p_cap[3] - (1 - m) * data.max_capacities[3]) 
    # Destilator
    model.addConstr(z_dest <= p_cap[4])         
    model.addConstr(z_dest <= data.max_capacities[4] * m)     
    model.addConstr(z_dest >= p_cap[4] - (1 - m) * data.max_capacities[4])

    # Meet demand
    balance = model.addConstr( data.Demand_m <= m_MeOH_dis, name='demand')

# OBJECTIVE FUNCTION
    objective = p_cap @ data.cost_cap + e_buy @ data.Prices + m_CO2_buy @ data.CO2_prices
     
    model.setObjective(objective, GRB.MINIMIZE)

    # Start timer
    start_time = time.time()

    # Write the equations
    model.update()

    model.write('equations.lp')

    # Optimize the Gurobi model
    #model.setParam('BarConvTol', 1e-4) #(to reduce Crossover time)
    model.setParam('Method', 2)  # Barrier method
    model.setParam('Crossover', 0)
    model.setParam('MIPGap', 0.01)  # Allow a 0.01% gap
    model.optimize()

    # End timer
    end_time = time.time()

    # Calculate elapsed time
    elapsed_time = end_time - start_time
    
    # Print optimal objective value
    optimal_objective = model.objVal

# Extract results

    results = {
        "Optimal_obj": optimal_objective,
        "Capacities": p_cap.x,   
        "H2_production": m_h2_elec.x,
        "raw_MeOH_production": m_rawMeOH_reac.x,
        "MeOH_production": m_MeOH_dest.x, 
        "Wind production": data.Wind*p_cap[0].x,
        "Solar production": data.Solar*p_cap[1].x,
        "Electricity bought": e_buy.x,
        "CO2 bought": m_CO2_buy.x,
        "Battery SOC": e_bat_SOC.x,
        "H2 SOC": m_h2_sto_SOC.x,
        "CO2 SOC": m_CO2_sto_SOC.x,
        "raw_MeOH SOC": m_rawMeOH_sto_SOC.x,
        "MeOH SOC": m_MeOH_SOC.x,
        "E Electrolyzer": e_elec.x,
        "E Reactor": e_reac.x,
        "E Compressor": e_compr.x,
        "E Desalinator": e_des.x,
        "E Heat": q_reac.x + q_dest.x,
        "E spilled":e_sell.x,
        "CAPEX":p_cap.x * data.cost_cap,
        "E bought":e_buy.x @ data.Prices,
        "CO2 cost":m_CO2_buy.x @ data.CO2_prices
        }

    return optimal_objective, elapsed_time, results