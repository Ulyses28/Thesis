from load_data_one_year import *
import numpy as np


class X1_Results:
    def __init__(self, 
                 # Power and grid
                 p_flow, max_grid,
                 
                 # Electrolyzer operation
                 elec_max, elec_min, elec_op, elec_H2O, elec_ramp_up, elec_ramp_down,
                 
                 # Battery
                 init_bat, fin_bat, bat_sto, bat_min, bat_max, C_bat_1, C_bat_2,
                 
                 # Desalinator
                 des_H2O,
                 
                 # Hydrogen storage
                 init_h2, fin_h2, h2_storage, h2_min, h2_max, C_rate_h2_1, C_rate_h2_2,
                 
                 # H2 Compressor
                 compressor_h2,
                 
                 # CO2 Storage
                 init_CO2, fin_CO2, CO2_sto, CO2_min, CO2_max, CO2_str_min, CO2_str_max, CO2_comp,
                 
                 # Feed syngas Compressor
                 compressor_e,
                 
                 # Methanol reactor operation
                 reac_max, reac_min, reac_H2, reac_CO2, reac_el, reac_heat, Reac_ramp_up, Reac_ramp_down,
                 
                 # Methanol compressor
                 compr_MeOH,
                 
                 # Raw Methanol storage
                 init_MeOH_sto, fin_MeOH_sto, rawMeOH_sto, raw_MeOH_min, raw_MeOH_max, C_raw_MeOH_1, C_raw_MeOH_2,
                 
                 # Methanol distillator
                 dest_max, dest_min, dest_rawMeOH, dest_water, dest_heat, dest_sto, dest_ramp_up, dest_ramp_down,
                 
                 # Pure Methanol storage
                 init_pure, fin_pure, pure_sto, MeOH_min, MeOH_max, C_rate_MeOH_1, C_rate_MeOH_2,
                 
                 # Meet demand
                 balance):
        
        # Power and grid
        self.p_flow = p_flow
        self.max_grid = max_grid
        
        # Electrolyzer operation
        self.elec_max = elec_max
        self.elec_min = elec_min
        self.elec_op = elec_op
        self.elec_H2O = elec_H2O
        self.elec_ramp_up = elec_ramp_up
        self.elec_ramp_down = elec_ramp_down
        
        # Battery
        self.init_bat = init_bat
        self.fin_bat = fin_bat
        self.bat_sto = bat_sto
        self.bat_min = bat_min
        self.bat_max = bat_max
        self.C_bat_1 = C_bat_1
        self.C_bat_2 = C_bat_2
        
        # Desalinator
        self.des_H2O = des_H2O
        
        # Hydrogen storage
        self.init_h2 = init_h2
        self.fin_h2 = fin_h2
        self.h2_storage = h2_storage
        self.h2_min = h2_min
        self.h2_max = h2_max
        self.C_rate_h2_1 = C_rate_h2_1
        self.C_rate_h2_2 = C_rate_h2_2
        
        # H2 Compressor
        self.compressor_h2 = compressor_h2
        
        # CO2 Storage
        self.init_CO2 = init_CO2
        self.fin_CO2 = fin_CO2
        self.CO2_sto = CO2_sto
        self.CO2_min = CO2_min
        self.CO2_max = CO2_max
        self.CO2_str_min = CO2_str_min
        self.CO2_str_max = CO2_str_max
        self.CO2_comp = CO2_comp
        
        # Feed syngas Compressor
        self.compressor_e = compressor_e
        
        # Methanol reactor operation
        self.reac_max = reac_max
        self.reac_min = reac_min
        self.reac_H2 = reac_H2
        self.reac_CO2 = reac_CO2
        self.reac_el = reac_el
        self.reac_heat = reac_heat
        self.Reac_ramp_up = Reac_ramp_up
        self.Reac_ramp_down = Reac_ramp_down
        
        # Methanol compressor
        self.compr_MeOH = compr_MeOH
        
        # Raw Methanol storage
        self.init_MeOH_sto = init_MeOH_sto
        self.fin_MeOH_sto = fin_MeOH_sto
        self.rawMeOH_sto = rawMeOH_sto
        self.raw_MeOH_min = raw_MeOH_min
        self.raw_MeOH_max = raw_MeOH_max
        self.C_raw_MeOH_1 = C_raw_MeOH_1
        self.C_raw_MeOH_2 = C_raw_MeOH_2
        
        # Methanol distillator
        self.dest_max = dest_max
        self.dest_min = dest_min
        self.dest_rawMeOH = dest_rawMeOH
        self.dest_water = dest_water
        self.dest_heat = dest_heat
        self.dest_sto = dest_sto
        self.dest_ramp_up = dest_ramp_up
        self.dest_ramp_down = dest_ramp_down
        
        # Pure Methanol storage
        self.init_pure = init_pure
        self.fin_pure = fin_pure
        self.pure_sto = pure_sto
        self.MeOH_min = MeOH_min
        self.MeOH_max = MeOH_max
        self.C_rate_MeOH_1 = C_rate_MeOH_1
        self.C_rate_MeOH_2 = C_rate_MeOH_2
        
        # Meet demand
        self.balance = balance

data = load_data_one()

X1_init = X1_Results(
    # Power and grid
    p_flow=np.zeros(data.hours), 
    max_grid=np.zeros(data.hours),
    
    # Electrolyzer operation
    elec_max=np.zeros(data.hours), 
    elec_min=np.zeros(data.hours), 
    elec_op=np.zeros(data.hours), 
    elec_H2O=np.zeros(data.hours), 
    elec_ramp_up=np.zeros(data.hours - 1), 
    elec_ramp_down=np.zeros(data.hours - 1),
    
    # Battery
    init_bat=np.zeros(1), 
    fin_bat=np.zeros(1), 
    bat_sto=np.zeros(data.hours - 1), 
    bat_min=np.zeros(data.hours), 
    bat_max=np.zeros(data.hours), 
    C_bat_1=np.zeros(data.hours), 
    C_bat_2=np.zeros(data.hours),
    
    # Desalinator
    des_H2O=np.zeros(data.hours),
    
    # Hydrogen storage
    init_h2=np.zeros(1), 
    fin_h2=np.zeros(1), 
    h2_storage=np.zeros(data.hours - 1), 
    h2_min=np.zeros(data.hours), 
    h2_max=np.zeros(data.hours), 
    C_rate_h2_1=np.zeros(data.hours), 
    C_rate_h2_2=np.zeros(data.hours),
    
    # H2 Compressor
    compressor_h2=np.zeros(data.hours),
    
    # CO2 Storage
    init_CO2=np.zeros(1), 
    fin_CO2=np.zeros(1), 
    CO2_sto=np.zeros(data.hours - 1), 
    CO2_min=np.zeros(data.hours), 
    CO2_max=np.zeros(data.hours), 
    CO2_str_min=np.zeros(data.hours), 
    CO2_str_max=np.zeros(data.hours), 
    CO2_comp=np.zeros(data.hours),
    
    # Feed syngas Compressor
    compressor_e=np.zeros(data.hours),
    
    # Methanol reactor operation
    reac_max=np.zeros(data.hours), 
    reac_min=np.zeros(data.hours), 
    reac_H2=np.zeros(data.hours), 
    reac_CO2=np.zeros(data.hours), 
    reac_el=np.zeros(data.hours), 
    reac_heat=np.zeros(data.hours), 
    Reac_ramp_up=np.zeros(data.hours - 1), 
    Reac_ramp_down=np.zeros(data.hours - 1),
    
    # Methanol compressor
    compr_MeOH=np.zeros(data.hours),
    
    # Raw Methanol storage
    init_MeOH_sto=np.zeros(1), 
    fin_MeOH_sto=np.zeros(1), 
    rawMeOH_sto=np.zeros(data.hours - 1), 
    raw_MeOH_min=np.zeros(data.hours), 
    raw_MeOH_max=np.zeros(data.hours), 
    C_raw_MeOH_1=np.zeros(data.hours), 
    C_raw_MeOH_2=np.zeros(data.hours),
    
    # Methanol distillator
    dest_max=np.zeros(data.hours), 
    dest_min=np.zeros(data.hours), 
    dest_rawMeOH=np.zeros(data.hours), 
    dest_water=np.zeros(data.hours), 
    dest_heat=np.zeros(data.hours), 
    dest_sto=np.zeros(data.hours), 
    dest_ramp_up=np.zeros(data.hours - 1), 
    dest_ramp_down=np.zeros(data.hours - 1),
    
    # Pure Methanol storage
    init_pure=np.zeros(1), 
    fin_pure=np.zeros(1), 
    pure_sto=np.zeros(data.hours - 1), 
    MeOH_min=np.zeros(data.hours), 
    MeOH_max=np.zeros(data.hours), 
    C_rate_MeOH_1=np.zeros(data.hours), 
    C_rate_MeOH_2=np.zeros(data.hours),
    
    # Meet demand
    balance=np.zeros(data.hours)
)
