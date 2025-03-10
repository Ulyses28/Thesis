import numpy as np

class Sub_Results_hour:
    def __init__(self, hours):
        self.hours = hours

        # Initialize all variables as numpy arrays
        self.p_flow = np.zeros(hours)
        self.max_grid = np.zeros(hours)

        # Battery variables
        self.init_bat = None
        self.fin_bat = None
        self.bat_sto = np.zeros(hours-1)
        self.bat_min = np.zeros(hours)
        self.bat_max = np.zeros(hours)
        self.C_bat_1 = np.zeros(hours)
        self.C_bat_2 = np.zeros(hours)

        self.des_H2O = np.zeros(hours)

        # Electrolyzer variables
        self.elec_min = np.zeros(hours)
        self.elec_max = np.zeros(hours)
        self.elec_op = np.zeros(hours)
        self.elec_H2O = np.zeros(hours)
        self.elec_ramp_up = np.zeros(hours-1)
        self.elec_ramp_down = np.zeros(hours-1)

        # Hydrogen storage variables
        self.init_h2 = None
        self.fin_h2 = None
        self.h2_storage = np.zeros(hours-1)
        self.h2_min = np.zeros(hours)
        self.h2_max = np.zeros(hours)
        self.C_rate_h2_1 = np.zeros(hours)
        self.C_rate_h2_2 = np.zeros(hours)

        self.compressor_h2 = np.zeros(hours)

        # CO2 storage variables
        self.init_CO2 = None
        self.fin_CO2 = None
        self.CO2_sto = np.zeros(hours-1)
        self.CO2_min = np.zeros(hours)
        self.CO2_max = np.zeros(hours)
        self.CO2_str_min = np.zeros(hours)
        self.CO2_str_max = np.zeros(hours)
        self.CO2_comp = np.zeros(hours)

        self.compressor_e = np.zeros(hours)

        # Reactor variables
        self.reac_min = np.zeros(hours)
        self.reac_max = np.zeros(hours)
        self.reac_H2 = np.zeros(hours)
        self.reac_CO2 = np.zeros(hours)
        self.reac_el = np.zeros(hours)
        self.reac_heat = np.zeros(hours)
        self.Reac_ramp_up = np.zeros(hours-1)
        self.Reac_ramp_down = np.zeros(hours-1)

        self.compr_MeOH = np.zeros(hours)

        # Raw methanol storage variables
        self.init_MeOH_sto = None
        self.fin_MeOH_sto = None
        self.rawMeOH_sto = np.zeros(hours-1)
        self.raw_MeOH_min = np.zeros(hours)
        self.raw_MeOH_max = np.zeros(hours)
        self.C_raw_MeOH_1 = np.zeros(hours)
        self.C_raw_MeOH_2 = np.zeros(hours)

        # Distillator variables
        self.dest_min = np.zeros(hours)
        self.dest_max = np.zeros(hours)
        self.dest_rawMeOH = np.zeros(hours)
        self.dest_water = np.zeros(hours)
        self.dest_heat = np.zeros(hours)
        self.dest_sto = np.zeros(hours)
        self.dest_ramp_up = np.zeros(hours-1)
        self.dest_ramp_down = np.zeros(hours-1)

        # Pure methanol storage variables
        self.init_pure = None
        self.fin_pure = None
        self.pure_sto = np.zeros(hours-1)
        self.MeOH_min = np.zeros(hours)
        self.MeOH_max = np.zeros(hours)
        self.C_rate_MeOH_1 = np.zeros(hours)
        self.C_rate_MeOH_2 = np.zeros(hours)

        # Demand balance and capacity
        self.balance = np.zeros(hours)
