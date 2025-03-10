class DW_duals:
    def __init__(self, 
                 # Battery
                 e_bat_SOC_0, e_bat_SOC,e_bat_SOC_h,
                 
                 # Electrolyzer
                 e_elec_0, e_elec, e_elec_h,
                 
                 # Hydrogen storage
                 m_h2_sto_SOC_0, m_h2_sto_SOC, m_h2_sto_SOC_h,

                 # CO2 Storage
                 m_CO2_sto_SOC_0, m_CO2_sto_SOC, m_CO2_sto_SOC_h,

                 # Methanol reactor operation
                 m_rawMeOH_reac_0, m_rawMeOH_reac, m_rawMeOH_reac_h,

                 # Raw Methanol storage
                 m_rawMeOH_sto_SOC_0, m_rawMeOH_sto_SOC, m_rawMeOH_sto_SOC_h,
                 
                 # Methanol distillator
                 m_MeOH_dest_0, m_MeOH_dest, m_MeOH_dest_h,
                 
                 # Pure Methanol storage
                 m_MeOH_SOC_0, m_MeOH_SOC, m_MeOH_SOC_h):
        
        # Battery
        self.e_bat_SOC_0 = e_bat_SOC_0
        self.e_bat_SOC = e_bat_SOC
        self.e_bat_SOC_h = e_bat_SOC_h
        
        # Electrolyzer operation
        self.e_elec_0 = e_elec_0
        self.e_elec = e_elec
        self.e_elec_h = e_elec_h
        
        # Hydrogen storage
        self.m_h2_sto_SOC_0 = m_h2_sto_SOC_0
        self.m_h2_sto_SOC = m_h2_sto_SOC
        self.m_h2_sto_SOC_h = m_h2_sto_SOC_h
        
        # CO2 Storage
        self.m_CO2_sto_SOC_0 = m_CO2_sto_SOC_0
        self.m_CO2_sto_SOC = m_CO2_sto_SOC
        self.m_CO2_sto_SOC_h = m_CO2_sto_SOC_h
        
        # Methanol reactor operation
        self.m_rawMeOH_reac_0 = m_rawMeOH_reac_0
        self.m_rawMeOH_reac = m_rawMeOH_reac
        self.m_rawMeOH_reac_h = m_rawMeOH_reac_h
        
        # Raw Methanol storage
        self.m_rawMeOH_sto_SOC_0 = m_rawMeOH_sto_SOC_0
        self.m_rawMeOH_sto_SOC = m_rawMeOH_sto_SOC
        self.m_rawMeOH_sto_SOC_h = m_rawMeOH_sto_SOC_h
        
        # Methanol distillator
        self.m_MeOH_dest_0 = m_MeOH_dest_0
        self.m_MeOH_dest = m_MeOH_dest
        self.m_MeOH_dest_h = m_MeOH_dest_h
        
        # Pure Methanol storage
        self.m_MeOH_SOC_0 = m_MeOH_SOC_0
        self.m_MeOH_SOC = m_MeOH_SOC
        self.m_MeOH_SOC_h = m_MeOH_SOC_h