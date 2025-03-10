class DW_duals:
    def __init__(self, 
                 # Battery
                 e_bat_SOC,
                 
                 # Electrolyzer
                 e_elec,
                 
                 # Hydrogen storage
                 m_h2_sto_SOC, 

                 # CO2 Storage
                 m_CO2_sto_SOC, 

                 # Methanol reactor operation
                 m_rawMeOH_reac, 

                 # Raw Methanol storage
                 m_rawMeOH_sto_SOC, 
                 
                 # Methanol distillator
                 m_MeOH_dest,
                 
                 # Pure Methanol storage
                  m_MeOH_SOC):
        
        # Battery
        self.e_bat_SOC = e_bat_SOC
        
        # Electrolyzer operation
        self.e_elec = e_elec
        
        # Hydrogen storage
        self.m_h2_sto_SOC = m_h2_sto_SOC
        
        # CO2 Storage
        self.m_CO2_sto_SOC = m_CO2_sto_SOC
        
        # Methanol reactor operation
        self.m_rawMeOH_reac = m_rawMeOH_reac
        
        # Raw Methanol storage
        self.m_rawMeOH_sto_SOC = m_rawMeOH_sto_SOC
        
        # Methanol distillator
        self.m_MeOH_dest = m_MeOH_dest
        
        # Pure Methanol storage
        self.m_MeOH_SOC = m_MeOH_SOC
