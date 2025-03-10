import numpy as np
from load_data_one_year import *

from benders_sub import *

from Sub import *
from Ray import *

#from energy_model_c import *

input_data = load_data_one()

P_cap = np.array([ 150,  100.,  100.,    10.,    7.,  50.,   5., 600., 600., 600.])
Master_obj = P_cap @ input_data.cost_cap

from Benders import *

elapsed_time, results, bounds, duals = Benders(input_data)