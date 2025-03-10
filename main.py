from load_data import *
from load_data_one_year import *
from energy_model import *
from energy_model_maint import *
from extract_results import extract_results

Long_model = True
Maintenance = True
Grid = 150 # Maximum value for the grid connection


if Long_model:
    input_file = load_data()
else:
    input_file = load_data_one()

if Maintenance == False:
    opt_obj, timed, results = energy_model(input_file,Grid)
else:
    opt_obj, timed, results = energy_model_maint(input_file,Grid)

extract_results(results)
