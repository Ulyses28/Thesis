from load_data_one_year import *

from benders_sub import *

from DW_X1_results import *

from DW_sub_1000 import *

from DW_master_1000 import *

from Sub import *

from DW_sub_results_hour import *

data = load_data_one()

P_cap = np.array([ 50,  50.,  100.,    10.,    7.,  10.,   3., 10., 10., 15.])

X_columns = [X1_init]  # List to store all X values, starting with X1_init
num_iterations = 50  # Number of iterations (you can adjust this)
i = 0
sub_res = Sub_Results_hour(data.hours)

# # NOTE Full model without hour separation at 4380

for i in range(1, num_iterations + 1):
    # Call DW_master with the current list of X values
    print(f"-------------------------MAS (Iteration {i})-------------------------")
    mas_obj, my_pi, my_kappa = DW_master(data, P_cap, X_columns)

    sub_res = Sub_Results_hour(data.hours)

    # Call DW_sub to get the next X value
    print(f"------------------------SUB 1 (Iteration {i})--------------------------")
    sub_val, sub_res = DW_sub_1(data, P_cap, my_pi, my_kappa,sub_res)

    # Call DW_sub to get the next X value
    print(f"------------------------SUB 2 (Iteration {i})--------------------------")
    sub_val_2, sub_res = DW_sub_2(data, P_cap, my_pi, my_kappa,sub_res)

    # Append the new X value to the list
    X_columns.append(sub_res)

# a,b,c,d = sub(data, P_cap)