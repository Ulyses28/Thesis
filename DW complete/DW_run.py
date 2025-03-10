from load_data_one_year import *

from benders_sub import *

from DW_X1_results import *

from DW_sub import *

from DW_master import *

# from DW_simp_sub import *

# from DW_simp_master import *

from DW_sub_hour import *

from DW_sub_results_hour import *

data = load_data_one()

P_cap = np.array([ 50,  50.,  100.,    10.,    7.,  10.,   3., 10., 10., 15.])

X_columns = [X1_init]  # List to store all X values, starting with X1_init
num_iterations = 5  # Number of iterations (you can adjust this)

# NOTE Full model without hour separation

for i in range(1, num_iterations + 1):
    # Call DW_master with the current list of X values
    mas_obj, my_pi, my_kappa = DW_master(data, P_cap, X_columns)
    print(f"-------------------------MAS (Iteration {i})-------------------------", mas_obj)

    # Call DW_sub to get the next X value
    sub_val, X_next = DW_sub(data, P_cap, my_pi, my_kappa)
    print(f"------------------------SUB (Iteration {i})--------------------------", sub_val)

    # Append the new X value to the list
    X_columns.append(X_next)

# NOTE Simple model
# for i in range(1, num_iterations + 1):

#     mas_obj, my_pi, my_kappa = DW_master_simp(data, P_cap, X_columns)
#     print(f"-------------------------MAS (Iteration {i})-------------------------", mas_obj)

#     sub_val, X_next = DW_sub_simp(data, P_cap, my_pi, my_kappa)
#     print(f"------------------------SUB (Iteration {i})--------------------------", sub_val)

#     X_columns.append(X_next)

# mas_obj, my_pi, my_kappa = DW_master(data, P_cap, X_columns)

# import time
# start = time.time()

# sub_res = Sub_Results_hour(data.hours)

# DW_sub_hour_0(data,P_cap,my_pi,my_kappa,sub_res)

# for h in range(8757):
#     DW_sub_hour(data,P_cap,my_pi,my_kappa,sub_res,h+1)

# DW_sub_hour_h_1(data,P_cap,my_pi,my_kappa,sub_res)

# DW_sub_hour_h(data,P_cap,my_pi,my_kappa,sub_res)

# end = time.time()
# print(end-start)