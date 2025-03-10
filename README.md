# Thesis
Whole model of my master thesis which consisted in a investment model to determine the optimal sizing of a methanol plant


ABOUT THE MODEL

main: contains the running script. In oder to run it is necessary to download both the input and results excel files. Additionally, all the files that are imported in main.py are also necessary, such as the load_data classes are the energy_model.

BENDERS C: contains Benders decomposition with the simplified model. The one where the heating and electricity demand of both the reactor and the distillator is removed.

BENDERS Divided: is an attempt of a simplified model. Is an intermediate version

BENDERS Simple: is a even simplier version, just to track mistakes. This is also not an useful version to use.

DW Complete: here there is the full Dantzig-Wolfe decomposition applied. There is hourly division where all the small subproblems are exposed.

DW Simple: an initial version where there is only a big subproblem

DW 1000: Dantzig-Wolfe where the subproblems are divided in stacks of 1000 hours.

ABOUT THE EXCEL FILES: in order to run the models it is necessary to make changes in the excel files. Due to the fact that GitHub only admit files under 25MB, the files have been converted to binary workbooks. To run the models it may be necessary to undo the changes
