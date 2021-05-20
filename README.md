# Gibbs Ensemble Monte Carlo Simulation

Parallel Implementation of the Gibbs Ensemble Monte Carlo Simulation used to predict vapor-gas coexistence. 

Run File [MD_GEMC.py](MD_GEMC.py) and get output in an excel file [Data.csv](Data.csv).

To get a visual output from the plot, please run [Plot_Mult.gnuplot](Plot_Mult.gnuplot)

The Molecular Dynamics Function used are in [MD_MC_LJ.py](MD_MC_LJ.py) and Parallel Functionality from [Parallel_Functions.py](Parallel_Functions.py)

Requires numpy,matplotlib,mpi4py,gnuplot to run.

Originally Meant as a course project for AM5080 High Performance Computing course offered in IIT Madras.
