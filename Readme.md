We can run the code using:

python main.py
The data is saved in a variable named mass_dtemp_dt. To run the code faster, we can comment out the data generation step and use the precomputed data.
Otherwise, to generate the data from the self-heating network, uncomment this line mass_dtemp_dt = evolution()

all the computation is done in compute.py file 
call_self_heat.py is used to call the self_heat network and  get the E_nuclear,E_neutrino and C_V for a given density 
