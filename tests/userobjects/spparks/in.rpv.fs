# SPPARKS application RPV input  
seed	        23810 

# define if 1NN or 2NN interaction 
app_style       rpv 1 1 2 

lattice         bcc 1.0
region	        box block 0 8 0 8 0 4  

create_box	box
create_sites	box value i1 1
set             i2 value 1 
set             i2 value 2 if i1 = 1 fraction 0.2  # vacancy

sector          yes nstop 1 
solve_style     tree 

# define 1NN and 2NN (optional) bond energy 
# ebond1 for 1NN; ebond2 for 2NN
# in the order of: 11 12 ... 1N; 21 22 ... 2N; ...; NN 

ebond1          5 -0.611 -0.163 -0.480 -0.651 -0.446     0.126 -0.102 -0.213 -0.038    -0.414 -0.540 -0.365    -0.626 -0.366    -0.271    
ebond2          5 -0.611 -0.163 -0.571 -0.596 -0.631    -0.014 -0.180 -0.193 -0.203    -0.611 -0.576 -0.621    -0.611 -0.496    -0.611 

#               Fe        Cu        Ni        Mn 
migbarrier      1 0.65    3 0.56    4 0.70    5 1.03

# ballistic       10 5.0 10.0 2 1 

# temperature in units of Kelvin  
temperature     873 

diag_style      rpv stats yes list energy vac  
stats           1000  
#dump            1 text 1000000 *.dump id i2 x y z d5 

reset_time      0
#run             10 
