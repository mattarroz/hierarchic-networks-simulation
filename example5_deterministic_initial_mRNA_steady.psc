# Stochastic Simulation Algorithm input file
# --> mRNA --> 

# Reactions
R1:
    $pool > S1
    Ksyn

R2:
    S1 > $pool
    Kmdeg*S1
    
R3:
	S1 > S1 + S2
	Ktransl*S1

R4:
	S2 > $pool
	Kpdeg*S2
 
# Fixed species

 
# Parameters
Ksyn = 1.0
Kmdeg = 0.5
Ktransl = 20.0
Kpdeg = 0.5

 
# Variable species
S1 = 2
S2 = 0
