
import sys
sys.path.append('C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/functions')
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import heatmachines as hm 

def cycleClaude(Tinit = 30,Tp = 10,Pinit = 3000, Pcond = 1500, 
                eta_t = 0.85, dmdt = 100):
    """Compute Claude
    Pressures   : [kPa]
    Temperature : [°C]
    dmdt : mass flow [kg/s]
    """

    #
    # Initialization
    #
    #arrays
    Ns = 7
    h  = np.zeros(Ns+1)
    P  = np.zeros(Ns+1)
    T  = np.zeros(Ns+1)
    X  = np.zeros(Ns+1)   
    #
    # fluid
    fluid        = ct.Water()
    fluid.basis  = 'mass'
    
    #
    # Cycle
    #
    # Stage 1
    T[1] = hm.C2K(Tinit)  
    X[1] = 0
    # compute
    fluid.TX    = T[1],X[1]  
    h[1]        = fluid.h
    
    # stage 2
    # data

    #
    P[2]        = Pinit
    h[2]        = h[1]
    fluid.HP    = h[2],P[2]
    X[2]        = fluid.X
    
    # stage 3
    # data
    P[3]        = P[2]
    X[3]        = 1
    #compute
    fluid.PX    = P[3], X[3]
    h[3]        = fluid.h
    density     = fluid.density
    
    # stage 4
    #
    P[4]        = P[3]
    X[4]        = 0
    #compute
    fluid.PX    = P[4],X[4]    
    h[4]        = fluid.h
    
    #stage 5
    
    P[5]        = Pcond
    #Reset to step №3
    fluid.PX    = P[3],X[3]
    fluid       = hm.turbine(fluid, P[5], eta_t)
    #compute
    h[5]        = fluid.h  
    
    
    # stage 6
    # data
    T[6]        = hm.C2K(Tp)
    X[6]        = 0
    fluid.TX    = T[6], X[5]

    # compute
    h[6]        = fluid.h  

    # stage 7
    # data
    P[7]        = P[5]
    X[7]        = 0
    fluid.PX    = P[7],X[7]
    #compute
    h[7]        = fluid.h
    
    #
    # mass flow 
    dmdt             = np.zeros(Ns+1)
    dmdt[1]          = 100 #[kg/s]
    dmdt[2]          = dmdt[1]
    dmdt[3]          = dmdt[2]*X[2]
    dmdt[4]          = dmdt[2]*(1-X[2])
    dmdt[5]          = dmdt[3]
    dmdt[6]          = dmdt[5]*(h[7]-h[5])/(h[6]-h[7])
    dmdt[7]          = dmdt[6] + dmdt[5]
    
    #
    #density
    dvdt3            = dmdt[3]/density 
    ###
    #
    ###
    # W and Q exchanges
    P_turb      = ( h[5] - h[3] ) * X[2]
    P_lost      = ( h[4] - h[1] ) * ( 1 - X[2] )
    #
    R           =  (P_turb / P_lost)*100
    #
    
    print('-----------------------------------------')
    print('Rankine cycle Claude')
    print('-----------------------------------------')
    print('-')
    print('Author : COSTRITA Alexandru')
    print('Date   : 02/02/2021')
    print('-') 
    #
    print('Question 1 : ')
    print('-')
    print('Debit masique a lentree : dmdt3 = ',dmdt[3], 'kg/s')
    print('Debit Volumique a lentree : dvdt3 = ',dvdt3,'m^3/s')
    print('-------------------')
    print('Question 2 :')
    print('-')
    print('Puissance de la turbine = ', P_turb, 'W')
    print('Rendement thermique de la centrale= ',R,'%')
    print('-------------------')
    print('Question 3 :')
    print('-')
    print('Dmdt6 eau froide en profondeur = ',dmdt[6], 'kg/s')

cycleClaude()
    
    