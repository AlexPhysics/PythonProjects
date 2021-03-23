
import sys
sys.path.append('C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/functions')
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import heatmachines as hm 

def cycleRankineCogen02(Pboil=9000, Phigh=1600, Plow=10, PCond=5, Tboil=400,
                         dmdt=1, y = 0.35):
    """Compute Rankine Cycle + Cogen >2
    Pressures   : [kPa]
    Temperature : [°C]
    dmdt : mass flow [kg/s]
    """

    #
    # Initialization
    #
    #arrays
    Ns = 9
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
    # Stage 6
    P[6] = hm.kPa2Pa(Pboil)    
    T[6] = hm.C2K(Tboil)
    # compute
    fluid.TP    = T[6],P[6]  
    h[6]        = fluid.h
    
    # stage 7
    # data

    #
    P[7]        = hm.kPa2Pa(Phigh)
    fluid       = hm.turbine(fluid, P[6], 1)
    #compute
    h[7]        = fluid.h   
    
    # stage 8
    # data
    fluid.TP    = T[6],P[6]
    P[8]        = hm.kPa2Pa(Plow)  
    #compute
    fluid       = hm.turbine(fluid, P[8],1)
    h[8]        = fluid.h
    
    # stage 9 
    #
    P[9]        = P[7]
    X[9]        = 0
    #compute
    fluid.PX    = P[9],X[9]    
    h[9]        = fluid.h
    
    #stage 1
    
    P[1]        = hm.kPa2Pa(Plow)  
    X[1]        = 0
    fluid.PX    = P[1],X[1]
    
    #compute
    h[1]        = fluid.h  
    
    
    # stage 2
    # data
    P[2]        = hm.kPa2Pa(Phigh)
    fluid       = hm.pump_comp(fluid, P[2], 1)

    # compute
    h[2]        = fluid.h  

    # stage 3
    # data
    P[3]        = P[7]
    X[3]        = 0
    fluid.PX    = P[3],X[3]
    #compute
    h[3]        = fluid.h
    
    # Process heater and condenser flow 
    y        = 0.35                    
     
    z = (1-y)*((h[3]-h[2])/(h[7]-h[3]))   
    #
    # mass flow 
    dmdt             = np.zeros(Ns+1)
    dmdt[6]          = 1 #[kg/s]
    dmdt[7]          = dmdt[6]*y
    dmdt[8]          = dmdt[6]*(1-y)
    dmdt[9]          = dmdt[7]*(y-z)
    dmdt[1]          = dmdt[8]
    dmdt[2]          = dmdt[1]
    dmdt[3]          = dmdt[2] + dmdt[7]*z
    dmdt[4]          = dmdt[9] + dmdt[3]
    dmdt[5]          = dmdt[4]
    
    # Stage-4
    P[4]        = hm.kPa2Pa(Phigh)                     # [Pa]
    h[4]        = ((y-z) * h[9]+(1-y+z)*h[3] )  # [J/kg]
    fluid.HP    = h[4],P[4]
    
    # Stage-5
    P[5]        = hm.kPa2Pa(Pboil)                  # [Pa]
    #compute
    fluid       = hm.pump_comp(fluid, P[5], 1)
    h[5]        = fluid.h                           # [J/kg]
    
    ###
    #
    ###
    # W and Q exchanges
    w_pump1     = (1-y)*(h[2] - h[1])
    q_FWH       = (1-y+z)*h[3] - (1-y)*h[2] - z*h[7]
    q_mix       = h[4] - h[3] * (1-y+z) - (y-z)*h[9]
    w_pump2     = h[5] - h[4]
    q_boil      = h[6] - h[5]
    w_turb      = y*h[7] + h[8] * (1-y) - h[6]
    q_process   = (y-z)*(h[9] - h[7])
    q_cond      = (1-y)*(h[1] - h[8])
    
    # check energy bilan
    
    check    = w_pump1 + q_FWH + q_mix + w_pump2 + q_boil + w_turb + q_process + q_cond
    hm.checkFirstLaw(check)
    
    #
    P_boiler = h[6]-h[5]
    
    # we want  dmdt[6]*(h[6]-h[5]) = 25E6
    
    dmdt25 = (25*10**6) / P_boiler
    eta_cogen2 = -(w_turb + q_process) / (q_boil + w_pump1 + w_pump2)
    eta_therm2 = -(q_process) / (q_boil + w_pump1 + w_pump2)
    eta_elec2 = -(w_turb) / (q_boil + w_pump1 + w_pump2)
    return z, dmdt25, eta_cogen2, eta_therm2, eta_elec2


if __name__ == "__main__":    
    
   
    print('-----------------------------------------')
    print('Rankine - Cogen02')
    print('-----------------------------------------')
    print('-')
    print('Author : COSTRITA Alexandru')
    print('Date   : 01/02/2021')
    print('-')  
    
    z1, flux25, eta_cycle, eta_q, eta_w = cycleRankineCogen02()

    print('z         = {:.4e} [-]'.format(z1))
    print('flow to get 25MW    = {:.4f} [kg/s]'.format(flux25))
    print('elec efficiency     = {:.4f} [-]'.format(eta_w))
    print('therm efficiency    = {:.4f} [-]'.format(eta_q))
    print('global efficiency   = {:.4f} [-]'.format(eta_cycle))
    print('-----------------------------------------')
    print('END')
    print('-----------------------------------------')
    
    # Question 2&3
    
    Na = 100
    extractPressure  = np.linspace(10,3000,Na)
    massFlowBoiler   = np.zeros(Na)
    
    for i in range(len(extractPressure)):
        z1, massFlowBoiler[i], eta_cycle, eta_q, eta_w = cycleRankineCogen02(Phigh=extractPressure[i])
    
    plt.figure(1) # setting eta_th to temperature figure
    plt.plot(extractPressure,massFlowBoiler,'--k')
    plt.xlabel('exractPressure [kPa]',fontsize=13)
    plt.ylabel('m_flow in boiler [kg/s]',fontsize=13)
    plt.title(' Mass flow in boiler to extract pressure for Rankine Cogen 2',fontsize=11)
    plt.grid(True)

    plt.tight_layout()

    plt.savefig('C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/mod4_rankine/Rankine Cycle Cogeneration 2/Graphs/mass flow to extractPressure.pdf')

    
cycleRankineCogen02()