import sys
sys.path.append('C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/functions')
import cantera as ct
import numpy as np
import heatmachines as hm 
import matplotlib.pyplot as plt
import scipy.interpolate as interpo


"""
- Etude sur 10  ans (degradation efficacite turbine)
- 

"""

def effTurbineYear(dec_years):
    """
    return efficcite turbine interpolee
    """
    yearTab=[0,0.44,0.74,1.00,1.22,1.38,1.53,1.68,1.79,1.92,2.00,2.12,2.29,2.62,2.93,3.44,3.94,4.46,4.91,5.41,06.01,6.54,6.90,7.15,7.41,7.63,7.94,8.43,8.83,9.33,9.64,10]
    effTab=[0.85,0.85,0.84,0.82,0.79,0.77,0.74,0.71,0.67,0.64,0.62,0.60,0.57,0.55,0.53,0.52,0.51,0.51,0.52,0.52,0.51,0.47,0.43,0.40,0.37,0.36,0.36,0.36,0.36,0.36,0.36,0.36]
    f = interpo.splrep(yearTab,effTab, s=0)
    return(interpo.splev(dec_years, f, der=0))

    


def cycleRankineBase(Pturb=100,Tturb=520,Pcond=0.08,
                     eta_p=0.47,eta_t=0.84,
                     PlantPower=400,
                     Tin=20,Tout=35):
    """Compute basic Rankine Cycle - Module 4 - 020
       Pressures   : bar
       Temperature : °C
       mfr : water mass flow rate kg/s
    """
    
    #
    # Initialization
    #  
    Ns = 4
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
    #stage 2
    #data
    P[2]        = hm.bar2Pa(Pturb)
    T[2]        = hm.C2K(Tturb)    
    fluid.TP    = T[2],P[2]   
    #compute
    h[2]        = fluid.h
    #
    #stage 3
    #data
    P[3]        = hm.bar2Pa(Pcond)
    fluid       = hm.turbine(fluid, P[3], eta_t)
    #compute
    h[3]        = fluid.h    
    #
    #stage 4
    #data
    P[4]        = P[3]
    X[4]        = 0
    fluid.PX    = P[4],X[4]
    #compute
    h[4]        = fluid.h
    #
    #stage 1
    #data
    P[1]        = P[2] 
    fluid       = hm.pump_comp(fluid, P[2], eta_p)   
    #compute
    h[1]        = fluid.h    
    
    #
    # energy
    #
    #
    # w and q exchanges
    #
    w_pump  = h[1] - h[4]
    q_boil  = h[2] - h[1]
    w_turb  = h[3] - h[2]
    q_cond  = h[4] - h[3]
    #
    # net specific work
    w_net   = w_pump + w_turb
    #check energy bilan
    check   = w_pump + w_turb + q_boil + q_cond    
    hm.checkFirstLaw(check)
    #
    # efficiency
    eta_th  = - w_net / q_boil
    #
    # bwr
    bwr     = -w_pump/w_turb    

    #
    # Cooling water 
    #  
    Patm        = ct.one_atm 
    # in cw (cooling water)
    fluid.TP    = hm.C2K(Tin),Patm
    hin         = fluid.h   
    # out cw
    fluid.TP    = hm.C2K(Tout),Patm
    hout        = fluid.h    
    #
    # Energy bilan 
    #  
    mfr = -hm.MW2W(PlantPower)/w_net

    mfr_cw = -mfr*(h[4]-h[3])/(hout-hin)    
    #
    return eta_th,mfr,mfr_cw 


if __name__ == "__main__":    
    
    
    print('-----------------------------------------')
    print('----  Rankine base  cycle resolution ----')
    print('-----------------------------------------')
    print('-')
    print('Author : COSTRITA Alexandru')
    print('Date   : 09/02/2021')
    print('-')  
    
    
    
    Nt = 100
    effTab  = np.zeros(Nt)
    mfrTab  = np.zeros(Nt)
    years   = np.linspace(0,10,Nt)
    print(years)
    eta_turbineTab = effTurbineYear(years) #interpolation efficacite turbine
    for i in range(len(years)):
        effTab[i],mfrTab[i],mfrcw = cycleRankineBase(eta_t=eta_turbineTab[i])
    
    print(effTab)
    
    plt.plot(years,effTab)
    plt.xlabel('Duration')
    plt.ylabel('Efficiency')
    plt.title('Efficiency to duration graph of a rankine cycle')
    plt.grid(True)
    plt.savefig("C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/mod4_rankine/Cycle Rankine  Simple/Graphs/efficiencytoDuration.png")
    plt.show()


    print('-----------------------------------------')
    print('----              The End            ----')
    print('-----------------------------------------')
    
    
cycleRankineBase()