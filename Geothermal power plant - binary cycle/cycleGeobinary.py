import sys
sys.path.append('C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/functions')
import cantera as ct
import numpy as np
import heatmachines as hm
import scipy.interpolate as interpo
import matplotlib.pyplot as plt


"""
Centrale Géothermique – cycle binaire &
Geothermal power plant - binary cycle

"""



def cycleGeobinary(T0 = 25, Ti = 150, P6 = 1, dm = 150,
                   Tech = 90, P3 = 3.25, T3 = 145, T4 = 80, Tout = 70, Twater = 10,
                   P4 = 0.6, eta_pump = 0.9 ):
    
    """
    T in celsius (C)
    dm in kg/s 
    P in MPa 
    """

    

    #START
    
    Ns = 11
    m  = np.zeros(Ns+1)
    h  = np.zeros(Ns+1)
    P  = np.zeros(Ns+1)
    T  = np.zeros(Ns+1)
    X  = np.zeros(Ns+1)
    V  = np.zeros(Ns+1)

    #
    # multiple fluid initialisations
    
    fluid        = ct.CarbonDioxide()
    fluid.basis  = 'mass'
    
    fluid1        = ct.Water()
    fluid1.basis  = 'mass'
    
    fluid2        = ct.Oxygen()
    fluid2.basis  = 'mass'
    
    
    #
    # Thermodynamic cycle - Air
    
    
    #Stage 8
    
    P[8]            = hm.MPa2Pa(P6)
    T[8]            = hm.C2K(T0)
    
    #compute
    fluid2.TP       = T[8], P[8]
    h[8]            = fluid2.h
    
    
    #Stage 9
    
    P[9]            = P[8]
    T[9]            = hm.C2K(Tout)        
    
    #compute
    fluid2.TP       = T[9], P[9]
    h[9]            = fluid2.h
    
    
    #
    # Thermodynamic cycle - CO2
    
    
    fluid2.HP       = h[9], P[9]
    
    
    #Stage 1
    
    P[1]            = hm.MPa2Pa(P4)
    X[1]            = 0
    
    #compute
    fluid.PX        = P[1], X[1]
    h[1]            = fluid.h
    

    #Stage 2
    
    P[2]            = hm.MPa2Pa(P3)
    fluid           = hm.pump_comp(fluid,P[2],eta_pump)   
    
    #compute
    h[2]            = fluid.h
       
    
    #Stage 3
    
    P[3]            = P[2]
    T[3]            = hm.C2K(T3)
    
    #compute
    fluid.TP        = T[3], P[3]
    V[3]            = fluid.v
    h[3]            = fluid.h
    
    
    # Stage 4
    # Initialisation fluid in pump 3
    
    fluid.HP        = h[3], P[3]
    
    T[4]            = hm.C2K(T4)   
    P[4]            = hm.MPa2Pa(P4)
    
    # compute
    fluid.TP        = T[4], P[4]
    h[4]            = fluid.h    
    
    
    #
    # Thermodynamic cycle - Water
    
    
    #Stage 5
    
    T[5]            = hm.C2K(Ti)
    X[5]            = 0
    
    #compute
    fluid1.TX       = T[5], X[5]
    h[5]            = fluid1.h


    # Stage 6
    
    P[6]            = hm.MPa2Pa(P6)                                                        
    
    #compute
    fluid1          = hm.pump_comp(fluid1,P[6],eta_pump)
    h[6]            = fluid1.h
    T[6]            = fluid1.T
    
    # Stage 7 
    T[7]            = hm.C2K(Tech)                         
    P[7]            = P[6]            
                      
    #compute
    fluid1.TP       = T[7], P[7]
    h[7]            = fluid1.h
    
    
    # Stage 0 
    T[0]            = hm.C2K(T0)                         
    X[0]            = 0 
                           
    #compute
    fluid1.TX       = T[0],X[0]
    h[0]            = fluid1.h
    
    
    #
    # Stage for the calculus of water mass flow rate
    #
    
    #Stage 10
    
    T[10]           = hm.C2K(Twater)  
    P[10]           = hm.MPa2Pa(P6)

    # compute
    fluid1.TP       = T[10],P[10]
    h[10]           = fluid1.h                                 
        
    # Stage 11 
    
    P[11]           = P[8] 
    X[11]           = 1           
                                                                
    # compute
    fluid1.PX       = P[11],X[11]
    h[11]           = fluid1.h                                  
    
    #
    # Mass flow rate - Water
    #
    
    m[6]            = dm  
    m[5]            = m[6]
    m[7]            = m[6]
    
    
    #
    # Mass flow rate - CO2
    #
    
    m[2]            = -m[6]*(h[7] - h[6])/(h[3] - h[2])
    m[1]            = m[2]
    m[3]            = m[2]
    m[4]            = m[3]
          
    #
    # Mass flow rate - CO2
    #

    m[8]            = -m[1]*(h[4]-h[1])/(h[8]-h[9])
    m[9]            = m[8]         
    
    m[10]           = -m[1]*(h[4]-h[1])/(h[10]-h[11])
    m[11]           = m[10]

    #Pnet and efficiency calculus
    
    q_in            = m[3]*(h[3]-h[2])
    
    Pturb1          = -(h[4]-h[3])*m[3]
    
    Pnet            = Pturb1- 0.15*Pturb1
    
    eta_tot         = Pnet / (q_in)

    
    return  m[10], m[8], m[3], eta_tot, Pnet
    
if __name__ == "__main__":  
    
    print('')
    print('Author : COSTRITA Alexandru')
    print('Date   :  02/04/2021')
    print('')
    
    m10, m8, m3, eta_tot, Pnet = cycleGeobinary()
    
    print('-')
    print('=====Geothermal power plant - binary cycle=====')     
    print('-')
    print('')
    print('Part I)')
    print('')    
    print('Mass flow rate of CO2                      = {:.4f} [kg/s]'.format(m3))
    print('Net Power  of the binary cycle             = {:.4f} [MW]'.format(Pnet/1000000))
    print('Efficiency of the binary cycle             = {:.4f} [%]'.format(eta_tot*100))
    print('Enforced mass flow rate of air             = {:.4f} [kg/s]'.format(m8))
    print('Enforced mass flow rate of water           = {:.4f} [kg/s]'.format(m10))
    

def cycleGeoflashbinary(Th2o=250, dm=150, Pch=0.5, Pcond=0.01, Text=25, eta_t=0.83, 
                              Tout=90, P8=3.25, T8=145, P9=0.6,
                              T9=80, eta_p=0.9):

    """
    T in celsius (C)
    dm in kg/s 
    P in MPa 
    
    """

    #START
    
    Ns = 11
    m  = np.zeros(Ns+1)
    h  = np.zeros(Ns+1)
    P  = np.zeros(Ns+1)
    T  = np.zeros(Ns+1)
    X  = np.zeros(Ns+1)

    #
    # multiple fluid initialisations
    
    air          = ct.CarbonDioxide()
    air.basis    = 'mass'
    
    fluid        = ct.Water()
    fluid.basis  = 'mass'
    
    
    #
    # Cycle resolution Water part
    #


    #Stage 0
    
    T[0]            = hm.C2K(Text)  
    X[0]            = 0
    
    #compute
    fluid.TX        = T[0], X[0]
    h[0]            = fluid.h
    
    
    #Stage 1
    
    T[1]            = hm.C2K(Th2o) 
    X[1]            = 0 
    
    #compute
    fluid.TX        = T[1], X[1]
    h[1]            = fluid.h
    
    
    #Stage 2
    
    h[2]            = h[1]
    P[2]            = hm.MPa2Pa(Pch) 
    
    #compute
    fluid.HP        = h[2], P[2]
    X[2]            = fluid.X
        
    
    #Stage 3
    
    P[3]            = P[2]
    X[3]            = 1
    
    #compute
    fluid.PX        = P[3], X[3]
    h[3]            = fluid.h
        
    
    #Stage 4
    
    P[4]            = hm.MPa2Pa(Pcond)
    
    #compute
    fluid           = hm.turbine(fluid, P[4], eta_t)
    h[4]            = fluid.h
        
    
    #Stage 5
    
    P[5]            = P[4]
    X[5]            = 0
    
    #compute
    fluid.PX        = P[5], X[5]
    h[5]            = fluid.h
        
    
    #Stage 6
    
    P[6]            = P[3]
    X[6]            = 0
    
    #compute
    fluid.PX        = P[6], X[6]
    h[6]            = fluid.h
        
    
    #Stage 7
    
    T[7]            = hm.C2K(Tout)   
    P[7]            = P[6]
    
    #compute
    fluid.TP        = T[7], P[7]
    h[7]            = fluid.h


    #
    # Cycle resolution CarbonDioxide part
    #
    

    #Stage 8
    
    T[8]            = hm.C2K(T8)   
    P[8]            = hm.MPa2Pa(P8)
    
    #compute
    air.TP          = T[8], P[8]
    h[8]            = air.h
    
    
    #Stage 9
    
    T[9]            = hm.C2K(T9)   
    P[9]            = hm.MPa2Pa(P9)
    
    #compute
    air.TP          = T[9], P[9]
    h[9]            = air.h
        
    
    #Stage 10
    
    P[10]           = P[9]
    X[10]           = 0
    
    #compute
    air.PX          = P[10], X[10]
    h[10]           = air.h
            
    
    #Stage 11
    
    P[11]           = P[8]
    
    #compute
    air             = hm.pump_comp(air, P[11], eta_p)
    h[11]           = air.h
    
       
    #
    # Mass flow rate calculus
    #
    
    m[1]            = dm                                                               
    m[2]            = m[1]                                          
    m[3]            = m[2]*X[2]                      
    m[4]            = m[3]                                                                       
    m[5]            = m[4]     
    m[6]            = m[1]*(1-X[2])                                          
    m[7]            = m[6]                                                            
    m[8]            = -m[6]*(h[7] - h[6])/(h[8] - h[11])
    m[9]            = m[8]
    m[10]           = m[9]
    m[11]           = m[10]

    # - dmdtC02 * ( h[8] - h[11] ) = dmdtWater * ( h[7] - h[6] )

    #W and Q calculus
    
    W_turb1         = m[8] * (h[9]  -h[8])
    W_pump          = m[10]* (h[11] -h[10])
    q_ot            = m[11]* (h[8]  -h[11])
    q_in            = m[1] * (h[1]  -h[0])
    W_turb2         = m[3] * (h[4]  -h[3])
    
    #Pnet and efficiency calculus  
    
    P_net2          = -W_turb2
    P_net1          = -((0.85)*W_turb1 + W_pump)
    
    eta_flash       = P_net2/q_in
    eta_binary      = P_net1/q_ot

    return m[8], P_net1, eta_binary, P_net2, eta_flash


if __name__ == "__main__":    
    
    airmassflow, Pnetbinary, eta_b, Pnetflash, eta_f = cycleGeoflashbinary()
    
    print('') 
    print('-')
    print('=====Geothermal power plant - combined flash and binary cycle=====')     
    print('-')
    print('') 
    print('Part II)')
    print('') 
    print('Mass flow rate of CO2                     = {:.4f} [kg/s]'.format(airmassflow))
    print('Net Power of the flash cycle              = {:.4f} [MW]'.format(Pnetflash/1000000))
    print('Net Power of the flash binary cycle       = {:.4f} [MW]'.format(Pnetbinary/1000000))
    print('Combined cycle efficiency                 = {:.4f} [%]'.format((eta_b+eta_f)*100))

cycleGeoflashbinary()
