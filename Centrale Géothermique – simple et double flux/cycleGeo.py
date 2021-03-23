import sys
sys.path.append('C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/functions')
import cantera as ct
import numpy as np
import heatmachines as hm
import scipy.interpolate as interpo
import matplotlib.pyplot as plt


"""
Centrale Géothermique – simple et double flux &
Geothermal power plant - single and double flow

"""



def cycleGeo(Ti = 200, dm = 100, Pce = 500, eta_turb = 0.83, Pcond = 10, Pturb = 10,Tinit = 25):
    
    """
    T en celsius (C)
    m en kg/s 
    P en kPa 
    """

    

    #START
    
    Ns = 6
    m  = np.zeros(Ns+1)
    h  = np.zeros(Ns+1)
    P  = np.zeros(Ns+1)
    T  = np.zeros(Ns+1)
    X  = np.zeros(Ns+1)
    
    #
    # fluid initialisation
    fluid        = ct.Water()
    fluid.basis  = 'mass'

    #Stage 0 
    
    T[0]        = hm.C2K(Tinit)
    X[0]        = 0
    
    #compute
    fluid.TX    = T[0], X[0]
    h[0]        = fluid.h
    
    #Stage 1 
    
    T[1]        = hm.C2K(Ti)
    X[1]        = 0
    
    #compute
    fluid.TX    = T[1],X[1]
    h[1]        = fluid.h
    
    #Stage 2
    
    h[2]        = h[1]
    P[2]        = hm.kPa2Pa(Pce)
    
    #compute
    fluid.HP    = h[2], P[2]
    T[2]        = fluid.T
    X[2]        = fluid.X

    #Stage 3

    P[3]        = P[2]
    X[3]        = 1
    
    #compute
    fluid.PX    = P[3], X[3]
    h[3]        = fluid.h
    
    #Stage 4
    
    P[4]        = hm.kPa2Pa(Pturb)
    
    #compute
    fluid       = hm.turbine(fluid, P[4], eta_turb)
    h[4]        = fluid.h

    #Stage 5
    
    P[5]        = P[4]
    X[5]        = 0
    
    #compute
    fluid.PX    = P[5], X[5]
    h[5]        = fluid.h
    
    
    #Stage 6
    
    X[6]        = 0
    P[6]        = P[2]
    
    #compute
    fluid.PX    = P[6], X[6]
    h[6]        = fluid.h
    
    #mass flow calculus
    
    m[1]        = dm
    m[2]        = m[1]
    m[3]        = m[2]*X[2]
    m[4]        = m[3]
    m[5]        = m[4]
    m[6]        = m[2]*(1-X[2])

    
    #Part b - energy&efficiency calculus
    q_in        = m[1]*(h[1]-h[0])
    Pnet        = m[3]*(h[4]-h[3]) #J/s
    eta_tot1    = -Pnet/q_in

    
    return m[3], m[6], X[2], Pnet, eta_tot1

if __name__ == "__main__":    
    
    
    print('-----------------------------------------')
    print('------  Centrale Géothermique ------')
    print('-----------------------------------------')
    print('-')
    print('Author : COSTRITA Alexandru')
    print('Date   :  23/03/2021')
    print('')
    print('-')
    
    m3, m6, X2, Pnet, eta_tot1 = cycleGeo()
    
    print('=====Geothermal power plant - single flow=====')      
    print('-')
    print('')
    print('Partie a)')
    print('')
    print('m[3]                              = {:.4f} [kg/s]'.format(m3))
    print('m[6]                              = {:.4f} [kg/s]'.format(m6))
    print('X[2]                              = {:.4f} [-]'.format(X2))
    print('')
    print('Partie b)')
    print('')
    print('Net Power of single flow cycle     = {:.4f} [MW]'.format(Pnet/1000000))
    print('Efficacity of single flow cycle    = {:.4f} [%]'.format(eta_tot1*100))
    
    
    
    
def cycleGeo2flux(Ti = 200, dm = 100, Pce = 500, Pce2 = 150, eta_turb = 0.83, Pcond = 10, Pturb = 10,Tinit = 25):
    

    #START
    
    Ns = 10
    m  = np.zeros(Ns+1)
    h  = np.zeros(Ns+1)
    P  = np.zeros(Ns+1)
    T  = np.zeros(Ns+1)
    X  = np.zeros(Ns+1)
    
    #
    # fluid initialisation
    fluid        = ct.Water()
    fluid.basis  = 'mass'

    
    #Stage 0 
    
    T[0]        = hm.C2K(Tinit)
    X[0]        = 0
    
    #compute
    fluid.TX    = T[0], X[0]
    h[0]        = fluid.h
    
    
    #Stage 1 
    
    T[1]        = hm.C2K(Ti)
    X[1]        = 0
    
    #compute
    fluid.TX    = T[1],X[1]
    h[1]        = fluid.h
    
    #Stage 2
    
    h[2]        = h[1]
    P[2]        = hm.kPa2Pa(Pce)
    
    #compute
    fluid.HP    = h[2], P[2]
    T[2]        = fluid.T
    X[2]        = fluid.X

    #Stage 3

    P[3]        = P[2]
    X[3]        = 1
    
    #compute
    fluid.PX    = P[3], X[3]
    h[3]        = fluid.h
    
    #Stage 4
    
    P[4]        = hm.kPa2Pa(Pturb)
    
    #compute
    fluid       = hm.turbine(fluid, P[4], eta_turb)
    h[4]        = fluid.h

    #Stage 5
    
    P[5]        = P[4]
    X[5]        = 0
    
    #compute
    fluid.PX    = P[5], X[5]
    h[5]        = fluid.h
    
    
    #Stage 6
    
    X[6]        = 0
    P[6]        = P[2]
    
    #compute
    fluid.PX    = P[6], X[6]
    h[6]        = fluid.h
    
    #Stage 7
    
    P[7]        = hm.kPa2Pa(Pce2)
    h[7]        = h[6]
    
    #compute
    fluid.HP    = h[7], P[7]
    X[7]        = fluid.X
    
    #Stage 8
    
    P[8]        = P[7]
    X[8]        = 1
    
    #compute
    fluid.PX    = P[8], X[8]
    h[8]        = fluid.h
    
    #Stage 9 
    
    X[9]        = 1-X[8]
    P[9]        = P[7]
    
    #compute
    fluid.PX    = P[9], X[9]
    h[9]        = fluid.h
    
    
    #mass flow calculus
    
    m[1]        = dm
    m[2]        = m[1]
    m[3]        = m[2]*X[2]
    m[4]        = m[3]
    m[5]        = m[4]
    m[6]        = m[2]*(1-X[2])
    m[7]        = m[6]
    m[8]        = m[7]*X[7]
    m[9]        = m[7]*(1-X[7])
    
    #pnet and efficiency
    q_in       = m[1]*(h[1]-h[0])
    Pnet2      = m[8]*(h[4]-h[8]) + m[3]*(h[4]-h[3])
    eta_tot2   = -Pnet2/q_in
    
    return Pnet2, eta_tot2
    
if __name__ == "__main__":  

    Pnet2, eta_tot2 = cycleGeo2flux()
      
    print('')
    print('-')
    print('=====Geothermal power plant - double flow=====')     
    print('-')
    print('')
    print('Partie c)')
    print('')    
    print('Net Power  of double flow cycle    = {:.4f} [MW]'.format(Pnet2/1000000))
    print('Efficacity of double flow cycle    = {:.4f} [%]'.format(eta_tot2*100))
    
    
    """ A FINIR
    Na = 100
    Pnetgraph        = np.linspace(0.1,30,Na)
    temperature      = np.zeros(Na)
    T0               = 150
    
    for i in range(Na):
    Ti               = T0
    T0               = T0+1
    
    plt.figure(1) # setting eta_th to temperature figure
    plt.plot(Pnetgraph,temperature,'--k')
    plt.xlabel('exractPressure [kPa]',fontsize=13)
    plt.ylabel('m_flow in boiler [kg/s]',fontsize=13)
    plt.title(' Mass flow in boiler to extract pressure for Rankine Cogen 2',fontsize=11)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/mod4_rankine/Cycle Rankine  Cogénération 2/Graphs/mass flow to extractPressure.pdf')


    """
cycleGeo2flux()
