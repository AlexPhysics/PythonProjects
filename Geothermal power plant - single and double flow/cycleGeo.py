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
    
    return Pnet2, eta_tot2, X[7]
    
if __name__ == "__main__":  

    Pnet2, eta_tot2, X7 = cycleGeo2flux()
      
    print('')
    print('-')
    print('=====Geothermal power plant - double flow=====')     
    print('-')
    print('')
    print('Partie c)')
    print('')    
    print('Net Power  of double flow cycle    = {:.4f} [MW]'.format(Pnet2/1000000))
    print('Efficacity of double flow cycle    = {:.4f} [%]'.format(eta_tot2*100))
    
    

    Na = 150
    Pnetgraph        = np.zeros(Na)
    Pnetgraph2       = np.zeros(Na)
    temperature      = np.linspace(150,300,Na)
    X2               = np.zeros(Na)
    X7               = np.zeros(Na)
    
    for i in range(len(temperature)):
       m3, m6, X2[i], Pnetgraph[i], eta_tot1  = cycleGeo(Ti = temperature[i])
       Pnetgraph2[i], eta_tot2, X7[i]         = cycleGeo2flux(Ti = temperature[i])
    
    plt.figure(1)
    plt.plot(temperature,-Pnetgraph,'--',label = 'Pnet single flow')
    plt.plot(temperature,-Pnetgraph2,'--',label = 'Pnet double flow')
    plt.xlabel('Temperature [C]',fontsize=13)
    plt.ylabel('Net Power [MW]',fontsize=13)
    plt.title('The evolution of power plant output power in relation to temperature',fontsize=11)
    plt.grid(True)
    plt.legend(loc = 'lower right')
    plt.tight_layout()
    plt.savefig('C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/mod4_rankine/Geothermal power plant - single and double flow/Graphs/Pnet to Temp.pdf')

    plt.figure(2)
    plt.plot(temperature, X2,'--', label = 'X[2]')
    plt.plot(temperature, X7,'--', label = 'X[7]')
    plt.xlabel('X[2]',fontsize=13)
    plt.ylabel('X[7]',fontsize=13)
    plt.title('X[7] to X[2]',fontsize=11)
    plt.legend(loc = 'lower right')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/mod4_rankine/Geothermal power plant - single and double flow/Graphs/X2toX7.pdf')

cycleGeo2flux()
