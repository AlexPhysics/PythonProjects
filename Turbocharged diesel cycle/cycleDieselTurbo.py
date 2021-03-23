import sys
sys.path.append('C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/functions')
import cantera as ct
import numpy as np
import heatmachines as hm 
import matplotlib.pyplot as plt
import scipy.interpolate as interpo



def cycleDieselTurbo(Tamb = 20, eta_c=0.76, dTac = 10,
                     Patm = 1.01325, Ncyl = 4, bore = 0.103, stroke = 0.132, 
                     CR = 17.5, N = 50, beta = 3.27, eta_turb = 0.78, PR = 3.2):
    
    """Cycle Diesel turbochargé
       Pressures   : bar
       Temperature : C
       distance    : m
       N           : tour/s
    """
    
    V_dis = stroke*3.14*((bore**2)/4)
    
    V_TDC = V_dis / (CR - 1)
    
    V_BDC = V_TDC + V_dis
    
    Ns  = 9
    h   = np.zeros(Ns+1)
    P   = np.zeros(Ns+1)
    T   = np.zeros(Ns+1)
    D   = np.zeros(Ns+1)
    Pd  = np.zeros(Ns+1)
    Td  = np.zeros(Ns+1)
    Vd  = np.zeros(Ns+1)
    ud  = np.zeros(Ns+1)
    sd  = np.zeros(Ns+1)
    #
    # gas initialisation
    gas                     = ct.Solution('air.xml')
    gas.basis               = 'mass'
    
    #stage 1

    P[1]            = Patm
    T[1]            = hm.C2K(Tamb)
    #compute
    gas.TP          = T[1],P[1]   
    h[1]            = gas.h
    
    
    #stage 2

    P[2]            = P[1]*PR
    #compute
    gas             = hm.pump_comp(gas, P[2], 0.76)
    h[2]            = gas.h

    #stage 3
    
    T[3]            = T[1] + 10
    P[3]            = P[2]
    #compute
    gas.TP          = T[3],P[3]
    h[3]            = gas.h
    
    #stage 1 diesel - cycle diesel starting
    
    Td[1]           = T[3]
    Pd[1]           = P[3]
    #compute
    gas.TP          = Td[1],Pd[1]
    sd[1]           = gas.s
    ud[1]           = gas.u
    Vd[1]           = gas.v
    

    #stage 2 diesel
    
    md              = V_BDC/Vd[1]
    sd[2]           = sd[1]
    Vd[2]           = V_TDC/md
    #compute
    gas.SV          = sd[2],Vd[2]
    ud[2]           = gas.u
    Pd[2]           = gas.P
    
    
    #stage 3 diesel
    
    Pd[3]           = Pd[2]
    Vd[3]           = Vd[2]*beta
    D[3]            = 1/Vd[3]
    #compute
    gas.DP          = D[3],Pd[3]
    ud[3]           = gas.u
    sd[3]           = gas.s
    
    #stage 4 diesel  
    
    Vd[4]           = Vd[1]
    sd[4]           = sd[3]
    #compute
    gas.SV          = sd[4],Vd[4]
    Td[4]           = gas.T
    ud[4]           = gas.u
    
    #stage 4 - diesel finished
    
    P[4]            = P[3]
    T[4]            = 0.5*(Td[1] + Td[4])
    #compute
    gas.TP          = T[4],P[4]
    h[4]            = gas.h
    
    #stage 5 
    
    P[5]            = Patm
    #compute
    gas             = hm.turbine(gas, P[5], 0.78)
    h[5]            = gas.h
    
    #stage 6 - valve (consider the 1-f energy)
    
    h[6]            = h[4]
    P[6]            = Patm
    #compute
    gas.HP          = h[6],P[6]
  
    
   # f calculus
   
    f        = -(h[2] - h[1])/(h[5]-h[4])

    
   # W and Q exchanges
    w_comp   = h[2] - h[1]
    q_cooler = h[3] - h[2]
    w_diesel = h[4] - h[3]
    w_turb   = f*(h[5] - h[4])
    q_valve  = (1-f)*(h[6]-h[4])
    q_out    = (h[1] - h[5]*f - h[6]*(1-f))
    
    # check energy bilan
    check    = w_comp + q_cooler + w_diesel + w_turb + q_valve + q_out
    hm.checkFirstLaw(check)
    #

    
    wd_compr         = ud[2] - ud[1]
    wd_comb          = -Pd[2] * (Vd[3]-Vd[2])
    qd_comb          = ud[3] - ud[2] - wd_comb
    wd_out           = ud[4] - ud[3]
    q_loss           = ud[1] - ud[4]

    # check energy bilan
    check    = wd_compr + wd_comb + qd_comb + wd_out + q_loss
    hm.checkFirstLaw(check)
    #
    
    wd_net           = wd_out + wd_compr + wd_comb
    
    eta_th           = -wd_net/qd_comb
    
    Pnet             = -(wd_net*100)*md
    
    return eta_th, wd_net, Pnet, V_TDC, V_BDC, f, V_dis
    
    
if __name__ == "__main__":
    
    #Display parameters
    print('-')
    print('-----------------------------------------')
    print('------  Cycle Diesel Turbochargé ------')
    print('-----------------------------------------')
    print('-')
    print('Author : COSTRITA Alexandru')
    print('Date   :  09/03/2021')
    print('-')
    
    eta_t,w_net, P_ne, Vtdc, Vbdc, f1, vdis = cycleDieselTurbo()
    
    print('')
    print('Partie a)')
    print('')
    print('Swept Volume Vdis                 = {:.4f} [l]'.format(vdis*1000))
    print('Clearance volume Vtdc             = {:.4f} [l]'.format(Vtdc*1000))
    print('Bottom point volume Vbdc          = {:.4f} [l]'.format(Vbdc*1000))
    print('')
    print('Partie b)')
    print('')
    print('Net work                          = {:.4f} [MJ/kg]'.format(w_net/1000000))
    print('Net power                         = {:.4f} [MW]'.format(P_ne))
    print('f                                 = {:.4f} [-]'.format(f1))
    print('Efficacity of the diesel cycle    = {:.4f} [-]'.format(eta_t))

    print('-----------------------------------------')
    print('----              The End            ----')
    print('-----------------------------------------')
    

    f1Array      = np.zeros(100)
    PRArray      = np.linspace(3.2,20,100)
    w_netArray   = np.zeros(100)
    P_neArray    = np.zeros(100)
    eta_tArray   = np.linspace(0,25,100)
    

for i in range(len(eta_tArray)):
            eta_tArray[i], w_netArray[i], P_neArray[i], Vtdc, Vbdc, f1Array[i], vdis = cycleDieselTurbo(PR = PRArray[i])


plt.figure(1)        
plt.plot(PRArray,f1Array,'-',color = 'red') 
plt.xlabel('PR ',fontsize=14)
plt.ylabel('f',fontsize=14)
plt.title('Pressure Ratio to f graph')
plt.savefig("C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/mod4_rankine/Turbocharged diesel cycle/Graphs/ftoPR.png")
plt.grid(True)


plt.figure(2)
plt.plot(PRArray,eta_tArray,'-',color='green')
plt.xlabel('PR',fontsize=14)
plt.ylabel('Efficacity',fontsize=14)
plt.title('Cycle efficienct to Pressure Ratio graph')
plt.ylim()
plt.savefig("C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/mod4_rankine/Turbocharged diesel cycle/Graphs/EfficacitytoPR.png")
plt.grid(True)


plt.figure(3)
plt.plot(PRArray, P_neArray,'-',color='black')
plt.xlabel('PR',fontsize=14)
plt.ylabel('Net Power [MW]',fontsize=14)
plt.title('Net Power to Pressure Ratio graph')
plt.grid(True)
plt.savefig("C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/mod4_rankine/Turbocharged diesel cycle/Graphs/NetPowertoPR.png")
plt.show()


cycleDieselTurbo()