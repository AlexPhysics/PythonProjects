import sys
sys.path.append('C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/functions')
import cantera as ct
import numpy as np
import heatmachines as hm
import scipy.interpolate as interpo
import matplotlib.pyplot as plt


"""
Centrale solaire
- 

"""

def etaHours(hour):

    HourTab=[0,0.44,0.74,1.00,1.22,1.38,1.53,1.68,1.79,1.92,2.00,2.12,2.29,2.62,2.93,3.44,3.94,4.46,4.91,5.41,06.01,6.54,6.90,7.15,7.41,7.63,7.94,8.43,8.83,9.33,9.64,10]
    effTab=[0.85,0.85,0.84,0.82,0.79,0.77,0.74,0.71,0.67,0.64,0.62,0.60,0.57,0.55,0.53,0.52,0.51,0.51,0.52,0.52,0.51,0.47,0.43,0.40,0.37,0.36,0.36,0.36,0.36,0.36,0.36,0.36]
    f = interpo.splrep(HourTab,effTab, s=0)
    return(interpo.splev(hour, f, der=0))

def Hour(dec_hour):
    HourTab = [0,1,2,3,4,5,5.68,5.76,5.880000000000001,6.039999999999999,6.159999999999998,6.239999999999998,
6.359999999999999,6.639999999999999,6.879999999999999,7.16,7.4,7.879999999999999,
8.399999999999999,8.399999999999999,8.92,9.36,10.16,10.839999999999998,11.399999999999999,
12.119999999999997,13.079999999999998,13.679999999999996,14.479999999999997,15.279999999999994,
15.879999999999999,16.399999999999995,16.919999999999995,17.359999999999996,17.719999999999995,
17.959999999999997,18.159999999999997,18.239999999999995,18.359999999999996,19,20,21,22,23,24]
    rayTab = [ 0,0,0,0,0,0,54.34782608695684, 100.00000000000011, 141.30434782608722, 180.43478260869574
, 221.73913043478262, 245.6521739130435, 278.2608695652175, 321.7391304347825, 358.695652173913
, 389.13043478260863, 415.21739130434776, 452.17391304347814, 480.4347826086955, 480.4347826086955
, 493.47826086956513, 506.52173913043464, 515.2173913043476, 519.5652173913043, 521.7391304347824
,521.7391304347824, 519.5652173913043, 515.2173913043476, 506.52173913043464, 489.1304347826085
, 467.391304347826, 439.1304347826085, 389.13043478260863, 319.5652173913044, 243.47826086956536
, 180.43478260869574, 143.47826086956536, 97.82608695652198, 52.17391304347848,0,0,0,0,0,0]
    return(np.interp(dec_hour, HourTab, rayTab))





def cycleSolaire(Thfin= 450,Tc=30,
                     DTreheat=10,DTcond=5,DTcfwh=2,DTboil=15,Tamb = 20,
                     Pb=10,Pext1=1.1,Pext2=0.25,Pext3=0.1,
                     eta_p1=0.65,eta_p2 = 0.67,eta_p3 = 0.69,eta_t1=0.87,eta_t2 = 0.9,eta_t3 =0.92,eta_t4 = 0.93,
                     Acol = 40000,Psol =1000,Lc = 0.35,ray = 0.2):

    

    #START
    
    Ns = 14
    h  = np.zeros(Ns+1)
    P  = np.zeros(Ns+1)
    T  = np.zeros(Ns+1)
    X  = np.zeros(Ns+1)
    #
    # fluid initialisation
    fluid        = ct.Water()
    fluid.basis  = 'mass'


    #stage 2

    P[2]            = hm.MPa2Pa(Pb)
    T[2]            = hm.C2K(Thfin)-DTboil
    #compute
    fluid.TP        = T[2],P[2]   
    h[2]            = fluid.h


    #stage 3

    P[3]            = hm.MPa2Pa(Pext1)
    #compute
    fluid           = hm.turbine(fluid, P[3], eta_t1)
    h[3]            = fluid.h   
    
    
    #stage 4

    P[4]            = hm.MPa2Pa(Pext2)
    #compute
    fluid           = hm.turbine(fluid, P[4], eta_t2)
    h[4]            = fluid.h     
     
    
    #stage 5

    P[5]            = hm.MPa2Pa(Pext2)
    T[5]            = hm.C2K(Thfin)-DTreheat  
    #compute    
    fluid.TP        = T[5],P[5]
    h[5]            = fluid.h       


    #stage 6

    P[6]            = hm.MPa2Pa(Pext3)
    #compute
    fluid           = hm.turbine(fluid, P[6], eta_t3)
    h[6]            = fluid.h
    
    
    #stage 8

    T[8]            = hm.C2K(Tc)+ DTcond
    X[8]            = 0
    #compute
    fluid.TX        = T[8],X[8]
    h[8]            = fluid.h
    P[8]            = fluid.P
    
    
    #stage 7

    fluid.HP        = h[6],P[6]
    P[7]            = P[8]
    #compute
    fluid           =  hm.turbine(fluid, P[7], eta_t4)
    h[7]            = fluid.h 
    
    
    #stage 9
    
    fluid.TX        = T[8], X[8]
    P[9]            = hm.MPa2Pa(Pext3)
    #compute
    fluid           = hm.pump_comp(fluid, P[9], eta_p1)
    h[9]            = fluid.h
    
    
    #stage 10

    P[10]           = hm.MPa2Pa(Pext3)
    X[10]           = 0
    #compute
    fluid.PX        = P[10],X[10]
    h[10]           = fluid.h 
    
    
    #stage 11

    P[11]           = P[2]
    #compute
    fluid           = hm.pump_comp(fluid, P[2], eta_p2)
    h[11]           = fluid.h   
    
    
    #stage 13

    P[13]           = P[3]
    X[13]           = 0
    #compute
    fluid.PX        = P[13],X[13]   
    h[13]           = fluid.h  
    T[13]           = fluid.T
    
    
     #stage 12
     
    P[12]           =  P[2]
    T[12]           = T[13] - DTcfwh
    #compute
    fluid.TP        = T[12],P[12]   
    h[12]           = fluid.h  

       
    #stage 14

    fluid.HP        = h[13],P[13]
    P[14]           = P[2]
    #compute
    fluid           = hm.pump_comp(fluid, P[14], eta_p3)
    h[14]           = fluid.h
    
    
    #stage 8
        
    # energy conservation : f1, f2, h1
    f1              = (h[11]-h[12])/(h[13]-h[3]-h[12]+h[11])
    f2              = (h[10]-h[9])/(h[6]-h[9])

    h[1]            = f1*h[14] + (1-f1)*h[12]
    
    
    #Stage 1
    
    P[1]            = P[2]
    #compute
    fluid.HP        = h[1],P[1]
    
    
    #Energy Calculus
    
    #
    w_pump          = (1-f2)*(h[9]-h[8])*(1-f1)+(h[11]-h[10])*(1-f1)+(h[14]-h[13])*f1     
    w_turb          = (h[3]-h[2])+(h[4]-h[3])*(1-f1)+(1-f1)*(h[6]-h[5])+((1-f1)*(1-f2))*(h[7]-h[6])        
    q_boil_rh       = (h[2]-h[1]) + (h[5]-h[4])*(1-f1)
    q_cond          = (1-f2)*(h[8]-h[7])*(1-f1)

    
    # Net Work
    w_net           = w_pump + w_turb
    
    # efficiency
    eta_th          = - w_net / q_boil_rh
    
    # bwr
    bwr             = w_pump/w_turb + q_cond   

    
    #Power Calculus
    
    Pabs            = Acol * Psol
    Ploss           = Lc*Acol*(Thfin-Tamb)
    Psoltoplant     = Pabs - Ploss 
    eta_pannel      = Psoltoplant/Pabs
    
    #Flow
    
    dm              = Psoltoplant/(q_boil_rh)
    Pnet            = -dm*w_net
    Eta_tot         = Pnet/Pabs
    
    
    if Pnet < 0:
        Pnet = 0
        rend_tot = 0
        rend     = 0
    if Pabs != 0 :
        rend_tot     = -Pnet/Pabs
        rend         = Psoltoplant/Pabs

    #Mass - Flow
    
    dmdt        = np.zeros(Ns+1)
    dmdt[1]     = dm
    dmdt[2]     = dmdt[1]
    dmdt[3]     = dmdt[1]
    dmdt[4]     = dmdt[1]*(1-f1)
    dmdt[5]     = dmdt[1]*(1-f1)
    dmdt[6]     = dmdt[1]*(1-f1)
    dmdt[7]     = dmdt[1]*(1-f1)*(1-f2)
    dmdt[8]     = dmdt[1]*(1-f1)*(1-f2)
    dmdt[9]     = dmdt[1]*(1-f1)*(1-f2)
    dmdt[10]    = dmdt[1]*(1-f1)
    dmdt[11]    = dmdt[1]*(1-f1)
    dmdt[12]    = dmdt[1]*(1-f1)
    dmdt[13]    = dmdt[1]*f1
    dmdt[14]    = dmdt[1]*f1
   
    #check function
    
    check       = w_pump + w_turb + q_cond + q_boil_rh
    hm.checkFirstLaw(check)
    
    return eta_th,bwr,-w_net,f1,f2,q_boil_rh,Pabs,Psoltoplant,eta_pannel,Pnet,dm,Eta_tot,-rend_tot
    
    
    
      

if __name__ == "__main__":    
    
    
    print('-----------------------------------------')
    print('------  Circuit Solaire ------')
    print('-----------------------------------------')
    print('-')
    print('Author : COSTRITA Alexandru')
    print('Date   :  03/03/2021')
    print('-')
    
    eta_th,bwr,wnet,fa,fb,q_boil,Pab,Psole,etapl,Pnet,dm,Eta_t,rdm = cycleSolaire()

    #
    print('Part - 1')
    #
    print('Net Work                 = {:.5e} [W]'.format(wnet))
    print('Injected energy (q_boil) = {:.4f} [J]'.format(q_boil))
    print('f1                       = {:.4f} [-]'.format(fa))
    print('f2                       = {:.4f} [-]'.format(fb))
    print('efficiency               = {:.4f} [-]'.format(eta_th)) 
    #
    print('Part - 2')
    #
    print('P Net                    = {:.5f} [MW]'.format(Pnet/1000000))
    print('P absorbed              = {:.5f} [MW]'.format(Pab/1000000))
    print('P received               = {:.5f} [MW]'.format(Psole/1000000))
    print('Pannel efficiency        = {:.4f} [-]'.format(etapl))
    print('Plant  efficiency        = {:.4f} [-]'.format(Eta_t))
    print('Mass flow rate           = {:.4f} [Kg.s-1]'.format(dm))
    #
    print('Part - 3 : ')
    #
    print('-')

    #Graph 1 - efficiency

    TempAr = np.linspace(400,1000,100) # 400 - min T ; 1000 - max T [K]
    EffAr  = np.zeros(100)   #table of 100
    QAr    = np.zeros(100) 
    
    for i in range(len(TempAr)):
            eta_th,bwr,wnet,fa,fb,q_boil,Pab,Psole,etapl,Pnet,dm,Eta_t,EffAr[i]  = cycleSolaire(Thfin=TempAr[i])   
    
   #Graph 2 - max eff to temp
    Valmax       = None
    Imax         = None
    for i, num in enumerate(EffAr):
        if (Valmax is None or num > Valmax):
            Valmax = num
            Imax  = i
    print('Max efficiency  = {:.5f} [-]'.format(Valmax),'Max Temp  = {:.5f} [K]'.format(TempAr[Imax]))
    print("Above this temperature the thermal capacity of the saturated heat transfer fluid does not allow a good efficiency.") 
    
  
    # For Graph 3 : 
        
    rd_tot  = np.zeros(100)
    Pnet    = np.zeros(100)
    effTab  = np.zeros(100)
    hour    = np.linspace(0,24,100)
    eta_p   = Hour(hour)

for i in range(len(eta_p)):
            eta_th,bwr,wnet,fa,fb,q_boil,Pab,Psole,etapl,Pnet[i],dm,etat_tot,rd_tot[i] = cycleSolaire(Psol=eta_p[i])


plt.figure(1)        
plt.plot(TempAr,EffAr,'-',label='Efficacity',) 
plt.legend(loc='best')
plt.xlabel('Fluid Temp [K] ',fontsize=11)
plt.ylabel('Efficiency',fontsize=11)
plt.title('Temp to efficiency graph')
plt.savefig("C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/mod4_rankine/Solar power plant/Graphs/temptoeff.png")
plt.grid(True)


plt.figure(2)
plt.plot(hour,eta_p,'-',label = 'average daily radiation',color='green',)
plt.xlabel('hours [h]',fontsize=13)
plt.ylabel('radiation [W/m2]',fontsize=13)
plt.title('radiation Antanan to hours')
plt.savefig("C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/mod4_rankine/Solar power plant/Graphs/radiation.png")
plt.grid(True)

   
plt.figure(3)
plt.plot(hour,Pnet,'-',label = 'Daily power',color='black',)
plt.xlabel('hours [h]',fontsize=11)
plt.ylabel('Net P [W]',fontsize=11)
plt.title('Net P to hours')
plt.grid(True)
plt.savefig("C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/mod4_rankine/Solar power plant/Graphs/netPtoh.png")
plt.show()


plt.figure(4)
plt.plot(hour,rd_tot,'k',label = 'daily efficiency',color='red',)
plt.xlabel('Time [h]',fontsize=13)
plt.ylabel('Total efficiency [-]',fontsize=13)
plt.title('Daily efficiency')
plt.grid(True)
plt.savefig("C:/Users/user/Desktop/COURS/UNIVERSITE/Semestre 5/Thermo TD/HeatMachinesCantera/mod4_rankine/Solar power plant/Graphs/efficiency.png")
plt.show()

cycleSolaire()