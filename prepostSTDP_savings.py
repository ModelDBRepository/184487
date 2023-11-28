#!/usr/bin/env python
'''
Long-term synaptic dynamics
Costa et al 2015
:: Figure 3 ::
:: Savings simulation (with two Guassian input profiles) ::
'''

import matplotlib as mpl
mpl.use('TkAgg')


from brian import *
from time import time
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.pylab as plab


exportOn = 1
plotOn = -1 #-1: none, 0: all, 1: only for input1, 2: only for input2
nruns = 1; #Number of runs

for nrun in range(0,nruns):

    close('all')

    #Reset
    #if 'neurons' in vars():
    #    reinit()

    realtime = 0
    stime = 1 * second
    stime2 = 50 * second


    extractFiringRatesOn = 1
    resolution_export = 10; # every x ms
    

    N = 100
    taum = 10 * ms
    Ee = 0 * mV
    taue = 2 * ms

    Fon = 50 * Hz
    Foff = 3 * Hz

    #s = 55.0000e-10
    s = 100.0000e-10
    Amax = 2.
    Amin = 0
    Ainit = 0.1
    Umax = 1.
    Umin = 0
    Uinit = 0.1

    dFBn = 0
    dFBp = 0
    dFFp = 0

    #Short-term plasticity params
    tau_u = 50 * ms
    tau_r = 200 * ms

    #prepostSTDP params: AFBn tau_FBn AFBp tau_FBp AFFp tau_FFp
    params = [0.1771,    0.0327,    0.1548,    0.2302,    0.0618,    0.0666];
    AFBn = params[0]
    tau_FBn = params[1]*1e3 * ms
    AFBp = params[2]
    tau_FBp = params[3]*1e3 * ms
    AFFp = params[4]
    tau_FFp = params[5]*1e3 * ms
    #etaU = 0.35
    etaU = 0.15
    etaA = 0.15
    #etaA = 0.35    

    defaultclock.dt = 1*ms

    # Adex Parameters
    C = 281 * pF
    gL = 30 * nS
    taum = C / gL
    EL = -70.6 * mV
    DeltaT = 2 * mV
    vti = -50.4 * mV
    #vtrest = vti + 5 * DeltaT
    vtrest = -45 * mV
    VTmax = 18 * mV
    tauvt = 50 * ms

    tauw, c, b, Vr = 144 * ms, 4 * nS, 0.0805 * nA, -70.6 * mV # Regular spiking (as in the paper)

    eqs_neuron = """
        dvm/dt=(gL*(EL-vm)+gL*DeltaT*exp((vm-vt)/DeltaT)+I-x)/C : volt
        dvt/dt=-(vt-vtrest)/tauvt : volt
        dx/dt=(c*(vm-EL)-x)/tauw : amp #In the standard formulation x is w
        I : amp
    """

    input1_pos = 25
    input2_pos = 75
    rad = 5

    #Define input 1
    F_input1 = ones(N)*Foff
    #F_input1[input1_pos-rad:input1_pos+rad] = Fon
    for i in range(0,N): #Define gaussian input
        F_input1[i] = exp(-((((i+1)-input1_pos)**2)/(2.0*rad**2)))*(Fon-Foff)+Foff; 

    #Define input 2
    F_input2 = ones(N)*Foff
    #F_input2[input2_pos-rad:input2_pos+rad] = Fon
    for i in range(0,N): #Define gaussian input
        F_input2[i] = exp(-((((i+1)-input2_pos)**2)/(2.0*rad**2)))*(Fon-Foff)+Foff; 


    input = PoissonGroup(N, rates=F_input1)
    neurons = NeuronGroup(1, model=eqs_neuron, threshold='vm>vt', reset="vm=Vr;x+=b;vt=VTmax", freeze = True)
    neurons.vt = vtrest
    neurons.vm = EL
    neurons.I = 0
    neurons.x = 0


    model='''w : 1
             FFp : 1
             FBp : 1
             FBn : 1
             R : 1
             u : 1
             U : 1
             A : 1         
             dFFp/dt=-FFp/tau_FFp : 1 (event-driven)
             dFBp/dt=-FBp/tau_FBp : 1 (event-driven)
             dFBn/dt=-FBn/tau_FBn : 1 (event-driven)
             dR/dt=(1-R)/tau_r : 1 (event-driven)
             du/dt=(U-u)/tau_u : 1 (event-driven)            
             '''

    syn = Synapses(input, neurons, model, pre='''I=s*A*R*u; 
                                                U=clip(U+etaU*(-AFBn*FBn*FBp + AFBp*FBp*FFp),Umin,Umax);
                                                w=U*A;
                                                FFp+=1; R-=R*u; u+=U*(1-u)''',
                                        post='''A=A+etaA*(AFFp*FFp*FBn);
                                        A=A-etaA*0.5*mean(AFFp*FFp*FBn);
                                        A=clip(A,Amin,Amax);
                                        w=U*A;
                                        FBp+=1.;FBn+=1.''')                                    


    #syn.connect_one_to_one(input, neurons)
    syn[:,:]=True
    syn.FBp=0
    syn.FBn=0
    syn.R=1
    #syn.U='rand()*Uinit'
    #syn.A='rand()*Ainit'
    #syn.U[:]=Umin
    #syn.U[:]=0.5
    
    for i in range(0,size(syn.U[:])): #Define gaussian input
        syn.U[i] = exp(-((((i+1)-input1_pos)**2)/(2.0*(rad+0)**2)))*(Umax-Umin)+Umin;
        syn.A[i] = exp(-((((i+1)-input1_pos)**2)/(2.0*(rad+3)**2)))*(Amax-Amin)+Amin;        
    
#    syn.U[input1_pos-rad:input1_pos+rad]=Umax
    #syn.A[:]=Amin
    #syn.A[input1_pos-rad-3:input1_pos+rad+3]=Amax

    Si = SpikeMonitor(input)
    So = SpikeMonitor(neurons)
    Mpost = MultiStateMonitor(neurons, record=True)
    Mrate = PopulationRateMonitor(neurons,bin=100*ms)
    Msyn = MultiStateMonitor(syn, record=True)
    synU = StateMonitor(syn, 'U', record=True)
    synA = StateMonitor(syn, 'A', record=True)
    #Mstdp_post = MultiStateMonitor(stdp.post_group, record=True)
    #Msyn = StateMonitor(synapses, 'W', record=True)


    ref = 200

    if(realtime):
        ion()
        subplot(411)
        raster_plot(Si, refresh=ref*ms, showlast=stime2*3+stime)
        subplot(412)
        raster_plot(So, refresh=ref*ms, showlast=stime2*3+stime)
        plab.ylim([-0.5,0.5])
        subplot(413)
        synU.plot([input1_pos, 50, input2_pos], refresh=ref*ms, showlast=stime2*3+stime)
        plab.ylim([-0.05, 1.05])
        subplot(414)
        synA.plot([input1_pos, 50, input2_pos], refresh=ref*ms, showlast=stime2*3+stime)
        plab.ylim([-0.05, 3.05])    
        show()
    
    start_time = time()

    run(stime)
    #input1.rate = 0.1
    input.rate = F_input2
    run(stime2)
    input.rate = F_input1
    run(stime2)
    input.rate = F_input2
    run(stime2)

    print "Simulation time:", time() - start_time


    #G = NeuronGroup(...)
    #spikemon = SpikeMonitor(G)
    #statemon = StateMonitor(G, 'V', record=range(5))
    #subplot(211)
    #raster_plot(spikemon, refresh=10*ms, showlast=200*ms)
    #subplot(212)
    #statemon.plot(refresh=10*ms, showlast=200*ms)
    #run(1*second)

    if plotOn>=0:
        plt.figure()
        ion()
        subplot(411)
        raster_plot(Si)
        subplot(412)
        raster_plot(So)
        subplot(413)
        plot(syn.U[:], '.')
        subplot(414)
        plot(syn.A[:], '.')
        show()

        plt.figure()
        i1 = 25
        pFBn, = plt.plot(Msyn['FBn'].times, Msyn['FBn'][i1,:])
        pFBp, = plt.plot(Msyn['FBp'].times, Msyn['FBp'][i1,:])
        pPreLTD, = plt.plot(Msyn['FBn'].times, AFBn*Msyn['FBn'][i1,:]*Msyn['FBp'][i1,:])
        pPreLTP, = plt.plot(Msyn['FBp'].times, AFBp*Msyn['FBp'][i1,:]*Msyn['FFp'][i1,:])
        pFFp, = plt.plot(Msyn['FFp'].times, Msyn['FFp'][i1,:])
        #pu, = plt.plot(Msyn['u'].times, Msyn['u'][i,:])
        #pR, = plt.plot(Msyn['R'].times, Msyn['R'][i,:])
        plt.legend([pFFp, pFBn, pFBp, pPreLTD, pPreLTP], ['FFp', 'FBn', 'FBp', 'preLTP', 'preLTD'])
        ion()
        show()


        plt.figure()
        i1 = input1_pos
        i2 = input2_pos
        pu1, = plt.plot(Msyn['U'].times, Msyn['U'][i1,:])
        pA1, = plt.plot(Msyn['A'].times, Msyn['A'][i1,:])
        pu2, = plt.plot(Msyn['U'].times, Msyn['U'][i2,:])
        pA2, = plt.plot(Msyn['A'].times, Msyn['A'][i2,:])
        plt.legend([pu1, pA1, pu2, pA2], ['U1', 'A1', 'U2', 'A2'])
        ion()
        show()


        plt.figure()
        pv, = plt.plot(Mpost['vm'].times, Mpost['vm'][0,:])
        pvt, = plt.plot(Mpost['vt'].times, Mpost['vt'][0,:])
        #pge, = plt.plot(Mpost['ge'].times, Mpost['ge'][0,:])
        plt.legend([pv, pvt], ['Vm', 'vt'])
        plt.show()

        #Plot firing rate
        rates = Mrate.smooth_rate(width=1000*ms,filter='gaussian')    
        if extractFiringRatesOn==0:
            plt.figure()
            pv, = plt.plot(Mrate.times, rates)
            plt.show()


    if exportOn:
        #Export results to be plotted in matlab
        import os as os

        path = 'fromBrian/'

        if os.path.exists(path + 'outParams.br'):
            os.remove(path + 'outParams.br')

        if os.path.exists(path + 'outU_run' + str(nrun) + '.br'):
            os.remove(path + 'outU_run' + str(nrun) + '.br')

        if os.path.exists(path + 'outA_run' + str(nrun) + '.br'):
            os.remove(path + 'outA_run' + str(nrun) + '.br')     

        #f_handle = file(filename, 'a')
        savetxt(path + 'outParams.br', [N, resolution_export, stime, stime2, input1_pos, input2_pos, rad, nruns, Amax], fmt='%f', newline='\n') # Number of postsynaptic neurons
        savetxt(path + 'outU_run' + str(nrun) + '.br', Msyn['U'][:,::resolution_export], fmt='%f', newline='\n') # Save Us
        savetxt(path + 'outA_run' + str(nrun) + '.br', Msyn['A'][:,::resolution_export], fmt='%f', newline='\n') # Save As
        #savetxt(f_handle, Msyn['U'][:,0:2], fmt='%f', newline='\n') # Save Us
        #savetxt(f_handle, Msyn['A'][:,0:2], fmt='%f', newline='\n') # Save As
        #f_handle.close()
    



    if extractFiringRatesOn:    
        #Extract postsynaptic firing rate for input1 and input2
    
        Usim = Msyn['U'][:,:];
        Asim = Msyn['A'][:,:];
    
        post_nspikes1 = So.nspikes
        forget(syn)
        reinit()
    
        neurons.vt = vtrest  
        neurons.vm = EL
        neurons.I = 0
        neurons.x = 0
    
        modelAfter='''R : 1
                 u : 1                        
                 w : 1
                 U : 1            
                 A : 1                       
                 dR/dt=(1-R)/tau_r : 1 (event-driven)
                 du/dt=(U-u)/tau_u : 1 (event-driven)            
                 '''                                  
        synAfter=Synapses(input, neurons, modelAfter, pre='''I=s*A*R*u;
                                                w=U*A;
                                                R-=R*u; u+=U*(1-u)''',
                                                post='''w=U*A''')                                            
             
        synAfter[:,:]=True
        synAfter.R=1
        synAfter.u=Uinit

        synAfter.U=Umin
        synAfter.U[input1_pos-rad:input1_pos+rad]=Umax
        synAfter.A[:]=Amin
        synAfter.A[input1_pos-rad:input1_pos+rad]=Amax
    
        Uaux = synAfter.U[:]
        Uaux[:] = Umin
        Uaux[input2_pos-rad:input2_pos+rad]=Umax
        Aaux = synAfter.A[:]
        Aaux[:] = Amin
        Aaux[input2_pos-rad:input2_pos+rad]=Amax

        @network_operation
        def loadUandA(clock): 
            synAfter.U[:] = Usim[:,clock.t/clock.dt]
            synAfter.A[:] = Asim[:,clock.t/clock.dt]
            #synAfter.U[:] = Uaux[:]
            #synAfter.A[:] = Aaux[:]


        #INPUT 1    
        MsynAfter = MultiStateMonitor(synAfter, record=True)
        Mpost = MultiStateMonitor(neurons, record=True)    

        start_time = time()
        F_input1n = F_input1
        F_input1n[F_input1n<(Fon/2)] = 0
        input.rate = F_input1n
        #run(stime+stime2)
    
        run(stime+stime2*3)
    
        rates_1 = Mrate.smooth_rate(width=1000*ms,filter='gaussian')
        print "Simulation time (for input1 alone):", time() - start_time
    
    
    
    
    
        #INPUT 2
        post_nspikes2 = So.nspikes
        reinit()
    
        neurons.vt = vtrest  
        neurons.vm = EL
        neurons.I = 0
        neurons.x = 0
    
        Mpost = MultiStateMonitor(neurons, record=True)
    
        start_time = time()
        F_input2n = F_input2
        F_input2n[F_input2n<(Fon/2)] = 0
        input.rate = F_input2
    
        run(stime+stime2*3)    
    
        rates_2 = Mrate.smooth_rate(width=1000*ms,filter='gaussian')    
        print "Simulation time (for input2 alone):", time() - start_time
        post_nspikes3 = So.nspikes
    
    
    
        if plotOn>=1:
            '''
            plt.figure()
            ion()
            subplot(411)
            raster_plot(Si)
            subplot(412)
            raster_plot(So)
            subplot(413)
            plot(synAfter.U[:], '.')
            subplot(414)
            plot(synAfter.A[:], '.')
            show()
        
            plt.figure()
            i1 = 25
            i2 = 75
            pu1, = plt.plot(MsynAfter['U'].times, MsynAfter['U'][i1,:])
            pA1, = plt.plot(MsynAfter['A'].times, MsynAfter['A'][i1,:])
            pu2, = plt.plot(MsynAfter['U'].times, MsynAfter['U'][i2,:])
            pA2, = plt.plot(MsynAfter['A'].times, MsynAfter['A'][i2,:])
            plt.legend([pu1, pA1, pu2, pA2], ['U1', 'A1', 'U2', 'A2'])
            ion()
            show()
        
            plt.figure()
            pv, = plt.plot(Mpost['vm'].times, Mpost['vm'][0,:])
            pvt, = plt.plot(Mpost['vt'].times, Mpost['vt'][0,:])
            #pge, = plt.plot(Mpost['ge'].times, Mpost['ge'][0,:])
            plt.legend([pv, pvt], ['Vm', 'vt'])
            plt.show()
            '''
        
            #Plot firing rate
            plt.figure()        
            prates, = plt.plot(Mrate.times, rates)
            prates_1, = plt.plot(Mrate.times, rates_1)
            prates_2, = plt.plot(Mrate.times, rates_2)        
            plt.legend([prates, prates_1, prates_2], ['Learning', 'Input_1', 'Input_2'])
            plt.show()
        
        print "Nspikes before: ", post_nspikes1, "| Nspikes after (Input1): ", post_nspikes2, "| Nspikes after (Input2): ", post_nspikes3
    
        if exportOn:
            if os.path.exists(path + 'rateInput1_run' + str(nrun) + '.br'):
                os.remove(path + 'rateInput1_run' + str(nrun) + '.br')
            if os.path.exists(path + 'rateInput2_run' + str(nrun) + '.br'):
                os.remove(path + 'rateInput2_run' + str(nrun) + '.br')

            savetxt(path + 'rateInput1_run' + str(nrun) + '.br', rates_1, fmt='%f', newline='\n') # Save post spike times
            savetxt(path + 'rateInput2_run' + str(nrun) + '.br', rates_2, fmt='%f', newline='\n') # Save post spike times
        clear()    