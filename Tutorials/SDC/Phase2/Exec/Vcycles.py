# Script to do a timing test on parareal
import numpy as np
import matplotlib.pyplot as plt
import re, sys, os

Nnodes=3

exestr="./main2d.gnu.ex inputs_2d "
outstr=" >  foobar"
Nx_arr=[128,128,128,128]
#Nx_arr=[64,64,64,64]
#Nx_arr=[32,32,32,32]
Nsteps_arr=[2,4,8,16]
maxV_arr=[1,2,4,6]
maxR_arr=[1,3,6]
Iters=np.zeros([len(maxR_arr),len(maxV_arr),len(Nsteps_arr)])
for kmaxR in range(len(maxR_arr)):
    maxR=maxR_arr[kmaxR]

    Viters=np.zeros([len(maxV_arr),len(Nsteps_arr)])
    Riters=np.zeros([len(maxV_arr),len(Nsteps_arr)])
    Siters=np.zeros([len(maxV_arr),len(Nsteps_arr)])
    for kmaxV in range(len(maxV_arr)):
        maxV=maxV_arr[kmaxV]
        dat=np.empty([0,6])    
        for kNsteps in range(len(Nsteps_arr)):
            Nsteps=Nsteps_arr[kNsteps]
            Nx=Nx_arr[kNsteps]
            runstr=  exestr+ 'Nsteps=' + str(Nsteps) +  ' n_cell=' + str(Nx) 
            runstr+=  ' max_Vcycle=' + str(maxV)
            runstr+=  ' max_res_iter=' + str(maxR)
            runstr+=  ' Nnodes=' + str(Nnodes)
            runstr+= outstr
            print(runstr)
            os.system(runstr)
            file=open("foobar")
            for line in file:
                if re.search("dat",line):
                    dd=line[4:]
                    d=np.fromstring(dd,sep=',')
                    dat=np.vstack([dat,d])
        print(dat)
        Viters[kmaxV,:]=dat[:,4]/(Nnodes-1)
        Riters[kmaxV,:]=dat[:,3]/(Nnodes-1)
        Siters[kmaxV,:]=dat[:,2]

    mksz=12
    kstr='$'+str(int(maxR_arr[kmaxR]))
    plt.figure(4)
    Nrow=len(maxR_arr)
    plt.subplot(Nrow,3,kmaxR*3+3)
    for k in range(0,len(maxV_arr)):
        mstr=kstr+str(int(maxV_arr[k]))+'$'
        plt.plot(dat[:,1],Viters[k,:]/dat[:,1],marker=mstr,markersize=mksz)    
    if (kmaxR == Nrow-1):
        plt.xlabel('N_t')
    if (kmaxR < 1):
        plt.title('Ave V-cycles')


    plt.subplot(Nrow,3,kmaxR*3+2)
    for k in range(0,len(maxV_arr)):
        mstr=kstr+str(int(maxV_arr[k]))+'$'
        plt.plot(dat[:,1],Riters[k,:]/dat[:,1],marker=mstr,markersize=mksz)
    if (kmaxR == Nrow-1):
        plt.xlabel('N_t')
    if (kmaxR < 1):
        plt.title('Ave Res. Iters')

    plt.subplot(Nrow,3,kmaxR*3+1)
    for k in range(0,len(maxV_arr)):
        mstr=kstr+str(int(maxV_arr[k]))+'$'
        plt.plot(dat[:,1],Siters[k,:]/dat[:,1],marker=mstr,markersize=mksz)
    if (kmaxR == Nrow-1):
        plt.xlabel('N_t')
    if (kmaxR < 1):
        plt.title('Ave SDC-sweeps')


#    amrex::Print() << "dat " << stop_time  << ", "<< Nsteps << ", " << tot_SDC_sweep <<", " << tot_res_iter <<", "<< tot_Vcycle <<", "<< phi_error << std::endl;
#plt.figure(1)
#plt.semilogy(np.log2(dat[:,1]),dat[:,5],'-*')
#plt.xlabel('log_2(N_t)')
#plt.ylabel('Error')
#plt.title('Error versus Number of Steps')

#plt.figure(2)
#plt.plot(dat[:,1],dat[:,2]/dat[:,1],'-*')    
#plt.plot(dat[:,1],dat[:,3]/dat[:,1],'-*')    
#plt.plot(dat[:,1],dat[:,4]/dat[:,1],'-*')    
#plt.legend(['Sweeps','Resid','Vcycle'])
#plt.xlabel('log_2(N_t)')
#plt.ylabel('Total')
#plt.title('Total SDC Sweeps, Res, and Vcycles')

plt.show()
