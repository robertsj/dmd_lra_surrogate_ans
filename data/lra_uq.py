from detran import *

from lra_mesh import *
from lra_uq_db import get_input

import sys 
import time
import cPickle as pickle
import numpy as np
np.random.seed(1234)

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
except:
    rank = 0
    size = 1

def sample(option=0):
    """ Sample inputs for the LRA problems.
    
    Inputs:
        option : int
            Possible values are:
                < 0 use means
                = 1 sample initial temperature and power
                > 1 sample initial temperature and power and the 
                  values of alpha and gamma (which are then fixed)
     
    """
    
    T0 = 300
    P0 = 1e-6
    SigA20 = 0.08344
    alpha = 3.830e-11
    gamma = 3.034e-3
    if option >= 0:
        
        # T0 is uniformly distributed between 297 and 303
        T0 = T0 + np.random.rand()*6-3           
    
        # Initial power is normal, 1%
        P0 = np.random.normal(loc=P0, scale=0.01*P0)
                      
    if option > 1:
        # The thermal absorption cross section of the control blade
        # normally distributed with a stdev of 1%        
        SigA20 = np.random.normal(loc=SigA20, scale=0.01*SigA20)
                        
        # Alpha is normally distributed with a stdev of 1%
        alpha = np.random.normal(loc=alpha, scale=0.01*alpha)

        # Gamma is normally distributed with a stdev of 1%
        gamma = np.random.normal(loc=gamma, scale=0.01*gamma)
    
    return T0, P0, SigA20, alpha, gamma


def run(num_samples, option) :

    chunk = num_samples // size
    start, end = chunk*rank, chunk*(rank+1)
    if rank == size - 1:
        end = num_samples

    for j in range(start, end):
        np.random.seed(j)

        # Input, mesh, and base material
        inp = get_input()
        fmm = 2 # fine mesh per node
        mesh, cmm = get_mesh_2D(fmm)
        name_base = "diffusion{0}x{0}".format(fmm)
        transport = False
        mat = LRA.Create(mesh, transport, False)
        
        name = name_base + '_{}'.format(j)

        # perturb the base material
        T0, P0, SigA20, alpha, gamma = sample(option)
        as_lra(mat).perturb(T0, SigA20, alpha, gamma)  

        # STEADY STATE
        mat.update(0.0, 0, 1, False)
        manager = Eigen2D(inp, downcast(mat), mesh)
        manager.solve()
        ic = manager.state()
        mat.set_eigenvalue(ic.eigenvalue())
        mat.update(0, 0, 1, False)
        # normalize state.
        F = 0.0;
        matmap = mesh.mesh_map("MATERIAL");
        for i in range(0, mesh.number_cells()) :
            m = matmap[i]
            F += mesh.volume(0) * \
                 (ic.phi(0)[i] * mat.sigma_f(m, 0) + \
                  ic.phi(1)[i] * mat.sigma_f(m, 1))
        F *=  as_lra(mat).d_KAPPA / 17550.0;
        ic.scale(P0/F);
        
        #---------------------------------------------------------------------------#
        # TIME STEP
        #---------------------------------------------------------------------------#

        ts = Time2D(inp, mat, mesh, True)

        # Initialize containers for saving
        times = []
        powers = []
        temps = []
        maxtemps = []
        asspower = []
        asstemp  = []
        phi0=[]
        phi1=[]
        meshpower=[]
        meshtemp=[]
        V = np.zeros(mesh.number_cells())
        mat_map = mesh.mesh_map("COARSEMESH")
        for i in range(0, mesh.number_cells()) :
            if mat_map[i] != 4 : 
                V[i] = mesh.volume(i)  
        area = 17550.0
        rates = ReactionRates(downcast(mat), mesh, ts.state())

        #---------------------------------------------------------------------------#
        def monitor(ts, step, t, dt, it, conv) :

            T    = np.asarray(as_lra(mat).physics().variable(0))
            Tavg = np.sum(T * V) / area
            Tmax = np.max(T)
            state = ts.state()
            totP = 0.0;
            matmap = mesh.mesh_map("MATERIAL");
            mpower = np.zeros(mesh.number_cells())
            for i in range(0, mesh.number_cells()) :
                m = matmap[i]
                tmpP = mesh.volume(0) * \
                       (state.phi(0)[i] * mat.sigma_f(m, 0) + \
                        state.phi(1)[i] * mat.sigma_f(m, 1))
                mpower[i] = tmpP
                totP += tmpP
            totP *=  as_lra(mat).d_KAPPA


            P = totP / 17550.0
            if (conv) : 
                phi0.append(np.asarray(state.phi(0)))
                phi1.append(np.asarray(state.phi(1)))
                meshpower.append(mpower)
                meshtemp.append(1.0*T) # copy
                times.append(t)
                powers.append(P)
                temps.append(Tavg)
                maxtemps.append(Tmax)
                regP = np.asarray(rates.region_power("ASSEMBLY", 78.0))[0:78]
                #print len(regP), np.asarray(regP)
                asspower.append( np.asarray(regP) )
                asstemp.append(np.asarray(rates.edit("ASSEMBLY", as_lra(mat).physics().variable(0), True)) )
            if step == 0 and it  == 1:
                print "rank j|step it t P Pmax Tag Tmax"
            if step % 25 == 0:
                print "{} {} | {} {} {:.3e} {:.3e} {:.3e} {:.3e} {:.3e}".format(rank, j, step, it, t, P, max(powers), Tavg, Tmax)
        # end monitor
        ts.set_monitor(monitor)
        set_lra_physics(ts, as_lra(mat).physics(), mat)
        ts.solve(ic)

        # Save the data 
        data = {}
        data['keff']=ic.eigenvalue()
        data['T0'] = T0
        data['P0'] = P0
        data['SigA20'] = SigA20
        data['alpha'] = alpha
        data['gamma'] = gamma
        data['meshpower']=meshpower
        data['meshtemp']=meshtemp
        data['phi0']=phi0
        data['phi1']=phi1
        data['times'] = times
        data['power'] = powers
        data['temp'] = temps
        data['maxtemp'] = maxtemps
        data['asspower'] = asspower
        data['asstemp'] = asstemp

        pickle.dump(data, open(name+'.p', 'wb'), protocol=2)

if __name__ == "__main__":
    solvers.Manager.initialize(sys.argv)
    run(num_samples=50, option=0)
    solvers.Manager.finalize()
