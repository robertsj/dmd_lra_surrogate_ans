"""
This script is a driver for the LRA benchmark problem.  Currently,
the LRA problem is built into Detran.  
"""

from detran import *

from lra_mesh import get_mesh_2D
from lra_uq_db import get_input

import sys 
import time
import cPickle as pickle
import numpy as np


def run(fmm) :


    # Input, mesh, and base material
    inp = get_input()
   # fmm = 15 # fine mesh per node
    mesh, cmm = get_mesh_2D(fmm)
    name_base = "diffusion{0}x{0}".format(fmm)
    transport = False
    mat = LRA.Create(mesh, transport, False)
        
    name = name_base# + '_{}'.format(j)
        
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
    P0 = 1e-6 # W/cm^3
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

    #-------------------------------------------------------------------------#
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
        P = totP / area
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
            asstemp.append(np.asarray(rates.edit("ASSEMBLY", 
                                                 as_lra(mat).physics().variable(0), True)) )
        if step == 0 and it  == 1:
            print "rank j|step it t P Pmax Tag Tmax"
        if step % 25 == 0:
            tmpl = " {} {} {:.3e} {:.3e} {:.3e} {:.3e} {:.3e}"
            print tmpl.format(step, it, t, P, max(powers), Tavg, Tmax)      
    # end monitor
    
    ts.set_monitor(monitor)
    set_lra_physics(ts, as_lra(mat).physics(), mat)
    ts.solve(ic)

    # Save the data 
    data = {}
    data['keff']=ic.eigenvalue()
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
    run(fmm=2)
    solvers.Manager.finalize()
