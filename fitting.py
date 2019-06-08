#!/usr/bin/env python

from ase import io
import lmfit
import numpy as np


name = {'2.0':'-0.482801',
        '2.2':'-0.503056',
        '2.4':'-0.397483',
        '2.5':'-0.325519',
        '2.6':'-0.245901',
        '2.8':'-0.176194',
        '3.0':'-0.371022',
        '3.2':'-0.310139',
        '3.4':'-0.095139',
        '3.6':'-0.5556',
        '3.8':'-0.03996',
        '4.0':'-0.03097',
        '4.2':'-0.032061',
        '4.4':'-0.017275',
        '4.5':'-0.019503',
        '5.0':'-0.004825'
        }
length = len(name)

Vreal = np.zeros(length)

p = lmfit.Parameters()
p.add_many(('D0_O', 0.51899099, True, 0.4, 0.7),('D0_H',0.02683872, True, 0.000, 0.04),('alpha_O',2.710708,True, 2.3, 3.0),('alpha_H', 0.94496310, True, 0.8, 1.2),('r0_O',2.03852937, True, 1.9, 2.4), ('r0_H',3.36387944, True, 3.2, 3.8), ('R0_O',3.2,True, 2.5, 3.2), ('R0_H',0,False), ('beta_O', 1.23675076, True, 0.9, 1.3), ('beta_H',0, False), ('gamma_O', 0.02196, True,0.01, 0.04), ('gamma_H',0, False))

def residul_func(p):
    Vmodel =np.zeros(length)
    v = p.valuesdict()
    i = 0
    for n in name:
        atoms = io.read('./'+n+'/qn.traj',index=-1)
        Vreal[i] = name[n]

        atomO = atoms[-3]
        atomH1 = atoms[-1]
        atomH2 = atoms[-2]
        for atom in atoms:
            if atom.symbol == 'Ni':
                rho_O = np.sqrt((atom.position[0]- atomO.position[0])**2 + (atom.position[1]- atomO.position[1])**2)
                frho_O = np.exp(-v['gamma_O'] * rho_O * rho_O)
                rho_H1 = np.sqrt((atom.position[0]- atomH1.position[0])**2 + (atom.position[1]- atomH1.position[1])**2)
                frho_H1 = np.exp(-v['gamma_H'] * rho_H1 * rho_H1)
                rho_H2 = np.sqrt((atom.position[0]- atomH2.position[0])**2 + (atom.position[1]- atomH2.position[1])**2)
                frho_H2 = np.exp(-v['gamma_H'] * rho_H2 * rho_H2)
                rO = np.linalg.norm(atom.position - atomO.position)
                rH1 = np.linalg.norm(atom.position - atomH1.position)
                rH2 = np.linalg.norm(atom.position - atomH2.position)
                
                Vmodel[i] += v['D0_O']*(np.exp(-2*v['alpha_O']*(rO -v['r0_O'])) - 2*np.exp(-v['alpha_O']*(rO -v['r0_O']))) * frho_O + v['R0_O'] * (np.exp(-v['beta_O']*(rO - v['r0_O'])))*(1-frho_O)
                Vmodel[i] += v['D0_H']*(np.exp(-2*v['alpha_H']*(rH1-v['r0_H'])) - 2*np.exp(-v['alpha_H']*(rH1-v['r0_H']))) * frho_H1 + v['R0_H'] * (np.exp(-v['beta_H']*(rH1-v['r0_H'])))*(1-frho_H1)
                Vmodel[i] += v['D0_H']*(np.exp(-2*v['alpha_H']*(rH2-v['r0_H'])) - 2*np.exp(-v['alpha_H']*(rH2-v['r0_H']))) * frho_H2 + v['R0_H'] * (np.exp(-v['beta_H']*(rH2-v['r0_H'])))*(1-frho_H2)
        i = i + 1
    print Vreal - Vmodel
    print np.sum((Vreal-Vmodel)**2)
    return np.abs(Vreal - Vmodel)
    
mi = lmfit.minimize(residul_func, p, method = 'cg')
lmfit.printfuncs.report_fit(mi.params, min_correl=.5)
mi.params.pretty_print()
