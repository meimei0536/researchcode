#!/usr/bin/env python

from amp.model.neuralnetwork import NeuralNetwork
from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model import LossFunction
from ase import io
from ase.io.trajectory import Trajectory

def make_symmetry_functions(elements):
    G = {}
    for element0 in elements:
        #etas = [0.05, 4., 20., 80.]                                                                               
        etas = [0.003214*42.25, 0.035711*42.25, 0.071421*42.25, 0.124987*42.25, 0.214264*42.25, 0.357106*42.25, 0.714213*42.25, 1.428426*42.25]
        _G = [{'type': 'G2', 'element': element, 'eta': eta}
              for eta in etas
              for element in elements]

        etas = [0.000357*42.25, 0.028569*42.25, 0.089277*42.25] 
        #etas = [0.005]
        zetas = [1., 2., 4.]  
        #zetas = [1., 4.]
        gammas = [+1., -1.]
        for eta in etas:
            for zeta in zetas:
                for gamma in gammas:
                    for i1, el1 in enumerate(elements):
                        for el2 in elements[i1:]:
                            els = sorted([el1, el2])
                            _G.append({'type': 'G4',
                                       'elements': els,
                                       'eta': eta,
                                       'gamma': gamma,
                                       'zeta': zeta})
        G[element0] = _G
    return G


elements =['O','Ru']
Gs = make_symmetry_functions(elements)

filehandle = open('nodes','r')
linelist = filehandle.readlines()
filehandle.close()
cores = {}
for i in range(len(linelist)):
    node =  linelist[i][0:5]
    cores[node]=32

train_images = io.read('/work/common/hxin_lab/jiamin/non_adiabatic/Langevin/Training_3nd/trajs/fingerprints/identified.traj',index=':')

print len(train_images)
#calc = Amp(descriptor=Gaussian(Gs=Gs, elements=['O','Ru']), model=NeuralNetwork(hiddenlayers=(70,70)), cores = 32)
convergence = {'energy_rmse':0.001,
               'energy_maxresid':None,
               'force_rmse':0.005,
               'force_maxresid':None}
lossfunction = LossFunction(convergence=convergence)
calc = Amp.load('./amp-checkpoint.amp',cores = cores)
#calc=Amp.load('/work/common/hxin_lab/jiamin/non_adiabatic/Langevin/Training_2nd/sym70_2L70/amp-checkpoint.amp', cores=cores)
calc.model.lossfunction=lossfunction
calc.train(images = train_images,)

