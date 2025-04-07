#!/usr/bin/env python

import os
import argparse

import sys
sys.path.append('./conf')

use_numpy = True

parser = argparse.ArgumentParser(description='Generate LCIO::MCParticles with specified parameters')
parser.add_argument('output', metavar='FILE_OUT.slcio', help='Output LCIO file')
parser.add_argument('-s', '--seed', metavar='seed', type=int, help='Seed to use for random generator', default=12345)
parser.add_argument('-c', '--comment', metavar='TEXT',  help='Comment to be added to the run header', type=str)
parser.add_argument('-e', '--events', metavar='N', type=int, default=1,  help='# of events to generate (default: 1)')
parser.add_argument('-p', '--particles', metavar='N', type=int, default=1,  help='# of particles/event to generate (default: 1)')
parser.add_argument('-o', '--overwrite', action='store_true',  help='Overwrite existing output file')
parser.add_argument('--pdg', metavar='ID', type=int, default=[13], nargs='+',  help='PdgIds of the allowed particles (default: [13])')
parser.add_argument('--pt', metavar='V', type=float, default=1,  help='Minimum ranverse momentum [GeV] (default: 1)')
parser.add_argument('--theta', metavar='A', type=float, default=[10,170], nargs='+',  help='Polar angle [deg] (default: 10 170)')
parser.add_argument('--phi', metavar='A', type=float, default=[0,360], nargs='+',  help='Azimuthal angle [deg] (default: 0 360)')
parser.add_argument('--dz', metavar='V', type=float, default=0,  help='Beam spread along Z [mm] (default: 0)')
parser.add_argument('--dt', metavar='V', type=float, default=0,  help='Time spread [ns] (default: 0)')
#parser.add_argument('--dz', metavar='V', type=float, nargs='*', default=0,  help='Vertex position along Z [mm] (default: 0)')
#parser.add_argument('--d0', metavar='V', type=float, nargs='*', default=0,  help='Vertex position along R [mm] (default: 0)')
 
args = parser.parse_args()


from pyLCIO import UTIL, EVENT, IMPL, IO, IOIMPL
from pdgs import PDG_PROPS
from array import array
if use_numpy:
    import numpy as np
else:
    import random as rng
import math

# Validating the arguments
if not args.overwrite and os.path.isfile(args.output):
    raise FileExistsError(f'Output file already exists: {args.output:s}')
for pdg in args.pdg:
    if pdg not in PDG_PROPS:
        raise RuntimeError(f'Particle properties not defined for pdgId: {pdg}')

if use_numpy:
    rng = np.random.default_rng(args.seed)
else:
    rng.seed(args.seed)
print(f'Generation seed = {args.seed}')


configs = {
    'invpt': [-1./args.pt, 1./args.pt],
    'dt': args.dt,
    'dz': args.dz,
    #'d0': args.d0,
    # Converting degrees to radians
    'theta': [math.radians(a) for a in args.theta],
    'phi': [math.radians(a) for a in args.phi]
}


# Opening the output file
wrt = IOIMPL.LCFactory.getInstance().createLCWriter( )
wrt.open( args.output , EVENT.LCIO.WRITE_NEW )
print(f'Opening output file: {args.output}')

# Writing the run headers
run = IMPL.LCRunHeaderImpl()
run.setRunNumber(0)
run.parameters().setValue('pdgIds', str(args.pdg))
run.parameters().setValue('events', args.events)
run.parameters().setValue('particles/event', args.particles)
if args.comment:
    run.parameters().setValue('comment', args.comment)
for name, values in configs.items():
    header = str(values) if isinstance(values, list) else values
    run.parameters().setValue(name, header)
wrt.writeRunHeader(run)

# Setting counters
n_events = 0
n_particles = 0

# Choosing pdgId of each particle randomly if # of pdgIds is different from # of particles/event
n_pdgs = len(args.pdg)
choose_random_pdg = True if args.particles != n_pdgs else False


# Generate events
for ievent in range(args.events):

    col = IMPL.LCCollectionVec(EVENT.LCIO.MCPARTICLE)
    evt = IMPL.LCEventImpl()
    evt.setEventNumber(ievent)
    evt.setRunNumber(run.getRunNumber())
    evt.addCollection(col, "MCParticle")
    
    # Generate particles
    for ipart in range(args.particles):

        # particle momentum
        invpt = rng.uniform(configs['invpt'][0],configs['invpt'][1])
        phi = rng.uniform(configs['phi'][0],configs['phi'][1])

        eta_min = -math.log(math.tan(0.5*configs['theta'][1]))
        eta_max = -math.log(math.tan(0.5*configs['theta'][0]))
        eta = rng.uniform(eta_min, eta_max)
        theta = 2.*math.atan(math.exp(-eta))

        px = math.cos(phi) / math.fabs(invpt)
        py = math.sin(phi) / math.fabs(invpt)
        pz = 1. / ( math.tan(theta)*math.fabs(invpt) )

        momentum = array('f', [px, py, pz])   

        # particle production vertex
        vx = 0.
        vy = 0.
        vz = rng.normal(0., configs['dz']) if use_numpy else rng.gauss(0., configs['dz'])

        vertex = array('d',[ vx, vy, vz ] )

        # particle production time
        time_prod = rng.normal(0., configs['dt']) if use_numpy else rng.gauss(0., configs['dt'])

        
        # particle pdg code
        pdg_idx = ipart
        if choose_random_pdg:
            pdg_idx = rng.choice(n_pdgs, 1)[0]
            pdg = args.pdg[pdg_idx]

        if invpt < 0:
            pdg = -pdg

        
        # create the MCParticle
        mcp = IMPL.MCParticleImpl()
        mcp.setGeneratorStatus(1)
        mcp.setMass(PDG_PROPS[pdg][1])
        mcp.setCharge(PDG_PROPS[pdg][0])
        mcp.setPDG(pdg)
        mcp.setMomentum(momentum)
        mcp.setVertex(vertex)
        mcp.setTime(time_prod)

        # add the particle to the event
        col.addElement(mcp)
        n_particles += 1

        #print (ipart, pdg, PDG_PROPS[pdg][0], px, py, pz, vx, vy, vz)

    # Writing the event
    wrt.writeEvent(evt)
    n_events += 1
    if n_events % (args.events / 10) == 0:
        print(f'Wrote event {n_events}/{args.events}')

# Closing the output file
wrt.close()
print(f'Wrote {n_particles} partiles in {n_events} events to file: {args.output}')
