import numpy as np 
import math
import os 
import subprocess 
from algo_sbm import *


N = 1000 # number of datapoints to generate

# COMMUNITY PROPERTIES
nc = 4 # number of communities
npc = 5 # nodes per community
d = 3 # hyper-edge degree
k_avg = 6 # average node degree
mu = 0.00 # mixing parameter
beta = 0.55 # inverse temperature: controls noise level

n = nc * npc; n_in = nc * math.comb(npc, d); n_abc = npc**3 * math.comb(nc, d)

# CALCULATE CONNECTION PROBABILITIES
p = n * k_avg * (1 - mu) / d / n_in
q = n * k_avg * mu / d / n_abc

if p > 1 or q > 1:
	raise ValueError(f'p: {p}, q: {q}')

# GENERATE HYPER-GRAPH
in_edges, out_edges, deg_dist, mu_real = hgsbm(nc, npc, p, q, d) # generate graph
in_ops = [np.sum(np.power(2,edge)) for edge in in_edges] # convert internal edges to spin operators
out_ops = [np.sum(np.power(2,edge)) for edge in out_edges] # convert external edges to spin operators
all_ops = np.array(sorted(in_ops + out_ops)) # combine into single model
model = all_ops

# SAVE GRAPH PROPERTIES
g_id = str(0).zfill(2) # specify identifier (if generating multiple graphs)
np.savetxt(f'./graphs/HG{g_id}_k{k_avg}_mu{mu:.2f}_({n},{nc},{npc})_degrees.dat', deg_dist, fmt='%s')
np.savetxt(f'./graphs/HG{g_id}_k{k_avg}_mu{mu:.2f}_({n},{nc},{npc})_interactions.dat', all_ops, fmt='%s')

# CONVERT INTERACTIONS TO METROPOLIS-READABLE FILE
pars = [beta for _ in range(len(model))] # initialize parameters
interactions = np.zeros((len(model),2))
interactions[:,0] = model
interactions[:,1] = pars
np.savetxt(f'./temp/interactions_{k_avg}_{beta}.dat', interactions, delimiter=';', fmt=['%d','%.5f'])

fname = f'HG{g_id}_k{k_avg}_mu{mu:.2f}_({n},{nc},{npc})_B{beta}_N{N}'
metro = subprocess.Popen(f'./metropolis_data.exe {N} ./temp/interactions_{k_avg}_{beta}.dat ./generated_data/{fname}')
metro.communicate()

# remove temporary interaction file 
os.remove(f'./temp/interactions_{k_avg}_{beta}.dat')











