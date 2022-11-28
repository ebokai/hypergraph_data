import numpy as np 
import math
import os 
import subprocess 
from algo_sbm import *
from tools import * 
import networkx as nx
import networkx.algorithms.community as nx_comm


tuples = [(i,j) for i in range(n) for j in range(i+1,n)]
ops = [2**i + 2**j for (i,j) in tuples]
stups = sort_by(ops, tuples)

nc = 4; npc = 5; d = 3; N = 1000;
n = nc * npc; n_in = nc * math.comb(npc, d); n_abc = npc**3 * math.comb(nc, d);

for k_avg in [6]:

	for mu in np.linspace(0,1,11):

		# CALCULATE CONNECTION PROBABILITIES
		p = n * k_avg * (1 - mu) / d / n_in
		q = n * k_avg * mu / d / n_abc

		if p > 1 or q > 1:
			raise ValueError(f'p: {p}, q: {q}')

		# GENERATE 20 GRAPHS FOR THIS MU

		for graph in range(20):

			gn = str(graph).zfill(2)

			in_edges, out_edges, deg_dist, mu_real = hgsbm(nc, npc, p, q, d)
			in_ops = [np.sum(np.power(2,edge)) for edge in in_edges]
			out_ops = [np.sum(np.power(2,edge)) for edge in out_edges]
			all_ops = np.array(sorted(in_ops + out_ops))

			model = all_ops

			np.savetxt(f'./graphs/HG{gn}_k{k_avg}_mu{mu:.2f}_({n},{nc},{npc})_degrees.dat', deg_dist, fmt='%s')
			np.savetxt(f'./graphs/HG{gn}_k{k_avg}_mu{mu:.2f}_({n},{nc},{npc})_interactions.dat', all_ops, fmt='%s')
			
			for beta in [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75]:

				pars = [beta for _ in range(len(model))]

				interactions = np.zeros((len(model),2))
				interactions[:,0] = model
				interactions[:,1] = pars
				np.savetxt(f'./temp/interactions_{k_avg}_{beta}.dat', interactions, delimiter=';', fmt=['%d','%.5f'])


				for dataset in range(20):

					metro_args = (f'./bin/metropolis_data_n{n}.out', f'{N}', f'./temp/interactions_{k_avg}_{beta}.dat', f'./temp/temp_data_{k_avg}')

					metro = subprocess.Popen(metro_args)
					metro.communicate()

					data_str = np.loadtxt(f'./temp/temp_data_{k_avg}.dat', dtype=str)
					data = [int(s,2) for s in data_str]
					data_x = np.array([[1 - 2 * int(s) for s in state] for state in data_str])

					cij = np.corrcoef(data_x, rowvar=False)				
					os.remove(f'./temp/temp_data_{k_avg}.dat')

					fname = f'HG{gn}_k{k_avg}_mu{mu:.2f}_({n},{nc},{npc})_B{beta}_N{N}_{str(dataset).zfill(2)}'

					np.savetxt(f'./data/{fname}.dat', data_str, fmt='%s')

					rise_args = ('./bin/main_RISE.out', f'{fname}')
					rise = subprocess.Popen(rise_args)
					rise.communicate()

					mcm_args = ('./bin/mcm_greedy.out', f'{fname}')

					mcm = subprocess.Popen(mcm_args)
					mcm.communicate()					

					jij = np.loadtxt(f'./jij/{fname}_jij_RISE.dat', delimiter='\t')
					
					k = 0
					GJ = nx.Graph()
					for k in range(190):
						i,j = stups[k]
						GJ.add_edge(i,j,weight=jij[k,1])
	
					comms_jij = nx_comm.louvain_communities(GJ)

					np.savetxt(f'./louvain/{fname}_louvain_comms_jij.dat', comms_jij, fmt='%s')

					GC = nx.Graph()
					for i in range(n):
						for j in range(i+1,n):
							GC.add_edge(n-1-i,n-1-j,weight=cij[i,j])

					comms_cij = nx_comm.louvain_communities(GC)

					np.savetxt(f'./louvain/{fname}_louvain_comms_cij.dat', comms_cij, fmt='%s')

				os.remove(f'./temp/interactions_{k_avg}_{beta}.dat')













tuples = [(i,j) for i in range(n) for j in range(i+1,n)]
ops = [2**i + 2**j for (i,j) in tuples]
stups = sort_by(ops, tuples)

k = 0
GJ = nx.Graph()
for k in range(190):
	i,j = stups[k]
	GJ.add_edge(i,j,weight=jij[k,1])