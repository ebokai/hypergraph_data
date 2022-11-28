"""
SPIN TOOLS

------------------------------
author: Ebo Peerbooms
contact: e.peerbooms@uva.nl 
------------------------------

algo_sbm.py 

hypergraph stochastic block model algorithms 

"""

import numpy as np 
from itertools import combinations, permutations
from scipy.special import gamma

def n_edges(n, d):

	return gamma(n + 1) / gamma(d + 1) / gamma(n - d + 1)

def hgsbm(nc, npc, p, q, d):

	n = nc * npc 

	print(f'generating hypergraph with {n} nodes, {nc} communities of {npc} nodes')

	n_tot = n_edges(n, d)
	n_in = nc * n_edges(npc, d)
	n_out = n_tot - n_in 
	n_abc = npc**3 * n_edges(nc, d)
	in_edges = []
	out_edges = []

	deg_dist = np.zeros((n,2))

	print(f'attempting to generate {int(n_in*p)} internal links of degree {d}')

	while (len(in_edges) < n_in * p):

		node = np.random.randint(n)

		cn = node // npc 

		edge = sorted(list(np.random.choice(np.arange(cn * npc, (cn+1)*npc), d, replace = False)))

		if edge not in in_edges:
			in_edges.append(edge)
			deg_dist[edge,0] += 1
	
	
	print(f'attempting to generate {int(n_abc*q)} [distinct] external links of degree {d}')

	while (len(out_edges) < n_abc * q):

		c = np.random.choice(np.arange(nc), d, replace = False)

		edge = c * npc + np.random.choice(np.arange(npc), d)

		edge = sorted(list(edge))

		if edge not in out_edges:
			out_edges.append(edge)
			deg_dist[edge,1] += 1

	mean_deg = np.mean(deg_dist, axis=0)
	std_deg = np.std(deg_dist, axis=0)

	mixing = deg_dist[:,1]/(deg_dist[:,0]+deg_dist[:,1])



	print(f'average internal degree: {mean_deg[0]:.1f} +/- {std_deg[0]:.1f}')
	print(f'average external degree: {mean_deg[1]:.1f} +/- {std_deg[1]:.1f}')
	print(f'average mixing: {np.mean(mixing):.2f} +/- {np.std(mixing):.2f}')

	return in_edges, out_edges, deg_dist, np.mean(mixing)








		

		



