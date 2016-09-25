#!/usr/bin/env python

import numpy
import os, sys
import pdb
import argparse
import simulateGenomes

"""
main call function: simulate_genomes. Params: n, k, s, tau, m, mu_gain, mu_loss, verbose

this script simulates a bunch of genomes according to various parameters such that these can be run
through PANINI and then, using a seperate RScript, turned into the points on a graph

i.e. if we want to compare the effect of increasing migration rates, we can set n, k, s, tau, mu as constants and
run this over a series of m values. The parameters are recorded in the file name as so:
genomes.n.k.s.tau.m.mu_gain.mu_loss.repeat.Rtab

The membership of genomes to the original lineages is regorded in the genome name:
name: lxgy -> x: lineage number, y: genome number (within lineage)

"""

def explore_pan_genome_size():
	#### explore the effect of pan genome size across a number of lineages, with no migration and c. 100 members / lineage
	nvec = [ 2, 3]
	m = 0
	s = 100
	kvec = [ 100, 1000 ]
	mu_gain, mu_loss = 2, 2
	rvec = [ 0 ] ## number of repeats
	tau = 100
	experiment = "panGenomeSize"
	for n in nvec:
		for k in kvec:
			for r in rvec:
				fname = "genomes/{}.{}.{}.{}.{}.{}.{}.{}.{}.Rtab".format(experiment, n, k, s, tau, m, mu_gain, mu_loss, r)
				print fname

				while True:
					try:
						genomes = simulateGenomes.simulate_genomes(n, k, s, tau, m, mu_gain, mu_loss, False)
					except AssertionError:
						pass # run through the while loop again
					else:
						break

				simulateGenomes.write_genomes_to_Rtab(fname, genomes)


if __name__ == "__main__":
	explore_pan_genome_size()




	# genomes = simulateGenomes.simulate_genomes(3, 200, 50, 10, 0, 5, 5, True)
	# simulateGenomes.write_genomes_to_Rtab('STDOUT', genomes)

