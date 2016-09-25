#!/usr/bin/env python

import numpy
import os, sys
import pdb
import argparse

### PARSE OPTIONS, INCLUDING ALL PARAMS ###
def get_user_options():
	parser = argparse.ArgumentParser()
	parser.description='Pan-genome simulation, for t-SNE panini'
	parser.epilog="..."
	parser.add_argument('-n', required=False, action='store', dest="n", help="number of lineages", default=2, metavar="INT")
	parser.add_argument('-k', required=False, action='store', dest="k", help="(max) number of genes in genomes", default=20, metavar="INT")
	parser.add_argument('-s', required=False, action='store', dest="s", help="starting population size per lineage", default=5, metavar="INT")
	parser.add_argument('-t', required=False, action='store', dest="tau", help="number of generations", default=3, metavar="INT")
	parser.add_argument('-g', '--gain', required=False, action='store', dest="mu_gain", help="(percentage) of genes to gain each generation", default=10, metavar="FLOAT")
	parser.add_argument('-l', '--loss', required=False, action='store', dest="mu_loss", help="(percentage) of genes to gain each generation", default=10, metavar="FLOAT")
	parser.add_argument('-m', required=False, action='store', dest="m", help="(percentage) migration rate per generation", default=5, metavar="FLOAT")
	parser.add_argument('-o', required=False, action='store', dest="o", help="output file name (.Rtab) [default: print to screen]", default='STDOUT', metavar="STRING")

	return parser.parse_args()




def sample_gene_freq_params(k):
	"""returns ndarray of length k"""
	return numpy.random.beta(.5, .5, size=(k))

def generate_starting_genomes(K, s):
	"""
	returns a numpy.ndarray of size s * length(k)
	binary entries represent pres/abs of genes in genomes
	"""
	ret = numpy.zeros((s,len(K)), dtype=int)
	for x in xrange(ret.shape[0]): # loop through genomes
		for y in xrange(ret.shape[1]): # loop through genes
			ret[x][y] = numpy.random.binomial(1, K[y])
	return ret

def perturb_pop_size(pop_old):
	"""
	samples with replacement from pop (ndarray shape(s,x))
	returns new pop with size drawn from poission(s)
	"""
	s_old = len(pop_old)
	if s_old == 0:
		return pop_old
	s_new = numpy.random.poisson(lam=s_old, size=None)
	pop_new = numpy.empty((s_new,len(pop_old[0])), dtype=int)
	haplos_to_use = numpy.random.choice(range(0,s_old), size=s_new, replace=True)

	for i in xrange(0,s_new):
		pop_new[i] = numpy.copy(pop_old[haplos_to_use[i]])
	return pop_new

def mu_rate_to_int(k, mu_gain, mu_loss):
	return (
		int(numpy.floor(float(mu_gain)/100 * int(k))),
		int(numpy.floor(float(mu_loss)/100 * int(k)))
	)

def gain_lose_genes(pop, mu_gain, mu_loss, verbose):
	"""
	in: a population (nd array, s x #genes), gain rate, loss rate
	NB rates are integers (lambda for poisson)
	"""
	for i in xrange(0, len(pop)):
		tmp = sum(pop[i])
		tmpp = numpy.copy(pop[i])
		## the following indexes are relative to the haplotypes
		idxs_present = [x for x in xrange(0,len(pop[i])) if pop[i][x]]
		idxs_missing = [x for x in xrange(0,len(pop[i])) if not pop[i][x]]
		n_to_lose = numpy.random.poisson(lam=mu_loss, size=None)
		n_to_gain = numpy.random.poisson(lam=mu_gain, size=None)
		try:
			assert(n_to_lose < len(idxs_present))
			assert(n_to_gain < len(idxs_missing))
		except AssertionError:
			if verbose:
				print "\t\tSkipping genome #{} as genes to add / lose are more than avaiable!".format(i)
			return pop

		#### GAIN
		if n_to_gain: # i.e. not 0
			idxs_to_gain = numpy.random.choice(idxs_missing, size=n_to_gain, replace=False)
			for j in idxs_to_gain:
				pop[i][j] = 1

		#### LOSS
		if n_to_lose: # i.e. not 0
			idxs_to_drop = numpy.random.choice(idxs_present, size=n_to_lose, replace=False)
			for j in idxs_to_drop:
				pop[i][j] = 0

		if verbose:
			print "\t\tgenome #{} #genes {} -> {} (gain {} loss {})".format(i, tmp, sum(pop[i]), n_to_gain, n_to_lose)

	return pop

def migrate_genomes_from_i_to_j(genomes, i, j, m, verbose):
	"""
	this is a bit misleading, as the incoming genomes stay in the original lineage, and they also replace those in the new lineage... could be improved
	"""
	n_to_migrate = numpy.random.poisson(lam=m, size=None)
	try:
		assert(n_to_migrate < len(genomes[i]))
		assert(n_to_migrate < len(genomes[j]))
	except AssertionError:
		print "\tSkipping migration from {} to {}".format(i, j)
		return genomes
	idxs_to_migrate = numpy.random.choice(range(0,len(genomes[i])), size=n_to_migrate, replace=False)
	idxs_to_replace = numpy.random.choice(range(0,len(genomes[j])), size=n_to_migrate, replace=False)
	for k in xrange(0, n_to_migrate):
		genomes[j][idxs_to_replace[k]] = numpy.copy(genomes[i][idxs_to_migrate[k]])
	if verbose:
		print "\tmigrated {} haplotypes from lineage {} to {}".format(n_to_migrate, i, j)
	return genomes

def write_genomes_to_Rtab(fname, genomes):
	with open(fname,'w') as fh:
		header = ["gene"]
		for i in xrange(0, len(genomes)):
			header.append('\t'.join(['l{}g{}'.format(i,x) for x in xrange(0,len(genomes[i]))]))
		if fname == "STDOUT":
			print '\t'.join(header)
		else:
			fh.write('\t'.join(header) + "\n")

		for gene_idx in xrange(0,len(genomes[0][0])):
			line = ["gene{}".format(gene_idx)]
			for lineage in xrange(0, len(genomes)):
				for genome in xrange(0,len(genomes[lineage])):
					line.append(str(genomes[lineage][genome][gene_idx]))

			if fname == "STDOUT":
				print '\t'.join(line)
			else:
				fh.write('\t'.join(line) + "\n")

def simulate_genomes(n, k, s, tau, m, mu_gain, mu_loss, verbose):

	#### generate starting genomes
	genomes = [0] * int(n)
	for x in xrange(0,len(genomes)):
		K = sample_gene_freq_params(int(k))
		genomes[x] = generate_starting_genomes(K, int(s))

	##### simulate evolution
	if verbose:
		print "Simulating evolution...\n"
	for t in xrange(0, int(tau)):
		if verbose:
			print "\n*** Time step {} ***\n".format(t)

		### grow/shrink populations in each lineage:
		if verbose:
			print "Changing population sizes..."
		for lineage in xrange(0, len(genomes)):
			tmp = len(genomes[lineage])
			genomes[lineage] = perturb_pop_size(genomes[lineage])
			if verbose:
				print "\tlineage {} size {} -> {}".format(lineage, tmp, len(genomes[lineage]))
			### if population has crashed, throw?
			try:
				assert(len(genomes[lineage]) != 0)
			except AssertionError:
				if verbose:
					print "A lineage's population has crashed FATAL"
				### throw a general error here which can be caught elsewhere
				raise AssertionError

		### migrate genomes between lineages
		mrate = int(float(m) / 100 * int(k))
		if verbose:
			print "\nMigrating genomes (m_ij {})...".format(mrate)
		for i in xrange(0, len(genomes)):
			for j in [x for x in range(0, len(genomes)) if x != i]:
				migrate_genomes_from_i_to_j(genomes, i, j, mrate, verbose)

		### delete and gain genes within haplotypes
		rates = mu_rate_to_int(k, mu_gain, mu_loss)
		if verbose:
			print "\nLosing / Gaining genes (\lambda gain {} loss {})...".format(rates[0], rates[1])
		for lineage in xrange(0, len(genomes)):
			if verbose:
				print "\tlineage {}".format(lineage)
			genomes[lineage] = gain_lose_genes(genomes[lineage], rates[0], rates[1], verbose)

	return genomes



if __name__ == "__main__":
	options = get_user_options()

	genomes = simulate_genomes(options.n, options.k, options.s, options.tau, options.m, options.mu_gain, options.mu_loss, True)

	write_genomes_to_Rtab(options.o, genomes)

