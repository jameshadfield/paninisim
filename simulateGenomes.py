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
	parser.add_argument('-k', required=False, action='store', dest="k", help="(max) number of genes in genomes", default=10, metavar="INT")
	parser.add_argument('-s', required=False, action='store', dest="s", help="starting population size per lineage", default=2, metavar="INT")
	return parser.parse_args()


if __name__ == "__main__":
	options = get_user_options()
	pdb.set_trace()
