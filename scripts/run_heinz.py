#!/usr/bin/env python

#Load required modules
import sys, os
import common
import logging
import subprocess

logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
#handler.setFormatter(NiceFormatter())
root = logging.getLogger()
root.addHandler(handler)
root.setLevel(logging.INFO)

params = open(snakemake.input["bum_params"], 'r')
l = params.readline().strip('\n\r')
a = params.readline().strip('\n\r')

FDR = snakemake.wildcards.FDR
experiment = snakemake.wildcards.experiment
network = snakemake.wildcards.network

fullnetwork = 'networks/' + snakemake.config['networks'][network]

for i in range(1):

    subprocess_heinz = ["heinz",
                        "-a", a,
                        "-lambda", l,
                        "-FDR", FDR,
                        "-n", snakemake.input.pvals,
                        "-e", fullnetwork,
                        #"-e", "networks/HPRD_Release9_062910_Heinz.txt",
                        #"-e", "networks/9606.protein.links.v10_mapped_geneIDs_700_no_UBC.txt",
                        "-o", "scores/{}_{}_{}_module.txt".format(experiment, network, FDR),
                        #"-o", "scores/{}_{}_{}_module_{}.txt".format(experiment, network, FDR, i),
                        #"-r", "TP53",
                        "-v", '0']
                        # ">", "scores/{}_{}_module.dot".format(experiment, FDR)]

    logger.info('Running Heinz on {} with FDR {} and network {}, computing module {} ...'.format(experiment, FDR, network, i))
    #logger.info('Running Heinz on {} with FDR={}...'.format(experiment, FDR))
    print(" ".join(subprocess_heinz))
    #dotfile = open("scores/{}_{}_{}_module_{}.dot".format(experiment, network, FDR, i), 'w')
    dotfile = open("scores/{}_{}_{}_module.dot".format(experiment, network, FDR), 'w')
    #errfile = open("/tmp/heinz.err", 'w')
    #print(" ".join(subprocess_heinz))
    p = subprocess.run(subprocess_heinz, stdout=dotfile)
    #subprocess.call("grep -v NaN scores/{}_{}_{}_module_{}.txt | grep -v \"#\" > scores/{}_{}_{}_module_{}.mod".format(experiment, network, FDR, i, experiment, network, FDR, i), shell=True)
    subprocess.call("grep -v NaN scores/{}_{}_{}_module.txt | grep -v \"#\" > scores/{}_{}_{}_module.mod".format(experiment, network, FDR, experiment, network, FDR), shell=True)
    #subprocess.call("neato -Tpdf scores/{}_{}_{}_module_{}.dot -o scores/{}_{}_{}_module_{}.pdf".format(experiment, network, FDR, i, experiment, network, FDR, i), shell=True)
    subprocess.call("neato -Tpdf scores/{}_{}_{}_module.dot -o scores/{}_{}_{}_module.pdf".format(experiment, network, FDR, experiment, network, FDR), shell=True)
