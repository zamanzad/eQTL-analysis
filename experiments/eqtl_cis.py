import sys
import limix
sys.path.append('./..')
from CFG.settings import *
import include.data as DATA
from include.utils import smartDumpDictHdf5
from include.utils import dumpDictHdf5
from include.utils import getLambda
from include.preprocess import rankStandardizeNormal
import include.fdr as FDR
import include.qtl as QTL
from limix.qtl import st_scan
import limix.qtl as QTL
from limix.stats import linear_kinship
from include.lmm import test_qtl_lmm

import numpy as np
import scipy as SP
import scipy.linalg as LA
import os
import _pickle as cPickle
import pdb
import time
import h5py

if __name__ == '__main__':

    out_dir = '/Users/fatemehzamanzad/Desktop/limix/out/runs'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if 'debug' in sys.argv:
        nfolds = 1
        fold_j = 0
        out_dir = os.path.join(out_dir,'..')
        fname = 'debug.hdf5'
        pdb.set_trace()
    else:
        nfolds = 100
        fold_j = 10
        fname  = '%d_%.3d.hdf5'%(nfolds,fold_j)

    # load data and split in jobs
	data  = DATA.data()
	#Kpop  = data.getKpopExpr(normalize=True) 
	Kpop = None
    cov   = data.getCovariates()
    cov = np.reshape(cov, (-1, 1)).shape
	genes = data.getGeneIDs()
	n_genes = genes.shape[1]
	Icv = SP.floor(nfolds*SP.arange(n_genes)/n_genes)
	I = Icv==fold_j
	genes_t=np.transpose(genes)
	genes = list(genes_t[I])

    # create output file
    out_file = os.path.join(out_dir,fname)
    fout = h5py.File(out_file,'w')

    for gene in genes:
        
        print (".. protein %s"%gene)

        #1. get geno and pheno data
        Y = data.getGeneExpression(gene,standardize=True)
        #Yc = data.getGeneExpression('P41235',standardize=True)
        Y = rankStandardizeNormal(Y)
        #Yc = rankStandardizeNormal(Yc)
        try:
            Xc,geno_info = data.getGermlineExpr(gene,cis_window=1e6)
        except:
            continue		
        qtl = limix.qtl.st_scan(Xc, Y, 'normal', K=None, verbose=False)
		pv = qtl.variant_pvalues
		pv = pv.sortby(pv).to_dataframe()
		pv["-log10(pv)"] = -np.log10(pv["pv"])
		print(pv.head())
	    print(qtl.variant_effsizes.sel(candidate=pv.index).to_dataframe().head())		
		print(qtl.variant_effsizes.sel(candidate=pv.index).to_dataframe())
		beta = qtl.variant_effsizes.sel(candidate=pv.index).to_dataframe()


        RV = {}
        pv = np.array(pv)
        RV['pv'] = pv
        RV['qv'] = FDR.qvalues(pv,m = None, return_pi0 = False, lowmem = False, pi0 = None, fix_lambda = None)
        #RV['lambda'] = getLambda(pv)
        RV['beta'] = beta

        # add gene info
        for key in geno_info.keys():
            RV[key] = geno_info[key]

        # export
        gene_group = fout.create_group("gene")
        dumpDictHdf5(RV,gene_group)

    fout.close()
