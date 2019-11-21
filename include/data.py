import sys
sys.path.append('./..')
from CFG.settings import *

from include.utils import matchIDs
from include.normalization import *

import scipy as SP
import numpy as np
import h5py
import pdb
import copy
import warnings

class data():
    def __init__(self):
        self.load()

    def load(self):
        """load data file:
        cache_genotype: load genotypes fully into memory (False)
        """
        #self.f  = h5py.File(CFG['data']['data_file'],'r')
        # self.fp = self.f['drugRes']
        # self.fs = self.f['somatic']
        #self.fc = self.f['covariate']
        # self.fc = h5py.File(CFG['data']['covariate_file'],'r')
        self.fe = h5py.File(CFG['data']['expr_file'],'r')
        # self.fg = h5py.File(CFG['data']['germ_file_drugRes'],'r')
        self.fge = h5py.File(CFG['data']['germ_file_expr'],'r')
        # self.drugID = self.fp['drug_ID'][:]
        self.geneID = self.fe['gex']['gex_name'][:]

    def getGeneIDs(self):
        """ get geneIDs """
        _chr = self.fe['gex']['chrom'][:]
        Iin = (_chr!='X') or (_chr!='X')or(_chr!='MT')
        rv = self.geneID[Iin]
        return rv 

    def getKpopExpr(self,Isample=None,normalize=True):
        """
        get Kpop for expression data 
        """
        RV = self.fge['KpopZ'][:]
        if Isample!=None:
            RV = RV[Isample,:][:,Isample]
        if normalize:
            RV/=RV.diagonal().mean()
        return RV

    def getKpopExprDel(self,Isample=None,normalize=True):
        """
        get Kpop for expression data
        """
        RV = self.fge['matrix'][:]
        RV = np.delete(RV, 6, 1)
        # standardize SNP to mean 0 and unit variance
        RV -= RV.mean(0)
        RV /= RV.std(0)
        K = SP.dot(RV,RV.T)
        if Isample!=None:
            K = K[Isample,:][:,Isample]
        if normalize:
            K/=K.diagonal().mean()
        return K

    def getGeneExpression(self,geneID,standardize=False):
        """
        Get gene expression levels
        """
        self.geneID = np.transpose(self.geneID)
        self.geneID = list(self.geneID)
        
        #idx = self.geneID==geneID
        idx=self.geneID.index(geneID)
        Y = self.fe['gex']['gex'][idx,:].T
        if standardize:
            Y-=Y.mean(0)
            Y/=Y.std(0)
        return Y

    def getGenePos(self,geneID):
        """
        get position of the gene
        """
        self.geneID = np.transpose(self.geneID)
        self.geneID = list(self.geneID)
        #idx = SP.where(self.geneID==geneID)[0][0]
        idx=self.geneID.index(geneID)
        gene_chrom = float(self.fe['gex']['chrom'][idx])
        gene_start = float(self.fe['gex']['gene_start_bp'][idx])
        gene_end = float(self.fe['gex']['gene_stop_bp'][idx])
        rv = SP.array([gene_chrom,gene_start,gene_end])
        return rv
    
    def getCovariates(self):
        """
        get covariates
        """
        rv = self.fc['covariate'][:].T
        col = 1.*(rv.sum(1)==0)[:,SP.newaxis]
        rv = SP.concatenate([rv,col],1)
        #rv = SP.concatenate([rv],1)
        return rv

    def getGermlineExpr(self,geneID,cis_window=1e6,standardize=False,Is=None,debug=False):
        """
        get genotypes, chrom, pos
        """
        genePos = self.getGenePos(geneID)
        pos = self.fge['col_header']['pos'][:]
        chrom = self.fge['col_header']['chrom'][:]
        # Icis  = (chrom==genePos[0])
        # Icis *= (pos>genePos[1]-cis_window)
        # Icis *= (pos<genePos[2]+cis_window)
        # assert Icis.sum()>0, 'nothing here bro'
        # X = self.fge['matrixT'][:]
        X = self.fge['matrixT'][:].T
        # X = np.delete(X, 6, 0)
        # X = X.T
        info = {}
        for key in self.fge['col_header'].keys():
            info[key] = self.fge['col_header'][key][:] 
            # info[key] = np.delete(info[key],50) 
        return X, info
