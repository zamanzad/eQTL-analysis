import sys
sys.path.append('./..')
from CFG.settings import *
#sys.path.insert(0,CFG['limix_folder'])
import include.fdr as FDR
import include.data as DATA
from include.utils import smartAppend
from include.utils import smartDumpDictHdf5

import scipy as SP
import h5py
import os
import pdb
import glob
import _pickle as cPickle

fdr_thr = 0.01

def getRelPos(pos,gene_pos):
    if gene_pos[3]==1:
        rv = pos-gene_pos[1]
    else:
        rv = gene_pos[2]-pos
    return rv

data = DATA.data()

if __name__=='__main__':

    pdb.set_trace()

    out_dir = '/Users/fatemehzamanzad/Desktop/limix/out'
    out_file = os.path.join(out_dir,'summary.hdf5')

    if not os.path.exists(out_file) or 'recalc' in sys.argv:
        
        table = {}
    
        fname = os.path.join(out_dir,'runs','*.hdf5')
        files = glob.glob(fname)
        files = SP.sort(files)
        
        first = True
        for file in files:

            print file
            try:
                f = h5py.File(file,'r')
            except:
                print 'file corrupted'
                continue
            geneIDs = f.keys()
            for geneID in geneIDs:

                fgene = f[geneID]
                try:
                    temp = {}
                    temp['geneID'] = SP.array([str(geneID)])
                    temp['file']   = SP.array([str(file)])

                    # linear mixed model
                    idx = fgene['pv'][0,:].argmin()
                    temp['pv'] = fgene['pv'][:,idx]
                    temp['qv'] = fgene['qv'][:,idx]
                    temp['lambda'] = fgene['lambda'][:,0]
                    temp['beta'] = fgene['beta'][:,idx]
                    temp['qv_total'] = fgene['qv'][:]
                    temp['pv_total'] = fgene['pv'][:]
                    temp['beta_total'] = fgene['beta'][:]
                    # position info
                    gene_pos = data.getGenePos(geneID)
                    pos = fgene['pos'][[idx]]
                    temp['gene_pos'] = gene_pos
                    temp['relpos'] = pos-gene_pos[1:].mean() 
                    temp['pos'] = pos
                    temp['rs'] = fgene['rs'][[idx]]
                    temp['chrom'] = fgene['chrom'][[idx]]
                except:
                    print "poppo"

                #append the temp table into the big table
                for key in temp.keys():
                    smartAppend(table,key,temp[key])

            f.close()


        for key in table.keys():
            table[key] = SP.concatenate(table[key])

        # add corrected qvalues also across genes
        table['qv_all'] = FDR.qvalues(table['qv'])
        fout = h5py.File(out_file,'w')
        smartDumpDictHdf5(table,fout)
        fout.close()

    else:
        f = h5py.File(fout_name,'r')
        table = {}
        for key in f.keys():
            table[key] = f[key][:]
        f.close()
