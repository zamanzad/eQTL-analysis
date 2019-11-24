import os

CFG = {}

# data files
CFG['data'] = {}
CFG['data']['base_folder'] ='/Users/fatemehzamanzad/Desktop/limix/data'
CFG['data']['gpfs_folder'] = '/Users/fatemehzamanzad/Desktop/limix/data'
CFG['data']['expr_file']   = os.path.join(CFG['data']['base_folder'],'proteins.h5')
CFG['data']['germ_file_expr']  = os.path.join(CFG['data']['base_folder'],'somatic.h5')
CFG['data']['covariate_file'] = os.path.join(CFG['data']['base_folder'],'covariatenew.h5')

#settings
CFG['settings'] = {}
CFG['settings']['ld_cutoff'] = 0.4

# chrom length
CFG['length'] = {}
CFG['length']['chrom1'] = 249250621
CFG['length']['chrom2'] = 243199373
CFG['length']['chrom3'] = 198022430
CFG['length']['chrom4'] = 191154276
CFG['length']['chrom5'] = 180915260
CFG['length']['chrom6'] = 171115067
CFG['length']['chrom7'] = 159138663
CFG['length']['chrom8'] = 146364022
CFG['length']['chrom9'] = 141213431
CFG['length']['chrom10'] = 135534747
CFG['length']['chrom11'] = 135006516
CFG['length']['chrom12'] = 133851895
CFG['length']['chrom13'] = 115169878
CFG['length']['chrom14'] = 107349540
CFG['length']['chrom15'] = 102531392
CFG['length']['chrom16'] = 90354753
CFG['length']['chrom17'] = 81195210
CFG['length']['chrom18'] = 78077248
CFG['length']['chrom19'] = 59128983
CFG['length']['chrom20'] = 63025520
CFG['length']['chrom21'] = 48129895
CFG['length']['chrom22'] = 51304566