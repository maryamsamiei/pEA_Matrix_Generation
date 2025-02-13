#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 10:50:34 2025

@author: maryamsamieinasab
"""


#!/usr/bin/env python
import re
from itertools import product
import numpy as np
import pandas as pd
from pysam import VariantFile
import argparse
from joblib import Parallel, delayed
from tqdm import tqdm




def validate_EA(ea):
    """
    Checks for valid EA score
    Args:
        ea (str/float/None): EA score as string
    Returns:
        float: EA score between 0-100 if valid, otherwise returns NaN
    """
    try:
        ea = float(ea)
    except ValueError:
        if type(ea) == str and (ea == 'fs-indel' or 'STOP' in ea):
            ea = 100
        else:
            ea = np.nan
    except TypeError:
        ea = np.nan
    return ea


# def pEA(dmatrix, ea, gts, ft_name):

#     dmatrix[ft_name] *= (1 - ea/100) ** gts


def af_check(rec, af_field='AF', max_af=None, min_af=None):
    """
    Check if variant allele frequency passes filters
    Args:
        rec (VariantRecord)
        af_field (str): Name of INFO field containing allele frequency information
        max_af (float): Maximum allele frequency for variant
        min_af (float): Minimum allele frequency for variant
    Returns:
        bool: True of AF passes filters, otherwise False
    """
    if max_af is None and min_af is None:
        return True
    elif max_af is None:
        max_af = 1
    elif min_af is None:
        min_af = 0
    af = rec.info[af_field]
    if type(af) == tuple:
        af = af[0]
    return min_af < af < max_af


def convert_zygo(genotype):
    """
    Convert a genotype tuple to a zygosity integer
    Args:
        genotype (tuple): The genotype of a variant for a sample
    Returns:
        int: The zygosity of the variant (0/1/2)
    """
    if genotype in [(1, 0), (0, 1)]:
        zygo = 1
    elif genotype == (1, 1):
        zygo = 2
    else:
        zygo = 0
    return zygo
### functions for VEP annottated vcf
def fetch_EA_VEP(EA, canon_ensp, all_ensp, csq, EA_parser='canonical'):
    if 'stop_gained' in csq or 'frameshift_variant' in csq or 'stop_lost' in csq or 'splice_donor_variant' in csq or 'splice_acceptor_variant' in csq:
        return 100
    if EA_parser == 'canonical':
        try:
            canon_idx = all_ensp.index(canon_ensp)
        except ValueError:
            return np.nan
        else:
            return validate_EA(EA[canon_idx])
    else:
        newEA = []
        for score in EA:
            newEA.append(validate_EA(score))
        if np.isnan(newEA).all():
            return np.nan
        elif EA_parser == 'mean':
            return np.nanmean(newEA)
        elif EA_parser == 'max':
            return np.nanmax(newEA)
        else:
            return newEA

### pEA/sumEA matrix for wavelet and EPIMUTESTR and Sigma diff pipelines
def parse_VEP(vcf_fn, gene, gene_ref, samples, method, min_af=None, max_af=None, af_field='AF', EA_parser='canonical'):
    vcf = VariantFile(vcf_fn)
    vcf.subset_samples(samples)
    dmatrix = pd.DataFrame(np.ones((len(samples), 1)), index=samples, columns=[gene])
    for var in vcf:
        if re.search(r'chr', var.chrom):
            contig = 'chr'+str(gene_ref.chrom)
        else:
            contig = str(gene_ref.chrom)
        break
    def _fetch_anno(anno):
        # for fields that could return either direct value or tuple depending on header
        if type(anno) == tuple:
            return anno[0]
        else:
            return anno
    for rec in vcf.fetch(contig=contig, start=gene_ref.start, stop=gene_ref.end):
        all_ea = rec.info.get('EA', (None,))
        all_ensp = rec.info.get('Ensembl_proteinid', (rec.info['ENSP'][0],))
        canon_ensp = _fetch_anno(rec.info['ENSP'])
        csq = _fetch_anno(rec.info['Consequence'])
        rec_gene = _fetch_anno(rec.info['SYMBOL'])
        ea = fetch_EA_VEP(all_ea, canon_ensp, all_ensp, csq, EA_parser=EA_parser)
        pass_af_check = af_check(rec, af_field=af_field, max_af=max_af, min_af=min_af)
        if not np.isnan(ea).all() and gene == rec_gene and pass_af_check:
            gts = pd.Series([convert_zygo(rec.samples[sample]['GT']) for sample in samples], index=samples, dtype=int)
             
            if method == 'pEA':
                dmatrix[gene] *= (1 - ea/100) ** gts
                return 1 - dmatrix
                
            elif method == 'sumEA':
                dmatrix[gene] += ea*gts            
                return dmatrix
            


def parse_args():
    """
    Parses the Big pipeline arguments.
    """
    parser = argparse.ArgumentParser(description="Big Pipeline arguments")
    parser.add_argument('--VCF', nargs='?', default='./',help='Location of Cohort VCF')
    parser.add_argument('--savepath', nargs='?', default='./',help='save path for output')
    parser.add_argument('--samples', nargs='?', default='./samples.txt',help='samples file path')
    parser.add_argument('--maxaf', type=float, default=0.01, help='maximum allele frequency cutoff')
    parser.add_argument('--minaf', type=float, default=0, help='minimum allele frequency cutoff')
    parser.add_argument('--method', default='pEA', choices=('pEA', 'sumEA'), help='how to compute input matrix for the pipelines')
    parser.add_argument('--transcript', default='canonical', choices=('all', 'max', 'mean', 'canonical'),help='how to parse EA scores from different transcripts')
    parser.add_argument('--cores', type=int, default=1, help='number of CPUs to use for multiprocessing')
    parser.add_argument('--chrX', type=int,default=1, help='1 if there is sex chromosome in the VCF, 0 if there is no sex chromosome in the VCF')
    return parser.parse_args()


def main(args):
    if args.chrX==1:
        ref = pd.read_csv('./refs/ENSEMBL-lite_GRCh38.v94.txt', delimiter='\t', header=0, index_col='gene')
    elif args.chrX==0:
        ref = pd.read_csv('./refs/ENSEMBL-lite_GRCh38.v94.noX.txt', delimiter='\t', header=0, index_col='gene')
     
    samples = pd.read_csv(args.samples, index_col=0,header=None).index.tolist()

    gene_dfs = Parallel(n_jobs=args.cores)(delayed(parse_VEP)(vcf_fn=args.VCF, gene=gene,gene_ref=ref.loc[gene],samples=samples, method = args.method,  min_af=args.minaf, max_af=args.maxaf, af_field='AF', EA_parser=args.transcript) for gene in tqdm(ref.index.unique()))
    design_matrix = pd.concat(gene_dfs, axis=1) 
    design_matrix.to_csv(args.savepath+'input_matrix.csv', header=True, index=True)    
    
if __name__ == "__main__":
    args = parse_args()
    main(args)
