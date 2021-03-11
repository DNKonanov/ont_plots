import h5py
import matplotlib.pyplot
import numpy as np
import os
from scipy.stats import mannwhitneyu
import argparse
from methods import compare_samples, extract_motifs

parser = argparse.ArgumentParser()

parser.add_argument('-s1_dir', type=str, help='path to sample1 multi-fast5 files', required=True)
parser.add_argument('-s2_dir', type=str, help='path to sample2 multi-fast5 files', required=True)
parser.add_argument('-maxreads', default=8000, type=int, help='The maximum of observed reads count for each sample. The default value is 8000')
parser.add_argument('-savepath', default='Motifs', type=str, help='An output directory. The default value is "Motifes" in the current directory')
parser.add_argument('-mlen', default=6, type=int, help='Motifs length. The default (and recommended) value is 6')
parser.add_argument('-s1_name', default=None, type=str, help='Sample 1 name which will be used in the output plots. The default value is the same as the -s1_dir parameter')
parser.add_argument('-s2_name', default=None, type=str, help='Sample 2 name which will be used in the output plots. The default value is the same as the -s2_dir parameter')
parser.add_argument('-mmincount', default=100, type=int, help='The minimal count of motif appearences to compute statistics. The default value is 100')
parser.add_argument('-min_effectsize', default=0.3, type=float, help='The minimal Cohen\'s effect size to plot signal distributions.The default value is 0.5')
parser.add_argument('-noplot', action='store_false')
args = parser.parse_args()


if args.s1_name is None:
    name1 = args.s1_dir
else:
    name1 = args.s1_name

if args.s2_name is None:
    name2 = args.s2_dir
else:
    name2 = args.s2_name

print('Sample 1 extracting..')
motifs_1 = extract_motifs(
    args.s1_dir, 
    mlen=args.mlen, 
    max_reads=args.maxreads)


print('Sample 2 extracting..')
motifs_2 = extract_motifs(
    args.s2_dir, 
    mlen=args.mlen, 
    max_reads=args.maxreads)

compare_samples(
    motifs_1, 
    motifs_2,
    name1=name1,
    name2=name2,
    save_path=args.savepath,
    mmincount=args.mmincount,
    min_effect_size=args.min_effectsize,
    plot=args.noplot,
    )