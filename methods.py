import numpy as np
import os
import h5py
import matplotlib.pyplot as plt

LETTERS = set(['A','G','C','T'])
TRACE_PATH = '/Analyses/RawGenomeCorrected_000/BaseCalled_template/Events'

def create_motifs_set(mlen):
    motifs = LETTERS
    for i in range(mlen - 1):
        new_motifs = set()
        for motif in motifs:
            for l in LETTERS:
                new_motifs.add(motif + l)
        motifs = set(new_motifs)
    motifs_dir = {
        m:[] for m in motifs
    }
    return motifs_dir

def cohend(d1, d2):
    n1, n2 = len(d1), len(d2)
    s1, s2 = np.var(d1, ddof=1), np.var(d2, ddof=1)
    s = np.sqrt(((n1-1) * s1 + (n2-1) * s2) / (n1 + n2 - 2))
    u1, u2 = np.mean(d1), np.mean(d2)
    return (u1 - u2) / s


def extract_motifs(dirpath, mlen=6, max_reads=8000):
    
    sample_motifs = create_motifs_set(mlen)


    left_shift = int(np.floor(mlen/2))
    right_shift = int(np.ceil(mlen/2))
    
    reads_counter = 0

    fast5_files = os.listdir(dirpath)
    fast5_files = [f for f in fast5_files if '.fast5' in f]

    for j in range(len(fast5_files)):
        print('\tBatch', j, '...')

        with h5py.File('{}/{}'.format(dirpath, fast5_files[j]), 'r', rdcc_nbytes=1024**3) as f:
            for i in list(f.items()):
                print('\t\tRead {} from {}'.format(reads_counter, max_reads), end='')
                readname = i[0]
                try:
                    trace =f['/{}'.format(readname) + TRACE_PATH][:]
                except KeyError:
                    print('\r', end='')
                    continue

                seq = ''.join([t[4].decode() for t in trace])
                signal = [t[0] for t in trace]
                for mi in range(left_shift, len(seq) - right_shift):
                    sample_motifs[seq [mi-left_shift:mi+right_shift] ].append(signal[mi])

                reads_counter += 1
                print('\r', end='')
                if reads_counter >= max_reads:
                    break
            print()
                
        if reads_counter >= max_reads:
            break
    print()
    return sample_motifs

def compare_samples(
    motifs_1, 
    motifs_2, 
    name1='sample 1', 
    name2='sample 2',  
    save_path='Motifs', 
    mmincount=100,
    min_effect_size=0.5,
    plot=True,
    ):

    summary_list = []

    try:
        os.mkdir(save_path)
    except FileExistsError:
        print('Save path exists! Files will be rewritten!')
    print('Processing...')
    motif_counter = 1
    for motif in motifs_1:
        print('\tMotif {} from {}'.format(motif_counter, len(motifs_1)), end='')
        motif_counter += 1
        if min(len(motifs_1[motif]), len(motifs_2[motif])) < 100:
            print('\r', end='')
            continue
        
        effect_size = cohend(motifs_1[motif], motifs_2[motif])
        summary_list.append((abs(effect_size), motif))
        
        if abs(effect_size) < min_effect_size:
            print('\r', end='')
            continue
        if plot:
            plt.hist(motifs_1[motif], label=name1, bins=100, alpha=0.5, density=True)
            plt.hist(motifs_2[motif], label=name2, bins=100, alpha=0.5, density=True)
            plt.title('{}: {}'.format(motif, round(effect_size, 4)))
            plt.legend()

            plt.tight_layout()
            plt.savefig('{}/{}.pdf'.format(save_path, motif), format='pdf')
            plt.close()
            print('\r', end='')

    summary_list.sort(reverse=True)
    print('Completed.')
    print('Writing summary..')

    f = open('{}/summary.csv'.format(save_path), 'w')
    f.write('MOTIF\tEffectSize\tplotted\n')

    for s in summary_list:
        if abs(s[0]) < min_effect_size:
            plotted = 'no'
        else:
            plotted = 'yes'
        f.write('{}\t{}\t{}\n'.format(s[1], s[0], plotted))
    f.close()

    f = open('{}/sequences.fasta'.format(save_path), 'w')
    for s in summary_list:
        if abs(s[0]) < min_effect_size:
            continue
        f.write('>{}; Effect_size {}\n'.format(s[1],s[0]))
        f.write(s[1] + '\n')
    f.close()
    
    print('Done!')