# ONT Plots

### Dependencies
* Python >3.6
* `h5py`, `numpy`, `matplotlib` Python libraries
* [Tombo](https://github.com/nanoporetech/tombo) software
* [ont_fast5_api](https://github.com/nanoporetech/ont_fast5_api)
### Manual

First, you should resquiggle you fast5 files using Tombo. Next you should convert single-fast5 reads to multi-fast5 with `single_to_multi_fast5` tool.

```
usage: start_script.py [-h] -s1_dir S1_DIR -s2_dir S2_DIR [-maxreads MAXREADS] [-savepath SAVEPATH] [-mlen MLEN] [-s1_name S1_NAME] [-s2_name S2_NAME]` ```[-mmincount MMINCOUNT] [-min_effectsize MIN_EFFECTSIZE]

optional arguments:
  -h, --help            show this help message and exit
  -s1_dir S1_DIR        path to sample1 multi-fast5 files
  -s2_dir S2_DIR        path to sample2 multi-fast5 files
  -maxreads MAXREADS    The maximum of observed reads count for each sample. The default value is 8000
  -savepath SAVEPATH    An output directory. The default value is "Motifes" in the current directory
  -mlen MLEN            Motifs length. The default (and recommended) value is 6
  -s1_name S1_NAME      Sample 1 name which will be used in the output plots. The default value is the same as the -s1_dir parameter
  -s2_name S2_NAME      Sample 2 name which will be used in the output plots. The default value is the same as the -s2_dir parameter
  -mmincount MMINCOUNT  The minimal count of motif appearences to compute statistics. The default value is 100
  -min_effectsize MIN_EFFECTSIZE
                        The minimal Cohen's effect size to plot signal distributions.The default value is 0.5
```
