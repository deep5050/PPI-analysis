---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.4.2
  kernelspec:
    display_name: Python 3
    name: python3
---

<!-- #region id="view-in-github" colab_type="text" -->
<a href="https://colab.research.google.com/github/deep5050/PPI-analysis/blob/master/genomics_common_patterns.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>
<!-- #endregion -->

```python id="84C0ivcZRbBO" colab_type="code" outputId="1925138f-60fb-4fac-ed03-83be422a9590" colab={"base_uri": "https://localhost:8080/", "height": 34}
from google.colab import drive
drive.mount('/gdrive')
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
```

<!-- #region id="OKS8ZIl3ylIO" colab_type="text" -->
# NECESSARY FUNCTIONS
<!-- #endregion -->

```python id="X37Lc_tVRlYk" colab_type="code" colab={}
from collections import OrderedDict

def revs_complement(dna):
    """
    help function for orf_identifier:
    to transform a sequence to reverse complementary sequence
    """
    pairs = {"A": "T", "C": "G", "G": "C", "T": "A"} # complementary code
    c_dna = [pairs[s] for s in dna] # complementary replace
    return "".join(c_dna)[::-1].strip() # reverse

def find_repeats( dna, n):
    """
    This help function for repeats_identifier find and count repeats for 
    each dna sequence
    dna: sequence, string
    n: number of repeats, int
    """
    repeats = {}
    for i in range(0, len(dna)):
        repeat = dna[i:i+n] # generate possible repeats
        if len(repeat) == n:
            if repeat not in repeats:
                repeats [repeat] = 1 # initiate record
            else:
                # count repeated repeats
                repeats[repeat] = repeats.get(repeat) + 1
    return dict(sorted(repeats.items()))

def get_seq_str(content):
  seq_str = ''
  seq = content[1:] # avoiding the header started with '>'
  for seq_ in seq:
    seq_str += str(seq_).replace('\n', '').replace('\r', '').replace(' ','')
  return seq_str.upper().strip()

def get_seq_neum(seq_str):
  table = {'A':'1','T':'2','G':'3','C':'4'}
  # seq_neum = seq_str.replace('A','1').replace('T','2').replace('G','3').replace('C','4')
  seq_neum = [ table[s] for s in seq_str]
  return "".join(seq_neum).strip()

def get_unique_patterns(dna,n):
  repeats = []
  for i in range(0, len(dna)):
    repeat = dna[i:i+n] # generate possible repeats
    if len(repeat) == n:
         repeats.append(repeat)
  return sorted(list(dict.fromkeys(repeats))) #eliminate duplicates

def get_seq_yr(seq_neum):
  table = { '1':'R','3':'R','2':'Y','4':'Y'}
  seq_yr = [ table[s] for s in seq_neum]
  # seq_yr = seq_neum.replace('1','R').replace('2','Y').replace('3','R').replace('4','Y')
  return "".join(seq_yr).strip()

def get_seq_yr_from_atcg(seq):
  table = { 'A':'R','G':'R','C':'Y','T':'Y'}
  seq_yr = [ table[s] for s in seq]
  # seq_yr = seq_neum.replace('1','R').replace('2','Y').replace('3','R').replace('4','Y')
  return "".join(seq_yr).strip()

def get_repeats_with_loc(dna,n):
  # get a sorted list of unique patterns
  patterns = get_unique_patterns(dna,n)
  #create a dictionary with pattern and its loc in the list
  patterns_loc = {}
  # to hold actual loc
  locs = []
  for k in range(0,len(patterns)):
    locs.append([])

  for i in range(0,len(patterns)):
    patterns_loc[patterns[i]] = i #insert the loc

  print(patterns_loc)
  # now actually including all the locs
  for i in range(0,len(dna)):
    repeat = dna[i:i+n]
    if len(repeat)==n:
      # get the location (where to be inserted in the pattern_loc) of the pattern just found
      arr_loc = patterns_loc[repeat]
      # print("arrloc{}".format(arr_loc))
      # include the actual location of occurance in dna
      locs[arr_loc].append(i)

  #returns the loc of each unique pattern as an array of array
  return locs

%matplotlib inline

def plot_protein_distr_by_name(dna):
  # protein distribution
  proteins = translate(seq_str)
  proteins_repeats = find_repeats(proteins,1)
  #
  frq_protein = list(proteins_repeats.values())
  proteins_names = list(proteins_repeats.keys())
  fig = plt.figure(figsize=[10,5])

  
  plt.bar(proteins_names,frq_protein,width=0.7,color='y')
  plt.title("PRODUCED PROTEINS FREQUENCY - sorted: name")
  plt.show()

def plot_protein_distr_by_frq(dna):
  # protein distribution
  proteins = translate(seq_str)
  proteins_repeats = find_repeats(proteins,1)
  # sort by values
  proteins_repeats = {k: v for k, v in sorted(proteins_repeats.items(), key=lambda item: item[1])}
  frq_protein = list(proteins_repeats.values())
  proteins_names = list(proteins_repeats.keys())
  fig = plt.figure(figsize=[10,5])

  
  plt.bar(proteins_names,frq_protein,width=0.7,color='r')
  plt.title("PRODUCED PROTEINS FREQUENCY - sorted: frqncy")
  plt.show()



def plot_pattern_distr_by_name(dna):
  codon_repeats = find_repeats(dna,3)
  frq = list(codon_repeats.values())
  seq_x = list(codon_repeats.keys())
  
  fig = plt.figure(figsize=[4,4])

  plt.bar(seq_x,frq,width=0.5)
  plt.ylabel('frequency')
  plt.title("CODON DISTRIBUTION - sorted: name")
  plt.show()



def plot_pattern_distr_by_frq(dna):
  codon_repeats = find_repeats(dna,3)
  codon_repeats = {k: v for k, v in sorted(codon_repeats.items(), key=lambda item: item[1])}
  frq = list(codon_repeats.values())
  seq_x = list(codon_repeats.keys())
  
  fig = plt.figure(figsize=[4,4])

  plt.bar(seq_x,frq,width=0.5,color='g')
  plt.ylabel('frequency')
  plt.title("CODON DISTRIBUTION - sorted: frqncy")
  plt.show() 








def percent_of_occurence(seq):
  repeats_ = find_repeats(seq,3)
  total_occur = 0
  percent_dist = {}
  for rpt in repeats_:
    total_occur += repeats_[rpt]
  # print("total",total_occur)
  for rpt in repeats_:
    percent_dist[rpt] = (repeats_[rpt]/total_occur)*100
  percents = [ per for per in percent_dist.values() ]
  return percent_dist,percents




def print_pattern_wise_locs(patterns,locs):
  patterns = patterns
  locs = locs
  for i in range(0,len(patterns)):
    print("Patterns:    {} | length: {} | loc: {}".format(patterns[i],len(locs[i]),locs[i]))



def translate(seq): 
       
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 
    protein ="" 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            protein+= table[codon] 
    return protein 
```

```python id="0msQmXSsRy9K" colab_type="code" colab={}

from os import listdir
from os.path import isfile, join

# trying to find common patterns of length 3 among all the genes
def common_patterns_diff_lengths_among_all(my_path_,length_,MODE='ATCG'):
  my_path = my_path_
  only_files = [f for f in listdir(my_path) if isfile(join(my_path, f))]
  all_seq = []
  file_names = []
  for file_ in only_files:

    file_path = my_path+'/'+file_
    file_names.append(file_.split('.txt')[0])

    with open(file_path,'r') as f:
      content = f.readlines()
      seq_str = get_seq_str(content)
      if MODE == 'YR':
        seq_str = get_seq_yr_from_atcg(seq_str)
      all_seq.append(seq_str)

    uniq_patterns = []
    for seq in all_seq:
      patt = get_unique_patterns(seq,length_)
      uniq_patterns.append(patt)

  # print(len(uniq_patterns))

  # print(uniq_patterns[0])
  # print(uniq_patterns[1])
  # xx = set(uniq_patterns[0]) & set(uniq_patterns[1])
  # print(sorted(xx))

  common_pattern_among_all = set(uniq_patterns[0])
  for i in range(0,len(uniq_patterns)):
    # common_pattern_among_all = common_pattern_among_all.intersection(set(uniq_patterns[i]))
    common_pattern_among_all = common_pattern_among_all & set(uniq_patterns[i])

  # for patt in uniq_patterns:
    
    # print(len(patt),":",patt)

  # print('-'*50)
  print(len(common_pattern_among_all),":",sorted(common_pattern_among_all))

```

<!-- #region id="k-obFVyAyDHL" colab_type="text" -->
# FINDING PATTERNS FOR DJ1
<!-- #endregion -->

<!-- #region id="SV4ZEU4F0AaL" colab_type="text" -->
## ATCG DISTRBUTION
<!-- #endregion -->

```python id="c0bs_1DHR4yM" colab_type="code" outputId="1e4cce03-f226-4e27-e66b-7d12d0b216b0" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',3)
```

```python id="ko_yHlgLR8aw" colab_type="code" outputId="87294686-dee4-40c1-ac6d-eebbb07323b7" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',4)
```

```python id="tIBp3j0uSfYi" colab_type="code" outputId="0334f8c7-369f-4409-fa46-a0e35a25c3d2" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',5)
```

```python id="iAdvkxJ-SnQy" colab_type="code" outputId="5c578561-fe46-4b49-a604-0edf241fe623" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',6)
```

```python id="ZADqdOuzSrSH" colab_type="code" outputId="a7894130-c787-4aaf-c09e-4af79657c4e6" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',7)
```

```python id="ZUlUJoqQSw8w" colab_type="code" outputId="421402dd-0be4-4895-ba64-af26b0db9674" colab={"base_uri": "https://localhost:8080/", "height": 34}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',8)
```

```python id="7p3xhnmUTEBJ" colab_type="code" outputId="b9570618-835b-4ecf-d9eb-9e1c28cd22fb" colab={"base_uri": "https://localhost:8080/", "height": 34}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',9)
```

<!-- #region id="Pp62ekDJ0IsI" colab_type="text" -->
## YR DISTRIBUTION
<!-- #endregion -->

```python id="KQIf4zPm-CIl" colab_type="code" outputId="a29665c2-b764-4061-9327-006c0e96db3f" colab={"base_uri": "https://localhost:8080/", "height": 34}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',3,MODE='YR')
```

```python id="EcVS9VDf-HHw" colab_type="code" outputId="b847a739-f2ae-47ab-cad6-75118795e275" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',4,MODE='YR')
```

```python id="dEFlQCia_99F" colab_type="code" outputId="4434b87b-aa40-4c4a-9a82-90c7f15ee6d9" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',5,MODE='YR')
```

```python id="FI_S17H_ACMJ" colab_type="code" outputId="44386685-6af3-452d-b52f-0e00f1de2c6f" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',6,MODE='YR')
```

```python id="7z-EKTAeAGIP" colab_type="code" outputId="aaeaa4a9-2846-431b-8028-5f830cbaff05" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',7,MODE='YR')
```

```python id="iIVCp3-cAJcs" colab_type="code" outputId="eff84784-9a6f-4c83-b023-4bea8504f8cb" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',8,MODE='YR')
```

```python id="LLdwtzvHANKL" colab_type="code" outputId="51239bd2-3280-43c8-aa36-1226fd015750" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',9,MODE='YR')
```

```python id="yZhpqcX7ARgC" colab_type="code" outputId="07e02750-f18d-48ea-f9c3-519c8e6a3c15" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',10,MODE='YR')
```

```python id="MmW8RNzCAXFQ" colab_type="code" outputId="e4591552-6f6d-4b8c-d458-bf18f1cafc15" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',11,MODE='YR')
```

```python id="IbQno-E7AbNE" colab_type="code" outputId="c105e43e-ae52-4859-dc38-80a9ede522a6" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',12,MODE='YR')
```

```python id="a4ss47uuAlvH" colab_type="code" outputId="88410c07-a0c8-4b3f-f683-c3de5dd82588" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',13,MODE='YR')
```

```python id="OT5S8HS8AqDE" colab_type="code" outputId="1b345675-69bc-447d-8b6e-fcb634e7b84e" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',14,MODE='YR')
```

```python id="UiBEQJipAuuX" colab_type="code" outputId="fa87bb42-c2df-4c2b-dce8-c99514606f81" colab={"base_uri": "https://localhost:8080/", "height": 34}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',15,MODE='YR')
```

```python id="azKFYmOXAyya" colab_type="code" outputId="ac0a314b-27b4-4c3a-d617-cc54c87eb005" colab={"base_uri": "https://localhost:8080/", "height": 34}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/DJ1',16,MODE='YR')
```

<!-- #region id="p4kAiQ8fy455" colab_type="text" -->
## RESULT

<!-- #endregion -->

```python id="qZZM_e6_BS1X" colab_type="code" outputId="2882d600-901c-4150-c2e8-7fcf3a2bf978" colab={"base_uri": "https://localhost:8080/", "height": 441}

# YR distribution plot
yr_data = np.array([8,16,32,64,128,255,480,736,672,319,97,13,0])
yr_x_axis = range(3,16)

atcg_data = [64,253,620,342,13,0,0]
atcg_x_axis = range(3,10)

fig,ax = plt.subplots(nrows=1,ncols=2,figsize=(10,6))

ax[0].bar(atcg_x_axis,atcg_data,color='m')
ax[0].plot(atcg_x_axis,atcg_data,color='k')
ax[0].set_title('Common ATCG distribution over the genes')
ax[0].set_xlabel("Length of pattern")
ax[0].set_ylabel("Count of pattern")


ax[1].bar(yr_x_axis,yr_data,color='r')
ax[1].set_title('Common YR distribution over the genes')
ax[1].plot(yr_x_axis,yr_data,color='k')
ax[1].set_xlabel("Length of pattern")
ax[1].set_ylabel("Count of pattern")
fig.tight_layout()
```

<!-- #region id="YZlwlY2xsYSd" colab_type="text" -->
# FINDING PATTERNS FOR PARKIN
<!-- #endregion -->

<!-- #region id="ClEfnYsrzt-N" colab_type="text" -->
## ATCG DISTRIBUTION
<!-- #endregion -->

```python id="bCstkikCrFPJ" colab_type="code" outputId="f483f254-3196-4abe-fef3-731b7a7bd197" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',3)
```

```python id="qGdn5rCqrRQn" colab_type="code" outputId="c6d3b4ab-e29e-44d7-b575-45d92f2a5a3a" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',4)
```

```python id="3qJvjtv1rfS1" colab_type="code" outputId="e9f697ce-d788-40fa-a495-3802d6d3eef7" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',5)
```

```python id="bVwalpSirn6q" colab_type="code" outputId="578dfc2a-7326-4edb-c36c-94565511081c" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',6)
```

```python id="hXu3u7Wjrwuw" colab_type="code" outputId="d7612f2d-68d3-45cc-fa33-438d0ca49dfc" colab={"base_uri": "https://localhost:8080/", "height": 34}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',7)
```

```python id="iEk24t9Tr69Z" colab_type="code" outputId="f49ec862-37f2-4423-e1ac-5e03a7a6c7e1" colab={"base_uri": "https://localhost:8080/", "height": 34}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',8)
```

<!-- #region id="3hfgUu4nzmbz" colab_type="text" -->
## YR DISTRIBUTION
<!-- #endregion -->

```python id="3K3EJTtHsxlp" colab_type="code" outputId="447871f6-0f7c-4b60-ddcd-c8d0ba903301" colab={"base_uri": "https://localhost:8080/", "height": 34}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',3,MODE='YR')
```

```python id="QP8o6CcLs4GF" colab_type="code" outputId="51946cb6-36d0-4e71-b90f-20ac2a3ace3c" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',4,MODE='YR')
```

```python id="7Aawv-zOtCHN" colab_type="code" outputId="03540d65-062d-40c5-cf04-81bee03a42a6" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',5,MODE='YR')
```

```python id="kfZ0Yd0itE29" colab_type="code" outputId="2b8b4198-b701-4dd3-e5e8-76f38762d53c" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',6,MODE='YR')
```

```python id="LVBWsA-Hthz_" colab_type="code" outputId="f6a93385-4907-4527-834c-3ebe284b20c7" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',7,MODE='YR')
```

```python id="HiaXwilLt2um" colab_type="code" outputId="c69ca1bd-41df-45e8-c382-1578d325afd2" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',8,MODE='YR')
```

```python id="gFe9HcnduBvu" colab_type="code" outputId="971f68d8-5aab-4452-98c1-56579abc0bda" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',9,MODE='YR')
```

```python id="IgBowm2BuMSA" colab_type="code" outputId="e7e80a3e-a358-4a3a-edc5-2b84f3f2db1c" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',10,MODE='YR')
```

```python id="4-tjLO7MuWL6" colab_type="code" outputId="21e888b0-44f4-48e6-9c3d-c68fd3c8bd54" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',11,MODE='YR')
```

```python id="ckbOOrFluf0k" colab_type="code" outputId="3703325f-a237-43a9-832f-6aa189827395" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',12,MODE='YR')
```

```python id="jBWmW0GIuqAf" colab_type="code" outputId="15065f45-3ad8-412c-8e8c-26e3809d8da5" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',13,MODE='YR')
```

```python id="OZx_CPpWu1Hk" colab_type="code" outputId="7537a8c5-7404-4f2d-8bf7-991801da6c53" colab={"base_uri": "https://localhost:8080/", "height": 34}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',14,MODE='YR')
```

```python id="uwxm6O8KvAy-" colab_type="code" outputId="c5ffea26-928e-4196-99ca-ddbed427ee12" colab={"base_uri": "https://localhost:8080/", "height": 34}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PARKIN',15,MODE='YR')
```

<!-- #region id="dSmYyT7CzFJS" colab_type="text" -->
## RESULT
<!-- #endregion -->

```python id="X2oPTN6TvPtr" colab_type="code" outputId="a7db9db5-ed0e-4716-aaa0-51d0647be0fd" colab={"base_uri": "https://localhost:8080/", "height": 441}
# YR distribution plot
yr_data = np.array([8,16,32,64,128,254,448,534,344,110,24,2,0])
yr_x_axis = range(3,16)

atcg_data = np.array([64,241,402,100,2,0])
atcg_x_axis = range(3,9)

fig,ax = plt.subplots(nrows=1,ncols=2,figsize=(10,6))

ax[0].bar(atcg_x_axis,atcg_data,color='m')
ax[0].plot(atcg_x_axis,atcg_data,color='k')
ax[0].set_title('Common ATCG distribution over the genes')
ax[0].set_xlabel("Length of pattern")
ax[0].set_ylabel("Count of pattern")


ax[1].bar(yr_x_axis,yr_data,color='r')
ax[1].set_title('Common YR distribution over the genes')
ax[1].plot(yr_x_axis,yr_data,color='k')
ax[1].set_xlabel("Length of pattern")
ax[1].set_ylabel("Count of pattern")
fig.tight_layout()
```

```python id="eg01hBnZwf6x" colab_type="code" colab={}

```

<!-- #region id="-CSamhJg0aAN" colab_type="text" -->
# FINDING PATTERNS FOR PINK1
<!-- #endregion -->

<!-- #region id="8j3SXeZI0jhJ" colab_type="text" -->
## ATCG DISTRIBUTION
<!-- #endregion -->

```python id="f5opG6ge0mDO" colab_type="code" outputId="3487f72c-48f0-4288-dbe3-dd2ee5542169" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',3)
```

```python id="ixTfDgeS0oN1" colab_type="code" outputId="aa280aa5-7204-4757-bddd-77958aee47d2" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',4)
```

```python id="zUn4jw8_0tl4" colab_type="code" outputId="3fe57174-2658-4eb7-daf7-da52561d78ae" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',5)
```

```python id="aC79-3ma0wxn" colab_type="code" outputId="8605c52d-7d6a-4fe5-ecd9-164d8293daae" colab={"base_uri": "https://localhost:8080/", "height": 54}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',6)
```

```python id="5vqW3ST20zvd" colab_type="code" outputId="48810d31-7d4b-4e00-edf4-047794956d2b" colab={"base_uri": "https://localhost:8080/", "height": 34}
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',7)
```

```python id="6vpQhBsN027c" colab_type="code" colab={"base_uri": "https://localhost:8080/", "height": 34} outputId="0d7120f4-f84d-4271-cc45-596f4fe1ef0a"
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',8)
```

<!-- #region id="xp_iZmrJ0-zw" colab_type="text" -->
## YR DISTRIBUTION
<!-- #endregion -->

```python id="iYj6kmEE1BKH" colab_type="code" colab={"base_uri": "https://localhost:8080/", "height": 34} outputId="618919bd-762c-4f53-b6a9-0d1b8502bf77"
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',3,MODE='YR')
```

```python id="LqIr94NS1GLc" colab_type="code" colab={"base_uri": "https://localhost:8080/", "height": 54} outputId="231723b2-f6f2-4c45-b8c0-8d4d7f81c38e"
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',4,MODE='YR')
```

```python id="B-4TYKun1JoT" colab_type="code" colab={"base_uri": "https://localhost:8080/", "height": 54} outputId="e2f6b080-7173-4ff5-8f8d-ab9fcfe1d614"
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',5,MODE='YR')
```

```python id="vMIxurAJ1RU-" colab_type="code" colab={"base_uri": "https://localhost:8080/", "height": 54} outputId="57a87d87-9c39-4347-c303-b0d8810050bf"
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',6,MODE='YR')
```

```python id="HC57-ZQZ1Ujt" colab_type="code" colab={"base_uri": "https://localhost:8080/", "height": 54} outputId="3c7ef0d9-0126-4570-b206-c6727c485496"
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',7,MODE='YR')
```

```python id="Crdo2uE_1YAS" colab_type="code" colab={"base_uri": "https://localhost:8080/", "height": 54} outputId="a052f927-66bd-475f-e000-7787ac47dad4"
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',8,MODE='YR')
```

```python id="11YzeT8M1bSp" colab_type="code" colab={"base_uri": "https://localhost:8080/", "height": 54} outputId="d98a6766-8d9d-48e4-de35-57b8628af55c"
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',9,MODE='YR')
```

```python id="nIZM-7gU1fYV" colab_type="code" colab={"base_uri": "https://localhost:8080/", "height": 54} outputId="d1d6f5ec-3d7c-4dc5-8630-19cf3b931ed6"
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',10,MODE='YR')
```

```python id="6GXID3nB1jNz" colab_type="code" colab={"base_uri": "https://localhost:8080/", "height": 54} outputId="21af93a5-86c2-46f8-bb0d-e8d9d8a38f3b"
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',11,MODE='YR')
```

```python id="-AR1pXej1m9r" colab_type="code" colab={"base_uri": "https://localhost:8080/", "height": 54} outputId="f49496d9-34f7-4096-8ddf-275313c1cc6e"
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',12,MODE='YR')
```

```python id="M3n7whIK10yZ" colab_type="code" colab={"base_uri": "https://localhost:8080/", "height": 34} outputId="dd220164-6a83-4f89-e8a8-4f671e2a2d47"
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',13,MODE='YR')
```

```python id="10YMTTVb14sW" colab_type="code" colab={"base_uri": "https://localhost:8080/", "height": 34} outputId="8557050b-d1e4-4d09-9040-6e99f8d6d30f"
common_patterns_diff_lengths_among_all('/gdrive/My Drive/Genomics/job/data/PINK1',14,MODE='YR')
```

<!-- #region id="LN_Gy80LDMfL" colab_type="text" -->
## RESULT
<!-- #endregion -->

```python id="riICtxVa18g8" colab_type="code" colab={"base_uri": "https://localhost:8080/", "height": 441} outputId="0e05fcf0-9849-4a99-b1cd-13e663d8f385"
# YR distribution plot
yr_data = np.array([8,16,32,64,128,255,455,546,279,68,5,0])
yr_x_axis = range(3,15)

atcg_data = np.array([64,247,445,87,40,0])

atcg_x_axis = range(3,9)

fig,ax = plt.subplots(nrows=1,ncols=2,figsize=(10,6))

ax[0].bar(atcg_x_axis,atcg_data,color='m')
ax[0].plot(atcg_x_axis,atcg_data,color='k')
ax[0].set_title('Common ATCG distribution over the genes')
ax[0].set_xlabel("Length of pattern")
ax[0].set_ylabel("Count of pattern")


ax[1].bar(yr_x_axis,yr_data,color='r')
ax[1].set_title('Common YR distribution over the genes')
ax[1].plot(yr_x_axis,yr_data,color='k')
ax[1].set_xlabel("Length of pattern")
ax[1].set_ylabel("Count of pattern")
fig.tight_layout()
```

<!-- #region id="3CghU054Bmlj" colab_type="text" -->
# FINAL GRAPH
<!-- #endregion -->

```python id="MBJZDGp-Bp1u" colab_type="code" colab={"base_uri": "https://localhost:8080/", "height": 295} outputId="41a3e01a-6c0b-4460-8be9-94ad82261e0a"
dj1 = np.array([8,16,32,64,128,255,480,736,672,319,97,13,0])
parkin = np.array([8,16,32,64,128,254,448,534,344,110,24, 2,0])
pink1 = np.array([8,16,32,64,128,255,455,546,279,68, 5, 0, 0])


purple_patch = mpatches.Patch(color='purple', label='DJ1')
green_patch = mpatches.Patch(color='green',label ='PARKIN')
red_patch = mpatches.Patch(color="red",label="PINK1")
plt.legend(handles=[purple_patch,green_patch,red_patch])
x_axis = range(3,16)
plt.plot(x_axis,dj1,color='purple')
plt.plot(x_axis,parkin,color='green')
plt.plot(x_axis,pink1,color='red')
plt.title('YR pattern distributions over the 3 networks')
plt.xlabel("Pattern length")
plt.ylabel("Pattern count")
plt.show()
```
