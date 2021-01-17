# Bioinformatics

This repository contains scripts and programs that was used for the BIO314 course at IISER Pune. Some scripts were assignments whereas most of the other ones were written solely for the joy of writing it. Some important concepts covered are implemented were : 

## Biological Sequence Alignment 

1. Global Sequence Alignment 

- You can find the implementation of the **Needleman-Wunsch algorithm** to align two Protein or Nucleotides sequences at `global_dp.py`.
- We use a dynamic programming approach split into three subsections : {Initialization, Computation and Traceback} to find the best Sequence alignment.
- We also show an approach with **affine gaps** in `global_dp_affine.py`

2. Local Sequence Alignment

- You can find the implementation of the **Smith-Waterman algorithm** for the local alignment of two sequences at `local_dp.py`. 
- We follow a similar approach as (1) but show all possible local alignments that arise primarily due to multiple maxima while computing the scores matrix 
- We show an approach with **affine gaps** in `local_dp_affine.py`

Comparison between Smith-Waterman and Needleman-Wunsch : 

| Property | Smith-Waterman algorithm | Needleman-Wunsh algorithm | 
| -------- | ------------------------| ---------------------------|
| Initialization | First row and first column are set to 0 | First row and first column are subject to gap penalty (affine, linear etc) | 
| Scoring |	Negative score is set to 0 |	Score can be negative |
| Traceback | 	Begin with the highest score, end when 0 is encountered | Begin with the cell at the lower right of the matrix, end at top left cell |


## Next Generation Sequencing (NGS)

- The advent of NGS techniques has led to wonderful application of the [Burrows–Wheeler transformation](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform). In NGS, DNA is fragmented into small pieces, of which the first few bases are sequenced, yielding several millions of "reads", each 30 to 500 base pairs ("DNA characters") long
- In an effort to reduce the memory requirement for sequence alignment, we use BWT as a data-compression-algorithm.

## Hidden Markov Models (HMM) 

- Evaluation (Forward Algorithm)
- Decoding (Vitterbi Algorithm)
- Learning (Baum-Welch Algorithm or the Forward-Backward algorithm) (TBA)

## How to run these scripts locally ? 


1. Clone this Github repository using 

```git clone https://github.com/Anantha-Rao12/Bioinformatics-BIO314```


2. Create a Python (>=3.2) virtual environemnt and call it 'bioinfo-BIO314'.

    - Linux/Mac: `python3 -m venv bioinfo-BIO314`
    - Windows: `python -m venv bioinfo-BIO314`
    - A new directory called "bioinfo-BIO314" will be created.

3. Activate the Virtual Environment by running the following.

    - Linux/Mac: `source bioinfo-BIO314/bin/activate`
    - Windows: `.\bioinfo-BIO314\Scripts\activate`

4. In the new virtual environemnt , run `pip3 install -r requirements.txt` to install all dependencies. On Windows, pip3 will be replaced by pip.

5. Run `python3 global_dp.py` for Global Alignment or run `python3 local_dp.py` for Local Alignment. 
    


# References : 

1. Sean, R. (2004). EDDY: What is dynamic programming. Nature biotechnology, 22(7), 909-910.
2. Durbin, R. (1998). Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids (Illustrated ed.). Cambridge University Press. 
