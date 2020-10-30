# Bioinformatics

This repository contains scripts and programs that is being used for the BIO314 course at IISER Pune. Some scripts were assignments where as most of the other ones were written solely with the joy of writing it. Some important concepts covered ans implemented were : 

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




# References : 

1. https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
2. Durbin, R. (1998). Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids (Illustrated ed.). Cambridge University Press. 
