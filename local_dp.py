import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from global_dp import *

def local_dp(seq1,seq2,gap_penalty,DNA=True):

    m,n = len(seq1),len(seq2)

    if DNA == True : 
        Sub_score = pd.read_csv('DNA_substitution_scores.csv',index_col='0')
    else : 
        Sub_score = pd.read_csv('blosum62.csv',index_col='0')

    M = np.zeros((m+1,n+1))  
    M[:,0] = 0*np.arange(0,m+1)
    M[0,:] = 0*np.arange(0,n+1)

    tracer = np.zeros((M.shape[0],M.shape[1]))

    for i in range(1,m+1):
        for j in range(1,n+1):
            arr = np.array([M[i-1,j-1] + (Sub_score[seq1[i-1]][seq2[j-1]]), (M[i-1,j]+gap_penalty), (M[i,j-1]+gap_penalty),0])
            M[i,j] = arr.max()
            tracer[i,j] = np.argmax(arr)+1

    return M, tracer

def local_tracer(M,tracer,seq1,seq2):

### Obtain Local traceback route

    if np.count_nonzero(M==M.max()) == 1 :  
        l_traceback = np.zeros((M.shape)) # Initialize route matrix
        tm,tn = np.unravel_index(M.argmax(), M.shape) # find the cell with highest score
        l_traceback[tm,tn] = 1
        current_cell = (tm,tn)

        while M[tm,tn] > 0 : 
            if tracer[current_cell] == 1: 
                current_cell = (tm-1,tn-1)
                l_traceback[current_cell] = 1
                tm -= 1 ; tn -= 1
                
            elif tracer[current_cell] == 2 :
                current_cell = (tm-1,tn)
                l_traceback[current_cell] = 1
                tm -= 1

            elif tracer[current_cell] == 3:
                current_cell = (tm,tn-1)
                l_traceback[current_cell] = 1
                tn -= 1 
    else:
        print('NOTE : Multiple Local alignments possible') 
        
        l_traceback = np.zeros((M.shape)) # Initialize route matrix

    return l_traceback 

def local_traceback(M,tracer,seq1,seq2):

    """ Obtain Final score and Global ALignment """
    
    tm,tn = np.unravel_index(M.argmax(), M.shape) 
    score = [M[tm,tn]]
    align1 = []
    align2 = []

    current_cell = (tm,tn)

    while M[tm,tn] > 0 : 
        if tracer[current_cell] == 1: 
            current_cell = (tm-1,tn-1)
            score.append(M[current_cell])
            align1.append(seq1[tm-1])
            align2.append(seq2[tn-1])

            tm -= 1 ; tn -= 1
            
        elif tracer[current_cell] == 2 :
            current_cell = (tm-1,tn)
            score.append(M[current_cell])
            align2.append('-')
            align1.append(seq1[tm-1])
            tm -= 1

        elif tracer[current_cell] == 3:
            current_cell = (tm,tn-1)
            score.append(M[current_cell])
            align2.append(seq2[tn-1])
            align1.append('-')
            tn -= 1 

    alg1 = conv_list2str(align1[::-1])
    alg2 = conv_list2str(align2[::-1])
    alignment = alg1 + '\n' + alg2

    return score[0],alignment

if __name__ == '__main__':
    
    input('Welcome to Local Sequence Alignment \n Press any key to continue')
    dna_or_protein = str(input('Do you want to align a DNA or Protein sequence ? \n Press 1 for DNA \n Press 2 for Protein:\n' ))

    if dna_or_protein == '1':
        seq1 = str(input('Enter DNA sequence 1:'))
        seq2 = str(input('Enter DNA sequence 2:'))
        M,tracer = local_dp(seq1,seq2,gap_penalty=-2,DNA=True)
        T = local_tracer(M,tracer,seq1,seq2)
        
        score, algn = local_traceback(M,tracer,seq1,seq2)
        print('Score',score)
        print('Local Alignment for the sequences:\n')
        print(algn)
        
        plot_align_matrix(M,seq1,seq2,Global=False) 
        plot_tracer(T,seq1,seq2,Global=False) 
        
        plt.show()
        

    else : 
        seq1 = str(input('Enter Protein sequence 1:'))
        seq2 = str(input('Enter Protein sequence 2:'))

        M,tracer = local_dp(seq1,seq2,gap_penalty=-2,DNA=False)
        T = local_tracer(M,tracer,seq1,seq2)
        
        score, algn = local_traceback(M,tracer,seq1,seq2)
        print('Score',score)
        print('Local Alignment for the sequences:\n')
        print(algn)
        
        plot_align_matrix(M,seq1,seq2,Global=False) 
        plot_tracer(T,seq1,seq2,Global=False) 
        
        plt.show()
    
    
