import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from helper_functions import *

def local_dp(seq1,seq2,gap_penalty,sub_matrix,DNA=True):

    m,n = len(seq1),len(seq2)

    M = np.zeros((m+1,n+1))  
    M[:,0] = 0*np.arange(0,m+1)
    M[0,:] = 0*np.arange(0,n+1)

    sub_matrix = pd.read_csv(sub_matrix,index_col='0')
    tracer = np.zeros((M.shape[0],M.shape[1]))

    for i in range(1,m+1):
        for j in range(1,n+1):
            arr = np.array([M[i-1,j-1] + (sub_matrix[seq1[i-1]][seq2[j-1]]), (M[i-1,j]+gap_penalty), (M[i,j-1]+gap_penalty),0])
            M[i,j] = arr.max()
            tracer[i,j] = np.argmax(arr)+1

    return M, tracer

def local_tracer(M,tracer):
    
    ### Get a list of indices with top scores
    ### Returns a tensor of traceback routes

    h_rowid , h_colid = np.where(M==M.max())
    h_idx = list(zip(h_rowid,h_colid))

    local_traceback_tensor = np.zeros((M.shape[0],M.shape[1],len(h_idx))) # Initialize a zeros tensor

    for counter,idx in enumerate(h_idx) : 
        l_traceback = np.zeros((M.shape)) # Initialize route matrix

        tm,tn = idx # find the cell with highest score
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
        local_traceback_tensor[:,:,counter] = l_traceback

    return local_traceback_tensor

def local_traceback(M,tracer,seq1,seq2):

    ### Get a list of indices with top scores
    h_rowid , h_colid = np.where(M==M.max())
    h_idx = list(zip(h_rowid,h_colid))
    
    if len(h_idx) == 1:
        print('Found One Best Local alignment')
    else:
        print('NOTE : Multiple Local alignments possible. Suspect maximum of ', len(h_idx)) 
   
    local_alignments = []  # list to store all possible local alignments

    for counter, idx in enumerate(h_idx):
        tm,tn = idx
        score = [M[tm,tn]]
        align1,align2 = [],[]
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

        # if M[tm,tn] == 0 : 
        #     align1.append(seq1[tm-1])
        #     align2.append(seq2[tn-1])
        # else:
        #     pass

        alg1 = conv_list2str(align1[::-1])
        alg2 = conv_list2str(align2[::-1])
        alignment = alg1 + '\n' + alg2

        local_alignments.append((score,alignment))

    return local_alignments

def main(M,tracer,seq1,seq2,plot):
    T = local_tracer(M,tracer)
    local_alignments = local_traceback(M,tracer,seq1,seq2)
    print('\n')

    for counter,alignment in enumerate(local_alignments) : 
        print('Alignment no',counter+1)
        print('Score:',sum(alignment[0]))
        print(alignment[1])
        print('\n')
        
        plot_tracer(T[:,:,counter],seq1,seq2,Global=False) if plot == 1 else 0

    if plot == 1 : 
        input('Press any key to look at the Matrix scores and the Traceback')
        plot_align_matrix(M,seq1,seq2,Global=False) 
        plt.show()

    else :
        pass


if __name__ == '__main__':
    
    input('Welcome to Local Sequence Alignment \n Press any key to continue')
    dna_or_protein = str(input('Do you want to align a DNA or Protein sequence ? \n Press 1 for DNA \n Press 2 for Protein:\n' ))

    if dna_or_protein == '1':
        seq1 = str(input('Enter DNA sequence 1:'))
        seq2 = str(input('Enter DNA sequence 2:'))
        g = int(input('Enter the Gap Penalty (preferably an integer):'))
        plot = int(input('Plot the scores and traceback ? \n Press 1 for Yes \n Press 2 for No\n'))
        print()


        M,tracer = local_dp(seq1,seq2,g,sub_matrix='DNA_identity_scores.csv' ,DNA=True)
        main(M,tracer,seq1,seq2,plot)   

    else : 
        seq1 = str(input('Enter Protein sequence 1:'))
        seq2 = str(input('Enter Protein sequence 2:'))
        g = int(input('Enter the Gap Penalty (preferably an integer):'))
        sub_matrix = int(input('Which Substitution matrix do you want to use?\n Press 1 for blosum62 \n Press 2 for blosum50 \n Press 3 for PAM100 \n Press 4 for PAM250 \n'))
        subs_matrix = select_substitution_matrix(sub_matrix)
        plot = int(input('Plot the scores and traceback ? \n Press 1 for Yes \n Press 2 for No \n'))
        print()

        M,tracer = local_dp(seq1,seq2,g,subs_matrix,DNA=False)
        main(M,tracer,seq1,seq2,plot)
            
