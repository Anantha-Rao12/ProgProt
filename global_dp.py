import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from helper_functions import *

def global_dp(seq1,seq2,gap_penalty,sub_matrix,DNA=True):

    m,n = len(seq1),len(seq2)

    sub_score = pd.read_csv(sub_matrix,index_col='0')

    M = np.zeros((m+1,n+1))  
    M[:,0] = gap_penalty*np.arange(0,m+1)
    M[0,:] = gap_penalty*np.arange(0,n+1)

    tracer = np.zeros((M.shape[0],M.shape[1]))

    for i in range(1,m+1):
        for j in range(1,n+1):
            arr = np.array([M[i-1,j-1] + (sub_score[seq1[i-1]][seq2[j-1]]), (M[i-1,j]+gap_penalty), (M[i,j-1]+gap_penalty)])
            M[i,j] = arr.max()
            tracer[i,j] = np.argmax(arr)+1

    return M, tracer

def global_tracer(M,tracer):

    ### Obtain Global traceback route

    global_traceback = np.zeros((M.shape))
    global_traceback[-1,-1] = 1
    tm , tn = M.shape[0] -1 , M.shape[1] -1 
    current_cell = (tm,tn)

    while tm > 0 and tn > 0 : 
        if tracer[current_cell] == 1: 
            current_cell = (tm-1,tn-1)
            global_traceback[current_cell] = 1
            tm -= 1 ; tn -= 1
            
        elif tracer[current_cell] == 2 :
            current_cell = (tm-1,tn)
            global_traceback[current_cell] = 1
            tm -= 1

        elif tracer[current_cell] == 3:
            current_cell = (tm,tn-1)
            global_traceback[current_cell] = 1
            tn -= 1 

    if tm != 0 and tn == 0 : 
        global_traceback[:tm,0] = 1

    if tm == 0 and tn != 0 : 
        global_traceback[0,:tn] = 1


    return global_traceback 


def global_traceback(M,tracer,seq1,seq2):

    """ Obtain Final score and Global ALignment """

    tm , tn = M.shape[0] -1 , M.shape[1] -1 
    score = [M[tm,tn]]
    align1 = []
    align2 = []

    current_cell = (tm,tn)

    while tm > 0 and tn > 0 : 
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

    if tm != 0 and tn == 0 : 
        align1.append(seq1[:tm])
        align2.append('-'*tm)


    if tm == 0 and tn != 0 : 
        align1.append('-'*tn)
        align2.append(seq2[:tn])

    alg1 = conv_list2str(align1[::-1])
    alg2 = conv_list2str(align2[::-1])
    alignment = alg1 + '\n' + alg2


    return score[0],alignment

def main(M,tracer,seq1,seq2,plot):
    T = global_tracer(M,tracer)
    score,g_alignment = global_traceback(M,tracer,seq1,seq2)
    print('\n')

    print('Best Alignment')
    print('Score:',score)
    print(g_alignment)
    print('\n')
    
    if plot == 1:
        input('Press any key to look at the Matrix scores and the Traceback')
        plot_tracer(T,seq1,seq2,Global=True)
        plot_align_matrix(M,seq1,seq2,Global=True) 
        plt.show()
    else : 
        pass


if __name__ == '__main__':
    
    input('Welcome to Global Sequence Alignment \n Press any key to continue')
    dna_or_protein = str(input('Do you want to align a DNA or Protein sequence ? \n Press 1 for DNA \n Press 2 for Protein:\n' ))

    if dna_or_protein == '1':
        seq1 = str(input('Enter DNA sequence 1:'))
        seq2 = str(input('Enter DNA sequence 2:'))
        g = int(input('Enter the Gap Penalty (preferably an integer):'))
        plot = int(input('Plot the scores and traceback ? \n Press 1 for Yes \n Pres 2 for No'))
        print()
        
        M,tracer = global_dp(seq1,seq2,gap_penalty=g,sub_matrix='DNA_substitution_scores.csv',DNA=True)
        main(M,tracer,seq1,seq2,plot)


    else : 
        seq1 = str(input('Enter Protein sequence 1:'))
        seq2 = str(input('Enter Protein sequence 2:'))
        g = int(input('Enter the Gap Penalty (preferably an integer):'))
        sub_matrix = int(input('Which Substitution matrix do you want to use?\n Press 1 for blosum62 \n Press 2 for blosum50 \n Press 3 for PAM100 \n Press 4 for PAM250 \n'))
        subs_matrix = select_substitution_matrix(sub_matrix)
        plot = int(input('Plot the scores and traceback ? \n Press 1 for Yes \n Press 2 for No\n'))
        print()

        M,tracer = global_dp(seq1,seq2,gap_penalty=g,sub_matrix=subs_matrix,DNA=False)
        main(M,tracer,seq1,seq2,plot)
        print()
 
