import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def global_dp(seq1,seq2,gap_penalty):

    m,n = len(seq1),len(seq2)
    Sub_score = pd.read_csv('DNA_substitution_scores.csv',index_col='0')

    M = np.zeros((m+1,n+1))  
    M[:,0] = gap_penalty*np.arange(0,m+1)
    M[0,:] = gap_penalty*np.arange(0,n+1)

    tracer = np.zeros((M.shape[0],M.shape[1]))

    for i in range(1,m+1):
        for j in range(1,n+1):
            arr = np.array([M[i-1,j-1] + (S[seq1[i-1]][seq2[j-1]]), (M[i-1,j]+gap_penalty), (M[i,j-1]+gap_penalty)])
            M[i,j] = arr.max()
            tracer[i,j] = np.argmax(arr)+1

    return M, tracer

def conv_list2str(list_obj) :
    string = ''
    return string.join(list_obj)

def global_tracer(M,tracer,seq1,seq2):

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


def traceback(M,tracer,seq1,seq2):

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


    return sum(score),alignment

def plot_align_scores(M,seq1,seq2):

    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8])

    sns.heatmap(M,linewidth=0.9,ax=ax,annot=True)
    ax.set_title('Global Alignment Scores',fontsize=14)
    ax.set_xticklabels(['Gap']+[letter for letter in seq2],fontsize=10)
    ax.xaxis.tick_top()
    ax.set_yticklabels(['Gap']+[letter for letter in seq1],fontsize=10, rotation=45)

    return ax 

def plot_global_tracer(T,seq1,seq2):
    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8])

    sns.heatmap(T,linewidth=0.9,ax=ax)
    ax.set_title('Traceback')
    ax.set_xticklabels(['Gap']+[letter for letter in seq2],fontsize=10)
    ax.xaxis.tick_top()
    ax.set_yticklabels(['Gap']+[letter for letter in seq1],fontsize=10, rotation=45)

    return ax

