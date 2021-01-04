import numpy as np

class HMM_evaluate : 
    """
    Given an input Sequence and a HMM ie the set containing Transition probabilitites,
    Emission probabilitites, No of Hidden States and Emitted states, evaluate the probability 
    of the sequence
    """
    def __init__(self,Transition_prob,Emission_prob,N_states,N_emitted_states):
        self.Transition_prob = Transition_prob
        self.Emission_prob = Emission_prob
        self.N_states = N_states
        self.N_emitted_states = N_emitted_states
        
    
    def forward_algorithm_matrix(self,Sequence,initial_proba):
        
        Tp = self.Transition_prob
        Ep = self.Emission_prob
        NHstates = self.N_states
        NEstates = self.N_emitted_states

        state_M = np.zeros((len(NHstates),len(Sequence[0])+1))

        # Initialization 
        state_M[:,:1] = initial_proba.reshape(len(NHstates),1)

        # Forward Algorithm

        for i in range(len(Sequence.T)):   
            state_M[:,i+1] = (Tp.T * Ep[:,Sequence.T[i]]) @ state_M[:,i]
        
        # Termination 
        return state_M
    
    def evaluate_sequence(self,Sequence,initial_proba):
    
        computed_M = self.forward_algorithm_matrix(Sequence,initial_proba)
        
        return np.sum(computed_M[:,-1])
    