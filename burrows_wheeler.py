import numpy as np

class BWT:
    """ Implement the Burrow Wheeler Transform for a given string """

    def __init__(self,string:str) -> str:

        self.string = string
        self.processed_string = string + '$'
        self.length = len(self.string)
        self.processed_length = len(self.processed_string)

    def permute_matrix(self,lexical_sort) -> np.ndarray:
        """Return the permuted matrix for the given string"""
        m = self.processed_length
        char_list = [i for i in self.processed_string]
        p_array = np.array(np.zeros([m,m]),dtype=object)

        for i in range(m):
            new_items = char_list[1:] + char_list[:1]
            char_list = new_items
            p_array[i,:] = new_items
        
        p_words = np.sum(p_array,axis=1)

        if lexical_sort == False:
            p_words = p_words.reshape(m,1)
        
        else :
            p_vec = p_words.reshape(m,)
            p_vec.sort()
            p_words = p_vec.reshape(m,1)

        return p_words

    def bwt(self):

        p_mat = self.permute_matrix(lexical_sort=True)
        bwt_char = [word[0][-1] for word in p_mat]
        bwt = ''.join(bwt_char)

        return bwt

    def ibwt_matrix(self):

        m = self.length
        char_list = np.array([i for i in self.string])
        sorted_bwt = np.array(sorted(char_list))

        for i in range(m-1):
            arr = np.array(np.vstack((char_list,sorted_bwt)),dtype=object)
            arr_sum = np.sum(arr,axis=0)
            sorted_bwt = np.array(sorted(arr_sum))
        
        ibwt_matrix = sorted_bwt.reshape(m,1)
        
        return ibwt_matrix

    def ibwt(self):

        ibwt_mat = self.ibwt_matrix()
        last_col = [i[0][-1] for i in ibwt_mat]
        carot_index = last_col.index('$')

        ibwt = ibwt_mat[carot_index][0]

        return ibwt