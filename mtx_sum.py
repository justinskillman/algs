import numpy as np
import math as math
import time as time

def get_combins(mtx):
    binrows = int((math.factorial(mtx.shape[0])) / (math.factorial(mtx.shape[0] - mtx.shape[1])))
    bin = np.zeros((binrows, mtx.shape[1]))
    bin[:,0] = np.repeat(mtx[:, 0], binrows/mtx.shape[0])
    for jc in range(mtx.shape[0]):
        rows = int((math.factorial(mtx.shape[0]-1) / math.factorial(mtx.shape[1]- mtx.shape[1])))
        beg_row = jc*rows
        end_row = beg_row + rows
        mtx_temp = np.delete(mtx, jc, axis=0)
        mtx_temp = np.delete(mtx_temp, 0, axis=1)
        if mtx_temp.shape[1]>1:
            bin[beg_row:end_row,1:mtx.shape[1]] = get_combins(mtx_temp)
        else:
            bin[beg_row:end_row, 1:mtx.shape[1]] = mtx_temp
    return bin

def main():
    for jt in range(2,11):
        t0 = time.clock()
        llim = 1
        hlim = 1000
        mtx_size = (jt,jt)
        mtx = np.random.uniform(llim, hlim, mtx_size).round()
        combins = get_combins(mtx)
        comb_sums = combins.sum(axis=1)
        max_vec = combins[comb_sums.argmax(),:]
        t1 = time.clock()
        print(jt," by ", jt, "matrix takes: ", t1-t0, " seconds")




