from asa159 import rcont2
import numpy as np



def get_random_matrix(c_in):
    #```GNU Lesser General Public License v3.0 code from https://github.com/maclandrol/FisherExact```
    # f2py -c -m asa159 asa159.f90
    #c = array
    # remove empty columns
    empty_cols = (np.where(~c_in.any(axis=0))[0])
    empty_rows = (np.where(~c_in.any(axis=1))[0])
    c = np.delete(c_in, empty_cols, axis=1)
    c = np.delete(c, empty_rows, axis=0)

    key = np.array([False], dtype=bool)
    ierror = np.array([0], dtype=np.int32)
    sr, sc = c.sum(axis=1).astype(np.int32), c.sum(axis=0).astype(np.int32)
    nr, nc = len(sr), len(sc)
    #print(sr, sc)
    n = np.sum(sr)
    replicate=10000
    results = np.zeros(replicate)

    seed=None
    #seed=np.random.seed(123456789)
    # test to see if we can increase wkslimit for neutral sims!!!!
    #wkslimit=5000
    wkslimit=50000
    DFAULT_MAX_TOT = 5000
    # set default maxtot to wkslimit
    if wkslimit < DFAULT_MAX_TOT:
        wkslimit = 50000
    if seed is None:
        try:
            seed = random.SystemRandom().randint(1, 100000)
            seed = np.array([seed], dtype=np.int32)
        except:
            try:
                import time
                seed = int(time.time())
                seed = np.array([seed], dtype=np.int32)
            except:
                seed = 12345
                seed = np.array([seed], dtype=np.int32)

    if n < wkslimit:
        # we can just set the limit  to the table sum
        wkslimit = n
        pass
    else:
        # throw error immediately
        raise ValueError(
            "Limit of %d on the table sum exceded (%d), please increase workspace !" % (DFAULT_MAX_TOT, n))

    maxtot = np.array([wkslimit], dtype=np.int32)
    fact = np.zeros(wkslimit + 1, dtype=np.float, order='F')
    observed = np.zeros((nr, nc), dtype=np.int32, order='F')


    rcont2(nrow=nr, ncol=nc, nrowt=sr, ncolt=sc, maxtot=maxtot,
           key=key, seed=seed, fact=fact, matrix=observed, ierror=ierror)

    # if we do not have an error, make spcial action
    #ans = 0.
    tmp_observed = observed.ravel()
    if ierror[0] in [1, 2]:
        raise ValueError(
            "Error in rcont2 (fortran) : row or column input size is less than 2!")
    elif ierror[0] in [3, 4]:
        raise ValueError(
            "Error in rcont2 (fortran) : Negative values in table !")
    elif ierror[0] == 6:
        # this shouldn't happen with the previous check
        raise ValueError(
            "Error in rcont2 (fortran) : Limit on the table sum (%d) exceded, please increase workspace !" % DFAULT_MAX_TOT)
    else:

        #for empty_column in empty_c:
        #    np.insert(tmp_observed, empty_column, nr, axis=1)
        rndm_matrix = np.reshape(tmp_observed, (nr,nc))
        for empty_column in empty_cols:
            rndm_matrix = np.insert(rndm_matrix, empty_column, 0, axis=1)
        for empty_row in empty_rows:
            rndm_matrix = np.insert(rndm_matrix, empty_row, 0, axis=0)

        return rndm_matrix
        #return np.reshape(tmp_observed, (nr,nc))
