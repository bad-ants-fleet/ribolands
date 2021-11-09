import numpy as np

def get_p8_detbal(A):
    """ Calculate the equilibrium distribution vector p8.

    Given a rate matrix A, calculate the equilibrium distribution using
    detailed balance.  (This is a copy of the treekin routine
    MxEqDistrFromLinSys.)

    """
    dim = len(A)
    p8 = np.zeros(dim, dtype=np.float64)
    p8[0] = 1.

    i, count = 0, 1  # current index, number of solved states
    done = np.zeros(dim, dtype=int)
    while (count != dim):
        # Rewrite entries of p8
        for j in range(dim):
            if not np.isclose(p8[j], 0):
                assert p8[j] > 0
                continue
            if not np.isclose(A[i][j], 0) or not np.isclose(A[j][i], 0):
                assert i == j or A[j][i] > 0
                p8[j] = A[j][i] * p8[i] / A[i][j]
                assert 0 <= p8[j] <= 1 or np.isclose(p8[j], 0) or np.isclose(p8[j], 1)

                if not np.isclose(p8[j], 0):
                    count += 1
        done[i] += 1
        if done[i] > 2:
            # NOTE: This is different from original implementation!!!!
            break
 
        # Determine next i for solving
        for c, p in enumerate(p8):
            if not done[c] and not np.isclose(p, 0):
                i = c
                break
        assert not (i == dim and count < dim), "non-ergodic chain!"
    return p8 / sum(p8) # make the vector sum = 1.0

def symmetrize(A, p8):
    """
    U = _sqrPI * A * sprPI_
    """
    dim = len(A)
    sqrPI_ = np.diag(np.sqrt(p8, dtype = np.float64))
    _sqrPI = np.diag(1/np.sqrt(p8, dtype = np.float64))
    U = np.matmul(np.matmul(_sqrPI, A), sqrPI_)

    # correct for numerical errors
    err = 0
    for (i, j) in zip(*np.triu_indices(dim)):
        if i == j:
            continue
        err += abs((U[i][j] + U[j][i]) / 2 - U[i][j])
        U[i][j] = (U[i][j] + U[j][i]) / 2
        U[j][i] = U[i][j]

    #print(f'Corrected numerical error: {err} ({err/(dim*dim-dim)} per number)')
    return _sqrPI, U, sqrPI_

def decompose_sym(U):
    # Check if matrix has issues
    evals, evecs = np.linalg.eig(U)

    # Sort them by largest eigenvalue
    idx = evals.argsort()[::-1]
    evals = evals[idx]
    evecs = evecs[:,idx]

    return evecs, evals, evecs.T # eves.T == np.linalg.inv(evecs) due to symmetry

def mxprint(A):
    return '\n'.join(' '.join([f'{x:10.4f}' for x in row]) for row in A)

def get_pt_from_file(input_filename, p0string, output_times):

    A = np.loadtxt(input_filename).T
    np.fill_diagonal(A, -np.einsum('ij->j', A))
    dim = len(A)

    p0 = np.zeros(dim)
    for term in p0string:
        p, o = term.split('=')
        p0[int(p) - 1] = float(o)
    assert np.isclose(sum(p0), 1)

    p8 = get_p8_detbal(A)
    _P, U, P_= symmetrize(A, p8)
    _S, L, S_ = decompose_sym(U)
    CL = np.matmul(P_, _S)
    CR = np.matmul(S_, _P)
    assert np.allclose(np.matmul(CL, CR), np.identity(dim))
    tV = CR.dot(p0)
    assert np.allclose(tV, S_.dot(_P.dot(p0)))
    assert isinstance(tV[0], np.float64)

    for t in output_times:
        eLt = np.diag(np.exp(L*t))
        pt = CL.dot(eLt.dot(tV))
        if not np.isclose(sum(pt), 1., atol = 1e-4):
            print(f'# Unstable simulation: {sum(pt)=}')
            break
        yield pt

        #print(f"{t:22.20e} {' '.join([f'{x:8.6e}' for x in pt])}")
        #if np.allclose(pt, p8, rtol = 1e-4, atol = 1e-4):
        #    print(f'# Reached equilibrium at time {t=}.')
        #    break


