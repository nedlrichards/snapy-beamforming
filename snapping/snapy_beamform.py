import numpy as np
import sparse_beam

def snapy_beamform(sparse_cols, sparse_vals, threshold, beam_taus,
        num_bmax=200000, num_lab_max=5000):
    """Python wrapper to fortran snapy beamforming routine

    N will be used to indicate number of data channels input to beamformer

    sparse_cols: list of length N, column index of each values
       Each list entry is a seperate channel
       Each channel is a numpy array and will be cast to int32
    sparse_values: list of length N, values
       Each list entry is a seperate channel
       Each channel is a numpy array and will be cast to double

    sparse_cols and sparse_values should have same length, and each entry
    should be the same size

    Beamformer specifications:

    threshold: double, beamformer threshold
    beam_taus: M x N or L x M x N matrix of delay values

    Preallocation variables:
    num_bmax: number of beams exceding the threshold before clustering
    num_lab_max: number of labeled clusters
    """

    # make delays into contigous array
    if len(sparse_cols) != len(sparse_vals):
        raise(ValueError('length of indicies and value list do not match'))

    all_inds = np.array(np.hstack(sparse_cols) + 1, dtype=np.int32)
    sparse_vals = np.array(np.hstack(sparse_vals), dtype=np.float64)

    if all_inds.size != sparse_vals.size:
        raise(ValueError('number of indicies and values do not match'))

    chan_ind = [a.size for a in sparse_cols]
    chan_ind = np.cumsum(chan_ind) + 1
    chan_ind = np.array(np.hstack([1, chan_ind[: -1]]), dtype=np.int32)

    if len(beam_taus.shape) == 2:
        ncols = 1
        nrows = beam_taus.shape[0]
    elif len(beam_taus.shape) == 3:
        (ncols, nrows, _) = beam_taus.shape
        beam_taus = np.reshape(beam_taus, (ncols * nrows, -1))
    else:
        raise(ValueError('beamformer delays can have 2 or 3 dimensions'))

    # transpose for Fortran order
    beam_taus = beam_taus.astype(np.int32).T

    labels, label_val, nlabels = sparse_beam.snapy(all_inds,
                                                   sparse_vals,
                                                   chan_ind,
                                                   threshold,
                                                   beam_taus,
                                                   num_bmax,
                                                   num_lab_max,
                                                   nrows,
                                                   ncols)

    import ipdb; ipdb.set_trace()
    if nlabels == -1:
        raise(ValueError('Preallocate more beams'))
    if nlabels >= num_lab_max:
        raise(ValueError('Preallocate more labels, at least %i labels needed'%
                         nlabels))

    label_val = label_val[: nlabels]
    labels = labels[:, : nlabels]
    # remove 1 indexing
    labels[[0, 1], :] -= 1

    return labels.astype(np.int_, order='C'), label_val.astype(np.float_)
