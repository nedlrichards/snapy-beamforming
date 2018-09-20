import numpy as np
import sparse_beam

def snap_beamform(sparse_cols, sparse_vals, threshold, beam_taus,
        num_bmax=200000, num_lab_max=5000):


    # make delays into contigous array
    all_inds = np.array(np.hstack(sparse_cols) + 1, dtype=np.int32)
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

    label_val = label_val[: nlabels]
    labels = labels[:, : nlabels]
    # remove 1 indexing
    labels[[0, 1], :] -= 1

    return labels.astype(np.int_, order='C'), label_val.astype(np.float_)
