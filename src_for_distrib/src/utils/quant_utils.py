from functools import partial
from pprint import pprint
import numpy as np
import scipy.stats
import statsmodels.api as sm
import os, sys
import multiprocessing
import hdf_utils

import matplotlib
# Allow saving fig without x forwarding
matplotlib.use("Agg")
from matplotlib import pyplot as plt


def supplement_imputed_vals(data_arr, rep_idx, samp_idx, var_vec, mean_vec):
    '''Supplements the imputed data values to the approriate indices
    of data_arr.

    Args:
    -----
    data_arr : 3d np.array
        Array of shape (R,G,T), R=number of replicates, G=number of 
        genome positions, T=number of sample types. Contains all zeros
        (missing data in this case) at index (missing_idx,:,samp_idx).
        Those data will be replaced with 
        mean_vec[g] + np.random.normal(np.sqrt(var_vec[g])), where g
        is a given genome position.
    rep_idx : int
        Replicate index with missing data of sample type "samp_idx"
    samp_idx : int
        Sample type index with missing data at replicate "rep_idx"
    var_vec : 1d np.array
        Vector containing estimates of variance due to sequencing and
        library prep derived from bootstrap replicates from the sample
        type of interest here.
    mean_vec : 1d np.array
        Vector containing the mean values for the actual data of the
        sample type of interest here.

    Modifies:
    ---------
    data_arr : 3d np.array
        Modified in place to include the imputed values.
    '''

    # sample once from random normal with loc=0, sd=np.sqrt(bs_var[g])
    #   at each genome position
    imputed_noise = np.random.normal(
        0,
        np.sqrt(var_vec),
        len(var_vec),
    )

    imputed_vals = mean_vec + imputed_noise
    data_arr[rep_idx, :, samp_idx] = imputed_vals


def get_jackknife_repweights(data_arr, missing_arr, paired):
    '''Returns array of shape (J, R, T) containing weights to apply
    to each replicate (R) and sample type (T) for each jackknife
    replication (J).

    Args:
    -----
    data_arr : 3d np.array
        Array of shape (R,G,T) containing count data.
    missing_arr : 2d np.array
        Boolean array of shape (R,T) containing information on whether
        each replicate index / sample type index is missing.
    paired : bool
        If True, we have paired replicates. If False, we don't.

    Returns:
    --------
    weights_arr : 3d np.array
        Array of shape (J,R,T) containing replicate weights to be applied
        to all genome positions for a given replicate(R)/sample_type(T) pair
        in jackknife replication J.
    jack_coefs : 1d np.array
        Array of shape (J,) containing jackknife coefficients for each
        jackknife replication.
    '''

    rep_num,type_num = missing_arr.shape

    if paired:
        # here I'm going based on SAS docs for PROC SURVEYFREQ
        #   Even though we only have a single jackknife coefficient for
        #   paired data like this, it helps to repeat it rep_num times
        #   so that we can easily use it in calculating variance of our
        #   mean estimate of our statistic of interest later.
        jack_coefs = np.repeat(
            np.array([float(rep_num-1) / rep_num]),
            rep_num,
        )
        w_orig = 1.0 / rep_num
        # Grab just the first value of jack_coefs, since they're all the
        #   same for the paired case anyway
        include_weight = w_orig / jack_coefs[0]

        # initially, we set up the array of appropriate shape
        #   to be entirely populated by include_weight.
        # then we go through it to re-set excluded replicates
        #   from each jackknife replication to zero
        weights_arr = np.zeros( (rep_num, rep_num, type_num) )
        weights_arr[...] = include_weight

        mask_arr = np.zeros(rep_num, dtype='bool')

        for j in range(rep_num):
            # set j-th index of mask_arr to True
            mask_arr[j] = True
            # set appropriate jackknife replication's replicate indices,
            #   for all sample types, to zero
            weights_arr[j,mask_arr,:] = 0.0
            # re-set the mask to entirely True
            mask_arr[...] = False

    else:
        # the unpaired case is a bit more complicated, since some replicate/
        #   sample_type pairs may still be missing by design. 
        # Overall, we want the still missing to always have weights = 0, no
        #   matter which jackknife replication we're in.
        num_missing = np.sum(missing_arr, axis=0)
        extant_arr = ~missing_arr
        num_extant = np.sum(extant_arr, axis=0)
        # The number of jackknife replications we'll perform is just
        #   the total number of extant replicates.
        jack_rep_count = np.sum(num_extant)

        # calculate the jackknife coefficient and original weights
        #   for each sample type
        stratum_coefs = (num_extant.astype(float) - 1) / num_extant
        # initialize jack_coefs to store each sample type's jackknife
        #   coefficient in a way that will enable simple calculation
        #   of variance of our jackknife-based estimate of our test
        #   statistic later on.
        jack_coefs = np.zeros((jack_rep_count,))
        non_donor_weights = 1.0 / num_extant
        donor_weights = non_donor_weights / stratum_coefs

        extant_rep_idxs,extant_type_idxs = np.where(extant_arr)
        # allocate array of zeros to start with.
        weights_arr = np.zeros( (jack_rep_count, rep_num, type_num) )

        for j in range(jack_rep_count):
            # the j-th idx replicate will be left out
            # To which R/T pair does this j correspond?
            j_rep_idx = extant_rep_idxs[j]
            j_type_idx = extant_type_idxs[j]

            jack_coefs[j] = stratum_coefs[j_type_idx]

            for t in range(type_num):
                # if this is the type of the jackknife replication
                # set all of this type in jackknife j to donor_weights
                #   of the appropriate type, then re-set the removed
                #   replicate of this type's weights to 0
                if t == j_type_idx:
                    weights_arr[j, :, t] = donor_weights[t]
                    weights_arr[j, j_rep_idx, t] = 0.0
                else:
                    weights_arr[j, :, t] = non_donor_weights[t]
        
        missing_rep_idxs,missing_type_idxs = np.where(missing_arr)
        weights_arr[:,missing_rep_idxs,missing_type_idxs] = 0.0

    return weights_arr,jack_coefs


def calc_signed_log10p(vals):
    '''Given a set of values that can be considered z-score like,
    return the -log10 p value for each score in a standard 
    normal distribution
    '''

    # sf yields survival function, which is basically 1-cdf
    pvals_raw = scipy.stats.norm.sf(vals)
    pvals_log = np.log10(pvals_raw)
    pvals_signed = (-1 * pvals_log)
    return pvals_signed


def get_fn_over_axes(inp_mat, iter_axis, fn):
    '''Return the robust z score normalized version of a matrix,
    acting separately on each column.

    Args:
    -----
    inp_mat : 2d or 3d np.array
    iter_axis : array-like
        Which axis, or axes, to iterate over.
        For instance, passing iter_axis=0 will calculate
        for all values in each element of the 0-th axis of inp_mat.
        Passing iter_axis=[0,2] will iterate over each of the 0-th and
        second axes, returning values from fn in each
        slice of the 1-th axis.
    fn : function
        Evaluate this function over desirec axis
    '''

    ax_arr = np.array(iter_axis)
    # if iter_axis was a scalar, we can make the array 1-dimensional here.
    if ax_arr.ndim == 0:
        ax_arr = np.expand_dims(ax_arr, -1)

    out_mat = np.zeros_like(inp_mat)

    # Get number of indices in each axis we're going to iterate over.
    iter_count_arr = np.zeros_like(ax_arr)
    for i,iter_ax in enumerate(ax_arr):
        iter_count_arr[i] = inp_mat.shape[iter_ax]

    total_iters = np.product(iter_count_arr)
    # set up dynamic slicing. Here we make a list of slice objects 
    #   the same length as the number of dimensions of inp_mat
    # slice(None) is the same as a slice like this: ':'
    slc = []
    for i in range(total_iters):
        slc.append([slice(None)] * inp_mat.ndim)

    loop_count = 0
    if len(ax_arr) == 2:
        for i in range(iter_count_arr[0]):
            for j in range(iter_count_arr[1]):
                slc[loop_count][ax_arr[0]] = slice(i, i+1, 1)
                slc[loop_count][ax_arr[1]] = slice(j, j+1, 1)
                loop_count += 1

    elif len(ax_arr) == 1:
        for i in range(total_iters):
            slc[i][ax_arr[0]] = slice(i, i+1, 1)

    else:
        print("ERROR in get_rzscores: function can only handle iter_axis of\
                length 1 or 2.")

    for slice_list in slc:
        # replace the 'iter_axis' element of slc list with a new slice obj
        out_mat[tuple(slice_list)] = fn(
            inp_mat[tuple(slice_list)]
        )

    return out_mat


def calc_rzscores(x):
    '''Return the robust z-score normalized version of the input values.

    Args:
    -----
    x : 1d np.array
        Contains original values to be robust z-transformed

    Returns:
    --------
    z : 1d np.array
        Robust z-score normalized values from x
    '''

    this_median = np.nanmedian(x)
    dev = x - this_median
    this_mad = 1.4826 * np.nanmedian( np.abs( dev ) )
    z = dev / this_mad
    return z

 
def get_weighted_mean_within_jackknife_reps(data_arr, weights_arr):
    '''Multiples data by weights, then sums over axis 1.
    Returns weighted mean.

    Args:
    -----
    data_arr : np.array
        Array of shape (R,G,T), where R is replicate number, G is 
        genome position number, and T is the sample type number.
        Array contains the statistic to which weights will be applied
        to perform jackknife sampling.
    weights_arr : np.array
        Array of shape (J,R,T), where J is the number of jackknife
        replicates, containing jackknife replicate weights.

    Returns:
    --------
    weighted_mean : np.array
        Array of shape (J,G,T) containing the jacknife replicate means
        for the statistic in data_arr.
    '''
    weighted_data = (
        # Insert axis to broadcast over jackknife replications
        np.expand_dims(data_arr, 0)
        # Insert axis to broadcast over genome positions
        * np.expand_dims(weights_arr, 2)
    )
    # summing over axis 1 gives us appropriately weighted mean
    #   of jackknife replicates. data_arr will now be of shape
    #   (J,G,T).
    weighted_mean = np.sum(weighted_data, axis=1)
    return weighted_mean


def calc_lograt_vs_input(data_arr, type_lut, weights_arr=None):
    '''Determines which index of the final axis of a jackknife replicate array
    contains the input samples. Calculates log2(everythong/input).

    Args:
    -----
    data_arr : 3d np.array
        Array of shape (R, G, T), containing the data for each replicate, R,
        genome position, G, and sample type, T.
    type_lut : dict
        Lookup table for determining which index contains input samples
    weights_arr : 3d np.array
        Array of shape (J, R, T), where J is the number of jackknife
        replications. Contains the weights to apply for 
        jackknife resampling. Default is None, in which case we're just
        directly calculating the log2 ratios without jackknifing.

    Returns:
    --------
    log2(all/input) : 3d np.array
        Array containing the log2-ratios of all sample types vs input. 
        Here we leave the result of dividing input by input, even though
        it's all 1, since it will preserve the indexing scheme present in
        type_lut. If we provided replicate weights for jackknifing, then the
        array is of shape (J,G,T). If we did not provide weights, then we're
        not doing jackknifing, and the array is of shape (R,G,T).
    '''
   
    if weights_arr is not None:
        data_arr = get_weighted_mean_within_jackknife_reps(data_arr, weights_arr)

    if "inp" in type_lut:
        input_idx = type_lut['inp']['idx']
    elif "input" in type_lut:
        input_idx = type_lut['input']['idx']
    else:
        sys.exit("ERROR: your input sample MUST be named either 'inp' or 'input', but neither was found")
    # grab input data and append axis to maintain ability to broadcast
    #   over data_arr.
    input_arr = np.expand_dims(data_arr[:,:,input_idx], -1)
    log2_rat = np.log2( data_arr / input_arr )

    return log2_rat


def median_norm(data_arr, targetval=100.0, offset=0.25):
    '''Median normalize data so that within the given vector, the
    new median is equal to targetval.

    Args:
    -----
    data_arr : 3d np.array
        Array of shape (R,G,T). We'll take median across
        genome positions (axis 1) to get each replicate/sample type's
        median coverage.
    targetval : float
        Value to which the normalized median will be set.
    offset : float
        Small pseudocount offset

    Modifies:
    ---------
    data_arr : 3d np.array
       Array's values are now median-normalized. 
    '''

    # calculate medians and insert new axis in middle to make
    #   the median array broadcastable with data_arr
    curr_medians = np.expand_dims(np.median(data_arr, axis=1), 1)
    data_arr *= ((targetval-offset) / curr_medians)
    data_arr += offset

def impute_missing_hdf(data_arr, missing_arr, type_lut,
                    bs_num, paired, force=False, spike_name=None):
    '''Supplements missing data in data_arr with mean + noise 
    of extant replicates.

    Args:
    -----
    data_arr : 3d np.array
        Array of shape (R,G,T) containing this experiments original 
        coverage values. Missing data is a G-sized slice entirely of zeros.
    missing_arr : 2d np.array
        Boolean array indicating which R/T pairs are missing.
    type_lut : dict
        Lookup table with information about each sample type, including
        which replicate indices are missing.
    missing_arr : 2d np.array
        Boolean array of shape (R,T) indicating which replicate/sample types
        are missing.
    bs_num : int
        Number of boostrap samples taken at bootstrapping step. Used to set
        up appropriate size array to store bootstraps.
    paired : bool
        True if we have paired data, False otherwise.
    force : bool
        True if we are forcing imputation of data when unpaired and number
        of extant replicates is only 1. Setting to true is not recommended
        unless you really understand what you're doing. You really must
        go back to do more biological replicates.
    spike_name : str
        Name of the spike-in chromosome.

    Modifies:
    ---------
    data_arr
        Is filled with imputed values as appropriate, depending on whether
        we have paired data or not.
    missing_arr
        Appropriate True falues are set to False if they were imputed.
        All will be False if data are paired, but some True may remain
        if unpaired.
    '''

    for samp_type,samp_info in list(type_lut.items()):

        samp_idx = samp_info['idx']
        # grab True/False values from this sample's column
        these_missing = missing_arr[:,samp_idx]
        all_missing = np.where(these_missing)[0]
        these_extant = ~these_missing
        num_missing = np.sum(these_missing)
        num_extant = np.sum(these_extant)

        # skip sample if no missing replicates
        if num_missing == 0:
            continue

        # If we have paired data, get which replicates of this 
        #   sample type to impute
        if paired:
            impute_reps = all_missing
        # if we have unpaired data and more than 1 extant
        #   replicate of this sample type, skip this sample type.
        # If, however, we have unpaired data and only one extant replicate,
        #   set num_missing to 1 and identify the single R index to 
        #   be imputed. Also throw a big ole warning at the user 
        #   IF THEY HAVE SET THE I KNOW IT'S WRONG BUT DO IT ANYWAY
        #   flag set. If the user doesn't have that flag set, exit the
        #   program with a big ole error message.
        else:
            if num_extant > 1:
                continue
            else:
                if force:
                    print("WARNING:==============================")
                    print("YOU HAVE CHOSEN TO IMPUTE MISSING DATA FROM A SINGLE REPLICATE!! THIS WILL YIELD RESULTS WHICH ONLY ACCOUNT FOR TECHNICAL VARIATION, WHICH IS A VERY SMALL FRACTION OF TOTAL VARIATION! YOU MUST GO BACK TO DO MORE BIOLOGICAL REPLICATES. USE THE RESULTING NUMBERS AT YOUR PERIL!")
                    print("--------------------------------------")
                    num_missing = 1
                    impute_reps = np.array([all_missing[0]])
                else:
                    sys.exit("ERROR: YOU HAVE AT LEAST ONE SAMPLE TYPE WITH ONLY ONE BIOLOGICAL REPLICATE!! THE VALUES OUTPUT FROM THIS PROGRAM WILL ONLY ACCOUNT FOR TECHNICAL VARIATION, WHICH IS A VERY SMALL FRACTION OF TOTAL VARIATION!! GO BACK AND DO MORE BIOLOGICAL REPLICATES!! EXITING THE PROGRAM WITHOUT SAVING ANY OUTPUTS.")

        # set up array containing extant replicate indices
        extant_rep_idxs = np.where(these_extant)[0]

        # set up list of filenames with bootstrap replicates available
        hdf_names = [
            fname
            for rep_idx,fname in samp_info['rep_idx_fname_lut'].items()
        ]
        # now read in the extant replicates' bootstrapped coverages
        # allocate array of zeros to store bootstraps, then replace with data
        bs_array = np.zeros((len(hdf_names), data_arr.shape[1], bs_num))
        
        for i,fname in enumerate(hdf_names):
            concat_data = hdf_utils.concatenate_contig_data(
                fname,
                sample_num = bs_num,
                #NOTE: add the following as parameter to this function
                #  and as cfg option in conf file
                dset_basename = "bs_qnorm",
                positions = data_arr.shape[1],
                #spikein_name = spike_name,
            )
            # Place concatenated contig data into bs_array.
            bs_array[i, :, :] = concat_data
    
        # calculate variance within each extant replicate's bootstraps,
        #   then get mean variance accross bootstraps
        bs_vars = np.var(bs_array, axis=2)
        bs_var = np.mean(bs_vars, axis=0)
    
        extant_data_mean = np.mean(
            data_arr[extant_rep_idxs,:,samp_idx],
            axis = 0
        )

        # bear in mind that this will only replace one of potentially
        #   multiple missing replicates if we have unpaired data.
        # We handle this later.
        for missing_idx in impute_reps:
        
            supplement_imputed_vals(
                data_arr,
                missing_idx,
                samp_idx,
                bs_var,
                extant_data_mean,
            )
            missing_arr[missing_idx, samp_idx] = False


def calc_jackknife_cl(jackrep_stat_arr, jack_coefs, alpha=[0.95]):
    '''Estimate confidence limits using jackknife samples. Using method
    described in documentation of SAS' PROQ SURVEYFREQ (https://documentation.sas.com/?cdcId=pgmsascdc&cdcVersion=9.4_3.4&docsetId=statug&docsetTarget=statug_surveyfreq_details31.htm&locale=en).
    
    Args:
    -----
    jackrep_stat_arr : 3d np.array
        Array of shape (J,G,T), where J is the number of jackknife replications,
        G is the number of genome positions considered, and T is the number of
        sample types. Each slice of the J axis contains the test statistic for
        each position, G, over the entire genome, for each sample type, 
        for a given jackknife replication.
    jack_coefs : 1d np.array
        Array of shape (J,). Contains jackknife coefficient for
        each sample type.
    alpha : list
        List of confidence interval widths to return. If 95% CI is desired, 
        use [0.95] (the default). If you want 80% and 90% CIs, use
        [0.8, 0.9].
    
    Returns:
    --------
    jacked_mean : 2d np.array
        Shape (G,T); contains jackknife-based estimate of mean statistic.
    cl_lo : 3d np.array
        Shape (G,T,A); contains jackknife-based lower confidence limit for each
        genome position (G), sample type (T), and confidence level in alpha (A).
    cl_hi : 3d np.array
        Same as cl_lo above, but for upper confidence limit.
    '''

    # caluclate mean over jackknife axis
    jacked_mean = np.mean(jackrep_stat_arr, axis=0)
    # get squared deviation of each jack rep from mean estimate
    sq_dev = np.power(
        jackrep_stat_arr - np.expand_dims(jacked_mean, 0),
        2
    )
    # multiply squared dviation by each jack rep's jackknife coefficient
    weighted_sq_dev = jack_coefs[:,np.newaxis,np.newaxis] * sq_dev
    # sum the weighted squared deviations over all jackknife replicates.
    #   The result is the jackknife-based variance of our estimate
    #   of the statistic of interest.
    jack_var = np.sum(weighted_sq_dev, axis=0)
    jack_se = np.sqrt(jack_var)

    lower = [(1-a)/2 for a in alpha]
    upper = [1-l for l in lower]
    ci = [scipy.stats.norm.ppf(u) * jack_se for u in upper]
    # cl_up and cl_lo will be 3d arrays after stacking. If there's only
    #   one alpha value the final axis will be of length 1, but if there
    #   are more alpha values, the final axis will be the length of
    #   the list supplied to the alpha argument.
    cl_up = [jacked_mean + c for c in ci]
    cl_up = np.stack(cl_up, axis=-1)
    cl_lo = [jacked_mean - c for c in ci]
    cl_lo = np.stack(cl_lo, axis=-1)

    return jacked_mean,cl_lo,cl_up


def get_chipsub_arglist(log_rats, type_lut, numerator_list,
                        plot_diagnostics, chipsub_percentile = 98,
                        slope_threshold = 0.95):
    '''Provides a list of arguments for each chipsub operation to be
    performed. This enables using multiprocessing.Pool.starmap.

    Args:
    -----
    log_rats : np.array
        Array of shape (R,G,T), where R is either sample replicate count
        or jacknife replicate count, G is genome position number, and T
        is the number of sample types in this experiment.
    type_lut : dict
        Dictionary with information on each sample type. Used here for
        proper indexing of log_rats array.
    numerator_list : list
        List of strings, each string corresponds to a sample type that
        will have the chip trend subtracted from it.
    plot_diagnostics : bool
        If True, we plot useful plots to diagnose potential issues in
        chip subtraction.
    chipsub_percentile : int
        Between 0 and 100, sets the percentile of chip data below which
        none will be considered in linear regression to fit the
        contribution of chip signal to ipod signal.
    slope_threshold : float
        Between 0 and 1, sets the point at which we stop incrementing
        the ipod ~ chip slope upward. i.e., once the fraction of the
        data indicated by this argument falls below the fitted line,
        stop increasing the slope.

    Returns:
    --------
    chipsub_arg_list : list
        List of tuples. Each tuple contains all the arguments needed to run
        calc_chipsub.
    '''

    chipsub_arg_list = []
    # Loop over replicates
    for rep_idx in range(log_rats.shape[0]):
        # Grab data for this replicate
        replicate_data = log_rats[rep_idx,...]
        # Get chip data for this replicate
        chip_lograt_vec = replicate_data[:,type_lut['chip']['idx']]
        # Loop over numerators for chipsub calculation

        for numer_name in numerator_list:
            # Grab numerator data for this replicate
            numer_idx = type_lut[numer_name]['idx']
            numer_lograt_vec = replicate_data[:,numer_idx]
            # place appropriate stuff into tuple, appended to arg list
            chipsub_arg_list.append((
                numer_lograt_vec,
                chip_lograt_vec,
                rep_idx,
                numer_idx,
                numer_name,
                plot_diagnostics,
                chipsub_percentile,
                slope_threshold
            ))

    return chipsub_arg_list


def do_linear_interpolation(xvals, slope, padding):
    '''Return the predicted value for each point in xvals based on a linear
    regression with slope slope and intercept 0, plus padding. We return
    0 for negative values of x.
    '''

    predvals = xvals * slope + padding

    return (predvals > 0) * predvals


def do_linear_fit(xdat, ydat, minval, thresh=0.95):
    '''Do linear regression on the selected data where xdat > minval
    and file and return the minimal slope for which 95% of the fitted
    data are below the line.

    Args:
    -----
    xdat : np.array
        1d numpy array containing log2(Chip/input) values.
    ydat : np.array
        1d numpy array containing log2(ipod/input) values.
    minval : float
        Minimum log2(ChIP/input) value to consider for regression.
    thresh : float
        Sets the threshold for when to stop increasing the slope
        Fit to the data for chip subtraction. This fraction of
        the data are below the fit line, stop increasing the slope.
        Must be between 0 and 1. Default is 0.95.

    Returns:
    --------
    v : np.array
        The slope of the line to subtract from log2(ipod/input) data
        to remove the RNAP contribution to the protein occupancy signal.
    '''
    
    # filter data, keeping values for which x > minval
    goodflags = (xdat > minval)
    X = xdat[goodflags]
    y = ydat[goodflags]

    # Use statsmodels' OLS for linear fit.
    #   sm.OLS will not fit an intercept here, i.e., intercept = 0,
    #   because x is just a 1d array. sm.OLS expects x to be a design
    #   matrix, so if the first column of x is entirely 1, it is used
    #   to fit the "constant", or sm.OLS-speak for intercept.
    res = sm.OLS(endog=y, exog=X).fit()
    slope_init = res.params[0]
    print("==============================================")
    print("Initial chipsub slope: {}".format(slope_init))
    print("BASED ON THIS VALUE, CONSIDER WHETHER IT IS APPROPRIATE TO PERFORM CHIP SUBTRACTION. LOOK AT THE DIAGNOSTIC PLOTS.")
    # increase initial slope by 1e-5 until > 95% of y values and under line
    while np.sum(y > (slope_init * X)) > ((1-thresh) * len(y)):
        slope_init += 1e-4
    
    v = slope_init
    print("Final chipsub slope: {}".format(v))
    return np.asarray(v)


def calc_chipsub(numer_lograt_vec, chip_lograt_vec,
                  rep_idx, numer_idx,
                  numer_name, plot_diagnostics,
                  minperc = 98, slope_threshold=0.95):

    print("-----")
    print("Subtracting RNAP contribution to IPOD signal from {}, replicate {}.".format(
        numer_name.upper(),
        rep_idx + 1)
    )

    minval = scipy.stats.scoreatpercentile(chip_lograt_vec, minperc)

    fit_slope = do_linear_fit(
        chip_lograt_vec,
        numer_lograt_vec,
        minval,
        thresh = slope_threshold,
    )

    predvals = do_linear_interpolation(
        chip_lograt_vec,
        fit_slope,
        padding = 0,
    )

    newvals = (
        numer_lograt_vec # the log2(sample type/input) data
        - np.fmax(predvals, 0.0) # subtract the max of zero or predval
    )

    print("Finished subtracting RNAP contribution to IPOD signal from {}, replicate {}.".format(numer_name.upper(), rep_idx + 1))
    print("======")

    if plot_diagnostics and (rep_idx == 0):
        save_diagnostics(
            numer_lograt_vec,
            chip_lograt_vec,
            numer_name,
            predvals,
            newvals
        )

    return (rep_idx, numer_name, newvals)


def do_chipsub(log_arr, 
                 type_lut, chipsub_lut, numerator_list=['ipod'],
                 plot_diagnostics=True, chipsub_percentile=98,
                slope_threshold=0.95,
                nproc=1):
    '''For each replicate, performs regression of log(ipod/input) ~
    log(chip/input) for the top 2 percent highest log(chip/input) signal.
    Uses result of that regression to subtract the RNAP occupance contribution
    to the ipod signal.

    Args:
    -----
    log_arr : 3d np.array
        Array of shape (R, G, T) containing log2(T/input) numbers. R could
        be paired replicate count for calculating empirical chipsub, or
        could represent the number of jackknife replicates.
    type_lut : dict
        Dictionary. Keys are sample type names, values are themselves
        dictionaries with information about each sample type. Used
        in this function for indexing log_arr to get correct numerator
        and denominator for log2(sample/input) - log2(ChIP/input).
    chipsub_lut : dict
        Dictionary with info to help with indexing chipsub data in out_mat.
    numerator_list : list
        List of strings denoting which sample types to consider numerator
        in these calculations. If a string is passed, it's converted to
    plot_diagnostics : bool
        If True, some diagnostic plots are saved.
    chipsub_percentile : int
        Must be between 0 and 100. This will set the initial amount
        of chip data that will be considered for ipod ~ chip regression
        for chip subtraction.
    slope_threshold : float
        Must be between 0 and 1. Will keep incrementing ipod ~ chip slope
        upward until this fraction of the data is below the line.
    nproc : int
        Number of cores to use for multiprocessing. Default is 1.

    Returns:
    --------
    out_mat : 3d np.array
        Array of shape (R, G, len(numerator_list)). If input is sample reps,
        we're just calculating empirical chipsub using paired data. In
        that case, R is the number of paired replicates. If above input is 
        jackknife reps, R represents the numer of jackknife replications
        we performed. out_mat contains values from subtraction of 
        RNAP chip contribution to the data of interest
        listed in numerator_list.
    '''

    out_mat = np.zeros((
        log_arr.shape[0],
        log_arr.shape[1],
        len(numerator_list)
    ))

    #NOTE: I really should put minval as argument and conf option
    chipsub_arg_list = get_chipsub_arglist(
        log_arr,
        type_lut,
        numerator_list,
        plot_diagnostics,
        chipsub_percentile,
        slope_threshold,
    )

    proc_pool = multiprocessing.Pool(nproc)
    #chipsub_list = []
    #for stuff in chipsub_arg_list:
    #    chipsub_list.append(calc_chipsub(*stuff))
    chipsub_list = proc_pool.starmap(
        calc_chipsub,
        chipsub_arg_list,
    )
    proc_pool.close()
    proc_pool.join()

    for res in chipsub_list:
        samp_idx = chipsub_lut[res[1]]['idx']
        rep_idx = res[0]
        out_mat[rep_idx, :, samp_idx] = res[2]

    return out_mat


def save_diagnostics(numer_vec, chip_vec, numer_name,
                     predvals, newvals):
    '''Saves plots demonstrating the effect of subtraction of
    the trend in association between log2(numerator/input) vs
    log2(ChIP/input).
    '''

    plt.figure()
    plt.rcParams.update({'font.size': 22})
    plt.hexbin(
        chip_vec,
        numer_vec,
        bins='log',
        cmap=plt.get_cmap("Blues"),
    )
    plt.xlabel('Log$_2$ ChIP vs. input')
    plt.ylabel('Log$_2$ {} vs. input'.format(numer_name.upper()))
    plt.tight_layout()
    plt.savefig("{}_noline.png".format(numer_name.upper()))

    chip_sort_idxs = np.argsort(chip_vec)
    plt.plot(
        chip_vec[chip_sort_idxs],
        predvals[chip_sort_idxs],
        'r--',
        linewidth=2.0,
    )
    plt.tight_layout()
    plt.savefig("{}_main.png".format(numer_name.upper()))

    plt.figure()
    plt.hexbin(
        chip_vec,
        newvals,
        bins = 'log',
        cmap = plt.get_cmap("Blues"),
    )
    plt.colorbar()
    plt.xlabel('Log$_2$ ChIP vs. input')
    plt.ylabel('ChIP-subtracted {}-HR'.format(numer_name.upper()))
    plt.tight_layout()
    plt.savefig("{}_chipsub.png".format(numer_name.upper()))


def get_chipsub_lut(type_lut, numerator_list=['ipod']):
    '''Make and return a lookuptable similar to type_lut, but only
    for chipsub numerators.
    '''

    if isinstance(numerator_list, str):
        numerator_list = [numerator_list]

    chipsub_lut = {}
    for numer_out_idx,numerator in enumerate(numerator_list):
        # add numerator name as key, index as val
        chipsub_lut[numerator] = {'idx' : numer_out_idx}
        chipsub_lut[numerator]['orig_idx'] = type_lut[numerator]['idx']
        chipsub_lut[numerator]['rep_idx_fname_lut'] = (
            type_lut[numerator]['rep_idx_fname_lut']
        )
    return chipsub_lut

def gather_norm_data(norm_lut, paired):
    '''Gathers important information from a dictionary containing
    data and log-ratio info for the two normalization types, quantile
    and spike-in.

    Args:
    -----
    norm_lut : dict
        Contains a bunch of information for each normalization type.
        Normalization types are keys, values are information for that
        given normalization type.
    chipsub_numerators : list

    Returns:
    --------
    type_lut : dict
    log_rats : np.array
    jacked_lograts : np.array
    weights_arr : np.array
    jack_coefs : np.array
    '''

    type_lut = {}
    type_idx = 0
    lograts = []
    jacked_lograts = []
    
    for norm_type,norm_info in norm_lut.items():

        if not norm_info['type_lut']:
            continue

        ctg_lut = norm_info['ctg_lut']
        rev_ctg_lut = norm_info['rev_ctg_lut']
        res = norm_info['res']
        weights_arr = norm_info['weights_arr']
        jack_coefs = norm_info['jack_coefs']

        this_type_lut = norm_info['type_lut']
        if paired:
            these_lograts = norm_info['log_rats']
        these_jacked_rats = norm_info['jacked_log_rats']

        for sample_type,type_info in this_type_lut.items():

            if sample_type not in type_lut:

                # get this norm type / sample types info
                type_lut[sample_type] = type_info
                # we'll use this index to grab data for this sample type
                this_idx = type_info['idx']

                # set appropriate index for gathered type lookup table
                type_lut[sample_type]['idx'] = type_idx
                type_idx += 1

                # grab data of interest for this sample type/norm type
                if paired:
                    lograts.append(these_lograts[:,:,this_idx])
                jacked_lograts.append(these_jacked_rats[:,:,this_idx])
                
    if paired:
        log_rats = np.stack(lograts, -1)
    jacked_log_rats = np.stack(jacked_lograts, -1)

    if paired:
        ret_vals = (
            type_lut,
            log_rats,
            jacked_log_rats,
            weights_arr,
            jack_coefs,
            ctg_lut,
            rev_ctg_lut,
            res,
        )
    else:
        ret_vals = (
            type_lut,
            jacked_log_rats,
            weights_arr,
            jack_coefs,
            ctg_lut,
            rev_ctg_lut,
            res,
        )

    return ret_vals            

def set_up_data_from_hdf2(norm_lut, conf_dict, bs_dir, pat):
    
    for norm_method,norm_info in norm_lut.items():

        type_lut = norm_info['type_lut']
        if not type_lut:
            continue
        dset_base = norm_info['dset']
        spike_name = norm_info['spikein_name']

        # type_lut is modifiecd in place
        data_arr,missing_arr,ctg_lut,res = set_up_data_from_hdf(
            type_lut,
            conf_dict,
            bs_dir,
            pat,
            dset_base,
            #spike_name,
        )

        # remove the spike-in chromosome data here since we filter
        #  them out in the set_up_data_from_hdf step
        #if spike_name != "None":
        #    del ctg_lut[spike_name]
        norm_info['data_arr'] = data_arr
        norm_info['missing_arr'] = missing_arr
        norm_info['ctg_lut'] = ctg_lut
        norm_info['res'] = res

def set_up_data_from_hdf(type_lut, conf_dict, bs_dir, pat,
                         norm_dset_base, spike_name=None):
    '''This function generates a 3d array containing the data for this
    experiment. The array is of shape (R,G,T), where R is the number
    of replicates, G is the number of genome positions, and T is the
    number of distinct sample types, i.e., for an experiment with input,
    IPOD, and ChIP, T=3. 

    Args:
    -----
    type_lut : dict
        Dictionary for mapping sample type names to information about samples
        of that type. Keys are sample type names, values are dictionaries.
    conf_dict : dict
        Dictionary containing configuration information for this experiment
    bs_dir : str
        Name of directory containing bootstrap and original data.
    pat : compiled regex pattern object
        Regular expression with a group to return the rep number in the
        resulting match object
    norm_dset_base : str
        basename of normalized coverage dataset

    Returns:
    --------
    data_arr : 3d np.array
        Array of shape (R,G,T), see above for details on dimension meanings.
        Any sample type that is missing a replicate will have entirely
        zero values, and its row/type index is recorded in missings.
    missing_arr : 2d np.array
        A boolean array of shape (R,T), where R is the number of replicates
        and T is the number of sample types. If a given sample type has all
        replicates, all values in its column will be False. If, for a given
        sample type, replicate 1 (rep_idx 0) is missing, then the 0th index
        of that type's column will be True.
    ctg_lut : dict
        Dictionary for looking up contigs.
    resolution : int
        Resolution for this experiment.

    Modifies:
    ---------
    type_lut : dict
        A dictionary containing information about the sample types in this
        experiment.
    '''

    rep_num = 0
    # compile regex to search for the replicate number 
    fname_dict = {}
    type_count = len(type_lut)

    # loop over each sample type, i.e., input, IPOD, ChIP, etc.
    #   modify rep_num if a sample type has a higher replicate id than
    #   we've previously seen. rep_num only really matters here for
    #   paired data.
    for samp_type,samp_info in list(type_lut.items()):

        these_reps = conf_dict[samp_type]["sample_names"]
        samp_dir = conf_dict[samp_type]["directory"]

        # Look through reps in this sample type.
        #   We append sample index and rep index to their
        #   respective lists, to later set the appropriate
        #   indices of pairs_arr to True
        samp_info['rep_idx_fname_lut'] = {}
        for rep_name in these_reps:

            mo = pat.search(rep_name)
            if mo is None:
                sys.exit("ERROR: your replicate files names must have the string \"rep<N>\" in them somewhere so that this program can find your replicates. Here, substitute and integer for <N>.")
            # get replicate id and adjust rep_num if it's higher
            #  than what's currently recorded
            rep_id = int(mo.groups()[0])
            if rep_id > rep_num:
                rep_num = rep_id
            rep_idx = rep_id - 1
            
            # Here, we're adding the appropriate tuple of
            #   (replicate index, sample type index) as a key to fname_dict,
            #   the value for which is the appropriate file name.
            # This will enable simple index-based file reading when acutally
            #   loading in data.
            hdf_name = os.path.join(
                os.getcwd(),
                samp_dir,
                bs_dir,
                rep_name + ".hdf5"
            )
            # rep_idx as key and fname as val
            samp_info['rep_idx_fname_lut'][rep_idx] = hdf_name

    # get number of genome positions in supercontig from single hdf5 file
    position_count = hdf_utils.calc_supercontig_posnum(hdf_name)#, spike_name)
    
    # now that we know the max rep_num for any given sample type, we
    #   can allocate our data array to the appropriate size
    data_arr = np.zeros((rep_num, position_count, type_count))

    for samp_type,samp_info in list(type_lut.items()):
        type_idx = samp_info['idx']
        rep_idx_fname_lut = samp_info['rep_idx_fname_lut']

        for rep_idx,rep_fname in list(rep_idx_fname_lut.items()):
            concat_data = hdf_utils.concatenate_contig_data(
                rep_fname,
                norm_dset_base,
                positions = position_count,
                #spikein_name = spike_name,
            )
            # Place concatenated contig data into data_arr.
            # Note that concatenate_contig_data function returns a "long"
            #   array of shape (G,1), where G is genome position number,
            #   so here we slice away that trailing axis.
            data_arr[ rep_idx,:,type_idx ] = concat_data[:,0]

    # get the rep/type indices which are missing
    miss = np.where(np.all(data_arr == 0, axis=1))
    # allocate array of False to store missingness info
    missing_arr = np.zeros((rep_num, type_count), dtype='bool')
    missing_arr[miss] = True
    ctg_lut = hdf_utils.get_ctg_lut(hdf_name)
    resolution = hdf_utils.get_resolution(hdf_name)

    return data_arr,missing_arr,ctg_lut,resolution


def make_name_arrs(dat_arr, type_lut, out_pref, info_str, pat = None):
    '''Given shape of dat_arr and info in type_lut, create array of 
    strings with appropriate output file names for each replicate/type.

    Args:
    -----
    dat_arr : 2d or 3d np.array
        Shape (G,T) if 2d or (R,G,T) if 3d.
    type_lut : dict
        Dictionary with info about each sample type.
    out_pref : str
        String to prifix to all file names.
    info_str : str
        Component of file name with information about what this
        sample actually is. For a given filename, we organize
        its name in the following schema: <1>_<2>.<3>, where <1>
        is substitued for the out_pref above, <2> is substituted for
        the contents of info_str, and <3> is substituted for either .npy
        or .gr.
    pat : Compiled regex pattern object
        If None (the default), not saving each replicate's info.
        Pattern to search if we're saving paired replicates. Will find
        appropriate number for each replicate.

    Returns:
    --------
    fname_arr : 2d np.array
        Sape (R,T). A file name of up to 200 characters for each 
        replicate/sample_type.
    '''

    # having a pattern indicates we're saving paired replicates' data
    if pat is not None:
        # initialize array to hold file names up to 200 characters long
        fname_arr = np.empty(dat_arr[:,0,:].shape, dtype='|S200')
        fname_arr[...] = ''
        dset_name_arr = np.empty(dat_arr[:,0,:].shape, dtype='|S200')
        dset_name_arr[...] = ''
        for rep_idx in range(dat_arr.shape[0]):
            for samp_type,samp_info in list(type_lut.items()):
                samp_idx = samp_info['idx']
                samp_fname_lut = samp_info['rep_idx_fname_lut']

                if not rep_idx in samp_fname_lut:
                    continue
                mo = pat.search(samp_fname_lut[rep_idx])
                rep_id = int(mo.groups()[0])

                out_info = info_str.format(samp_type.upper(), rep_id)
                # substitute suffix later
                fname = '{}_{}.{{}}'.format(out_pref, out_info)

                fname_arr[rep_idx, samp_idx] = fname
                dset_name_arr[rep_idx, samp_idx] = out_info

    # if pat is None, we're saving summaries
    else:
        fname_arr = np.empty(dat_arr.shape[-1], dtype='|S200')
        fname_arr[...] = ''
        dset_name_arr = np.empty(dat_arr.shape[-1], dtype='|S200')
        dset_name_arr[...] = ''
        for samp_type,samp_info in list(type_lut.items()):
            samp_idx = samp_info['idx']
            out_info = info_str.format(samp_type.upper())
            fname = '{}_{}.{{}}'.format(out_pref, out_info)
            fname_arr[samp_idx] = fname
            dset_name_arr[samp_idx] = out_info

    return fname_arr,dset_name_arr


def write_outs(in_arr, type_lut, out_prefix,# out_hdf_name,
                info_str, pat=None, spike_name=None):
    '''Writes information in in_arr to an hdf5 file and to a gff file.
    '''

    # set up lookup table for getting sample type from index
    reverse_type_lut = {
        type_lut[samp_type]['idx'] : samp_type
        for samp_type in type_lut.keys()
    }

    # get array of file names. Each still needs a suffix appended
    fname_arr,dset_name_arr = make_name_arrs(
        in_arr,
        type_lut,
        out_prefix,
        info_str,
        pat,
    )

    # a 2D fname array indicates we're doing replicate-level writes
    if fname_arr.ndim == 2:
        for i in range(fname_arr.shape[0]):
            for j in range(fname_arr.shape[1]):
                # get this sample's type, i.e., inp, chip, ipod...
                samp_type = reverse_type_lut[j]
                samp_info = type_lut[samp_type]
                # skip input
                if samp_type in ['inp','input']:
                    continue
                fname = fname_arr[i,j]
                # Skip missing data
                if fname == b'':
                    continue
                in_hdf_name = samp_info['rep_idx_fname_lut'][i]
                ctg_lut = hdf_utils.get_ctg_lut(in_hdf_name)
                res = hdf_utils.get_resolution(in_hdf_name)
                # grab every genome position's data for this rep/type
                out_data = in_arr[i,:,j]
                # get BEDGraphRecord object from data and hdf file
                print("===================")
                print("Writing to {}".format(fname.decode().format("bedgraph")))
                bg_record = hdf_utils.write_bedgraph(
                    out_data,
                    in_hdf_name,
                    fname.decode().format("bedgraph"),
                    #spikein_name = spike_name,
                )
                print("-------------------")

                # make that dataset name the base file name
                #   without the trailing '.'
                #dset_name = dset_name_arr[i,j].decode()
                #rep_id = dset_name.split("_")[-1]

                #grp_name = "contigs/{{}}/{}/replicates".format(
                #    samp_type,
                #)
                # Now write data to group as appropriately-named dataset
                #hdf_utils.decatenate_and_write_supercontig_data(
                #    out_hdf_name,
                #    out_data,
                #    dset_name,
                #    dtype = out_data.dtype,
                #    grp_fmt_str = grp_name,
                #)
    # a 1D fname array indicates we're doing summary-level writes
    elif fname_arr.ndim == 1:
        for i,fname in enumerate(fname_arr):
            # get this sample-s type
            samp_type = reverse_type_lut[i]
            samp_info = type_lut[samp_type]
            # skip input
            if samp_type == 'inp':
                continue
            # skip missing data
            if fname == b'':
                continue
            # NOTE: grabbing the 0-th index here is fine. Just using
            #   The hdf file here for getting ctg_lut and res.
            in_hdf_name = [
                name for name in 
                samp_info['rep_idx_fname_lut'].values()
            ][0]
            ctg_lut = hdf_utils.get_ctg_lut(in_hdf_name)
            res = hdf_utils.get_resolution(in_hdf_name)

            out_data = in_arr[:,i]

            # get BEDRecord object from data and hdf file
            print("===================")
            print("Writing to {}".format(fname.decode().format("bedgraph")))
            bg_record = hdf_utils.write_bedgraph(
                out_data,
                in_hdf_name,
                fname.decode().format("bedgraph"),
                #spikein_name = spike_name,
            )
            print("-------------------")

            # make the dataset name the base file name
            #   without the training '.'
            #dset_name = dset_name_arr[i].decode()

            #grp_name = "contigs/{{}}/{}".format(
            #    samp_type,
            #)
            ## Now write data to group as appropriately-named dataset
            #hdf_utils.decatenate_and_write_supercontig_data(
            #    out_hdf_name,
            #    out_data,
            #    dset_name,
            #    dtype = out_data.dtype,
            #    grp_fmt_str = grp_name,
            #)

    else:
        print("ERROR in fname creation!!")


