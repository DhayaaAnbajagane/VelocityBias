import numpy as np
from scipy import interpolate
import os

avail_sims = ['TNG', 'MGTM', 'BM', 'The300']

avail_params = ['mean', 'slope', 'scatter']

parameter_to_index = {'mean': 2, 'slope': 9, 'scatter': 16}

cwd_path = os.path.dirname(__file__)

Om = {'TNG' : 0.3089, 'MGTM' : 0.2726, 'BM' : 0.3175, 'The300' : 0.307}
h  = {'TNG' : 0.6774, 'MGTM' : 0.7040, 'BM' : 0.6711, 'The300' : 0.678}

def _Hofz(sim, z):

    return h[sim]*100*np.sqrt(Om[sim]*(1 + z)**3 + (1 - Om[sim]))

def _hofz(sim, z):

    return h[sim]*np.sqrt(Om[sim]*(1 + z)**3 + (1 - Om[sim]))

def sigmaDM(M200c, z, parameter, sim):

    '''
    Convenience function that interpolates from a precomputed table
    to provide scaling parameters of sigmaDM with halo mass, M200c,
    at redshifts 0 <= z <= 1.

    The interpolation is done linearly over log(M200c), and log(a), where
    a = 1/(1 + z) is the scale factor.


    ---------
    Params
    ---------

    M200c:  (float, int) or (list, numpy array)
        The halo mass to extract parameters at. In units
        of Msun

    z: int, float
        The redshift that the parameters should be extracted at.
        When needed the function will interpolate between available
        data to estimate parameters at the exact input redshift.

    parameter: str
        The scaling relation parameter to extract. Can take
        values of 'mean' (in dex), 'slope', 'scatter' (in natural log).

    sim: str
        The specific simulation run the scaling relation should
        be extracted from. See Interpolator.avail_sims for the list of
        available simulations

    --------
    Output
    --------

    numpy array:

        Array of dimension (M200c.size, 7). The first column gives the
        mean value for the parameter while the columns after that provide
        the upper and lower bounds at 1/2/3 sigma. If a requested M200c or z value
        is outside the interpolation range, the corresponding row of entries
        in the output will contain np.NaN values.

    '''

    if sim not in avail_sims:

        raise ValueError("Requested sim, %s, is not available. Choose one of "%sim + str(avail_sims))

    if parameter not in avail_params:

        raise ValueError("Requested parameter, %s, is not available. Choosen one of "%parameter + str(avail_params))


    M200c = np.atleast_1d(M200c)[:, None]
    z     = np.ones_like(M200c) * z
    a     = 1/(1 + z)

    input  = np.concatenate([a, M200c], axis = 1)
    input  = np.log10(input)

    data       = np.loadtxt(cwd_path +'/../Data/KLLR_Params_sigmaDM_%s.txt'%sim)
    data[:, 0] = np.log10(1/(1 + data[:, 0]))

    index  = parameter_to_index[parameter]
    output = np.zeros([M200c.size, 7])

    for i in range(7):

        data_modified = data[:, index + i]

        interp       = interpolate.LinearNDInterpolator(data[:, [0, 1]], data_modified)
        output[:, i] = interp(input)

    return output


def sigmaSat(M200c, Mstarsat_Th, z, parameter, sim):

    '''
    Convenience function that interpolates from a precomputed table
    to provide scaling parameters of sigmaSat with halo mass, M200c,
    at redshifts 0 <= z <= 1.

    The interpolation is done linearly over log(M200c), and log(a), where
    a = 1/(1 + z) is the scale factor.


    ---------
    Params
    ---------


    M200c:  (float, int) or (list, numpy array)
        The halo mass to extract parameters at. In units
        of Msun

    MStarsat_Th:  (float, int)
        The galaxy stellar mass threshold of the sample that
        the parameters should be extracted from.

    z: int, float
        The redshift that the parameters should be extracted at.
        When needed the function will interpolate between available
        data to estimate parameters at the exact input redshift.

    parameter: str
        The scaling relation parameter to extract. Can take
        values of 'mean' (in dex), 'slope', 'scatter' (in natural log).
        The slope and scatter are assumed to be independent of mass.

    sim: str
        The specific simulation run the scaling relation should
        be extracted from. See Interpolator.avail_sims for the list of
        available simulations

    --------
    Output
    --------

    numpy array:

        Array of dimension (M200c.size, 7). The first column gives the
        mean value for the parameter while the columns after that provide
        the upper and lower bounds at 1/2/3 sigma. If a requested M200c or z value
        is outside the interpolation range, the corresponding row of entries
        in the output will contain np.NaN values.

    '''

    if sim not in avail_sims:

        raise ValueError("Requested sim, %s, is not available. Choose one of "%sim + str(avail_sims))

    if parameter not in avail_params:

        raise ValueError("Requested parameter, %s, is not available. Choosen one of "%parameter + str(avail_params))


    M200c = np.atleast_1d(M200c)[:, None]
    z     = np.ones_like(M200c) * z
    a     = 1/(1 + z)
    Mstarsat_Th = np.ones_like(M200c) * Mstarsat_Th

    input  = np.concatenate([a, Mstarsat_Th], axis = 1)
    input  = np.log10(input)

    data   = np.loadtxt(cwd_path +'/../Data/Params_sigmaSat_%s.txt'%sim)
    data[:, 0] = np.log10(1/(1 + data[:, 0]))

    index  = parameter_to_index[parameter]
    output = np.zeros([M200c.size, 7])

    if parameter in ['slope', 'scatter']:

        for i in range(7):

            data_modified = data[:, index + i]

            interp       = interpolate.LinearNDInterpolator(data[:, [0, 1]], data_modified)
            output[:, i] = interp(input)

    elif parameter in ['mean']:

        norm_mean = interpolate.LinearNDInterpolator(data[:, [0, 1]], data[:, 2])(input)
        norm_upp  = interpolate.LinearNDInterpolator(data[:, [0, 1]], data[:, 3])(input)
        norm_low  = interpolate.LinearNDInterpolator(data[:, [0, 1]], data[:, 4])(input)

        slope_mean = interpolate.LinearNDInterpolator(data[:, [0, 1]], data[:, 2 + 7])(input)
        slope_upp  = interpolate.LinearNDInterpolator(data[:, [0, 1]], data[:, 3 + 7])(input)
        slope_low  = interpolate.LinearNDInterpolator(data[:, [0, 1]], data[:, 4 + 7])(input)

        N = 10000

        norm  = norm_mean  + (norm_upp  - norm_low)/2  * np.random.randn(N, 1)
        slope = slope_mean + (slope_upp - slope_low)/2 * np.random.randn(N, 1)

        mean  = norm + slope*(np.log10(M200c).flatten() - 14 + np.log10(_hofz(sim, z)).flatten())

        output = np.quantile(mean, [0.5, 0.84, 0.16, 0.975, 0.025, 0.9985, 0.0015], axis = 0).T

    return output

def velocity_bias(M200c, Mstarsat_Th, z, sims):

    '''
    Convenience function that interpolates from a precomputed table
    to provide the bias as a function of halo mass, galaxy stellar mass threshold,
    and redshift.

    The interpolation is done linearly over log(M200c), log(Mstarsat), and log(a), where
    a = 1/(1 + z) is the scale factor.

    When more than one sim is provided in <sims>, the output is a theoretical prior
    estimated using the ensemble of simulations. The prior is represented as a gaussian
    with a mean and \sigma.


    ---------
    Params
    ---------


    M200c:  (float, int) or (list, numpy array)
        The halo mass to extract parameters at. In units
        of Msun

    MStarsat_Th:  (float, int)
        The galaxy stellar mass threshold of the sample that
        the parameters should be extracted from.

    z: (float, int)
        The redshift that the parameters should be extracted at.
        When needed the function will interpolate between available
        data to estimate parameters at the exact input redshift.

    sim: list
        The specific simulation run the scaling relation should
        be extracted from. See Interpolator.avail_sims for the list of
        available simulations

    --------
    Output
    --------

    numpy array:

        Array of dimension (M200c.size,) and contains the
        mean bias (linear, not log quantity) at requested
        values of M200c, Mstarsat, and z using all requested sims.

    numpy array:

        Array of dimension (M200c.size,) and contains the
        1\sigma absolute (not fractional) uncertainty on
        the bias at requested values of M200c, Mstarsat, and z
        using all requested sims.

    '''
    Store_bv, Store_dbv, Store_xline = {}, {}, {}
    for sim in sims:

        data = np.loadtxt(cwd_path + '/../Data/Velocity_bias_' + sim + '.txt')

        Store_bv[sim]     = data[:, 3]
        Store_dbv[sim]    = data[:, 4]
        Store_xline[sim]  = data[:, :3]

        Store_xline[sim][:, 0] = np.log10(1/(1 + Store_xline[sim][:, 0]))

    M200c = np.atleast_1d(M200c)
    z     = np.ones_like(M200c)*z
    Mstarsat_Th = np.ones_like(M200c)*Mstarsat_Th
    a = 1/(1 +z)

    input_params = np.concatenate([a[:, None], Mstarsat_Th[:, None], M200c[:, None]], axis = 1)
    input_params = np.log10(input_params)

    mean_bv = np.zeros([len(sims), input_params.shape[0]])
    std_bv  = np.zeros([len(sims), input_params.shape[0]])

    for i, sim in enumerate(sims):

        interp        = interpolate.LinearNDInterpolator(Store_xline[sim], Store_bv[sim])
        mean_bv[i, :] = interp(input_params)

        interp        = interpolate.LinearNDInterpolator(Store_xline[sim], Store_dbv[sim])
        std_bv[i, :]  = interp(input_params)

    output = np.concatenate([mean_bv[i, :] + np.random.randn(1000, 1)*std_bv[i, :] for i in range(len(sims))])

    # A bit of acrobatics to avoid runtime errors
    # for case where the whole vector is filled with np.NaNs
    # (due to extrapolating into space where no sims are available)
    mean, std = np.zeros([2, input_params.shape[0]]) + np.NaN
    Mask      = np.invert(np.all(np.isnan(mean_bv) | np.isnan(std_bv), axis = 0))

    mean[Mask], std[Mask] = np.nanmean(output[:, Mask], axis = 0), np.nanstd(output[:, Mask], axis = 0)

    return mean, std
