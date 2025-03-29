import scipy
import utils  # Import of utils.py
import numpy as np
import time as systime
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy.stats import bootstrap
from scipy.optimize import minimize
from matplotlib.gridspec import GridSpec
from astropy.cosmology import FlatLambdaCDM
from astropy.stats import sigma_clip, sigma_clipped_stats

# Utility Functions ===================================================================================================
def bootstrap_errs(x):
    try:
        return bootstrap([x], np.std).standard_error
    except:
        return 0
def get_bin_num(data: np.ndarray) -> int:
    # Using Freedman Diaconis Rule: bn = (max-min)/(2*IQR*n^(1/3))
    q1, q3 = np.percentile(data, 25), np.percentile(data, 75)
    iqr = q3 - q1
    return int((np.max(data) - np.min(data)) / (2 * iqr * (len(data)**(-1/3))))
def mass_step_calc(mu: np.ndarray, mu_err: np.ndarray, resid: np.ndarray, mass: np.ndarray, z: np.ndarray,
                   cut: float = 10.0) -> (dict, dict):
    """
    Calculates the mass step given arrays data
    :param mu: numpy.array(float);
    :param mu_err: numpy.array(float);
    :param resid: numpy.array(float);
    :param mass: numpy.array(float);
    :param z: numpy.array(float);
    :param cut: float = 10.0;
    :return dict; two dictionary of mass step and error & lower/upper weighted averages
    """
    if cut == 'median':
        cut = round(np.median(mass), 4)

    try:
        upper_resid = np.average(resid[mass > cut], weights=(1/(mu_err[mass > cut]**2)))
        lower_resid = np.average(resid[mass < cut], weights=(1/(mu_err[mass < cut]**2)))

        upper_resid_err = np.std(resid[mass > cut]) / np.sqrt(len(mu_err[mass > cut]))
        lower_resid_err = np.std(resid[mass < cut]) / np.sqrt(len(mu_err[mass < cut]))

        mass_step = np.abs(upper_resid - lower_resid)
        mass_step_err = np.sqrt((lower_resid_err**2) + (upper_resid_err**2))
    except ZeroDivisionError:
        return ({'value': 0.00, 'err': 0.00},
                {'lower_resid': {'value': 0.00, 'err': 0.00},
                 'upper_resid': {'value': 0.00, 'err': 0.00}})

    return ({'value': mass_step, 'err': mass_step_err},
            {'lower_resid': {'value': lower_resid, 'err': lower_resid_err},
             'upper_resid': {'value': upper_resid, 'err': upper_resid_err}})
def get_chi2(intercept, x, y, sigma, slope):
    b = intercept[0]  # Extract intercept
    model = slope * x + b
    return np.sum(((y - model) / sigma) ** 2)
def plot_binned_param(axis, path, param_name, bin_num, bin_bounds, p_label, p_color, label_bins: bool = True):
    # Open data
    tb = utils.default_open(path, True)
    param = np.array(tb[param_name])
    resid = np.array(tb['mu'] - utils.current_cosmo().distmod(tb['z_cmb']).value)

    ## Bin Dust & Residuals with STD of residuals
    param_bins = np.linspace(bin_bounds[0], bin_bounds[1], bin_num)
    binned_vals = scipy.stats.binned_statistic(param, resid, bins=param_bins, statistic=np.std).statistic
    binned_errs = scipy.stats.binned_statistic(param, resid, bins=param_bins, statistic=bootstrap_errs).statistic
    param_bins_adj = (param_bins[:-1] + param_bins[1:]) / 2.

    # Remove points with <1
    x_color, y_resid, y_resid_err = np.array([]), np.array([]), np.array([])
    for i in range(len(param_bins) - 1):
        if len(resid[(param > param_bins[i]) & (param < param_bins[i + 1])]) > 1:
            x_color = np.append(x_color, param_bins_adj[i])
            y_resid = np.append(y_resid, binned_vals[i])
            y_resid_err = np.append(y_resid_err, binned_errs[i])

    ## Plot binned dust w/ scatter
    axis.errorbar(x_color, y_resid, yerr=y_resid_err, fmt='o-', color=p_color,
                    label=p_label + f"{len(param)}")

    ## Label number of points in bin
    if label_bins:
        for i in range(len(param_bins) - 1):
            if np.isnan(binned_vals[i]): continue
            if len(resid[(param > param_bins[i]) & (param < param_bins[i + 1])]) <= 1: continue
            axis.text(param_bins_adj[i], binned_vals[i],
                        len(resid[(param > param_bins[i]) & (param < param_bins[i + 1])]),
                        ha='left', va='bottom', size='small')
    return

# WIP Functions =======================================================================================================
def resid_v_mass_dust(path_91bg: str = 'merged_params_cut.txt',
                      path_norm: str = 'aaronDo_salt2_params_cut.txt',
                      path_dust: str = 'global_dust_params.txt',
                      save_loc: str = '', label: bool = False):
    """
    Plots the Hubble Residual v. Mass
    """
    fig, axs = plt.subplots(2, 2, figsize=(18, 10), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)
    plt.style.use('tableau-colorblind10')
    all_resid, all_mass = [], []
    c_norm, c_norm_mass = 'C2', 'C3'
    c_91bg, c_91bg_mass = 'C8', 'C1'

    # Plot Normals
    # -----------------------------------------------------------------------------------------------------------------
    tb_norm = utils.default_open(path_norm, True)
    tb_norm = tb_norm[tb_norm['hostMass'] > 7.5]  # Removes the low mass Normals

    ## Calculate Hubble Residual
    tb_norm['resid_mu'] = tb_norm['mu'] - utils.current_cosmo().distmod(tb_norm['z_cmb']).value
    tb_norm['resid_mu'] -= np.average(tb_norm['resid_mu'][~np.isnan(tb_norm['resid_mu'])])  # Centering around average
    tb_norm['mu_err'] = np.sqrt(tb_norm['mu_err'] ** 2.0 + 0.1 ** 2.0)  # intrinsic dispersion added in quadrature
    tb_norm['resid_mu_err'] = np.copy(tb_norm['mu_err'])

    # Adding 0.1 mag in quadrature (taylor+11)
    tb_norm['hostMass_err'] = np.sqrt(
        tb_norm['hostMass_err'] ** 2.0 + 0.1 ** 2.0)  # intrinsic dispersion added in quadrature

    ## Scatter plot & histogram
    axs[0, 0].errorbar(x=tb_norm['hostMass'], y=tb_norm['resid_mu'], xerr=tb_norm['hostMass_err'],
                       yerr=tb_norm['resid_mu_err'],
                       marker='o', alpha=0.5, color=c_norm, fmt='o', ms=6, elinewidth=0.8)
    axs[0, 1].hist(tb_norm['resid_mu'], bins=20, orientation="horizontal", color=c_norm)
    # int((np.max(tb_norm['resid_mu']) - np.min(tb_norm['resid_mu'])) / 0.02)

    # Labels
    if label:
        for x, y, name in zip(tb_norm['hostMass'], tb_norm['resid_mu'], tb_norm['objname']):
            axs[0, 0].text(x, y, name, ha='left', va='top', size='xx-small')

    # Plot 10dex & Median Mass Lines
    tol = 1
    for cut, ls, cl in zip([10, 10.55], ['-', '--'], [c_norm_mass, c_norm_mass]):
        if cut == 'median': cut = np.median(tb_norm['hostMass'])
        lin_details = {'linestyle': ls, 'linewidth': 3, 'color': cl, 'zorder': 5}
        mass_step_dict, resid_dict = mass_step_calc(tb_norm['mu'], tb_norm['mu_err'], tb_norm['resid_mu'],
                                                    tb_norm['hostMass'], tb_norm['z_cmb'], cut=cut)
        if resid_dict['lower_resid']['value'] > resid_dict['upper_resid']['value']: mass_step_dict['value'] = \
        mass_step_dict['value'] * -1
        axs[0, 0].hlines(y=resid_dict['lower_resid']['value'], xmin=np.min(tb_norm['hostMass']) - tol, xmax=cut,
                         **lin_details)  # Left
        axs[0, 0].hlines(y=resid_dict['upper_resid']['value'], xmin=cut, xmax=np.max(tb_norm['hostMass']) + tol,
                         **lin_details)  # Right
        axs[0, 0].axvline(cut, alpha=0.75, **lin_details,
                          label="$\gamma (M_{split}=$" + f"{round(cut, 2)}) = " +
                                f"${round(mass_step_dict['value'], 3)} \pm {round(mass_step_dict['err'], 3)}$")

    # Plot 91bg-like
    # -----------------------------------------------------------------------------------------------------------------
    tb_91bg = utils.default_open(path_91bg, True)
    tb_dust = utils.default_open(path_dust, True)

    ## Calculate Hubble Residual
    tb_91bg['resid_mu'] = tb_91bg['mu'] - utils.current_cosmo().distmod(tb_91bg['z_cmb']).value
    tb_91bg['resid_mu'] -= np.average(
        tb_91bg['resid_mu'][~np.isnan(tb_91bg['resid_mu'])])  # Centering around average
    tb_91bg['mu_err'] = np.sqrt(tb_91bg['mu_err'] ** 2.0 + 0.1 ** 2.0)  # intrinsic dispersion added in quadrature
    tb_91bg['resid_mu_err'] = np.copy(tb_91bg['mu_err'])

    # Adding 0.1 mag in quadrature (taylor+11)
    tb_91bg['hostMass_err'] = np.sqrt(
        tb_91bg['hostMass_err'] ** 2.0 + 0.1 ** 2.0)  # intrinsic dispersion added in quadrature

    # Get associated dust values
    dust_91bg = np.array([])
    for name in tb_91bg['objname']:
        if name in list(tb_dust['objname']):
            dust_91bg = np.append(dust_91bg, tb_dust['av_50'][tb_dust['objname'] == name].value[0])
        else:
            tb_91bg.remove_row(list(tb_91bg['objname']).index(name))  # Removes if no dust value found

    ## Scatter plot & histogram
    median_dust = np.median(dust_91bg)
    axs[1, 0].errorbar(x=tb_91bg['hostMass'][dust_91bg < median_dust],
                       y=tb_91bg['resid_mu'][dust_91bg < median_dust],
                       xerr=tb_91bg['hostMass_err'][dust_91bg < median_dust],
                       yerr=tb_91bg['resid_mu_err'][dust_91bg < median_dust],
                       marker='s', alpha=1, color='C8', fmt='o', ms=6, elinewidth=0.8,
                       label='$A_{V=50} < Median$ ('+f'{round(median_dust,2)})')
    axs[1, 0].errorbar(x=tb_91bg['hostMass'][dust_91bg > median_dust],
                       y=tb_91bg['resid_mu'][dust_91bg > median_dust],
                       xerr=tb_91bg['hostMass_err'][dust_91bg > median_dust],
                       yerr=tb_91bg['resid_mu_err'][dust_91bg > median_dust],
                       marker='s', alpha=1, color='C9', fmt='o', ms=6, elinewidth=0.8,
                       label='$A_{V=50} > Median$ ('+f'{round(median_dust,2)})')
    axs[1, 1].hist(tb_91bg['resid_mu'], bins=20, orientation="horizontal", color=c_91bg)

    # Labels
    if label:
        for x, y, name in zip(tb_91bg['hostMass'], tb_91bg['resid_mu'], tb_91bg['objname']):
            axs[1, 0].text(x, y, name, ha='left', va='top', size='xx-small')

    # # Plot 10dex & Median Mass Lines
    tol = 1
    for cut, ls, cl in zip([10, 10.55], ['-', '--'], [c_91bg_mass, c_91bg_mass]):
        if cut == 'median': cut = np.median(tb_91bg['hostMass'])
        lin_details = {'linestyle': ls, 'linewidth': 3, 'color': cl, 'zorder': 5}
        mass_step_dict, resid_dict = mass_step_calc(tb_91bg['mu'], tb_91bg['mu_err'], tb_91bg['resid_mu'],
                                                    tb_91bg['hostMass'], tb_91bg['z_cmb'], cut=cut)
        if resid_dict['lower_resid']['value'] > resid_dict['upper_resid']['value']: mass_step_dict['value'] = \
        mass_step_dict['value'] * -1
        axs[1, 0].hlines(y=resid_dict['lower_resid']['value'], xmin=np.min(tb_91bg['hostMass']) - tol, xmax=cut,
                         **lin_details)  # Left
        axs[1, 0].hlines(y=resid_dict['upper_resid']['value'], xmin=cut, xmax=np.max(tb_91bg['hostMass']) + tol,
                         **lin_details)  # Right
        axs[1, 0].axvline(cut, alpha=0.75, **lin_details,
                          label="$\gamma (M_{split}=$" + f"{round(cut, 2)}) = " +
                                f"${round(mass_step_dict['value'], 3)} \pm {round(mass_step_dict['err'], 3)}$")

    # Brount, Scolnic 2021 Dust Prediction
    axs[1, 0].hlines(y=np.average(tb_91bg['resid_mu'][tb_91bg['hostMass'] < 10]) - 0.25,
                     xmin=10, xmax=np.max(tb_91bg['hostMass']) + tol,
                     label='Brout et al. 2021 (c = 0.2)', linestyle=':', linewidth=3, color='C0', zorder=5)

    # Formatting
    # -----------------------------------------------------------------------------------------------------------------
    ## Label number of SNe and Scatter
    axs[0, 0].text(0.04, 0.96,
                   "Normal SNe Ia\n" +
                   "$N_{SNe}$ = " + f"{len(tb_norm)}\n" +
                   "$\sigma$ = " + f"{round(np.std(tb_norm['resid_mu']), 3)} mag",
                   transform=axs[0, 0].transAxes, ha='left', va='top', fontsize=12)
    axs[1, 0].text(0.04, 0.96,
                   "1991bg-like SNe Ia\n" +
                   "$N_{SNe}$ = " + f"{len(tb_91bg)}\n" +
                   "$\sigma$ = " + f"{round(np.std(tb_91bg['resid_mu']), 3)} mag",
                   transform=axs[1, 0].transAxes, ha='left', va='top', fontsize=12)

    ## Adjust Axises
    tol = 0.1
    x_min = np.min(np.hstack([tb_norm['hostMass'], tb_91bg['hostMass']])) - tol
    x_max = np.max(np.hstack([tb_norm['hostMass'], tb_91bg['hostMass']])) + tol
    y_min = np.min(np.hstack([tb_norm['resid_mu'], tb_91bg['resid_mu']])) - tol
    y_max = np.max(np.hstack([tb_norm['resid_mu'], tb_91bg['resid_mu']])) + tol
    axs[0, 0].set(xlim=(x_min, x_max), ylim=(y_min, y_max))
    axs[1, 0].set(xlim=(x_min, x_max), ylim=(y_min, y_max))
    axs[0, 1].set(ylim=(y_min, y_max))
    axs[1, 1].set(ylim=(y_min, y_max))
    axs[0, 0].tick_params(labelbottom=False)
    axs[0, 1].tick_params(labelleft=False, labelbottom=False)
    axs[1, 1].tick_params(labelleft=False, labelbottom=False)
    # cb = plt.colorbar(sc, ax=[axs[1, 0]], location='bottom')

    ## Labels
    axs[0, 0].set_ylabel('Hubble Residual (mag)', size=16)
    axs[1, 0].set_ylabel('Hubble Residual (mag)', size=16)
    axs[1, 0].set_xlabel("Host Stellar Mass ($\log M_{*}[M_{\odot}]$)", size=16)
    axs[0, 0].legend(loc='lower left')
    axs[1, 0].legend(loc='lower left')

    # Saving Figure
    if len(save_loc) != 0:
        print('Saved figure to... ', save_loc)
        plt.savefig(save_loc)
    plt.show()
    return
def dust_mass_step():
    fig, axs = plt.subplots(2, 2, figsize=(18, 10), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)
    plt.style.use('tableau-colorblind10')
    all_resid, all_mass = [], []
    c_norm, c_norm_mass = 'C2', 'C3'
    c_91bg, c_91bg_mass = 'C8', 'C1'

    tb_91bg = utils.default_open('results/old/merged_params_cut.txt', True)
    tb_norm = utils.default_open('results/old/aaronDo_salt2_params_cut.txt', True)
    tb_dust = utils.default_open('results/old/global_dust_params.txt', True)

    # Fix dust tables
    tb_91bg.add_column(name='av_50', col=np.nan)
    for n in tb_dust['objname']:
        if n in tb_91bg['objname']:
            tb_91bg['av_50'][tb_91bg['objname'] == n] = tb_dust['av_50'][tb_dust['objname'] == n]
    tb_91bg = tb_91bg[~np.isnan(tb_91bg['av_50'])]

    for tb, axis in zip([tb_91bg, tb_91bg], axs):
        resid, resid_err = tb['mu'] - utils.current_cosmo().distmod(tb['z_cmb']).value, tb['mu_err']
        mass, mass_err = tb['hostMass'], tb['hostMass_err']

        # axis[0].errorbar(mass[tb['av_50']], resid, xerr=mass_err, yerr=resid_err, fmt='o')




    ## Adjust Axises
    # tol = 0.1
    # x_min = np.min(np.hstack([tb_norm['hostMass'], tb_91bg['hostMass']])) - tol
    # x_max = np.max(np.hstack([tb_norm['hostMass'], tb_91bg['hostMass']])) + tol
    # y_min = np.min(np.hstack([tb_norm['resid_mu'], tb_91bg['resid_mu']])) - tol
    # y_max = np.max(np.hstack([tb_norm['resid_mu'], tb_91bg['resid_mu']])) + tol
    # axs[0,0].set(xlim=(x_min, x_max), ylim=(y_min, y_max))
    # axs[1,0].set(xlim=(x_min, x_max), ylim=(y_min, y_max))
    # axs[0,1].set(ylim=(y_min, y_max))
    # axs[1,1].set(ylim=(y_min, y_max))
    axs[0,0].tick_params(labelbottom=False)
    axs[0,1].tick_params(labelleft=False, labelbottom=False)
    axs[1,1].tick_params(labelleft=False, labelbottom=False)

    ## Labels
    axs[0,0].set_ylabel('Hubble Residual (mag)', size=16)
    axs[1,0].set_ylabel('Hubble Residual (mag)', size=16)
    axs[1,0].set_xlabel("Host Stellar Mass ($\log M_{*}[M_{\odot}]$)", size=16)
    axs[0,0].legend(loc='lower left')
    axs[1,0].legend(loc='lower left')

    # Saving Figure
    # if len(save_loc) != 0:
    #     print('Saved figure to... ', save_loc)
    #     plt.savefig(save_loc)
    plt.show()
    return

# Finished Functions ==================================================================================================
def resid_v_mass(path_91bg: str = 'merged_params_cut.txt',
                 path_norm: str = 'aaronDo_salt2_params_cut.txt',
                 save_loc: str = '', label: bool = False):
    """
    Plots the Hubble Residual v. Mass
    """
    fig, axs = plt.subplots(2, 2, figsize=(18, 10), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)
    plt.style.use('tableau-colorblind10')
    all_resid, all_mass = [], []
    c_norm, c_norm_mass = 'C2', 'C3'
    c_91bg, c_91bg_mass = 'C8', 'C1'

    # Get median of 91bg-like
    tb_91bg = utils.default_open(path_91bg, True, delimiter=', ')
    median_91bg_mass = round(np.median(tb_91bg['hostMass']), 3)

    # Plot Normals
    # -----------------------------------------------------------------------------------------------------------------
    tb_norm = utils.default_open(path_norm, True)
    tb_norm = tb_norm[tb_norm['hostMass']>7.5] # Removes the low mass Normals

    ## Calculate Hubble Residual
    tb_norm['resid_mu'] = tb_norm['mu'] - utils.current_cosmo().distmod(tb_norm['z_cmb']).value
    tb_norm['resid_mu'] -= np.average(tb_norm['resid_mu'][~np.isnan(tb_norm['resid_mu'])]) # Centering around average
    tb_norm['mu_err'] = np.sqrt(tb_norm['mu_err'] ** 2.0 + 0.1 ** 2.0)  # intrinsic dispersion added in quadrature
    tb_norm['resid_mu_err'] = np.copy(tb_norm['mu_err'])

    # Adding 0.1 mag in quadrature (taylor+11)
    tb_norm['hostMass_err'] = np.sqrt(tb_norm['hostMass_err'] ** 2.0 + 0.1 ** 2.0)  # intrinsic dispersion added in quadrature

    ## Scatter plot & histogram
    axs[0,0].errorbar(x=tb_norm['hostMass'], y=tb_norm['resid_mu'], xerr=tb_norm['hostMass_err'], yerr=tb_norm['resid_mu_err'],
                      marker='o', alpha=0.5, color=c_norm, fmt='o', ms=6, elinewidth=0.8)
    axs[0,1].hist(tb_norm['resid_mu'], bins=20, orientation="horizontal", color=c_norm)
    # int((np.max(tb_norm['resid_mu']) - np.min(tb_norm['resid_mu'])) / 0.02)

    # Labels
    if label:
        for x, y, name in zip(tb_norm['hostMass'], tb_norm['resid_mu'], tb_norm['objname']):
            axs[0,0].text(x, y, name, ha='left', va='top', size='xx-small')

    # Plot 10dex & Median Mass Lines
    tol = 1
    for cut, ls, cl in zip([10, median_91bg_mass], ['-', '--'], [c_norm_mass, c_norm_mass]):
        if cut == 'median': cut = np.median(tb_norm['hostMass'])
        lin_details = {'linestyle': ls, 'linewidth': 3, 'color': cl, 'zorder': 5}
        mass_step_dict, resid_dict = mass_step_calc(tb_norm['mu'], tb_norm['mu_err'], tb_norm['resid_mu'],
                                                    tb_norm['hostMass'], tb_norm['z_cmb'], cut=cut)
        if resid_dict['lower_resid']['value'] > resid_dict['upper_resid']['value']: mass_step_dict['value'] = mass_step_dict['value']*-1
        axs[0,0].hlines(y=resid_dict['lower_resid']['value'], xmin=np.min(tb_norm['hostMass']) - tol, xmax=cut, **lin_details)  # Left
        axs[0,0].hlines(y=resid_dict['upper_resid']['value'], xmin=cut, xmax=np.max(tb_norm['hostMass']) + tol, **lin_details)  # Right
        axs[0,0].axvline(cut, alpha=0.75, **lin_details,
                         label="$\gamma (M_{split}=$" + f"{round(cut, 2)}) = " +
                               f"${round(mass_step_dict['value'], 3)} \pm {round(mass_step_dict['err'], 3)}$ mag")

    # Plot 91bg-like
    # -----------------------------------------------------------------------------------------------------------------
    tb_91bg = utils.default_open(path_91bg, True, delimiter=', ')

    ## Calculate Hubble Residual
    tb_91bg['resid_mu'] = tb_91bg['mu'] - utils.current_cosmo().distmod(tb_91bg['z_cmb']).value

    # Adding 0.1 mag in quadrature for intrinsic dispersion
    tb_91bg['mu_err'] = np.sqrt(tb_91bg['mu_err'] ** 2.0 + 0.1 ** 2.0)
    tb_91bg['resid_mu_err'] = np.copy(tb_91bg['mu_err'])

    # Adding 0.1 mag in quadrature (taylor+11)
    tb_91bg['hostMass_err'] = np.sqrt(tb_91bg['hostMass_err'] ** 2.0 + 0.1 ** 2.0)  # intrinsic dispersion added in quadrature

    # Subtracting off Average Hubble Residual
    tb_91bg['resid_mu'] -= np.average(tb_91bg['resid_mu'][~np.isnan(tb_91bg['resid_mu'])])

    ## Scatter plot & histogram
    axs[1, 0].errorbar(x=tb_91bg['hostMass'][tb_91bg['algo'] == 'SNPY'],
                       y=tb_91bg['resid_mu'][tb_91bg['algo'] == 'SNPY'],
                       xerr=tb_91bg['hostMass_err'][tb_91bg['algo'] == 'SNPY'],
                       yerr=tb_91bg['resid_mu_err'][tb_91bg['algo'] == 'SNPY'],
                       marker='s', alpha=1, color=c_91bg, fmt='o', ms=6, elinewidth=0.8, label='SNooPy')
    axs[1, 0].errorbar(x=tb_91bg['hostMass'][tb_91bg['algo'] == 'SALT'],
                       y=tb_91bg['resid_mu'][tb_91bg['algo'] == 'SALT'],
                       xerr=tb_91bg['hostMass_err'][tb_91bg['algo'] == 'SALT'],
                       yerr=tb_91bg['resid_mu_err'][tb_91bg['algo'] == 'SALT'],
                       marker='^', alpha=1, color=c_91bg, fmt='o', ms=6, elinewidth=0.8, label='SALT3')
    axs[1, 1].hist(tb_91bg['resid_mu'], bins=20, orientation="horizontal", color=c_91bg)

    # Labels
    if label:
        for x, y, name in zip(tb_91bg['hostMass'], tb_91bg['resid_mu'], tb_91bg['objname']):
            axs[1, 0].text(x, y, name, ha='left', va='top', size='xx-small')

    # # Plot 10dex & Median Mass Lines
    tol = 1
    for cut, ls, cl in zip([10, median_91bg_mass], ['-', '--'], [c_91bg_mass, c_91bg_mass]):
        if cut == 'median': cut = np.median(tb_91bg['hostMass'])
        lin_details = {'linestyle': ls, 'linewidth': 3, 'color': cl, 'zorder': 5}
        mass_step_dict, resid_dict = mass_step_calc(tb_91bg['mu'], tb_91bg['mu_err'], tb_91bg['resid_mu'],
                                                    tb_91bg['hostMass'], tb_91bg['z_cmb'], cut=cut)
        if resid_dict['lower_resid']['value'] > resid_dict['upper_resid']['value']: mass_step_dict['value'] = mass_step_dict['value']*-1
        axs[1, 0].hlines(y=resid_dict['lower_resid']['value'], xmin=np.min(tb_91bg['hostMass']) - tol, xmax=cut, **lin_details)  # Left
        axs[1, 0].hlines(y=resid_dict['upper_resid']['value'], xmin=cut, xmax=np.max(tb_91bg['hostMass']) + tol, **lin_details)  # Right
        axs[1, 0].axvline(cut, alpha=0.75, **lin_details,
                          label="$\gamma (M_{split}=$" + f"{round(cut, 2)}) = " +
                                f"${round(mass_step_dict['value'], 3)} \pm {round(mass_step_dict['err'], 3)}$ mag")

    # Brount, Scolnic 2021 Dust Prediction
    axs[1, 0].hlines(y=np.average(tb_91bg['resid_mu'][tb_91bg['hostMass'] < 10]) - 0.4297,
                     xmin=10, xmax=np.max(tb_91bg['hostMass']) + tol,
                     label='Brout et al. 2021 (c = 0.2)', linestyle=':', linewidth=3, color='C0', zorder=5)

    # # Plot 10dex & Median Mass Lines -- with fill
    # tol = 0.3
    # for cut, ls, cl in zip([10, 'median'], ['-', '--'], [['C1', 'C5'], ['C4', 'C0']]):
    #     if cut == 'median': cut = np.median(tb_91bg['hostMass'])
    #     lin_details = {'linestyle': ls, 'linewidth': 1.5, 'color': cl[0], 'zorder': 5}
    #     fill_details = {'color': cl[1], 'alpha': 0.15}
    #     mass_step_dict, resid_dict = mass_step_calc(tb_91bg['mu'], tb_91bg['mu_err'], tb_91bg['resid_mu'],
    #                                                 tb_91bg['hostMass'], tb_91bg['z_cmb'], cut=cut)
    #     axs[1, 0].vlines(x=cut, ymin=resid_dict['lower_resid']['value'], ymax=resid_dict['upper_resid']['value'],
    #                      **lin_details)  # Vertical
    #     axs[1, 0].hlines(y=resid_dict['lower_resid']['value'], xmin=np.min(tb_91bg['hostMass']) - tol, xmax=cut,
    #                      **lin_details)  # Left
    #     axs[1, 0].hlines(y=resid_dict['upper_resid']['value'], xmin=cut, xmax=np.max(tb_91bg['hostMass']) + tol,
    #                      **lin_details)  # Right
    #     axs[1, 0].axvline(cut, alpha=0.75, **lin_details,
    #                       label="$\gamma (M_{split}=$" + f"{round(cut, 2)}) = " +
    #                             f"{round(mass_step_dict['value'], 3)} $\pm$ {round(mass_step_dict['err'], 3)}")
    #     axs[1, 0].fill_between([cut, np.max(tb_91bg['hostMass']) + tol],
    #                         resid_dict['upper_resid']['value'] - resid_dict['upper_resid']['err'],
    #                         resid_dict['upper_resid']['value'] + resid_dict['upper_resid']['err'],
    #                         **fill_details)  # Right
    #     axs[1, 0].fill_between([np.min(tb_91bg['hostMass']) - tol, cut],
    #                         resid_dict['lower_resid']['value'] - resid_dict['lower_resid']['err'],
    #                         resid_dict['lower_resid']['value'] + resid_dict['lower_resid']['err'],
    #                         **fill_details) # Left

    # # Over plotting mass lines
    # tol = 1
    # for l in range(2):
    #     for t, cl in zip([tb_norm, tb_91bg], ['C3', 'C1']):
    #         for cut, ls in zip([10, 10.55], ['-', '--']):
    #             if cut == 'median': cut = np.median(t['hostMass'])
    #             lin_details = {'linestyle': ls, 'linewidth': 3, 'color': cl, 'zorder': 5}
    #             mass_step_dict, resid_dict = mass_step_calc(t['mu'], t['mu_err'], t['resid_mu'],
    #                                                         t['hostMass'], t['z_cmb'], cut=cut)
    #             if resid_dict['lower_resid']['value'] > resid_dict['upper_resid']['value']: mass_step_dict['value'] = \
    #             mass_step_dict['value'] * -1
    #             axs[l, 0].hlines(y=resid_dict['lower_resid']['value'], xmin=np.min(t['hostMass']) - tol, xmax=cut,
    #                              **lin_details)  # Left
    #             axs[l, 0].hlines(y=resid_dict['upper_resid']['value'], xmin=cut, xmax=np.max(t['hostMass']) + tol,
    #                              **lin_details)  # Right
    #             axs[l, 0].axvline(cut, alpha=0.75, **lin_details,
    #                               label="$\gamma (M_{split}=$" + f"{round(cut, 2)}) = " +
    #                                     f"${round(mass_step_dict['value'], 3)} \pm {round(mass_step_dict['err'], 3)}$")

    # Formatting
    # -----------------------------------------------------------------------------------------------------------------
    ## Label number of SNe and Scatter
    norm_mn, norm_md, norm_std = sigma_clipped_stats(tb_norm['resid_mu'])
    sn91bg_mn, sn91bg_md, sn91bg_std = sigma_clipped_stats(tb_91bg['resid_mu'])
    axs[0,0].text(0.04, 0.96,
                  "Normal SNe Ia\n"+
                  "$N_{SNe}$ = "+f"{len(tb_norm)}\n"+
                  "$\sigma$ = "+f"{round(norm_std,3)} mag\n",
                  transform=axs[0, 0].transAxes, ha='left', va='top', fontsize=12)
    axs[1,0].text(0.04, 0.96,
                  "1991bg-like SNe Ia\n" +
                  "$N_{SNe}$ = " + f"{len(tb_91bg)}\n" +
                  "$\sigma$ = " + f"{round(sn91bg_std, 3)} mag\n",
                  transform=axs[1, 0].transAxes, ha='left', va='top', fontsize=12)

    ## Adjust Axises
    tol = 0.1
    x_min = np.min(np.hstack([tb_norm['hostMass'], tb_91bg['hostMass']])) - tol
    x_max = np.max(np.hstack([tb_norm['hostMass'], tb_91bg['hostMass']])) + tol
    y_min = np.min(np.hstack([tb_norm['resid_mu'], tb_91bg['resid_mu']])) - tol
    y_max = np.max(np.hstack([tb_norm['resid_mu'], tb_91bg['resid_mu']])) + tol
    axs[0,0].set(xlim=(x_min, x_max), ylim=(y_min, y_max))
    axs[1,0].set(xlim=(x_min, x_max), ylim=(y_min, y_max))
    axs[0,1].set(ylim=(y_min, y_max))
    axs[1,1].set(ylim=(y_min, y_max))
    axs[0,0].tick_params(labelbottom=False)
    axs[0,1].tick_params(labelleft=False, labelbottom=False)
    axs[1,1].tick_params(labelleft=False, labelbottom=False)

    ## Labels
    axs[0,0].set_ylabel('Hubble Residual (mag)', size=16)
    axs[1,0].set_ylabel('Hubble Residual (mag)', size=16)
    axs[1,0].set_xlabel("Host Stellar Mass ($\log M_{*}[M_{\odot}]$)", size=16)
    axs[0,0].legend(loc='lower left')
    axs[1,0].legend(loc='lower left')

    # Saving Figure
    if len(save_loc) != 0:
        print('Saved figure to... ', save_loc)
        plt.savefig(save_loc)
    plt.show()
    return
def mu_v_z(path_91bg: str = 'merged_params_cut.txt',
           path_norm: str = 'aaronDo_salt2_params_cut.txt',
           save_loc: str = '', label: bool = False):
    """
    Plots the Hubble Residual v. Redshift
    """
    fig = plt.figure(layout="constrained", figsize=(18, 8), constrained_layout=True)
    plt.style.use('tableau-colorblind10')
    gs = GridSpec(6, 9, figure=fig)
    ax1 = fig.add_subplot(gs[:4, :8])
    ax2 = fig.add_subplot(gs[4:, :8])
    ax3 = fig.add_subplot(gs[:4, 8:])
    ax4 = fig.add_subplot(gs[4:, 8:])
    all_resid = []
    c_norm, c_model = 'C2', 'C3'
    c_91bg = 'C8'

    # Plot Normals
    # -----------------------------------------------------------------------------------------------------------------
    hdr, data = utils.default_open(path_norm)
    names = data[:, hdr.index('objname')]
    z = data[:, hdr.index('z_cmb')].astype(float)
    mass, mass_err = data[:, hdr.index('hostMass')].astype(float), data[:, hdr.index('hostMass_err')].astype(float)
    mu, mu_err = data[:, hdr.index('mu')].astype(float), np.sqrt(data[:, hdr.index('mu_err')].astype(float) ** 2.0 + 0.1 ** 2.0)  # intrinsic dispersion added in quadrature

    ## Calculate Hubble Residual
    resid_mu = mu - utils.current_cosmo().distmod(z).value
    resid_mu -= np.average(resid_mu[~np.isnan(resid_mu)])  # Centering around average
    resid_mu_err = np.copy(mu_err)

    ## Scatter plot
    fmt_scatter_dict = {'marker': 'o', 'alpha': 0.2, 'color': c_norm, 'fmt': 'o', 'ms': 6, 'elinewidth': 0.8}
    ax1.errorbar(z, mu, yerr=mu_err, label='$Normal\\text{ }Ia\\text{ }SNe$', **fmt_scatter_dict)
    ax2.errorbar(z, resid_mu, yerr=resid_mu_err, **fmt_scatter_dict)

    ## Labels
    if label:
        for i in range(len(z)):
            ax1.text(z[i], mu[i], names[i], ha='left', va='top', size='xx-small')

    ## Make histogram
    fmt_hist_dict = {'orientation': "horizontal", 'color': c_norm,}
    ax3.hist(mu, bins=int((np.max(mu) - np.min(mu)) / 0.2), **fmt_hist_dict)
    ax4.hist(resid_mu, bins=20, **fmt_hist_dict)

    all_resid.append(resid_mu) # Save mu data

    # Label number of SNe and Scatter
    ax1.text(0.98, 0.20,
             "Normal SNe Ia\n" +
             "$N_{SNe}$ = " + f"{len(mu)}\n" +
             "$\sigma$ = " + f"{round(np.std(resid_mu), 3)} mag",
             transform=ax1.transAxes, ha='right', va='bottom', fontsize=12)

    # # Plot 91bg-like
    # # -----------------------------------------------------------------------------------------------------------------
    hdr, data = utils.default_open(path_91bg)
    names = data[:, hdr.index('objname')]
    origins = data[:, hdr.index('origin')]
    algo = data[:, hdr.index('algo')]
    z = data[:, hdr.index('z_cmb')].astype(float)
    mass, mass_err = data[:, hdr.index('hostMass')].astype(float), data[:, hdr.index('hostMass_err')].astype(float)
    mu, mu_err = data[:, hdr.index('mu')].astype(float), data[:, hdr.index('mu_err')].astype(float)
    mu_err = np.sqrt(mu_err ** 2.0 + 0.1 ** 2.0)  # intrinsic dispersion added in quadrature

    # Calculate Hubble Residual
    resid_mu = mu - utils.current_cosmo().distmod(z).value
    resid_mu -= np.average(resid_mu[~np.isnan(resid_mu)])  # Centering around average
    resid_mu_err = np.copy(mu_err)

    # Make main plot
    ax1.errorbar(x=z[algo == 'SNPY'], y=mu[algo == 'SNPY'], yerr=mu_err[algo == 'SNPY'],
                 marker='s', alpha=1, color=c_91bg, fmt='o', ms=6, elinewidth=0.8,
                 label='$1991bg\\text{-}like\\text{ }Ia\\text{ }SNe_{SNooPy}$')
    ax2.errorbar(x=z[algo == 'SNPY'], y=resid_mu[algo == 'SNPY'], yerr=resid_mu_err[algo == 'SNPY'],
                 marker='s', alpha=1, color=c_91bg, fmt='o', ms=6, elinewidth=0.8)
    ax1.errorbar(x=z[algo == 'SALT'], y=mu[algo == 'SALT'], yerr=mu_err[algo == 'SALT'],
                 marker='^', alpha=1, color=c_91bg, fmt='o', ms=6, elinewidth=0.8,
                 label='$1991bg\\text{-}like\\text{ }Ia\\text{ }SNe_{SALT3}$')
    ax2.errorbar(x=z[algo == 'SALT'], y=resid_mu[algo == 'SALT'], yerr=resid_mu_err[algo == 'SALT'],
                 marker='^', alpha=1, color=c_91bg, fmt='o', ms=6, elinewidth=0.8)

    # Labels
    if label:
        for i in range(len(z)):
            ax1.text(z[i], mu[i], names[i], ha='left', va='top', size='xx-small')

    # Make histogram
    ax3.hist(mu, bins=int((np.max(mu) - np.min(mu)) / 0.2), orientation="horizontal", color=c_91bg)
    ax4.hist(resid_mu, bins=20, orientation="horizontal", color=c_91bg)

    all_resid.append(resid_mu) # Save mu data

    # Label number of SNe and Scatter
    ax1.text(0.98, 0.02,
             "1991bg-like SNe Ia\n" +
             "$N_{SNe}$ = " + f"{len(mu)}\n" +
             "$\sigma$ = " + f"{round(np.std(resid_mu), 3)} mag",
             transform=ax1.transAxes, ha='right', va='bottom', fontsize=12)

    # Plot fit line
    # -----------------------------------------------------------------------------------------------------------------
    model_z = np.arange(0.015, 0.115, 0.001)
    ax1.plot(model_z, utils.current_cosmo().distmod(model_z).value,
             label='Model [$H_0 = 70$, $\Omega_m = 0.3$]', zorder=10, c=c_model)
    ax2.axhline(y=0, zorder=10, color=c_model)


    # axs[1, 0].text(0.04, 0.96,
    #                "1991bg-like SNe Ia\n" +
    #                "$N_{SNe}$ = " + f"{len(tb_91bg)}\n" +
    #                "$\sigma$ = " + f"{round(np.std(tb_91bg['resid_mu']), 3)} mag",
    #                transform=axs[1, 0].transAxes, ha='left', va='top', fontsize=12)

    # Formatting
    ax1.set_ylabel('$\mu$', size=16)
    ax2.set_ylabel('Residuals', size=16)
    ax2.set_xlabel('Host Galaxy CMB Redshift', size=16)
    ax1.legend(loc='best')
    ax1.tick_params(axis='x', labelbottom=False)
    ax3.tick_params(axis='x', labelbottom=False)
    ax3.tick_params(axis='y', labelleft=False)
    ax4.tick_params(axis='y', labelleft=False)

    # Saving Figure
    if len(save_loc) != 0:
        print('Saved figure to... ', save_loc)
        plt.savefig(save_loc)
    plt.show()
    return
def alpha_beta(path_91bg: str = 'salt_params_cov_cut.txt',
               path_norm: str = 'aaronDo_salt2_params_cut.txt',
               save_loc: str = '', label: bool = False):
    fig, ax = plt.subplots(1, 2, figsize=(21, 7), constrained_layout=True)
    plt.style.use('tableau-colorblind10')
    c_norm, c_norm_line = 'C2', 'C3'
    c_91bg, c_91bg_line = 'C8', 'C1'

    # Plot alpha & beta values
    CONSTANTS = utils.get_constants()
    alpha_91bg = -1*float(CONSTANTS['salt_alpha_91bg'])
    beta_91bg = float(CONSTANTS['salt_beta_91bg'])
    alpha_norm = -1*float(CONSTANTS['salt_alpha'])
    beta_norm = float(CONSTANTS['salt_beta'])

    # Scatter Plot for 91bg-like
    fmt_dict_91bg = {'fmt': 'o', 'marker': 's', 'alpha': 1.0, 'label': '$M_{1991bg\\text{-}like}$', 'color': c_91bg}
    hdr, data = utils.default_open(path_91bg)
    if len(data) == 0:
        print(f"[!!!] Missing SALT data! Read as empty...")
        return
    x0_91bg, x0_err_91bg = data[:, hdr.index('amplitude')].astype(float), data[:, hdr.index('amplitude_err')].astype(float)
    x1_91bg, x1_err_91bg = data[:, hdr.index('stretch')].astype(float), data[:, hdr.index('stretch_err')].astype(float)
    c_91bg, c_err_91bg = data[:, hdr.index('color')].astype(float), data[:, hdr.index('color_err')].astype(float)
    z_91bg = data[:, hdr.index('z_cmb')].astype(float)
    mu_91bg = utils.current_cosmo().distmod(z_91bg).value
    mB_91bg, mB_err_91bg = ((-2.5 * np.log10(x0_91bg)) + 10.635), np.sqrt((2.5 * (x0_err_91bg / (x0_91bg * np.log(10)))) ** 2.0 + 0.1 ** 2.0)
    absmB_91bg, absmB_err_91bg = (mB_91bg - mu_91bg), np.copy(mB_err_91bg)
    ax[0].errorbar(x=x1_91bg, y=absmB_91bg, xerr=x1_err_91bg, yerr=absmB_err_91bg, **fmt_dict_91bg)
    ax[1].errorbar(x=c_91bg, y=absmB_91bg, xerr=c_err_91bg, yerr=absmB_err_91bg, **fmt_dict_91bg)

    # Scatter Plot for Normals
    fmt_dict_norm = {'fmt': 'o', 'marker': 'o', 'alpha': 0.2, 'label': '$M_{Normal\\text{ }SNIa}$', 'color': c_norm}
    hdr, data = utils.default_open(path_norm)
    x0_norm, x0_err_norm = data[:, hdr.index('amplitude')].astype(float), data[:, hdr.index('amplitude_err')].astype(float)
    x1_norm, x1_err_norm = data[:, hdr.index('stretch')].astype(float), data[:, hdr.index('stretch_err')].astype(float)
    c_norm, c_err_norm = data[:, hdr.index('color')].astype(float), data[:, hdr.index('color_err')].astype(float)
    z_norm = data[:, hdr.index('z_cmb')].astype(float)
    mu_norm = utils.current_cosmo().distmod(z_norm).value
    mB_norm, mB_err_norm = ((-2.5 * np.log10(x0_norm)) + 10.635), np.sqrt((2.5 * (x0_err_norm / (x0_norm * np.log(10)))) ** 2.0 + 0.1 ** 2.0)
    absmB_norm, absmB_err_norm = (mB_norm - mu_norm), np.copy(mB_err_norm)
    ax[0].errorbar(x=x1_norm, y=absmB_norm, xerr=x1_err_norm, yerr=absmB_err_norm, **fmt_dict_norm)
    ax[1].errorbar(x=c_norm, y=absmB_norm, xerr=c_err_norm, yerr=absmB_err_norm, **fmt_dict_norm)

    # 91bg-like Fit Lines
    ax[0].axline((0, minimize(get_chi2, 0.00, args=(x1_91bg, absmB_91bg, absmB_err_91bg, alpha_91bg)).x[0]),
                 slope=alpha_91bg, color=c_91bg_line, label="$\\alpha_{1991bg\\text{-}like}" + f"={round(-1 * alpha_91bg, 2)}$", zorder=10)
    ax[0].axline((0, minimize(get_chi2, 0.00, args=(x1_norm, absmB_norm, absmB_err_norm, alpha_91bg)).x[0]),
                 slope=alpha_91bg, color=c_91bg_line, linestyle='--', zorder=10)
    ax[1].axline((0, minimize(get_chi2, 0.00, args=(c_91bg, absmB_91bg, absmB_err_91bg, beta_91bg)).x[0]),
                 slope=beta_91bg, color=c_91bg_line, label="$\\beta_{1991bg\\text{-}like}" + f"={round(beta_91bg, 2)}$", zorder=10)
    ax[1].axline((0, minimize(get_chi2, 0.00, args=(c_norm, absmB_norm, absmB_err_norm, beta_91bg)).x[0]),
                 slope=beta_91bg, color=c_91bg_line, linestyle='--', zorder=10)

    # Normal Fit Lines
    ax[0].axline((0, minimize(get_chi2, 0.00, args=(x1_norm, absmB_norm, absmB_err_norm, alpha_norm)).x[0]),
                 slope=alpha_norm, color=c_norm_line, label="$\\alpha_{Normal\\text{ }SNIa}" + f"={round(-1*alpha_norm, 2)}$", zorder=10)
    ax[0].axline((0, minimize(get_chi2, 0.00, args=(x1_91bg, absmB_91bg, absmB_err_91bg, alpha_norm)).x[0]),
                 slope=alpha_norm, color=c_norm_line, linestyle='--', zorder=10)
    ax[1].axline((0, minimize(get_chi2, 0.00, args=(c_norm, absmB_norm, absmB_err_norm, beta_norm)).x[0]),
                 slope=beta_norm, color=c_norm_line, label="$\\beta_{Normal\\text{ }SNIa}"+f"={round(beta_norm, 2)}$", zorder=10)
    ax[1].axline((0, minimize(get_chi2, 0.00, args=(c_91bg, absmB_91bg, absmB_err_91bg, beta_norm)).x[0]),
                 slope=beta_norm, color=c_norm_line, linestyle='--', zorder=10)

    # Formatting
    ax[0].set_xlabel('$x_1$', size=16)
    ax[0].set_ylabel('$m_{B} - \mu$', size=16)
    ax[1].set_xlabel('$c$', size=16)
    ax[0].invert_yaxis(); ax[1].invert_yaxis()
    ax[0].legend(); ax[1].legend()
    plt.subplots_adjust(wspace=0)
    plt.tick_params(labelleft=False)

    if len(save_loc) > 0:
        print(f"Saved figure to...  {save_loc}")
        plt.savefig(save_loc, dpi=300)
    plt.show()
    return
def param_hist(snpy_91bg_path: str, salt_91bg_path: str, snpy_norm_path: str, salt_norm_path: str,
               line_type: str = 'median', save_loc: str = ''):
    # Set colors
    c_norm, c_norm_line = 'C2', 'C6'
    c_91bg, c_91bg_line = 'C8', 'C1'

    # Open data
    tb_snpy_91bg = utils.default_open(snpy_91bg_path, True)
    tb_snpy_norm = utils.default_open(snpy_norm_path, True)
    tb_salt_91bg = utils.default_open(salt_91bg_path, True)
    tb_salt_norm = utils.default_open(salt_norm_path, True)

    # Print data ranges
    if False:
        # print('=====')
        print(f"s_BV: "
              f"{round(min(tb_snpy_91bg['stretch']), 3)} < "
              f"{round(np.median(tb_snpy_91bg['stretch']), 3)} < "
              f"{round(max(tb_snpy_91bg['stretch']), 3)}")
        print(f"E(B-V): "
              f"{round(min(tb_snpy_91bg['color']), 3)} < "
              f"{round(np.median(tb_snpy_91bg['color']), 3)} < "
              f"{round(max(tb_snpy_91bg['color']), 3)}")
        print(f"x_1: "
              f"{round(min(tb_salt_91bg['stretch']), 3)} < "
              f"{round(np.median(tb_salt_91bg['stretch']), 3)} < "
              f"{round(max(tb_salt_91bg['stretch']), 3)}")
        print(f"c: "
              f"{round(min(tb_salt_91bg['color']), 3)} < "
              f"{round(np.median(tb_salt_91bg['color']), 3)} < "
              f"{round(max(tb_salt_91bg['color']), 3)}")

    fig, axs = plt.subplots(2, 2, figsize=(20, 8), constrained_layout=True)
    plt.style.use('tableau-colorblind10')

    # Plot data
    bin_num = 20
    axs[0, 0].hist(tb_snpy_norm['stretch'], label="$s_{BV, Normal\\text{ }Ia\\text{ }SNe}$", color=c_norm,
                   bins=bin_num)
    axs[0, 0].hist(tb_snpy_91bg['stretch'], label="$s_{BV, 1991bg\\text{-}like\\text{ }Ia\\text{ }SNe}$", color=c_91bg, alpha=0.75,
                   bins=bin_num)

    axs[0, 1].hist(tb_snpy_norm['color'], label="$E(B-V)_{Normal\\text{ }Ia\\text{ }SNe}$", color=c_norm,
                   bins=bin_num)
    axs[0, 1].hist(tb_snpy_91bg['color'], label="$E(B-V)_{1991bg\\text{-}like\\text{ }Ia\\text{ }SNe}$", color=c_91bg, alpha=0.75,
                   bins=bin_num)

    axs[1, 0].hist(tb_salt_norm['stretch'], label="$x_{1, Normal\\text{ }Ia\\text{ }SNe}$", color=c_norm,
                   bins=bin_num)
    axs[1, 0].hist(tb_salt_91bg['stretch'], label="$x_{1, 1991bg\\text{-}like\\text{ }Ia\\text{ }SNe}$", color=c_91bg, alpha=0.75,
                   bins=bin_num)

    axs[1, 1].hist(tb_salt_norm['color'], label="$c_{Normal\\text{ }Ia\\text{ }SNe}$", color=c_norm,
                   bins=bin_num)
    axs[1, 1].hist(tb_salt_91bg['color'], label="$c_{1991bg\\text{-}like\\text{ }Ia\\text{ }SNe}$", color=c_91bg, alpha=0.75,
                   bins=bin_num)

    # Plot median/average lines
    if line_type == 'median':
        line_type = line_type[0].upper() + line_type[1:]

        axs[0, 0].axvline(np.median(tb_snpy_norm['stretch']), color=c_norm_line, linestyle='--', linewidth=3,
                          label=f"{line_type}"+
                                "$_{Normal\\text{ }Ia\\text{ }SNe}$"+
                                f" = ${round(np.median(tb_snpy_norm['stretch']), 3)}$")
        axs[0, 0].axvline(np.median(tb_snpy_91bg['stretch']), color=c_91bg_line, linestyle=':', linewidth=3,
                          label=f"{line_type}" +
                                "$_{1991bg\\text{-}like\\text{ }Ia\\text{ }SNe}$" +
                                f" = ${round(np.median(tb_snpy_91bg['stretch']), 3)}$")

        axs[0, 1].axvline(np.median(tb_snpy_norm['color']), color=c_norm_line, linestyle='--', linewidth=3,
                          label=f"{line_type}"+
                                "$_{Normal\\text{ }Ia\\text{ }SNe}$"+
                                f" = ${round(np.median(tb_snpy_norm['color']), 3)}$")
        axs[0, 1].axvline(np.median(tb_snpy_91bg['color']), color=c_91bg_line, linestyle=':', linewidth=3,
                          label=f"{line_type}" +
                                "$_{1991bg\\text{-}like\\text{ }Ia\\text{ }SNe}$" +
                                f" = ${round(np.median(tb_snpy_91bg['color']), 3)}$")

        axs[1, 0].axvline(np.median(tb_salt_norm['stretch']), color=c_norm_line, linestyle='--', linewidth=3,
                          label=f"{line_type}"+
                                "$_{Normal\\text{ }Ia\\text{ }SNe}$"+
                                f" = ${round(np.median(tb_salt_norm['stretch']), 3)}$")
        axs[1, 0].axvline(np.median(tb_salt_91bg['stretch']), color=c_91bg_line, linestyle=':', linewidth=3,
                          label=f"{line_type}" +
                                "$_{1991bg\\text{-}like\\text{ }Ia\\text{ }SNe}$" +
                                f" = ${round(np.median(tb_salt_91bg['stretch']), 3)}$")

        axs[1, 1].axvline(np.median(tb_salt_norm['color']), color=c_norm_line, linestyle='--', linewidth=3,
                          label=f"{line_type}"+
                                "$_{Normal\\text{ }Ia\\text{ }SNe}$"+
                                f" = ${round(np.median(tb_salt_norm['color']), 3)}$")
        axs[1, 1].axvline(np.median(tb_salt_91bg['color']), color=c_91bg_line, linestyle=':', linewidth=3,
                          label=f"{line_type}" +
                                "$_{1991bg\\text{-}like\\text{ }Ia\\text{ }SNe}$" +
                                f" = ${round(np.median(tb_salt_91bg['color']), 3)}$")
    elif line_type == 'average':
        line_type = line_type[0].upper() + line_type[1:]

        axs[0, 0].axvline(np.average(tb_snpy_norm['stretch']), color=c_norm_line, linestyle='--', linewidth=3,
                          label=f"{line_type}"+
                                "$_{Normal\\text{ }Ia\\text{ }SNe}$"+
                                f" = ${round(np.average(tb_snpy_norm['stretch']), 3)}$"+
                                f" $\pm {round(np.average(tb_snpy_norm['st_err']), 3)}$")
        axs[0, 0].axvline(np.average(tb_snpy_91bg['stretch']), color=c_91bg_line, linestyle=':', linewidth=3,
                          label=f"{line_type}" +
                                "$_{1991bg\\text{-}like\\text{ }Ia\\text{ }SNe}$" +
                                f" = ${round(np.average(tb_snpy_91bg['stretch']), 3)}$" +
                                f" $\pm {round(np.average(tb_snpy_91bg['st_err']), 3)}$")

        axs[0, 1].axvline(np.average(tb_snpy_norm['color']), color=c_norm_line, linestyle='--', linewidth=3,
                          label=f"{line_type}"+
                                "$_{Normal\\text{ }Ia\\text{ }SNe}$"+
                                f" = ${round(np.average(tb_snpy_norm['color']), 3)}$"+
                                f" $\pm {round(np.average(tb_snpy_norm['EBVhost_err']), 3)}$")
        axs[0, 1].axvline(np.average(tb_snpy_91bg['color']), color=c_91bg_line, linestyle=':', linewidth=3,
                          label=f"{line_type}" +
                                "$_{1991bg\\text{-}like\\text{ }Ia\\text{ }SNe}$" +
                                f" = ${round(np.average(tb_snpy_91bg['color']), 3)}$" +
                                f" $\pm {round(np.average(tb_snpy_91bg['EBVhost_err']), 3)}$")

        axs[1, 0].axvline(np.average(tb_salt_norm['stretch']), color=c_norm_line, linestyle='--', linewidth=3,
                          label=f"{line_type}"+
                                "$_{Normal\\text{ }Ia\\text{ }SNe}$"+
                                f" = ${round(np.average(tb_salt_norm['stretch']), 3)}$"+
                                f" $\pm {round(np.average(tb_salt_norm['x1_err']), 3)}$")
        axs[1, 0].axvline(np.average(tb_salt_91bg['stretch']), color=c_91bg_line, linestyle=':', linewidth=3,
                          label=f"{line_type}" +
                                "$_{1991bg\\text{-}like\\text{ }Ia\\text{ }SNe}$" +
                                f" = ${round(np.average(tb_salt_91bg['stretch']), 3)}$" +
                                f" $\pm {round(np.average(tb_salt_91bg['x1_err']), 3)}$")

        axs[1, 1].axvline(np.average(tb_salt_norm['color']), color=c_norm_line, linestyle='--', linewidth=3,
                          label=f"{line_type}"+
                                "$_{Normal\\text{ }Ia\\text{ }SNe}$"+
                                f" = ${round(np.average(tb_salt_norm['color']), 3)}$"+
                                f" $\pm {round(np.average(tb_salt_norm['c_err']), 3)}$")
        axs[1, 1].axvline(np.average(tb_salt_91bg['color']), color=c_91bg_line, linestyle=':', linewidth=3,
                          label=f"{line_type}" +
                                "$_{1991bg\\text{-}like\\text{ }Ia\\text{ }SNe}$" +
                                f" = ${round(np.average(tb_salt_91bg['color']), 3)}$" +
                                f" $\pm {round(np.average(tb_salt_91bg['c_err']), 3)}$")

    # Enable legends
    axs[0, 0].legend(loc='upper right')
    axs[0, 1].legend(loc='upper right')
    axs[1, 0].legend(loc='upper right')
    axs[1, 1].legend(loc='upper right')

    # Set labels
    axs[0, 0].set_ylabel('SNooPy\n$N_{SNe}$', size=16)
    axs[1, 0].set_ylabel('SALT3\n$N_{SNe}$', size=16)
    axs[1, 0].set_xlabel('Stretch', size=16)
    axs[1, 1].set_xlabel('Color', size=16)

    # Adjust formatting
    axs[0, 1].tick_params(labelleft=False, labelright=True)
    axs[1, 1].tick_params(labelleft=False, labelright=True)

    # Adjust bounds
    tol = 0.3
    axs[0, 0].set_xlim(-1 * np.max(np.hstack([tb_snpy_norm['stretch'], tb_snpy_91bg['stretch']])) + 1 - tol,
                       np.max(np.hstack([tb_snpy_norm['stretch'], tb_snpy_91bg['stretch']])) + 1 + tol)
    axs[1, 0].set_xlim(-1 * np.max(np.hstack([tb_salt_norm['stretch'], tb_salt_91bg['stretch']])) - tol,
                       np.max(np.hstack([tb_salt_norm['stretch'], tb_salt_91bg['stretch']])) + tol)
    axs[0, 1].set_xlim(-1*np.max(np.hstack([tb_snpy_norm['color'], tb_snpy_91bg['color']])) - tol,
                       np.max(np.hstack([tb_snpy_norm['color'], tb_snpy_91bg['color']])) + tol)
    axs[1, 1].set_xlim(-1*np.max(np.hstack([tb_salt_norm['color'], tb_salt_91bg['color']])) - tol,
                       np.max(np.hstack([tb_salt_norm['color'], tb_salt_91bg['color']])) + tol)

    if len(save_loc) > 0:
        print(f"Saved figure to...  {save_loc}")
        plt.savefig(save_loc, dpi=300)
    plt.show()
    return
def dust_hist(path_91bg: str = 'salt_params_cov_cut.txt',
              path_red_norm: str = 'redNormSNe_salt.txt',
              path_dust: str = 'global_dust_params.txt',
              save_loc: str = '', label: bool = False):
    fig, ax = plt.subplots(1, 1, figsize=(12, 6), constrained_layout=True)
    plt.style.use('tableau-colorblind10')

    # Set colors
    c_norm, c_norm_line = 'C2', 'C6'
    c_91bg, c_91bg_line = 'C8', 'C1'

    # Open data
    tb_91bg = utils.default_open(path_91bg, True)
    tb_red_norm = utils.default_open(path_red_norm, True)
    tb_dust = utils.default_open(path_dust, True)
    tb_combined = Table(
        names=('objname', 'source', 'av', 'av_upper', 'av_lower'),
        dtype=(str, str, float, float, float))
    for i, n in enumerate(tb_dust['objname']):
        # Get 91bg-like Data
        if len(tb_91bg[tb_91bg['objname'] == n]) > 0:
            source = '91bg'
        # Get Normal Data
        elif len(tb_red_norm[tb_red_norm['objname'] == n]) > 0:
            source = 'norm'
        else:
            continue  # Skip values with no dust value
        # Get Dust E(B-V)
        av = tb_dust[tb_dust['objname'] == n]['av_50'].value[0]
        av_upper = (tb_dust[tb_dust['objname'] == n]['av_84'].value[0] -
                    tb_dust[tb_dust['objname'] == n]['av_50'].value[0])
        av_lower = (tb_dust[tb_dust['objname'] == n]['av_50'].value[0] -
                    tb_dust[tb_dust['objname'] == n]['av_16'].value[0])

        # Add to new table
        tb_combined.add_row((n, source, av, av_upper, av_lower))

    # Plot histogram
    ax.hist(tb_combined['av'][tb_combined['source'] == 'norm'], color=c_norm, bins=20,
            label='$A_{V=50}$ ($c > 0.15$) Normal Ia SNe')
    ax.hist(tb_combined['av'][tb_combined['source'] == '91bg'], color=c_91bg, bins=20, alpha=0.75,
            label='$A_{V=50}$ ($c > 0.15$) 1991bg-like Ia SNe')

    # Median lines
    for s, cl, lb in zip(['norm', '91bg'],
                         [c_norm_line, c_91bg_line],
                         ["$Median_{Normal\\text{ }Ia\\text{ }SNe}$", "$Median_{1991bg\\text{-}like\\text{ }Ia\\text{ }SNe}$"]):
        n_median = np.median(tb_combined['av'][tb_combined['source'] == s])
        ax.axvline(n_median, color=cl, linestyle='--', linewidth=3,
                   label=lb + f" = {round(n_median, 3)}")

    # Enable legend
    ax.legend(loc='upper right')

    # Add labels
    ax.set_xlabel('$A_{V=50}$', size=16)
    ax.set_ylabel('$N_{SNe}$', size=16)

    if len(save_loc) > 0:
        print(f"Saved figure to...  {save_loc}")
        plt.savefig(save_loc, dpi=300)
    plt.show()
    return
def abs_mag_v_color(path_91bg: str = 'salt_params_cov_cut.txt',
                    path_norm: str = 'redNormSNe_salt.txt',
                    path_dust: str = 'global_dust_params.txt',
                    save_loc: str = '', label: bool = False):
    fig, ax = plt.subplots(1, 1, figsize=(14, 6), constrained_layout=True)
    plt.style.use('tableau-colorblind10')
    c_norm, c_norm_line = 'C2', 'C3'
    c_91bg, c_91bg_line = 'C8', 'C1'

    # Open data
    tb_91bg = utils.default_open(path_91bg, True)
    tb_norm = utils.default_open(path_norm, True)
    tb_dust = utils.default_open(path_dust, True)
    tb_combined = Table(
        names=('objname', 'source', 'av', 'av_upper', 'av_lower', 'mu', 'mu_err', 'absmB', 'absmB_err', 'c', 'c_err'),
        dtype=(str, str, float, float, float, float, float, float, float, float, float))
    for i, n in enumerate(tb_dust['objname']):
        # Get 91bg-like Data
        if len(tb_91bg[tb_91bg['objname'] == n]) > 0:
            source = '91bg'
            mu = tb_91bg[tb_91bg['objname'] == n]['mu'].value[0]
            mu_err = tb_91bg[tb_91bg['objname'] == n]['mu_err'].value[0]
            amplitude = tb_91bg[tb_91bg['objname'] == n]['amplitude'].value[0]
            amplitude_err = tb_91bg[tb_91bg['objname'] == n]['amplitude_err'].value[0]
            c = tb_91bg[tb_91bg['objname'] == n]['color'].value[0]
            c_err = tb_91bg[tb_91bg['objname'] == n]['color_err'].value[0]
        # Get Normal Data
        elif len(tb_norm[tb_norm['objname'] == n]) > 0:
            source = 'norm'
            mu = tb_norm[tb_norm['objname'] == n]['mu'].value[0]
            amplitude = tb_norm[tb_norm['objname'] == n]['amplitude'].value[0]
            amplitude_err = tb_norm[tb_norm['objname'] == n]['amplitude_err'].value[0]
            c = tb_norm[tb_norm['objname'] == n]['color'].value[0]
            c_err = tb_norm[tb_norm['objname'] == n]['color_err'].value[0]
        else:
            continue  # Skip values that don't have dust values

        # Get Dust E(B-V)
        av = tb_dust[tb_dust['objname'] == n]['av_50'].value[0]
        av_upper = (tb_dust[tb_dust['objname'] == n]['av_84'].value[0] -
                    tb_dust[tb_dust['objname'] == n]['av_50'].value[0])
        av_lower = (tb_dust[tb_dust['objname'] == n]['av_50'].value[0] -
                    tb_dust[tb_dust['objname'] == n]['av_16'].value[0])

        # Calculate Absolute Mag
        mB = ((-2.5 * np.log10(amplitude)) + 10.635)
        mB_err = np.sqrt((2.5 * (amplitude_err / (amplitude * np.log(10)))) ** 2.0 + 0.1 ** 2.0)
        absmB = mB - mu
        absmB_err = np.copy(mB_err)

        # Add to new table
        tb_combined.add_row((n, source, av, av_upper, av_lower, mu, mu_err, absmB, absmB_err, c, c_err))

    # Plot data
    fmt_dict_91bg = {'fmt': 'o', 'marker': 's', 'alpha': 1.0, 'color': c_91bg, 'label': '$M_{1991bg\\text{-}like}$'}
    fmt_dict_norm = {'fmt': 'o', 'marker': 'o', 'alpha': 1.0, 'color': c_norm, 'label': '$M_{Normal\\text{ }Ia\\text{ }SNe}$'}
    for s, fmt_dict, ln_cl in zip(['norm', '91bg'], [fmt_dict_norm, fmt_dict_91bg], [c_norm_line, c_91bg_line]):
        # Fix av_err
        av_err = []
        low = np.array(tb_combined['av_lower'][tb_combined['source'] == s])
        up = np.array(tb_combined['av_upper'][tb_combined['source'] == s])
        for i in range(len(low)):
            av_err.append(np.array([low[i], up[i]]))
        av_err = np.array(av_err).T

        # ax[0].errorbar(x=tb_combined['av'][tb_combined['source'] == s],
        #                y=tb_combined['absmB'][tb_combined['source'] == s],
        #                xerr=av_err,
        #                yerr=tb_combined['absmB_err'][tb_combined['source'] == s],
        #                **fmt_dict)
        ax.errorbar(x=tb_combined['c'][tb_combined['source'] == s],
                       y=tb_combined['absmB'][tb_combined['source'] == s],
                       xerr=tb_combined['c_err'][tb_combined['source'] == s],
                       yerr=tb_combined['absmB_err'][tb_combined['source'] == s],
                       **fmt_dict)

        # Fit Lines
        a, b = np.polyfit(tb_combined['c'][tb_combined['source'] == s],
                          tb_combined['absmB'][tb_combined['source'] == s], 1)
        ax.axline((0, b), slope=a, color=ln_cl, zorder=5) # label=f'{round(a, 2)}',

    # Formatting
    ax.set_ylabel('Absolute Magnitude, $M = m_{B} - \mu$ (mag)', size=16)
    # ax[0].set_xlabel('$A_{V=50}$', size=16)
    # ax[0].invert_yaxis()
    # ax[0].legend(loc='upper right')
    ax.legend(loc='upper right')
    ax.set_xlabel('SALT3 color, $c$', size=16)
    ax.invert_yaxis()
    ax.tick_params(labelleft=False)

    if len(save_loc) > 0:
        print(f"Saved figure to...  {save_loc}")
        plt.savefig(save_loc, dpi=300)
    plt.show()
    return
def color_v_scatter(path_snpy_91bg: str = 'results/combiend__snpy_params_cut.txt',
                    path_salt_91bg: str = 'results/combiend__salt_params_cut.txt',
                    path_snpy_norm: str = 'results/output/norm_snpy_params_cut.txt',
                    path_salt_norm: str = 'results/aaronDo_salt2_params_cut.txt',
                    bin_nums: list = [[40, 40], [40, 50]], bin_bounds: list = [[-1, 1], [-1, 1]], label: bool = False,
                    save_loc: str = ''):
    """
    :param path_snpy_91bg: File path to SNooPy 1991bg-like SNe data.
    :param path_salt_91bg: File path to SALT3 1991bg-like SNe data.
    :param path_snpy_norm: File path to SNooPy Normals SNe data.
    :param path_salt_norm: File path to SALT3 Normal SNe data.
    :param bin_nums: Number of bins of order... [[SNooPy-91bg, SNooPy-Normal], [SALT3-91bg, SALT3-Normal]]
    :param bin_bounds: Upper/Lower bounds of bin array... [[SNooPy-91bg, SNooPy-Normal], [SALT3-91bg, SALT3-Normal]].
    :param label: Whether or not to label amount in each bin on the plot.
    :param save_loc: File path to save PNG of plot.
    :return: None
    """
    fig, axs = plt.subplots(2, 1, figsize=(16, 10), constrained_layout=True)
    plt.style.use('tableau-colorblind10')
    c_91bg, c_norm = 'C8', 'C2'

    # Plot top SNooPy panel ===========================================================================================
    for path, bin_num, p_label, p_color in zip([path_snpy_91bg, path_snpy_norm],
                                               [bin_nums[0][0], bin_nums[0][1]],
                                               ["1991bg-like SNe Ia, $N_{SNe}$ = ", "Normal SNe Ia, $N_{SNe}$ = "],
                                               [c_91bg, c_norm]):
        ## Open data
        tb = utils.default_open(path, True)
        colors = np.array(tb['color'])
        resid = np.array(tb['mu'] - utils.current_cosmo().distmod(tb['z_cmb']).value)

        ## Bin Dust & Residuals with STD of residuals
        color_bins = np.linspace(bin_bounds[0][0], bin_bounds[0][1], bin_num)
        binned_vals = scipy.stats.binned_statistic(colors, resid, bins=color_bins, statistic=np.std).statistic
        binned_errs = scipy.stats.binned_statistic(colors, resid, bins=color_bins, statistic=bootstrap_errs).statistic
        color_bins_adj = (color_bins[:-1] + color_bins[1:]) / 2.

        # Remove points with <1
        x_color, y_resid, y_resid_err = np.array([]), np.array([]), np.array([])
        for i in range(len(color_bins) - 1):
            if len(resid[(colors > color_bins[i]) & (colors < color_bins[i + 1])]) > 1:
                x_color = np.append(x_color, color_bins_adj[i])
                y_resid = np.append(y_resid, binned_vals[i])
                y_resid_err = np.append(y_resid_err, binned_errs[i])

        ## Plot binned dust w/ scatter
        axs[0].errorbar(x_color, y_resid, yerr=y_resid_err, fmt='o-', color=p_color,
                        label=p_label + f"{len(colors)}")

        ## Label number of points in bin
        if label:
            for i in range(len(color_bins) - 1):
                if np.isnan(binned_vals[i]): continue
                if len(resid[(colors > color_bins[i]) & (colors < color_bins[i + 1])]) <= 1: continue
                axs[0].text(color_bins_adj[i], binned_vals[i],
                            len(resid[(colors > color_bins[i]) & (colors < color_bins[i + 1])]),
                            ha='left', va='bottom', size='small')

    # Plot bottom SALT panel ==========================================================================================
    for path, bin_num, p_label, p_color in zip([path_salt_91bg, path_salt_norm],
                                               [bin_nums[1][0], bin_nums[1][1]],
                                               ["1991bg-like SNe Ia, $N_{SNe}$ = ", "Normal SNe Ia, $N_{SNe}$ = "],
                                               [c_91bg, c_norm]):
        ## Open data
        tb = utils.default_open(path, True)
        colors = np.array(tb['color'])
        resid = np.array(tb['mu'] - utils.current_cosmo().distmod(tb['z_cmb']).value)

        ## Bin Dust & Residuals with STD of residuals
        color_bins = np.linspace(bin_bounds[1][0], bin_bounds[1][1], bin_num)
        binned_vals = scipy.stats.binned_statistic(colors, resid, bins=color_bins, statistic=np.std).statistic
        binned_errs = scipy.stats.binned_statistic(colors, resid, bins=color_bins, statistic=bootstrap_errs).statistic
        color_bins_adj = (color_bins[:-1] + color_bins[1:]) / 2.

        # Remove points with <1
        x_color, y_resid, y_resid_err = np.array([]), np.array([]), np.array([])
        for i in range(len(color_bins) - 1):
            if len(resid[(colors > color_bins[i]) & (colors < color_bins[i + 1])]) > 1:
                x_color = np.append(x_color, color_bins_adj[i])
                y_resid = np.append(y_resid, binned_vals[i])
                y_resid_err = np.append(y_resid_err, binned_errs[i])

        ## Plot binned dust w/ scatter
        axs[1].errorbar(x_color, y_resid, yerr=y_resid_err, fmt='o-', color=p_color,
                        label=p_label + f"{len(colors)}")

        ## Label number of points in bin
        if label:
            for i in range(len(color_bins) - 1):
                if np.isnan(binned_vals[i]): continue
                if len(resid[(colors > color_bins[i]) & (colors < color_bins[i + 1])]) <= 1: continue
                axs[1].text(color_bins_adj[i], binned_vals[i],
                            len(resid[(colors > color_bins[i]) & (colors < color_bins[i + 1])]),
                            ha='left', va='bottom', size='small')

    # Formatting ======================================================================================================
    axs[0].set_xlabel('Binned SNooPy Color, $E(B-V)_{host}$', size=16)
    axs[0].set_ylabel('Hubble Residual Scatter, $\sigma$', size=16)
    axs[0].legend(loc='upper left')
    axs[1].set_xlabel('Binned SALT3 Color, $c$', size=16)
    axs[1].set_ylabel('Hubble Residual Scatter, $\sigma$', size=16)
    axs[1].legend(loc='upper left')

    # Saving Figure
    if len(save_loc) != 0:
        print('Saved figure to... ', save_loc)
        plt.savefig(save_loc)
    plt.show()
    return
def params_v_scatter(path_snpy_91bg: str = 'results/combiend__snpy_params_cut.txt',
                     path_salt_91bg: str = 'results/combiend__salt_params_cut.txt',
                     path_snpy_norm: str = 'results/output/norm_snpy_params_cut.txt',
                     path_salt_norm: str = 'results/aaronDo_salt2_params_cut.txt',
                     bin_nums: list = [[40, 40], [40, 40], [40, 40], [40, 40]],
                     bin_bounds: list = [[-1, 1], [-1, 1], [-1, 1], [-1, 1]],
                     label: bool = False, save_loc: str = ''):
    fig, axs = plt.subplots(2, 2, figsize=(32, 10), constrained_layout=True)
    plt.style.use('tableau-colorblind10')
    c_91bg, c_norm = 'C8', 'C2'

    # Statisitcs
    tb_1 = utils.default_open(path_snpy_91bg, True)
    tb_2 = utils.default_open(path_salt_91bg, True)
    tb_3 = utils.default_open(path_snpy_norm, True)
    tb_4 = utils.default_open(path_salt_norm, True)
    print(f"sBV: 1991bg-like = [{round(np.min(tb_1['stretch'])-0.1, 3)}, {round(np.max(tb_1['stretch'])+0.1, 3)}], "
          f"{len(tb_1['stretch'])}, Avg. {np.average(tb_1['stretch'])}, Med. {np.median(tb_1['stretch'])} | "
          f"Normals = [{round(np.min(tb_3['stretch'])-0.1, 3)}, {round(np.max(tb_3['stretch'])+0.1, 3)}], "
          f"{len(tb_3['stretch'])}, Avg. {np.average(tb_3['stretch'])}, Med. {np.median(tb_3['stretch'])}")
    print(f"E(B-V): 1991bg-like = [{round(np.min(tb_1['color']) - 0.1, 3)}, {round(np.max(tb_1['color']) + 0.1, 3)}], "
          f"{len(tb_1['color'])}, Avg. {np.average(tb_1['color'])}, Med. {np.median(tb_1['color'])} | "
          f"Normals = [{round(np.min(tb_3['color']) - 0.1, 3)}, {round(np.max(tb_3['color']) + 0.1, 3)}], "
          f"{len(tb_3['color'])}, Avg. {np.average(tb_3['color'])}, Med. {np.median(tb_3['color'])}")
    print(f"x1: 1991bg-like = [{round(np.min(tb_2['stretch']) - 0.1, 3)}, {round(np.max(tb_2['stretch']) + 0.1, 3)}], "
          f"{len(tb_2['stretch'])}, Avg. {np.average(tb_2['stretch'])}, Med. {np.median(tb_2['stretch'])} | "
          f"Normals = [{round(np.min(tb_4['stretch']) - 0.1, 3)}, {round(np.max(tb_4['stretch']) + 0.1, 3)}], "
          f"{len(tb_4['stretch'])}, Avg. {np.average(tb_4['stretch'])}, Med. {np.median(tb_4['stretch'])}")
    print(f"c: 1991bg-like = [{round(np.min(tb_2['color']) - 0.1, 3)}, {round(np.max(tb_2['color']) + 0.1, 3)}], "
          f"{len(tb_2['color'])}, Avg. {np.average(tb_2['color'])}, Med. {np.median(tb_2['color'])} | "
          f"Normals = [{round(np.min(tb_4['color']) - 0.1, 3)}, {round(np.max(tb_4['color']) + 0.1, 3)}], "
          f"{len(tb_4['color'])}, Avg. {np.average(tb_4['color'])}, Med. {np.median(tb_4['color'])}")

    # Top Panels: SNooPy ===============================================================================================
    ## Top-Left s_BV ===================================================================================================
    for path, bin_num, p_label, p_color in zip([path_snpy_91bg, path_snpy_norm],
                                               [bin_nums[0][0], bin_nums[0][1]],
                                               ["1991bg-like SNe Ia, $N_{SNe}$ = ", "Normal SNe Ia, $N_{SNe}$ = "],
                                               [c_91bg, c_norm]):
        plot_binned_param(axis=axs[0, 0], path=path, param_name='stretch',
                          bin_num=bin_num, bin_bounds=bin_bounds[0],
                          p_label=p_label, p_color=p_color)

    ## Top-Right E(B-V) ================================================================================================
    for path, bin_num, p_label, p_color in zip([path_snpy_91bg, path_snpy_norm],
                                               [bin_nums[1][0], bin_nums[1][1]],
                                               ["1991bg-like SNe Ia, $N_{SNe}$ = ", "Normal SNe Ia, $N_{SNe}$ = "],
                                               [c_91bg, c_norm]):
        plot_binned_param(axis=axs[0, 1], path=path, param_name='color',
                          bin_num=bin_num, bin_bounds=bin_bounds[1],
                          p_label=p_label, p_color=p_color)
    # Bottom Panels: SALT3 =============================================================================================
    ## Bottom-Left x_1 =================================================================================================
    for path, bin_num, p_label, p_color in zip([path_salt_91bg, path_salt_norm],
                                               [bin_nums[2][0], bin_nums[2][1]],
                                               ["1991bg-like SNe Ia, $N_{SNe}$ = ", "Normal SNe Ia, $N_{SNe}$ = "],
                                               [c_91bg, c_norm]):
        plot_binned_param(axis=axs[1, 0], path=path, param_name='stretch',
                          bin_num=bin_num, bin_bounds=bin_bounds[2],
                          p_label=p_label, p_color=p_color)
    ## Bottom-Right color ==================================================================================================
    for path, bin_num, p_label, p_color in zip([path_salt_91bg, path_salt_norm],
                                               [bin_nums[3][0], bin_nums[3][1]],
                                               ["1991bg-like SNe Ia, $N_{SNe}$ = ", "Normal SNe Ia, $N_{SNe}$ = "],
                                               [c_91bg, c_norm]):
        plot_binned_param(axis=axs[1, 1], path=path, param_name='color',
                          bin_num=bin_num, bin_bounds=bin_bounds[3],
                          p_label=p_label, p_color=p_color)

    # Formatting ======================================================================================================
    axs[0, 0].set_ylabel('SNooPy\n Hubble Residual Scatter, $\sigma$', size=16)
    axs[0, 0].set_xlabel('Binned SNooPy $s_{BV}$', size=16)
    axs[0, 1].set_xlabel('Binned SNooPy $E(B-V)_{host}$', size=16)
    axs[1, 0].set_ylabel('SALT3\n Hubble Residual Scatter, $\sigma$', size=16)
    axs[1, 0].set_xlabel('Binned SALT3 $x_1$', size=16)
    axs[1, 1].set_xlabel('Binned SALT3 $c$', size=16)
    axs[0, 0].legend(loc='upper left')
    axs[0, 1].legend(loc='upper left')
    axs[1, 0].legend(loc='upper left')
    axs[1, 1].legend(loc='upper left')
    # axs[0, 0].set_ylim(0, 0.3)
    # axs[0, 1].set_ylim(0, 0.3)
    # axs[1, 0].set_ylim(0, 0.3)
    # axs[1, 1].set_ylim(0, 0.3)
    axs[0, 1].tick_params(labelleft=False)
    axs[1, 1].tick_params(labelleft=False)

    # Saving Figure
    if len(save_loc) != 0:
        print('Saved figure to... ', save_loc)
        plt.savefig(save_loc)
    plt.show()
    return
def dust_v_scatter(path_91bg: str = 'results/old/merged_params_cut.txt',
                   path_norm: str = 'results/old/redNormSNe.txt',
                   path_dust: str = 'results/old/global_dust_params.txt',
                   bin_num: int = 20, bin_bounds: list = [0.1, 6.3], hist_bins: int = 20,
                   label: bool = True, save_loc: str = ''):
    fig, axs = plt.subplots(2, 1, figsize=(16, 9), constrained_layout=True, height_ratios=[10, 2])
    plt.style.use('tableau-colorblind10')
    c_91bg, c_norm = 'C8', 'C2'

    # Open data
    tb_91bg = utils.default_open(path_91bg, True)
    tb_norm = utils.default_open(path_norm, True)
    tb_dust = utils.default_open(path_dust, True)
    tb_red_norm = tb_norm[tb_norm['color'] > 0.15]

    for tb, lb, cl in zip([tb_91bg, tb_red_norm],
                          ["1991bg-like SNe Ia", "Normal SNe (c > 0.15)"],
                          [c_91bg, c_norm]):
        # Get associated dust values
        dust = np.array([])
        for name in tb['objname']:
            if name in list(tb_dust['objname']):
                dust = np.append(dust, tb_dust['av_50'][tb_dust['objname'] == name].value[0])
            else:
                tb.remove_row(list(tb['objname']).index(name))  # Removes if no dust value found
        resid = np.array(tb['mu'] - utils.current_cosmo().distmod(tb['z_cmb']).value)

        # Bin Dust & Residuals with STD of residuals
        dust_bins = np.linspace(bin_bounds[0], bin_bounds[1], bin_num)
        binned_vals = scipy.stats.binned_statistic(dust, resid, bins=dust_bins, statistic=np.std).statistic
        binned_errs = scipy.stats.binned_statistic(dust, resid, bins=dust_bins, statistic=bootstrap_errs).statistic
        dust_bins_adj = (dust_bins[:-1] + dust_bins[1:]) / 2.

        # Remove points with <1
        y_resid, y_resid_err, x_dust = np.array([]), np.array([]), np.array([])
        for i in range(len(dust_bins) - 1):
            if len(resid[(dust > dust_bins[i]) & (dust < dust_bins[i + 1])]) > 1:
                x_dust = np.append(x_dust, dust_bins_adj[i])
                y_resid = np.append(y_resid, binned_vals[i])
                y_resid_err = np.append(y_resid_err, binned_errs[i])

        # Plot binned dust w/ scatter
        axs[0].errorbar(x_dust, y_resid, xerr=y_resid_err, fmt='o', color=cl,
                     label=lb+", $N_{SNe}$ = " + f"{len(dust)}")

        # Plot histogram to side
        axs[1].hist(dust, bins=hist_bins, color=cl)

        # Plot median line
        axs[0].axvline(np.median(dust), label=lb+f' Median: {round(np.median(dust), 3)}',
                       color=cl, linestyle='--')
        axs[1].axvline(np.median(dust), color=cl, linestyle='--')

        # Label number of points in bin
        if label:
            for i in range(len(dust_bins)-1):
                if np.isnan(binned_vals[i]): continue
                if len(resid[(dust > dust_bins[i]) & (dust < dust_bins[i + 1])]) <= 1: continue
                axs[0].text(dust_bins_adj[i], binned_vals[i],
                         len(resid[(dust > dust_bins[i]) & (dust < dust_bins[i+1])]),
                         ha='left', va='bottom', size='small')

    # STD of entire data set
    axs[0].axhline(np.std(resid),
                   label=f'Total Hubble Residual Scatter, $\sigma$ = {round(np.std(resid), 3)}',
                   color='C3')

    # Formatting
    axs[0].set_ylabel('Hubble Residual Scatter, $\sigma$', size=16)
    axs[1].set_xlabel('Binned Dust Parameter, $A_{V}$', size=16)
    axs[0].legend(loc='upper right')
    axs[0].tick_params(labelbottom=False)
    axs[0].set_xlim(np.min(x_dust)-0.1, np.max(x_dust)+0.1)
    axs[1].set_xlim(np.min(x_dust)-0.1, np.max(x_dust)+0.1)

    # Saving Figure
    if len(save_loc) != 0:
        print('Saved figure to... ', save_loc)
        plt.savefig(save_loc)
    plt.show()
    return


if __name__ == '__main__':
    start = systime.time()  # Runtime tracker
    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
