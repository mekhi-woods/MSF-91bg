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
    # tb_norm = tb_norm[tb_norm['hostMass']>7.5] # Removes the low mass Normals

    ## Calculate Hubble Residual
    tb_norm['resid_mu'] = tb_norm['mu'] - utils.current_cosmo().distmod(tb_norm['z_cmb']).value
    tb_norm['resid_mu_err'] = np.copy(tb_norm['mu_err'])

    # Subtracting off Average Hubble Residual
    tb_norm['resid_mu'] -= np.average(tb_norm['resid_mu'][~np.isnan(tb_norm['resid_mu'])]) # Centering around average

    ## Scatter plot & histogram
    axs[0,0].errorbar(x=tb_norm['hostMass'], y=tb_norm['resid_mu'], xerr=tb_norm['hostMass_err'], yerr=tb_norm['resid_mu_err'],
                      marker='o', alpha=0.25, color=c_norm, fmt='o', ms=6, elinewidth=0.8)
    axs[0,1].hist(tb_norm['resid_mu'], bins=20, orientation="horizontal", color=c_norm)

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
    tb_91bg['resid_mu_err'] = np.copy(tb_91bg['mu_err'])

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
        axs[1, 0].hlines(y=resid_dict['lower_resid']['value'], xmin=-7.5, xmax=cut, **lin_details)  # Left
        axs[1, 0].hlines(y=resid_dict['upper_resid']['value'], xmin=cut, xmax=13.5, **lin_details)  # Right
        axs[1, 0].axvline(cut, alpha=0.75, **lin_details,
                          label="$\gamma (M_{split}=$" + f"{round(cut, 2)}) = " +
                                f"${round(mass_step_dict['value'], 3)} \pm {round(mass_step_dict['err'], 3)}$ mag")

    # Brount, Scolnic 2021 Dust Prediction
    axs[1, 0].hlines(y=np.average(tb_91bg['resid_mu'][tb_91bg['hostMass'] < 10]) - 0.25,
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
    # norm_mn, norm_md, norm_std = sigma_clipped_stats(tb_norm['resid_mu'])
    # sn91bg_mn, sn91bg_md, sn91bg_std = sigma_clipped_stats(tb_91bg['resid_mu'])
    axs[0,0].text(0.04, 0.96,
                  "Normal SNe Ia\n"+
                  "$N_{SNe}$ = "+f"{len(tb_norm)}\n"+
                  "$\sigma$ = "+f"{np.std(tb_norm['resid_mu']):.3f} mag\n",
                  transform=axs[0, 0].transAxes, ha='left', va='top', fontsize=12)
    axs[1,0].text(0.04, 0.96,
                  "1991bg-like SNe Ia\n" +
                  "$N_{SNe}$ = " + f"{len(tb_91bg)}\n" +
                  "$\sigma$ = " + f"{np.std(tb_91bg['resid_mu']):.3f} mag\n",
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
    ax1.set_xlim(0.015, 0.096)
    ax2.set_xlim(0.015, 0.096)

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
                 slope=alpha_91bg, color=c_91bg_line, label="$\\alpha_{1991bg\\text{-}like}" + f"={-1*alpha_91bg:.2f}$", zorder=10)
    ax[0].axline((0, minimize(get_chi2, 0.00, args=(x1_norm, absmB_norm, absmB_err_norm, alpha_91bg)).x[0]),
                 slope=alpha_91bg, color=c_91bg_line, linestyle='--', zorder=10)
    ax[1].axline((0, minimize(get_chi2, 0.00, args=(c_91bg, absmB_91bg, absmB_err_91bg, beta_91bg)).x[0]),
                 slope=beta_91bg, color=c_91bg_line, label="$\\beta_{1991bg\\text{-}like}" + f"={beta_91bg:.2f}$", zorder=10)
    ax[1].axline((0, minimize(get_chi2, 0.00, args=(c_norm, absmB_norm, absmB_err_norm, beta_91bg)).x[0]),
                 slope=beta_91bg, color=c_91bg_line, linestyle='--', zorder=10)

    # Normal Fit Lines
    ax[0].axline((0, minimize(get_chi2, 0.00, args=(x1_norm, absmB_norm, absmB_err_norm, alpha_norm)).x[0]),
                 slope=alpha_norm, color=c_norm_line, label="$\\alpha_{Normal\\text{ }SNIa}" + f"={-1*alpha_norm:.2f}$", zorder=10)
    ax[0].axline((0, minimize(get_chi2, 0.00, args=(x1_91bg, absmB_91bg, absmB_err_91bg, alpha_norm)).x[0]),
                 slope=alpha_norm, color=c_norm_line, linestyle='--', zorder=10)
    ax[1].axline((0, minimize(get_chi2, 0.00, args=(c_norm, absmB_norm, absmB_err_norm, beta_norm)).x[0]),
                 slope=beta_norm, color=c_norm_line, label="$\\beta_{Normal\\text{ }SNIa}"+f"={beta_norm:.2f}$", zorder=10)
    ax[1].axline((0, minimize(get_chi2, 0.00, args=(c_91bg, absmB_91bg, absmB_err_91bg, beta_norm)).x[0]),
                 slope=beta_norm, color=c_norm_line, linestyle='--', zorder=10)

    # Formatting
    ax[0].set_xlabel('$x_1$', size=18)
    ax[0].set_ylabel('$m_{B} - \mu$', size=18)
    ax[1].set_xlabel('$c$', size=18)
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
               bins: list = [0.1, 0.8, 0.2, 0.15], hist_tol: list = [0.00, 0.00, 0.00, 0.00],
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
    ## sBV plot
    sbv_binwidth = bins[0]
    sbv_stack = np.hstack([tb_snpy_norm['stretch'], tb_snpy_91bg['stretch']])
    axs[0, 0].hist(tb_snpy_norm['stretch'],
                   color=c_norm,
                   bins=np.arange(min(sbv_stack), max(sbv_stack), sbv_binwidth))
    axs[0, 0].hist(tb_snpy_91bg['stretch'],
                   color=c_91bg, alpha=0.75,
                   bins=np.arange(min(sbv_stack), max(sbv_stack), sbv_binwidth))

    # x1 plot
    x1_binwidth = bins[1]
    x1_stack = np.hstack([tb_salt_norm['stretch'], tb_salt_91bg['stretch']])
    axs[1, 0].hist(tb_salt_norm['stretch'],
                   color=c_norm,
                   bins=np.arange(min(x1_stack), max(x1_stack), x1_binwidth))

    axs[1, 0].hist(tb_salt_91bg['stretch'],
                   color=c_91bg, alpha=0.75,
                   bins=np.arange(min(x1_stack), max(x1_stack), x1_binwidth))

    # EBVhost plot
    ebvhost_binwidth = bins[2]
    ebvhost_stack = np.hstack([tb_snpy_norm['color'], tb_snpy_91bg['color']])
    axs[0, 1].hist(tb_snpy_norm['color'],
                   color=c_norm,
                   bins=np.arange(min(ebvhost_stack), max(ebvhost_stack), ebvhost_binwidth))
    axs[0, 1].hist(tb_snpy_91bg['color'],
                   color=c_91bg, alpha=0.75,
                   bins=np.arange(min(ebvhost_stack), max(ebvhost_stack), ebvhost_binwidth))

    # c plot
    c_binwidth = bins[3]
    c_stack = np.hstack([tb_salt_norm['color'], tb_salt_91bg['color']])
    axs[1, 1].hist(tb_salt_norm['color'],
                   color=c_norm,
                   bins=np.arange(min(c_stack), max(c_stack), c_binwidth))
    axs[1, 1].hist(tb_salt_91bg['color'],
                   color=c_91bg, alpha=0.75,
                   bins=np.arange(min(c_stack), max(c_stack), c_binwidth))

    if line_type == 'median':
        line_type = line_type[0].upper() + line_type[1:]

        axs[0, 0].axvline(np.median(tb_snpy_norm['stretch']), color=c_norm_line, linestyle='--', linewidth=3,
                          label=f"{line_type}"+
                                "$_{Normal\\text{ }Ia\\text{ }SNe}$"+
                                f" = ${np.median(tb_snpy_norm['stretch']):.4f}$") # f" = ${round(np.median(tb_snpy_norm['stretch']), 3)}$")
        axs[0, 0].axvline(np.median(tb_snpy_91bg['stretch']), color=c_91bg_line, linestyle=':', linewidth=3,
                          label=f"{line_type}" +
                                "$_{1991bg\\text{-}like\\text{ }Ia\\text{ }SNe}$" +
                                f" = ${np.median(tb_snpy_91bg['stretch']):.4f}$")
        axs[0, 1].axvline(np.median(tb_snpy_norm['color']), color=c_norm_line, linestyle='--', linewidth=3,
                          label=f"{line_type}"+
                                "$_{Normal\\text{ }Ia\\text{ }SNe}$"+
                                f" = ${np.median(tb_snpy_norm['color']):.4f}$")
        axs[0, 1].axvline(np.median(tb_snpy_91bg['color']), color=c_91bg_line, linestyle=':', linewidth=3,
                          label=f"{line_type}" +
                                "$_{1991bg\\text{-}like\\text{ }Ia\\text{ }SNe}$" +
                                f" = ${np.median(tb_snpy_91bg['color']):.4f}$")
        axs[1, 0].axvline(np.median(tb_salt_norm['stretch']), color=c_norm_line, linestyle='--', linewidth=3,
                          label=f"{line_type}"+
                                "$_{Normal\\text{ }Ia\\text{ }SNe}$"+
                                f" = ${np.median(tb_salt_norm['stretch']):.4f}$")
        axs[1, 0].axvline(np.median(tb_salt_91bg['stretch']), color=c_91bg_line, linestyle=':', linewidth=3,
                          label=f"{line_type}" +
                                "$_{1991bg\\text{-}like\\text{ }Ia\\text{ }SNe}$" +
                                f" = ${np.median(tb_salt_91bg['stretch']):.4f}$")
        axs[1, 1].axvline(np.median(tb_salt_norm['color']), color=c_norm_line, linestyle='--', linewidth=3,
                          label=f"{line_type}"+
                                "$_{Normal\\text{ }Ia\\text{ }SNe}$"+
                                f" = ${np.median(tb_salt_norm['color']):.4f}$")
        axs[1, 1].axvline(np.median(tb_salt_91bg['color']), color=c_91bg_line, linestyle=':', linewidth=3,
                          label=f"{line_type}" +
                                "$_{1991bg\\text{-}like\\text{ }Ia\\text{ }SNe}$" +
                                f" = ${np.median(tb_salt_91bg['color']):.4f}$")
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
    axs[0, 1].legend(loc='upper left')
    axs[1, 0].legend(loc='upper right')
    axs[1, 1].legend(loc='upper left')

    # Set labels
    axs[0, 0].set_ylabel('SNooPy\n$N_{SNe}$', size=18)
    axs[1, 0].set_ylabel('SALT3\n$N_{SNe}$', size=18)
    axs[0, 0].set_xlabel('$s_{BV}$', size=18)
    axs[0, 1].set_xlabel('$E(B-V)_{host}$', size=18)
    axs[1, 0].set_xlabel('$x_1$', size=18)
    axs[1, 1].set_xlabel('$c$', size=18)

    # Adjust formatting
    axs[0, 1].tick_params(labelleft=False, labelright=True)
    axs[1, 1].tick_params(labelleft=False, labelright=True)

    # Adjust bounds
    ## sBV
    axs[0, 0].set_xlim(((np.std(tb_snpy_91bg['stretch'] - np.average(tb_snpy_91bg['stretch']))*-2.9) - hist_tol[0]) + 1,
                       ((np.std(tb_snpy_91bg['stretch'] - np.average(tb_snpy_91bg['stretch']))*2.9) + hist_tol[0]) + 1)
    ## x1
    axs[1, 0].set_xlim(np.std(tb_salt_91bg['stretch'] - np.average(tb_salt_91bg['stretch']))*-2.9 - hist_tol[1],
                       np.std(tb_salt_91bg['stretch'] - np.average(tb_salt_91bg['stretch']))*2.9 + hist_tol[1])
    ## EBVhost
    axs[0, 1].set_xlim(np.std(tb_snpy_91bg['color'] - np.average(tb_snpy_91bg['color']))*-2.9 - hist_tol[2],
                       np.std(tb_snpy_91bg['color'] - np.average(tb_snpy_91bg['color']))*2.9 + hist_tol[2])
    ## c
    axs[1, 1].set_xlim(np.std(tb_salt_91bg['color'][np.abs(tb_salt_91bg['color'])<2] - np.average(tb_salt_91bg['color'][np.abs(tb_salt_91bg['color'])<2]))*-2.9 - hist_tol[3],
                       np.std(tb_salt_91bg['color'][np.abs(tb_salt_91bg['color'])<2] - np.average(tb_salt_91bg['color'][np.abs(tb_salt_91bg['color'])<2]))*2.9 + hist_tol[3])

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
    axs[0, 0].set_ylabel('SNooPy\n Hubble Residual Scatter, $\sigma$', size=20)
    axs[0, 0].set_xlabel('Binned SNooPy $s_{BV}$', size=20)
    axs[0, 1].set_xlabel('Binned SNooPy $E(B-V)_{host}$', size=20)
    axs[1, 0].set_ylabel('SALT3\n Hubble Residual Scatter, $\sigma$', size=20)
    axs[1, 0].set_xlabel('Binned SALT3 $x_1$', size=20)
    axs[1, 1].set_xlabel('Binned SALT3 $c$', size=20)
    axs[0, 0].legend(loc='upper left')
    axs[0, 1].legend(loc='upper left')
    axs[1, 0].legend(loc='upper left')
    axs[1, 1].legend(loc='upper left')
    axs[0, 0].set_ylim(0, 0.5)
    axs[0, 1].set_ylim(0, 0.5)
    axs[1, 0].set_ylim(0, 0.5)
    axs[1, 1].set_ylim(0, 0.5)
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
                   path_dust: str = 'txts/global_dust_params.txt',
                   bin_width_91bg: float = 0.5, bin_width_norm: float = 0.1, binning_lim: int = 2,
                   label: bool = True, save_loc: str = ''):
    fig, axs = plt.subplots(1, 1, figsize=(18, 5), constrained_layout=True) # width_ratios=[10, 1]
    plt.style.use('tableau-colorblind10')
    # c_norm, c_norm_contrast, c_norm_mass = 'C7', 'C4', 'C0' # Blues
    c_norm, c_norm_contrast, c_norm_mass = 'C2', 'C6', 'C3' # Greys
    c_91bg, c_91bg_contrast, c_91bg_mass = 'C8', 'C5', 'C1' # Oranges

    # Open data
    tb_91bg = utils.default_open(path_91bg, True)
    tb_norm = utils.default_open(path_norm, True)
    tb_dust = utils.default_open(path_dust, True)
    # tb_norm = tb_norm[tb_norm['color'] > 0.15]


    # Fix dust tables
    # =================================================================================================================
    tb_91bg.add_column(name='av_50', col=np.nan)
    tb_norm.add_column(name='av_50', col=np.nan)
    for n in tb_dust['objname']:
        if n in tb_91bg['objname']:
            tb_91bg['av_50'][tb_91bg['objname'] == n] = tb_dust['av_50'][tb_dust['objname'] == n]
        if n in tb_norm['objname']:
            tb_norm['av_50'][tb_norm['objname'] == n] = tb_dust['av_50'][tb_dust['objname'] == n]

    tb_91bg = tb_91bg[~np.isnan(tb_91bg['av_50'])]
    tb_norm = tb_norm[~np.isnan(tb_norm['av_50'])]

    # Adjust and save resid
    tb_91bg['resid'], tb_91bg['resid_err'] = tb_91bg['mu'] - utils.current_cosmo().distmod(tb_91bg['z_cmb']).value, tb_91bg['mu_err']
    tb_91bg['resid'] = tb_91bg['resid'] - np.average(tb_91bg['resid'])
    tb_norm['resid'], tb_norm['resid_err'] = tb_norm['mu'] - utils.current_cosmo().distmod(tb_norm['z_cmb']).value, tb_norm['mu_err']
    tb_norm['resid'] = tb_norm['resid'] - np.average(tb_norm['resid'])

    # Bin scatter with dust
    # =================================================================================================================
    # dust_bins_91bg = np.arange(0.15, 6.0, bin_width_91bg)
    # binned_scatter_91bg, binned_dust_91bg, counts_91bg = [], [], []
    # for i in range(len(dust_bins_91bg) - 1):
    #     n_tb_91bg = tb_91bg[(tb_91bg['av_50'] > dust_bins_91bg[i]) & (tb_91bg['av_50'] < dust_bins_91bg[i + 1])]
    #
    #     if len(n_tb_91bg) > 2:
    #         binned_scatter_91bg.append(np.std(n_tb_91bg['resid']))
    #         binned_dust_91bg.append((dust_bins_91bg[i] + dust_bins_91bg[i + 1]) / 2)
    #         counts_91bg.append(len(n_tb_91bg))
    #     else:
    #         continue
    #
    # dust_bins_norm = np.arange(0.15, 6.0, bin_width_norm)
    # binned_scatter_norm, binned_dust_norm, counts_norm = [], [], []
    # for i in range(len(dust_bins_norm) - 1):
    #     n_tb_norm = tb_norm[(tb_norm['av_50'] > dust_bins_norm[i]) & (tb_norm['av_50'] < dust_bins_norm[i + 1])]
    #
    #     if len(n_tb_norm) > 2:
    #         binned_scatter_norm.append(np.std(n_tb_norm['resid']))
    #         binned_dust_norm.append((dust_bins_norm[i] + dust_bins_norm[i + 1]) / 2)
    #         counts_norm.append(len(n_tb_norm))
    #     else:
    #         continue

    # Make bins
    dust_bins_91bg = np.arange(0.15, 6.0, bin_width_91bg)
    binned_vals_91bg = scipy.stats.binned_statistic(tb_91bg['av_50'], tb_91bg['resid'], bins=dust_bins_91bg, statistic=np.std).statistic
    binned_errs_91bg = scipy.stats.binned_statistic(tb_91bg['av_50'], tb_91bg['resid'], bins=dust_bins_91bg, statistic=bootstrap_errs).statistic
    dust_bins_91bg_adj = (dust_bins_91bg[:-1] + dust_bins_91bg[1:]) / 2.

    dust_bins_norm = np.arange(0.15, 6.0, bin_width_norm)
    binned_vals_norm = scipy.stats.binned_statistic(tb_norm['av_50'], tb_norm['resid'], bins=dust_bins_norm, statistic=np.std).statistic
    binned_errs_norm = scipy.stats.binned_statistic(tb_norm['av_50'], tb_norm['resid'], bins=dust_bins_norm, statistic=bootstrap_errs).statistic
    dust_bins_norm_adj = (dust_bins_norm[:-1] + dust_bins_norm[1:]) / 2.

    # Remove points with less than 'binning_lim'
    binned_dust_91bg, binned_scatter_91bg, binned_scatter_err_91bg, counts_91bg = np.array([]), np.array([]), np.array([]), np.array([])
    for i in range(len(dust_bins_91bg) - 1):
        arr_len = len(tb_91bg['resid'][(tb_91bg['av_50'] > dust_bins_91bg[i]) & (tb_91bg['av_50'] < dust_bins_91bg[i + 1])])
        if arr_len > binning_lim:
            binned_dust_91bg = np.append(binned_dust_91bg, dust_bins_91bg_adj[i])
            binned_scatter_91bg = np.append(binned_scatter_91bg, binned_vals_91bg[i])
            binned_scatter_err_91bg = np.append(binned_scatter_err_91bg, binned_errs_91bg[i])
            counts_91bg = np.append(counts_91bg, arr_len)

    binned_dust_norm, binned_scatter_norm, binned_scatter_err_norm, counts_norm = np.array([]), np.array([]), np.array([]), np.array([])
    for i in range(len(dust_bins_norm) - 1):
        arr_len = len(tb_norm['resid'][(tb_norm['av_50'] > dust_bins_norm[i]) & (tb_norm['av_50'] < dust_bins_norm[i + 1])])
        if arr_len > binning_lim:
            binned_dust_norm = np.append(binned_dust_norm, dust_bins_norm_adj[i])
            binned_scatter_norm = np.append(binned_scatter_norm, binned_vals_norm[i])
            binned_scatter_err_norm = np.append(binned_scatter_err_norm, binned_errs_norm[i])
            counts_norm = np.append(counts_norm, arr_len)




    # Plot
    # =================================================================================================================
    plt.errorbar(x=binned_dust_91bg, y=binned_scatter_91bg,
                 yerr=binned_scatter_err_91bg,
                 fmt='o-', color=c_91bg, alpha=1.0, zorder=2,
                 label="1991bg-like SNe Ia\n"
                       "$N_{SNe}$ = "+f"{len(tb_91bg)}")
    plt.errorbar(x=binned_dust_norm, y=binned_scatter_norm,
                 fmt='o-', color=c_norm, alpha=1.0, zorder=2,
                 label="Normal SNe Ia\n"
                       "$N_{SNe}$ = "+f"{len(tb_norm)}")

    # Label
    if label:
        for x, y, c in zip(binned_dust_91bg, binned_scatter_91bg, counts_91bg.astype(int)):
            plt.text(x, y, c, ha='right', va='bottom', size='small')
        for x, y, c in zip(binned_dust_norm, binned_scatter_norm, counts_norm.astype(int)):
            plt.text(x, y, c, ha='right', va='bottom', size='small')

    # Formatting
    plt.ylabel('Hubble Residual Scatter, $\sigma$', size=16)
    plt.xlabel('Binned Dust Parameter, $A_{V}$', size=16)
    plt.legend(loc='upper right')

    # Saving Figure
    if len(save_loc) != 0:
        print('Saved figure to... ', save_loc)
        plt.savefig(save_loc)
    plt.show()
    return
def dust_mass_step(path_91bg: str, path_norm: str, path_dust: str, save_loc: str = '', label: bool = False):
    fig, axs = plt.subplots(2, 2, figsize=(18, 10), gridspec_kw={'width_ratios': [10, 1]},
                            constrained_layout=True)
    plt.style.use('tableau-colorblind10')
    c_norm, c_norm_contrast, c_norm_mass = 'C7', 'C4', 'C0' # Blues
    # c_norm, c_norm_contrast, c_norm_mass = 'C2', 'C6', 'C3' # Greys
    c_91bg, c_91bg_contrast, c_91bg_mass = 'C8', 'C5', 'C1' # Oranges

    tb_91bg = utils.default_open(path_91bg, True)
    tb_norm = utils.default_open(path_norm, True)
    tb_dust = utils.default_open(path_dust, True)

    # Fix dust tables
    # =================================================================================================================
    tb_91bg.add_column(name='av_50', col=np.nan)
    tb_norm.add_column(name='av_50', col=np.nan)
    for n in tb_dust['objname']:
        if n in tb_91bg['objname']:
            tb_91bg['av_50'][tb_91bg['objname'] == n] = tb_dust['av_50'][tb_dust['objname'] == n]
        if n in tb_norm['objname']:
            tb_norm['av_50'][tb_norm['objname'] == n] = tb_dust['av_50'][tb_dust['objname'] == n]

    with open('need_dust.txt', 'w') as f:
        print('objname,ra,dec', file=f)
        for n in tb_91bg[np.isnan(tb_91bg['av_50'])]:
            print(f"{n['objname']},{n['ra']},{n['dec']}", file=f)

    tb_91bg = tb_91bg[~np.isnan(tb_91bg['av_50'])]
    tb_norm = tb_norm[~np.isnan(tb_norm['av_50'])]

    # Adjust and save resid
    tb_91bg['resid'], tb_91bg['resid_err'] = tb_91bg['mu'] - utils.current_cosmo().distmod(tb_91bg['z_cmb']).value, tb_91bg['mu_err']
    tb_91bg['resid'] = tb_91bg['resid'] - np.average(tb_91bg['resid'])
    tb_norm['resid'], tb_norm['resid_err'] = tb_norm['mu'] - utils.current_cosmo().distmod(tb_norm['z_cmb']).value, tb_norm['mu_err']
    tb_norm['resid'] = tb_norm['resid'] - np.average(tb_norm['resid'])

    # Plot
    # =================================================================================================================
    # =================================================================================================================
    bin_width = 0.05

    # Plot 1991bg-like
    median_dust = np.median(tb_91bg['av_50'])
    ## Plot Above Median
    axs[1,0].scatter(tb_91bg[tb_91bg['av_50'] > median_dust]['hostMass'], tb_91bg[tb_91bg['av_50'] > median_dust]['resid'],
                   marker='o', s=55, alpha=1.0, zorder=2, color=c_91bg_contrast,
                   label='Above Median Dust ($A_{V,Med}=$'+f"{median_dust:.2f} mag)")
    # axs[1,0].errorbar(x=tb_91bg[tb_91bg['av_50'] > median_dust]['hostMass'],
    #                   y=tb_91bg[tb_91bg['av_50'] > median_dust]['resid'],
    #                   xerr=tb_91bg[tb_91bg['av_50'] > median_dust]['hostMass_err'],
    #                   yerr=tb_91bg[tb_91bg['av_50'] > median_dust]['resid_err'],
    #                   fmt='o', ms=1, ecolor=c_91bg_mass, alpha=1.0, zorder=1)
    axs[1,1].hist(tb_91bg[tb_91bg['av_50'] > median_dust]['resid'],
                bins=int((max(tb_91bg[tb_91bg['av_50'] > median_dust]['resid']) - min(tb_91bg[tb_91bg['av_50'] > median_dust]['resid'])) / bin_width),
                orientation="horizontal", color=c_91bg_contrast, alpha=1.0, zorder=1)
    ## Plot Below Median
    axs[1,0].scatter(tb_91bg[tb_91bg['av_50'] < median_dust]['hostMass'], tb_91bg[tb_91bg['av_50'] < median_dust]['resid'],
                   marker='o', s=55, alpha=1.0, zorder=2, color=c_91bg,
                   label=f'Below Median Dust') # 1991bg-like SNe Ia
    # axs[1,0].errorbar(tb_91bg[tb_91bg['av_50'] < median_dust]['hostMass'],
    #                   tb_91bg[tb_91bg['av_50'] < median_dust]['resid'],
    #                   xerr=tb_91bg[tb_91bg['av_50'] < median_dust]['hostMass_err'],
    #                   yerr=tb_91bg[tb_91bg['av_50'] < median_dust]['resid_err'],
    #                   fmt='o', ms=1, ecolor=c_91bg_mass, alpha=1.0, zorder=1)
    axs[1,1].hist(tb_91bg[tb_91bg['av_50'] < median_dust]['resid'],
                bins=int((max(tb_91bg[tb_91bg['av_50'] < median_dust]['resid']) - min(tb_91bg[tb_91bg['av_50'] < median_dust]['resid'])) / bin_width),
                orientation="horizontal", color=c_91bg, alpha=0.75, zorder=1)


    # Plot Mass Lines
    # =================================================================================================================
    tol = 1
    median_91bg_mass = np.median(tb_91bg['hostMass'])
    # Plot 10dex Mass Line
    cut, ls, cl = 10, '-', c_91bg_mass
    lin_details = {'linestyle': ls, 'linewidth': 3, 'color': cl, 'zorder': 5}
    mass_step_dict, resid_dict = mass_step_calc(tb_91bg['mu'], tb_91bg['mu_err'], tb_91bg['resid'],
                                                tb_91bg['hostMass'], tb_91bg['z_cmb'], cut=cut)
    if resid_dict['lower_resid']['value'] > resid_dict['upper_resid']['value']: mass_step_dict['value'] = mass_step_dict['value']*-1
    axs[1,0].hlines(y=resid_dict['lower_resid']['value'], xmin=np.min(tb_91bg['hostMass']) - tol, xmax=cut, **lin_details)  # Left
    axs[1,0].hlines(y=resid_dict['upper_resid']['value'], xmin=cut, xmax=np.max(tb_91bg['hostMass']) + tol, **lin_details)  # Right
    axs[1,0].axvline(cut, alpha=0.75, **lin_details,
                   label="$\gamma (M_{split}=$" + f"{round(cut, 2)}) = " +
                         f"${round(mass_step_dict['value'], 3)} \pm {round(mass_step_dict['err'], 3)}$ mag")
    # ## Plot Median Mass Line
    # cut, ls, cl = median_91bg_mass, '--', c_91bg_mass
    # lin_details = {'linestyle': ls, 'linewidth': 3, 'color': cl, 'zorder': 5}
    # mass_step_dict, resid_dict = mass_step_calc(tb_91bg['mu'], tb_91bg['mu_err'], tb_91bg['resid'],
    #                                             tb_91bg['hostMass'], tb_91bg['z_cmb'], cut=cut)
    # if resid_dict['lower_resid']['value'] > resid_dict['upper_resid']['value']: mass_step_dict['value'] = mass_step_dict['value']*-1
    # axs[1,0].hlines(y=resid_dict['lower_resid']['value'], xmin=np.min(tb_91bg['hostMass']) - tol, xmax=cut, **lin_details)  # Left
    # axs[1,0].hlines(y=resid_dict['upper_resid']['value'], xmin=cut, xmax=np.max(tb_91bg['hostMass']) + tol, **lin_details)  # Right
    # axs[1,0].axvline(cut, alpha=0.75, **lin_details,
    #                label="$\gamma (M_{split}=$" + f"{round(cut, 2)}) = " +
    #                      f"${round(mass_step_dict['value'], 3)} \pm {round(mass_step_dict['err'], 3)}$ mag")

    # Plot Normal
    # =================================================================================================================
    # =================================================================================================================
    bin_width = 0.05

    # Plot Normals
    median_dust = np.median(tb_norm['av_50'])
    ## Plot Above Median
    axs[0,0].scatter(tb_norm[tb_norm['av_50'] > median_dust]['hostMass'], tb_norm[tb_norm['av_50'] > median_dust]['resid'],
                   marker='o', s=55, alpha=1.0, zorder=1, color=c_norm_contrast,
                   label='Above Median Dust ($A_{V,Med}=$'+f"{median_dust:.2f} mag)")
    # axs[0,0].errorbar(tb_norm['hostMass'], tb_norm['resid'], xerr=tb_norm['hostMass_err'], yerr=tb_norm['resid_err'],
    #                 fmt='o', ms=1, ecolor=c_91bg_mass, alpha=1.0, zorder=3)
    axs[0,1].hist(tb_norm[tb_norm['av_50'] > median_dust]['resid'],
                bins=int((max(tb_norm[tb_norm['av_50'] > median_dust]['resid']) - min(tb_norm[tb_norm['av_50'] > median_dust]['resid'])) / bin_width),
                orientation="horizontal", color=c_norm_contrast, alpha=1.0, zorder=1)
    ## Plot Below Median
    axs[0,0].scatter(tb_norm[tb_norm['av_50'] < median_dust]['hostMass'], tb_norm[tb_norm['av_50'] < median_dust]['resid'],
                   marker='o', s=55, alpha=1.0, zorder=1, color=c_norm,
                   label=f'Below Median Dust') # 1991bg-like SNe Ia
    # axs[0,0].errorbar(tb_norm['hostMass'], tb_norm['resid'], xerr=tb_norm['hostMass_err'], yerr=tb_norm['resid_err'],
    #                 fmt='o', ms=1, ecolor=c_91bg_mass, alpha=1.0, zorder=3)
    axs[0,1].hist(tb_norm[tb_norm['av_50'] < median_dust]['resid'],
                bins=int((max(tb_norm[tb_norm['av_50'] < median_dust]['resid']) - min(tb_norm[tb_norm['av_50'] < median_dust]['resid'])) / bin_width),
                orientation="horizontal", color=c_norm, alpha=0.75, zorder=1)


    # Plot Mass Lines
    # =================================================================================================================
    tol = 1
    median_91bg_mass = np.median(tb_norm['hostMass'])
    # Plot 10dex Mass Line
    cut, ls, cl = 10, '-', c_norm_mass
    lin_details = {'linestyle': ls, 'linewidth': 3, 'color': cl, 'zorder': 5}
    mass_step_dict, resid_dict = mass_step_calc(tb_norm['mu'], tb_norm['mu_err'], tb_norm['resid'],
                                                tb_norm['hostMass'], tb_norm['z_cmb'], cut=cut)
    if resid_dict['lower_resid']['value'] > resid_dict['upper_resid']['value']: mass_step_dict['value'] = mass_step_dict['value']*-1
    axs[0,0].hlines(y=resid_dict['lower_resid']['value'], xmin=np.min(tb_norm['hostMass']) - tol, xmax=cut, **lin_details)  # Left
    axs[0,0].hlines(y=resid_dict['upper_resid']['value'], xmin=cut, xmax=np.max(tb_norm['hostMass']) + tol, **lin_details)  # Right
    axs[0,0].axvline(cut, alpha=0.75, **lin_details,
                   label="$\gamma (M_{split}=$" + f"{round(cut, 2)}) = " +
                         f"${round(mass_step_dict['value'], 3)} \pm {round(mass_step_dict['err'], 3)}$ mag")
    # ## Plot Median Mass Line
    # cut, ls, cl = median_91bg_mass, '--', c_91bg_mass
    # lin_details = {'linestyle': ls, 'linewidth': 3, 'color': cl, 'zorder': 5}
    # mass_step_dict, resid_dict = mass_step_calc(tb_norm['mu'], tb_norm['mu_err'], tb_norm['resid'],
    #                                             tb_norm['hostMass'], tb_norm['z_cmb'], cut=cut)
    # if resid_dict['lower_resid']['value'] > resid_dict['upper_resid']['value']: mass_step_dict['value'] = mass_step_dict['value']*-1
    # axs[0,0].hlines(y=resid_dict['lower_resid']['value'], xmin=np.min(tb_norm['hostMass']) - tol, xmax=cut, **lin_details)  # Left
    # axs[0,0].hlines(y=resid_dict['upper_resid']['value'], xmin=cut, xmax=np.max(tb_norm['hostMass']) + tol, **lin_details)  # Right
    # axs[0,0].axvline(cut, alpha=0.75, **lin_details,
    #                label="$\gamma (M_{split}=$" + f"{round(cut, 2)}) = " +
    #                      f"${round(mass_step_dict['value'], 3)} \pm {round(mass_step_dict['err'], 3)}$ mag")

    # Quadrant Labels
    # =================================================================================================================
    axs[0,0].text(0.04, 0.96,
                  "Normal SNe Ia\n"+
                  "$N_{SNe}$ = "+f"{len(tb_norm)}\n"+
                  "$\sigma$ = "+f"{np.std(tb_norm['resid']):.3f} mag\n",
                  transform=axs[0, 0].transAxes, ha='left', va='top', fontsize=12)
    axs[1,0].text(0.04, 0.96,
                  "1991bg-like SNe Ia\n" +
                  "$N_{SNe}$ = " + f"{len(tb_91bg)}\n" +
                  "$\sigma$ = " + f"{np.std(tb_91bg['resid']):.3f} mag\n",
                  transform=axs[1, 0].transAxes, ha='left', va='top', fontsize=12)

    # Adjust Axises
    # =================================================================================================================
    axs[1,0].set(xlim=(7.7, 11.4), ylim=(-0.7, 0.7))
    axs[1,1].set(ylim=(-0.7, 0.7))
    axs[0,0].set(xlim=(7.7, 11.4), ylim=(-0.7, 0.7))
    axs[0,1].set(ylim=(-0.7, 0.7))

    axs[0,1].tick_params(labelleft=False, labelbottom=False)
    axs[0,1].tick_params(labelleft=False)
    axs[1,1].tick_params(labelleft=False)
    axs[0,0].tick_params(labelbottom=False)

    # Labels
    axs[1,0].set_ylabel('Hubble Residual (mag)', size=16)
    axs[0,0].set_ylabel('Hubble Residual (mag)', size=16)
    axs[1,0].set_xlabel("Host Stellar Mass ($\log M_{*}[M_{\odot}]$)", size=16)
    axs[1,0].legend(loc='lower left')
    axs[0,0].legend(loc='lower left')

    # Saving Figure
    if len(save_loc) != 0:
        print('Saved figure to... ', save_loc)
        plt.savefig(save_loc)
    plt.show()
    return
def distance_modulus_residuals(path_snpy: str = 'results/norm_snpy_params_cut.txt',
                               path_salt: str = 'results/norm_salt_params_cut.txt'):

    tb_snpy = utils.default_open(path_snpy, True)
    tb_salt = utils.default_open(path_salt, True)

    mu_residuals = {'objname': np.array([]), 'hostMass': np.array([]),
                    'snpyMU': np.array([]), 'saltMU': np.array([]),
                    'dmu': np.array([])}
    for name in np.unique(list(tb_snpy["objname"])+list(tb_salt["objname"])):
        if (name in list(tb_snpy["objname"])) & (name in list(tb_salt["objname"])):
            row_snpy = tb_snpy[tb_snpy["objname"] == name]
            row_salt = tb_salt[tb_salt["objname"] == name]
            mu_residuals['objname'] = np.append(mu_residuals['objname'], name)
            mu_residuals['hostMass'] = np.append(mu_residuals['hostMass'], row_snpy['hostMass'])
            mu_residuals['snpyMU'] = np.append(mu_residuals['snpyMU'], row_snpy['mu'])
            mu_residuals['saltMU'] = np.append(mu_residuals['saltMU'], row_salt['mu'])
            mu_residuals['dmu'] = np.append(mu_residuals['dmu'], np.abs(row_snpy['mu']- row_salt['mu']))
            mu_residuals['dmu_s'] = np.append(mu_residuals['dmu'], row_snpy['mu'] - row_salt['mu'])

    # Correcting for average distance modulus residual
    avg_dmu = np.average(mu_residuals['dmu'])
    mu_residuals.update({'dmu_corr': np.array([])})
    mu_residuals.update({'saltMU_corr': np.array([])})
    for n_snpyMU, n_saltMU in zip(mu_residuals['snpyMU'], mu_residuals['saltMU']):
        if (n_snpyMU - n_saltMU) < 0:
            mu_residuals['saltMU_corr'] = np.append(mu_residuals['saltMU_corr'],
                                                    n_saltMU - avg_dmu)
            mu_residuals['dmu_corr'] = np.append(mu_residuals['dmu_corr'],
                                                 np.abs(n_snpyMU - (n_saltMU-avg_dmu)))
        elif (n_snpyMU - n_saltMU) > 0:
            mu_residuals['saltMU_corr'] = np.append(mu_residuals['saltMU_corr'],
                                                    n_saltMU + avg_dmu)
            mu_residuals['dmu_corr'] = np.append(mu_residuals['dmu_corr'],
                                                 np.abs(n_snpyMU - (n_saltMU + avg_dmu)))


    # # MU v. MU
    # fig = plt.figure(figsize=[8, 8])
    # plt.scatter(mu_residuals['snpyMU'], mu_residuals['saltMU'])
    # plt.title('SNooPy v. SALT3 Dist. Mod. Scatter')
    # plt.xlabel('SNooPy Dist. Mod.', size=16)
    # plt.ylabel('SALT3 Dist. Mod.', size=16)
    # plt.show()

    # # MU v. dMU
    # fig, axs = plt.subplots(1, 2, figsize=(16, 8), constrained_layout=True)
    # axs[0].scatter(mu_residuals['snpyMU'], mu_residuals['dmu'])
    # axs[0].set_xlabel('SNooPy Dist. Mod.', size=16)
    # axs[0].set_ylabel('Absoulute Residual Dist. Mod.', size=16)
    # axs[1].scatter(mu_residuals['saltMU'], mu_residuals['dmu'])
    # axs[1].set_xlabel('SALT3 Dist. Mod.', size=16)
    # axs[1].tick_params(labelleft=False, labelright=True)
    # plt.suptitle('SNooPy v. SALT3 Dist. Mod. Scatter')
    # plt.show()

    # MU v. MU, corrected for average distance modulus residual
    fig, axs = plt.subplots(2, 2, figsize=(16, 8), height_ratios=([10, 1]), constrained_layout=True)

    axs[0,0].scatter(mu_residuals['saltMU'], mu_residuals['snpyMU'])
    axs[1,0].scatter(mu_residuals['saltMU'], mu_residuals['dmu'], s=8)
    axs[0,0].text(0.04, 0.96,
               "Uncorrected Distance Moduli\n"+
               "$N_{SNe}$ = "+f"{len(mu_residuals['snpyMU'])}\n"+
               "$\sigma_{\mu}$ = "+f"{np.std(mu_residuals['dmu']):.3f} mag\n",
               transform=axs[0,0].transAxes, ha='left', va='top', fontsize=14)

    axs[0,1].scatter(mu_residuals['saltMU_corr'], mu_residuals['snpyMU'])
    axs[1,1].scatter(mu_residuals['saltMU_corr'], mu_residuals['dmu_corr'], s=8)
    axs[0,1].text(0.04, 0.96,
               f"Corrected Distance Moduli (Avg. dmu = {avg_dmu:.2f} mag)\n"+
               "$N_{SNe}$ = " + f"{len(mu_residuals['snpyMU'])}\n" +
               "$\sigma_{\mu}$ = "+f"{np.std(mu_residuals['dmu_corr']):.3f} mag\n",
               transform=axs[0,1].transAxes, ha='left', va='top', fontsize=14)

    axs[0,1].tick_params(labelleft=False, labelright=True)
    axs[0,0].tick_params(labelbottom=False)
    axs[0,1].tick_params(labelbottom=False)
    axs[1,1].tick_params(labelleft=False)

    mu_min = min([min(mu_residuals['snpyMU']), min(mu_residuals['saltMU'])]) + 0.1
    mu_max = max([max(mu_residuals['snpyMU']), max(mu_residuals['saltMU'])]) + 0.1
    dmu_max = max(mu_residuals['dmu']) + 0.1
    mu_min, mu_max, dmu_max = 34.25, 38, 0.5
    axs[0,0].set_xlim(mu_min, mu_max)
    axs[0,1].set_xlim(mu_min, mu_max)
    axs[1,0].set_xlim(mu_min, mu_max)
    axs[1,1].set_xlim(mu_min, mu_max)
    axs[0,0].set_ylim(mu_min, mu_max)
    axs[0,1].set_ylim(mu_min, mu_max)
    axs[1,0].set_ylim(0, dmu_max)
    axs[1,1].set_ylim(0, dmu_max)

    axs[0,0].set_ylabel('SNooPy Dist. Mod.', size=16)
    axs[1,0].set_ylabel('Residuals', size=16)
    axs[1,0].set_xlabel('SALT3 Dist. Mod.', size=16)
    axs[1,1].set_xlabel('Corrected SALT3 Dist. Mod.', size=16)


    plt.show()

    return

if __name__ == '__main__':
    start = systime.time()  # Runtime tracker
    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
