import matplotlib.pyplot as plt
import numpy as np
import utils
from scipy.stats import binned_statistic
from scipy.stats import bootstrap


def bootstrap_errs(x):
    try:
        return bootstrap([x], np.std).standard_error
    except:
        return 0
def plot_binned_param(axis, path, param_name, bin_num, bin_bounds, p_label, p_color, label_bins: bool = True):
    # Open data
    tb = utils.default_open(path, True)
    param = np.array(tb[param_name])
    resid = np.array(tb['mu'] - utils.current_cosmo().distmod(tb['z_cmb']).value)

    ## Bin Dust & Residuals with STD of residuals
    param_bins = np.linspace(bin_bounds[0], bin_bounds[1], bin_num)
    binned_vals = binned_statistic(param, resid, bins=param_bins, statistic=np.std).statistic
    binned_errs = binned_statistic(param, resid, bins=param_bins, statistic=bootstrap_errs).statistic
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
def params_v_scatter(path_snpy_91bg: str = '91bg_snpy_cut.txt',
                     path_salt_91bg: str = '91bg_salt_cut.txt',
                     path_snpy_norm: str = 'norm_snpy_cut.txt',
                     path_salt_norm: str = 'norm_salt_cut.txt', # aaronDo_salt2_params_cut.txt
                     bin_nums: list = [[15, 20], [11, 20], [10, 20], [10, 20]],
                     bin_bounds: list = [[0.076, 1.757], [-0.723, 1.078], [-3.289, 1.953], [-0.162, 0.7]],
                     label: bool = False, save_loc: str = 'params_v_scatter.png'):
    fig, axs = plt.subplots(2, 2, figsize=(32, 10), constrained_layout=True)
    plt.style.use('tableau-colorblind10')
    c_91bg, c_norm = 'C8', 'C2'

    # Statistics
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

# [[15, 20], [11, 20], [10, 20], [10, 20]]
# [[top_let_91bg, top_let_norms], [top_let_91bg, top_let_norms]]
# [[0.076, 1.757], [-0.723, 1.078], [-3.289, 1.953], [-0.162, 0.7]]


if __name__ == '__main__':
    params_v_scatter()
