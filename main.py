import dataManipulation
import plotter
import convert
import fitter
import utils
import os
import glob
import shutil

def run_fitter(subtype: str = '91bg', algo: str = 'snpy', rewrite: bool = True):

    sne_csp, sne_atlas, sne_ztf, sne_atlasztf = [], [], [], []
    sne_csp = fitter.fit(glob.glob(f'data/{subtype.upper()}_CSP_*.txt'), algo, rewrite)
    sne_atlas = fitter.fit(glob.glob(f'data/{subtype.upper()}_ATLAS_*.txt'), algo, rewrite)
    sne_ztf = fitter.fit(glob.glob(f'data/{subtype.upper()}_ZTF_*.txt'), algo, rewrite)
    sne_atlasztf = fitter.fit(glob.glob(f'data/{subtype.upper()}_ATLAS-ZTF_*.txt'), algo, rewrite)

    all_sne = sne_csp + sne_atlas + sne_ztf + sne_atlasztf

    dataManipulation.make_param_file(all_sne, '91bg', f'results/{subtype}_{algo}_params.txt')

    return
def run_lc_checker(subtype: str = '91bg', clear_old: bool = True):
    if subtype == '91bg':
        path = "results/91bg_snpy-salt_params_cut.txt"
    elif subtype == 'norm':
        path = "results/norm_snpy-salt_params_cut.txt"
    else:
        return

    if clear_old:
        if os.path.exists(f"fitting/final_lc/{subtype}/"):
            print(f"[~~~] Replacing files @ fitting/final_lc/{subtype}/...")
            shutil.rmtree(f"fitting/final_lc/{subtype}")
            os.mkdir(f"fitting/final_lc/{subtype}/")

    tb = utils.default_open(path, True)
    for name, algo, source in zip(tb['objname'], tb['algo'], tb['origin']):
        n_o = f"fitting/{algo.lower()}-plots/{name}_lc.png"
        n_f = f"fitting/final_lc/{subtype}/{name}_lc.png"
        shutil.copyfile(n_o, n_f)
        # print(f"[~~~] Copied {name}: '{n_o}' -> '{n_f}'...")
    print(f"[~~~] Copied '{path}' light curves to 'fitting/final_lc/{subtype}/'...")
    return
def run_plotter(final_dir: str = 'plots/'):
    pm_norm_salt_cut = 'results/norm_salt_params_cut.txt'
    pm_norm_salt_uncut = 'results/norm_salt_params.txt'
    pm_norm_snpy_cut = 'results/norm_snpy_params_cut.txt'
    pm_norm_snpy_uncut = 'results/norm_snpy_params.txt'
    pm_norm_merged_cut = 'results/norm_snpy-salt_params_cut.txt'
    pm_norm_merged_uncut = 'results/norm_snpy-salt_params.txt'

    pm_91bg_salt_cut = 'results/91bg_salt_params_cut.txt'
    pm_91bg_salt_uncut = 'results/91bg_salt_params.txt'
    pm_91bg_snpy_cut = 'results/91bg_snpy_params_cut.txt'
    pm_91bg_snpy_uncut = 'results/91bg_snpy_params.txt'
    pm_91bg_merged_cut = 'results/91bg_snpy-salt_params_cut.txt'
    pm_91bg_merged_uncut = 'results/91bg_snpy-salt_params.txt'

    # pm_dust = 'txts/global_dust_params.txt'
    pm_dust = 'txts/local_dust_params.txt'

    # # WORKING =========================================================================================================
    plotter.resid_v_mass(path_91bg=pm_91bg_merged_cut,
                         path_norm=pm_norm_merged_cut,
                         save_loc=final_dir + 'resid_v_mass.png',
                         label = True)

    # SALT3 Plots
    plotter.alpha_beta(path_91bg=pm_91bg_salt_cut,
                        path_norm=pm_norm_salt_cut,
                        save_loc=final_dir+'alpha_beta.png')

    # Paramater Histograms
    plotter.param_hist(snpy_91bg_path=pm_91bg_snpy_cut, snpy_norm_path=pm_norm_snpy_cut,
                       salt_91bg_path=pm_91bg_salt_cut, salt_norm_path=pm_norm_salt_cut,
                       bins = [0.08, 0.5, 0.05, 0.05], hist_tol = [0.40, 2.00, 0.00, 0.15],
                       save_loc=final_dir + 'param_hist_cut.png', line_type='median')
    plotter.param_hist(snpy_91bg_path=pm_91bg_snpy_uncut, snpy_norm_path=pm_norm_snpy_uncut,
                       salt_91bg_path=pm_91bg_salt_uncut, salt_norm_path=pm_norm_salt_uncut,
                       bins=[0.1, 0.8, 0.2, 0.15], hist_tol = [0.00, 0.00, 0.00, 0.00],
                       save_loc=final_dir + 'param_hist_uncut.png', line_type='median')

    # Brout+Scolnic 2021 style Paramaters (Color & Stretch) v. Scatter
    plotter.params_v_scatter(path_snpy_91bg=pm_91bg_snpy_cut, path_snpy_norm=pm_norm_snpy_cut,
                             path_salt_91bg=pm_91bg_salt_cut, path_salt_norm=pm_norm_salt_cut,
                             bin_nums=[[15, 20], [11, 20], [10, 20], [10, 20]],
                             bin_bounds=[[0.076, 1.757], [-0.723, 1.078], [-3.289, 1.953], [-0.162, 0.7]],
                             label=False, save_loc=final_dir + 'params_v_scatter.png')

    # Require Dust
    # =================================================================================================================
    # Absolute Mag v. Color
    plotter.abs_mag_v_color(path_91bg=pm_91bg_salt_cut, path_norm=pm_norm_salt_cut, path_dust=pm_dust,
                            save_loc=final_dir+'absMag_v_color.png')
    # Brout+Scolnic 2021 style Dust v. Scatter
    plotter.dust_v_scatter(path_91bg = pm_91bg_salt_cut, path_norm = pm_norm_salt_cut, path_dust = pm_dust,
                           bin_num = 50, bin_bounds = [0.1, 6.3], hist_bins = 60,
                           label = True, save_loc = final_dir+'dust_v_scatter.png')

    # WIP
    # =================================================================================================================
    plotter.dust_mass_step(path_91bg = pm_91bg_merged_cut, path_norm = pm_norm_merged_cut, path_dust = pm_dust,
                           save_loc = final_dir+'dust_mass_step.png', label = False)


    return


if __name__ == '__main__':
    # =================================================================================================================
    ###
    # 1991bg-like SNe Ia Process
    ##
    run_fitter('91bg', 'snpy', False)
    run_fitter('91bg', 'salt', False)
    dataManipulation.combine_snpy_salt(snpy_path='results/91bg_snpy_params.txt',
                                       salt_path='results/91bg_snpy_params.txt',
                                       save_loc='results/91bg_snpy-salt_params.txt')
    criteria_91bg = {
        'z_cmb':       [0.015, 999],
        'z_err':       [-999, 999],
        'Tmax_err':    [-999, 1.0],
        'mu_err':      [-999, 0.1],
        'chisquare':   [-999, 999999],
        'EBVhost':     [-0.3, 0.3],
        'EBVhost_err': [-999, 0.1],
        'st':          [0.00, 1.0],
        'st_err':      [-999, 0.1],
        'c':           [-0.6, 0.6],
        'c_err':       [-999, 0.1],
        'x1':          [-3.0, 3.0],
        'x1_err':      [-999, 0.2]
    }
    dataManipulation.selection_criteria(snpy_path='results/91bg_snpy_params.txt', salt_path='results/91bg_salt_params.txt',
                                        save_loc='results/91bg_snpy-salt_params_cut.txt', criteria=criteria_91bg)
    dataManipulation.selection_criteria(snpy_path='results/91bg_snpy_params.txt',
                                        save_loc='results/91bg_snpy_params_cut.txt', criteria=criteria_91bg)
    dataManipulation.selection_criteria(salt_path='results/91bg_salt_params.txt',
                                        save_loc='results/91bg_salt_params_cut.txt', criteria=criteria_91bg)
    run_lc_checker('91bg')
    # =================================================================================================================
    ###
    # Normal SNe Ia Process
    ###
    run_fitter('norm', 'snpy', False)
    run_fitter('norm', 'salt', False)
    dataManipulation.combine_snpy_salt('results/norm_snpy_params.txt',
                                       'results/norm_snpy_params.txt',
                                       'results/norm_snpy-salt_params.txt')
    criteria_norm = {
        'z_cmb':       [.015, 999], # hold
        'z_err':       [-999, 999],  # hold
        'Tmax_err':    [-999, 1.0], # hold, 2.0
        'mu_err':      [-999, 0.1], # hold
        'chisquare':   [-999, 110000],
        'EBVhost':     [-0.2, 0.2], #
        'EBVhost_err': [-999, 0.1], #
        'st':          [0.75, 1.18],# hold, given by David
        'st_err':      [-999, 0.2], # 0.2
        'c':           [-0.2, 0.3], #
        'c_err':       [-999, 999], # no cut
        'x1':          [-2.0, 2.0], #
        'x1_err':      [-999, 1.0]  # 1.0
    }
    dataManipulation.selection_criteria(snpy_path='results/norm_snpy_params.txt', salt_path='results/norm_salt_params.txt',
                                        save_loc='results/norm_snpy-salt_params_cut.txt', criteria=criteria_norm)
    dataManipulation.selection_criteria(snpy_path='results/norm_snpy_params.txt',
                                        save_loc='results/norm_snpy_params_cut.txt', criteria=criteria_norm)
    dataManipulation.selection_criteria(salt_path='results/norm_salt_params.txt',
                                        save_loc='results/norm_salt_params_cut.txt', criteria=criteria_norm)
    run_lc_checker('norm')
    # =================================================================================================================

    run_plotter()



