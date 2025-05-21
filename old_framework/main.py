# M. D. Woods
# 04/10/2025
import glob
import os
import shutil
import numpy as np
import time as systime

import utils              # Utility functions for various common actions accross project; Here to be initiated for later.
import queryCSP           # Retrieves "Photometric Data release 3" tarball from Krisciunas et al. 2017 hosted on
                          # https://csp.obs.carnegiescience.edu/data, then seperates data using target file.
import queryATLAS         # Retrieves data from ATLAS server and stores it in data\ -- requires 'atlas_key' in 'txts/api_keys.txt'.
import queryZTF           # Retrieves data from ZTF server sends out request to be send to an email.
import fitter             # Uses data in data\ to create SNooPy & SALT fits.
import plotter            # Plots data found in files in results\.
import dataManipulation

# RUN =================================================================================================================
def run_time(func):
  def wrapper():
    start = systime.time()  # Runtime tracker
    func()
    if (systime.time() - start) < 60:
        main_run_time = f"{round((systime.time() - start), 4)} seconds"
    else:
        main_run_time = f"{int((systime.time() - start) / 60)}m {round((systime.time() - start) % 60, 4)}s"
    print(f'|---------------------------|\n Run-time: {main_run_time}\n|---------------------------|')
  return wrapper
def run_queryCSP():
    queryCSP.download(save_loc='data/')  # location of CSP-norm & CSP-/91bg
    queryCSP.separate_data('91bg-like', data_loc='data/DR3/', save_loc='data/CSP-91bg/')  # Separate 1991bg-like out of all CSP-I data
    queryCSP.separate_data('normal', data_loc='data/DR3/', save_loc='data/CSP-norm/')  # Separate Normals out of all CSP-I data
    return
def run_queryATLAS():
    # # Get SN1991bg-like SNe
    # queryATLAS.download(tar_list="txts/target_files/sn1991bglike_tarlist.csv", save_loc='data/ATLAS-91bg/')
    # Get normal SNe
    # queryATLAS.download(tar_list="txts/target_files/old_normal_tarlist.csv", save_loc='data/ATLAS-norm/')
    return
def run_queryZTF():
    # Get SN1991bg-like SNe
    queryZTF.download(tar_list="txts/target_files/sn1991bglike_tarlist.csv", save_loc = 'data/ZTF-91bg/')
    return
def run_file_verification():
    # SN1991bg-like SNe -- Verify where the SNe from the target file are in the directories
    utils.verify_downloads(combined_tarlist_path="txts/target_files/sn1991bglike_tarlist.csv", source='CSP', subtype='91bg')
    utils.verify_downloads(combined_tarlist_path="txts/target_files/sn1991bglike_tarlist.csv", source='ATLAS', subtype='91bg')
    utils.verify_downloads(combined_tarlist_path="txts/target_files/sn1991bglike_tarlist.csv", source='ZTF', subtype='91bg')
    utils.verify_downloads(combined_tarlist_path="txts/target_files/sn1991bglike_tarlist.csv", source='ATLAS+ZTF', subtype='91bg')

    # normal SNe -- Verify where the SNe from the target file are in the directories
    utils.verify_downloads(combined_tarlist_path="txts/target_files/old_normal_tarlist.csv", source='CSP', subtype='norm')
    utils.verify_downloads(combined_tarlist_path="txts/target_files/old_normal_tarlist.csv", source='ATLAS', subtype='norm')
    return
def run_data_combination(clear_old: bool = True):
    dataManipulation.combine_like_data("txts/target_files/sn1991bglike_tarlist.csv",
                                       'data/combined_91bg/',
                                       '91bg', clear_old=clear_old)
    dataManipulation.check_combined_stats('91bg')
    return
def run_algo_combination(subtype: str = '91bg'):
    dataManipulation.combine_snpy_salt(snpy_path=f'results/{subtype}_snpy.txt',
                                       salt_path=f'results/{subtype}_salt.txt',
                                       save_loc=f'results/{subtype}_snpy-salt.txt')
    dataManipulation.selection_criteria(path=f'results/{subtype}_snpy-salt.txt',
                                        subtype=subtype,
                                        save_loc=f'results/{subtype}_snpy-salt_cut.txt')
    dataManipulation.visual_inspection(path=f'results/{subtype}_snpy-salt_cut.txt',
                                       save_loc=f"results/{subtype}_snpy-salt_cut.txt")
    dataManipulation.selection_criteria(path=f'results/{subtype}_snpy.txt',
                                        subtype=subtype,
                                        save_loc=f'results/{subtype}_snpy_cut.txt')
    dataManipulation.selection_criteria(path=f'results/{subtype}_salt.txt',
                                        subtype=subtype,
                                        save_loc=f'results/{subtype}_salt_cut.txt')
    # dataManipulation.seperate_cut_data(mixed_path=f'results/{subtype}_snpy-salt_cut.txt',
    #                                    snpy_path=f'results/{subtype}_snpy.txt',
    #                                    salt_path=f'results/{subtype}_salt.txt',
    #                                    snpy_path_new=f'results/{subtype}_snpy_cut.txt',
    #                                    salt_path_new=f'results/{subtype}_salt_cut.txt')
    # dataManipulation.seperate_cut_data(mixed_path=f'results/{subtype}_snpy-salt_cut.txt',
    #                                    snpy_path=f'results/{subtype}_snpy_cut.txt',
    #                                    salt_path=f'results/{subtype}_salt_cut.txt')
    return
def run_fitter(subtype: str = '91bg', algo: str = 'snpy', rewrite: bool = True):
    """
    :param algo: Algorithm to use for fitting
    :param rewrite: Whether to rewrite the fit data or not
    :return: None
    """
    # Run fitting
    sne_csp = fitter.fit(data_loc=f'data/combined_{subtype}/CSP*.txt', algo=algo, rewrite=rewrite)
    sne_atlas = fitter.fit(data_loc=f'data/combined_{subtype}/ATLAS*.txt', algo=algo, rewrite=rewrite)
    sne_ztf = fitter.fit(data_loc=f'data/combined_{subtype}/ZTF*.txt', algo=algo, rewrite=rewrite)
    sne_combined = fitter.fit(data_loc=f'data/combined_{subtype}/2*.txt', algo=algo, rewrite=rewrite)

    # Combined resultant object lists and make a txt file of parameters
    dataManipulation.make_param_file(sne=sne_csp+sne_atlas+sne_ztf+sne_combined, #
                                     save_loc=f'results/{subtype}_{algo}.txt')

    # # Apply selection criteria to parameter file
    # dataManipulation.selection_criteria(path=f'results/{subtype}_{algo}.txt',
    #                                     subtype=subtype,
    #                                     save_loc=f'results/{subtype}_{algo}_cut.txt')
    # Old logic
    # if subtype == '91bg':
    #     # Run 1991bg-like combined data
    #     sne_csp_91bg = fitter.fit(data_loc='data/combined_91bg/CSP*.txt', algo=algo, rewrite=rewrite)
    #     sne_atlas_91bg = fitter.fit(data_loc='data/combined_91bg/ATLAS*.txt', algo=algo, rewrite=rewrite)
    #     sne_ztf_91bg = fitter.fit(data_loc='data/combined_91bg/ZTF*.txt', algo=algo, rewrite=rewrite)
    #     sne_combined_91bg = fitter.fit(data_loc='data/combined_91bg/2*.txt', algo=algo, rewrite=rewrite)
    #     dataManipulation.make_param_file(sne=sne_csp_91bg+sne_atlas_91bg+sne_ztf_91bg+sne_combined_91bg,
    #                                      save_loc=f'results/SN1991bglikeParameters_{algo}.txt')
    #     dataManipulation.selection_criteria(f'results/SN1991bglikeParameters_{algo}.txt',
    #                                         f'results/SN1991bglikeParameters_{algo}_cut.txt')
    # elif subtype == 'norm':
    #     sne_csp_norm = fitter.fit(data_loc='data/combined_norm/CSP*.txt', algo=algo, rewrite=rewrite)
    #     sne_atlas_norm = fitter.fit(data_loc='data/combined_norm/ATLAS*.txt', algo=algo, rewrite=rewrite)
    #     sne_ztf_norm = fitter.fit(data_loc='data/combined_norm/ZTF*.txt', algo=algo, rewrite=rewrite)
    #     sne_combined_norm = fitter.fit(data_loc='data/combined_norm/2*.txt', algo=algo, rewrite=rewrite)
    #     # for a in ['snpy', 'salt']:
    #     #     norm_atlas_sne = fitter.fit(data_loc='data/ATLAS-norm/*.txt', algo=a, rewrite=True)
    #     #     dataManipulation.make_param_file(sne=norm_atlas_sne,
    #     #                                      save_loc=f'results/NormalParameters_{a}.txt')
    #     #     dataManipulation.selection_criteria(f'results/NormalParameters_{a}.txt',
    #     #                                         f'results/NormalParameters_{a}_cut.txt')
    #     pass
    return
def run_lc_checker(subtype: str = '91bg', clear_old: bool = True):
    if subtype == '91bg':
        path = "results/91bg_snpy-salt_cut.txt"
    elif subtype == 'norm':
        path = "results/norm_snpy-salt_cut.txt"

    if clear_old:
        if os.path.exists(f"results/final_lc/{subtype}/"):
            print(f"[~~~] Replacing files @ results/final_lc/{subtype}/...")
            shutil.rmtree(f"results/final_lc/{subtype}")
            os.mkdir(f"results/final_lc/{subtype}/")

    tb = utils.default_open(path, True)
    for name, algo, source, subtype in zip(tb['objname'], tb['algo'], tb['origin'], tb['subtype']):
        n_o = f"fitting/{algo}-plots/{name}_lc.png"
        n_f = f"results/final_lc/{subtype}/{name}_lc.png"
        shutil.copyfile(n_o, n_f)
        # print(f"[~~~] Copied {name}: '{n_o}' -> '{n_f}'...")
    return
def run_plotter(final_dir: str = 'plots/'):
    # OLDEST PARAMS
    # pm_norm_salt_cut = 'results/old/aaronDo_salt2_params_cut.txt'
    # pm_norm_salt_uncut = 'results/old/aaronDo_salt2_params.txt'
    # pm_norm_snpy_cut = 'results/old/dr3_params.txt'  # Needs to be cut?, no mu
    # pm_norm_snpy_uncut = 'results/old/dr3_params.txt'
    # pm_norm_merged_cut = 'results/old/aaronDo_salt2_params_cut.txt'  # only contains salt fitted
    # pm_norm_merged_uncut = 'results/old/aaronDo_salt2_params.txt'  # only contains salt fitted
    #
    # pm_91bg_salt_cut = 'results/old/combiend__salt_params_cut.txt'
    # pm_91bg_salt_uncut = 'results/old/combiend__salt_params.txt'
    # pm_91bg_snpy_cut = 'results/old/combiend__snpy_params_cut.txt'
    # pm_91bg_snpy_uncut = 'results/old/combiend__snpy_params.txt'
    # pm_91bg_merged_cut = 'results/old/merged_params_cut.txt'
    # pm_91bg_merged_uncut = 'results/old/merged_params.txt'

    # # OLD PARAMS
    # pm_norm_salt_cut = 'results/NormalParameters_salt_cut.txt'
    # pm_norm_salt_uncut = 'results/NormalParameters_salt.txt'
    # pm_norm_snpy_cut = 'results/NormalParameters_snpy_cut.txt'
    # pm_norm_snpy_uncut = 'results/NormalParameters_snpy.txt'
    # pm_norm_merged_cut = 'results/NormalParameters_snpy-salt_cut.txt'
    # pm_norm_merged_uncut = 'results/NormalParameters_snpy-salt.txt'
    #
    # pm_91bg_salt_cut = 'results/SN1991bglikeParameters_salt_cut.txt'
    # pm_91bg_salt_uncut = 'results/SN1991bglikeParameters_salt.txt'
    # pm_91bg_snpy_cut = 'results/SN1991bglikeParameters_snpy_cut.txt'
    # pm_91bg_snpy_uncut = 'results/SN1991bglikeParameters_snpy.txt'
    # pm_91bg_merged_cut = 'results/SN1991bglikeParameters_snpy-salt_cut.txt'
    # pm_91bg_merged_uncut = 'results/SN1991bglikeParameters_snpy-salt.txt'

    # NEW PARAMS
    pm_norm_salt_cut = 'results/norm_salt_cut.txt'
    pm_norm_salt_uncut = 'results/norm_salt.txt'
    pm_norm_snpy_cut = 'results/norm_snpy_cut.txt'
    pm_norm_snpy_uncut = 'results/norm_snpy.txt'
    pm_norm_merged_cut = 'results/norm_snpy-salt_cut.txt'
    pm_norm_merged_uncut = 'results/norm_snpy-salt.txt'

    pm_91bg_salt_cut = 'results/91bg_salt_cut.txt'
    pm_91bg_salt_uncut = 'results/91bg_salt.txt'
    pm_91bg_snpy_cut = 'results/91bg_snpy_cut.txt'
    pm_91bg_snpy_uncut = 'results/91bg_snpy.txt'
    pm_91bg_merged_cut = 'results/91bg_snpy-salt_cut.txt'
    pm_91bg_merged_uncut = 'results/91bg_snpy-salt.txt'

    pm_dust = 'txts/global_dust_params.txt'
    # pm_dust = 'txts/local_dust_params.txt'

    # # WORKING =========================================================================================================
    plotter.resid_v_mass(path_91bg=pm_91bg_merged_cut,
                         path_norm=pm_norm_merged_cut,
                         save_loc=final_dir + 'resid_v_mass.png',
                         label = True)

    # # # Brout+Scolnic 2021 style Dust v. Scatter
    # # plotter.dust_v_scatter(path_91bg = pm_91bg_salt_cut, path_norm = pm_norm_salt_cut, path_dust = pm_dust,
    # #                        bin_num = 50, bin_bounds = [0.1, 6.3], hist_bins = 60,
    # #                        label = True, save_loc = final_dir+'dust_v_scatter.png')
    #
    # # SALT3 Plots
    # plotter.alpha_beta(path_91bg=pm_91bg_salt_cut,
    #                     path_norm=pm_norm_salt_cut,
    #                     save_loc=final_dir+'alpha_beta.png')

    # # Absolute Mag v. Color
    # plotter.abs_mag_v_color(path_91bg=pm_91bg_salt_cut, path_norm=pm_norm_salt_cut, path_dust=pm_dust,
    #                         save_loc=final_dir+'absMag_v_color.png')
    #
    # # # Paramater Histograms
    # plotter.param_hist(snpy_91bg_path=pm_91bg_snpy_cut, snpy_norm_path=pm_norm_snpy_cut,
    #                    salt_91bg_path=pm_91bg_salt_cut, salt_norm_path=pm_norm_salt_cut,
    #                    bins = [0.08, 0.8, 0.05, 0.05], hist_tol = [0.40, 1.50, 0.00, 0.15],
    #                    save_loc=final_dir + 'param_hist_cut.png', line_type='median')
    # plotter.param_hist(snpy_91bg_path=pm_91bg_snpy_uncut, snpy_norm_path=pm_norm_snpy_uncut,
    #                    salt_91bg_path=pm_91bg_salt_uncut, salt_norm_path=pm_norm_salt_uncut,
    #                    bins=[0.1, 0.8, 0.2, 0.15], hist_tol = [0.00, 0.00, 0.00, 0.00],
    #                    save_loc=final_dir + 'param_hist_uncut.png', line_type='median')
    #
    # # Brout+Scolnic 2021 style Paramaters (Color & Stretch) v. Scatter
    # plotter.params_v_scatter(path_snpy_91bg=pm_91bg_snpy_cut, path_snpy_norm=pm_norm_snpy_cut,
    #                          path_salt_91bg=pm_91bg_salt_cut, path_salt_norm=pm_norm_salt_cut,
    #                          bin_nums=[[15, 20], [11, 20], [10, 20], [10, 20]],
    #                          bin_bounds=[[0.076, 1.757], [-0.723, 1.078], [-3.289, 1.953], [-0.162, 0.7]],
    #                          label=False, save_loc=final_dir + 'params_v_scatter.png')

    # # REDUNDANT =======================================================================================================
    # plotter.color_v_scatter(path_snpy_91bg=pm_91bg_snpy_cut, path_snpy_norm=pm_norm_snpy_cut,
    #                         path_salt_91bg=pm_91bg_salt_cut, path_salt_norm=pm_norm_salt_cut,
    #                         bin_nums=[[10, 10], [10, 10]], bin_bounds=[[-0.3, 0.6], [-0.3, 0.6]], label=True,
    #                         save_loc=final_dir + 'color_v_scatter.png')

    # plotter.dust_hist(path_91bg=pm_91bg_salt_cut, path_red_norm=pm_redNorms, path_dust=pm_dust,
    #                    save_loc=final_dir+'dust_params.png')
    # plotter.mu_v_z(path_91bg=pm_91bg_merged_cut,
    #                 path_norm=pm_norm_merged_cut,
    #                 save_loc=final_dir+'mu_v_z.png',
    #                 label = False)
    # BROKEN ==========================================================================================================
    # plotter.resid_v_mass_dust(path_91bg=pm_91bg_merged_cut, path_norm=pm_norm_merged_cut, path_dust=pm_dust,
    #                           save_loc=final_dir + 'dust_resid_v_mass.png')

    return
def run_stats(style: str = 'terminal'):
    if style == 'terminal':
        ### Data Statistics **Before** Fitting Selection Criteria
        ### Selection Criteria
        ### Data Statistics **After** Fitting Selection Criteria

        # Errors ======================================================================================================
        from astropy.table import Table
        report_loc = 'txts/error_report.txt'
        data = np.genfromtxt(report_loc, skip_header=7, dtype=str, delimiter=';')
        report_tbl = Table(names=data[0, :], data=data[1:, :], dtype=[str, list, list, list, list, list, list, list])
        print("-=-=-=- Error Table -=-=-=-\n"
              "99999: TNS Faliure! Needs manual TNS\n"
              "124: Not enough filters to fit!\n"
              "107: SALT3 couldn't fit with current parameter selection!\n"
              "434: SALT3: No data points with S/N > 5.0!\n"
              "222: GHOST failed to intiate!\n"
              "333: GHOST failed to find mass!\n"
              "404: Unknown fitting error!")
        print(report_tbl.colnames)
        for row in report_tbl: print(list(row))
    elif style == 'readme':
        # Errors ======================================================================================================
        from astropy.table import Table
        report_loc = 'txts/error_report.txt'
        data = np.genfromtxt(report_loc, skip_header=7, dtype=str, delimiter=';')
        report_tbl = Table(names=data[0, :], data=data[1:, :], dtype=[str, list, list, list, list, list, list, list])

        # Markdown Table Header
        print("| Source | "
              "99999: TNS Faliure! Needs manual TNS | "
              "124: Not enough filters to fit! | "
              "107: SALT3 couldn't fit with current parameter selection! | "
              "434: SALT3: No data points with S/N > 5.0! | "
              "222: GHOST failed to intiate! | "
              "333: GHOST failed to find mass! | "
              "404: Unknown fitting error! |")

        # Markdown Table Designation
        print("|"+str([":----------------------|"]*8)[1:-1].replace("'", "").replace(", ", ""))

        # Markdown Table Data
        for row in report_tbl:
            line = '| '
            for col in row:
                line += f"{str(col)} | "
            print(line)
    elif style == 'latex':
        # =============================================================================================================
        # Variable call
        # Ex.
        # % Variable Results
        # % SN~1991bg-like SNe Stats
        # \def \msfOurNumSNe {36}
        # \def \msfOurMassStepLogTen {0.014}
        # \def \msfOurMassStepLogTenErr {0.083}
        # \def \msfOurMassStepMedian {0.051}
        # \def \msfOurMassStepMedianErr {0.066}
        #
        # % Normal SNe Stats
        # \def \msfNormNumSNe {256}
        # \def \msfNormMassStepLogTen {-0.076}
        # \def \msfNormMassStepLogTenErr {0.022}
        # \def \msfNormMassStepMedian {-0.084}
        # \def \msfNormMassStepMedianErr {0.027}
        # tb_91bg = utils.default_open('results/91bg_snpy-salt_cut.txt', True)
        # tb_91bg['resid_mu'] = tb_91bg['mu'] - utils.current_cosmo().distmod(tb_91bg['z_cmb']).value
        # tb_91bg['resid_mu_err'] = np.copy(tb_91bg['mu_err'])
        # mass_step_10_91bg, resid_10_91bg = plotter.mass_step_calc(tb_91bg['mu'], tb_91bg['mu_err'],
        #                                                           tb_91bg['resid_mu'], tb_91bg['hostMass'],
        #                                                           tb_91bg['z_cmb'], cut = 10.0)
        # mass_step_med_91bg, resid_med_91bg = plotter.mass_step_calc(tb_91bg['mu'], tb_91bg['mu_err'],
        #                                                             tb_91bg['resid_mu'], tb_91bg['hostMass'],
        #                                                             tb_91bg['z_cmb'], cut=np.median(tb_91bg['hostMass']))
        # print("% Variable Results\n"
        #       "% SN~1991bg-like SNe Stats\n"
        #       "\def \msfOurNumSNe {"+f"{len(tb_91bg)}"+"}\n"
        #       "\def \msfOurMassStepLogTen {"+f"{round(mass_step_10_91bg['value'], 3)}"+"}\n"
        #       "\def \msfOurMassStepLogTenErr {"+f"{round(mass_step_10_91bg['err'], 3)}"+"}\n"
        #       "\def \msfOurMassStepMedian {"+f"{round(mass_step_med_91bg['value'], 3)}"+"}\n"
        #       "\def \msfOurMassStepMedianErr {"+f"{round(mass_step_med_91bg['err'], 3)}"+"}")
        # tb_norm = utils.default_open('results/norm_snpy-salt_cut.txt', True)
        # tb_norm['resid_mu'] = tb_norm['mu'] - utils.current_cosmo().distmod(tb_norm['z_cmb']).value
        # tb_norm['resid_mu_err'] = np.copy(tb_norm['mu_err'])
        # mass_step_10_norm, resid_10_norm = plotter.mass_step_calc(tb_norm['mu'], tb_norm['mu_err'],
        #                                                           tb_norm['resid_mu'], tb_norm['hostMass'],
        #                                                           tb_norm['z_cmb'], cut=10.0)
        # mass_step_med_norm, resid_med_norm = plotter.mass_step_calc(tb_norm['mu'], tb_norm['mu_err'],
        #                                                             tb_norm['resid_mu'], tb_norm['hostMass'],
        #                                                             tb_norm['z_cmb'], cut=np.median(tb_91bg['hostMass']))
        # print("% Normal SNe Stats\n"
        #       "\def \msfNormNumSNe {"+f"{len(tb_norm)}"+"}\n"
        #       "\def \msfNormMassStepLogTen {"+f"{round(mass_step_10_norm['value'], 3)}"+"}\n"
        #       "\def \msfNormMassStepLogTenErr {"+f"{round(mass_step_10_norm['err'], 3)}"+"}\n"
        #       "\def \msfNormMassStepMedian {"+f"{round(mass_step_med_norm['value'], 3)}"+"}\n"
        #       "\def \msfNormMassStepMedianErr {"+f"{round(mass_step_med_norm['err'], 3)}"+"}")

        # =============================================================================================================
        # # Table 1. Summary of SN 1991bg-like SNe Samples
        # # Ex.
        # # \begin{table*}
        # #     \centering
        # #     \caption{Summary of SN~1991bg-like SNe Samples}
        # #     \begin{tabular}{lcrrrrrr}
        # #         \hline
        # #          & Number of SNe & \multicolumn{2}{r}{Redshift Range}& \multicolumn{2}{r}{Declination Range} & Average Mag. & Average Mag. Err\\
        # #          &&min.&max.&min.&max.&$\pm5$~days from peak&$\pm5$~days from peak\\
        # #         \hline
        # #         CSP & $11$ & $0.0039$ & $0.0424$ & $-24.9441$ & $20.5262$ & $17.1799$ & $0.0259$\\
        # #         ATLAS & $10$ & $0.0149$& $0.0594$ & $-44.2625$ & $74.8300$& $17.3974$ & $0.0870$\\
        # #         ZTF & $20$ & $0.0160$ & $0.1002$ & $-22.4477$ & $65.4675$ & $18.8785$ & $0.1230$\\
        # #         ATLAS-ZTF$^{\ast}$ & $61$ & $0.0055$ & $0.0706$ & $-27.2149$ & $86.9325$ & $18.2714$ & $0.1098$\\
        # #         All Surveys &  $102$ & $0.0039$ & $0.1002$ & $-44.2625$ & $86.9325$ & $18.2427$ & $0.1049$\\
        # #         \hline
        # #     \end{tabular}
        # #
        # #     *ATLAS and ZTF have 61 SNe that overlap, so these data contains the \textit{c}, \textit{o}, \textit{g}, \textit{r}, and \textit{i}-bands.
        # #     \label{tab:sample_summary_uncut}
        # # \end{table*}
        tb_uncut = utils.default_open('results/91bg_snpy-salt.txt', True)

        all_data_paths = glob.glob('data/combined_91bg/*.txt')
        dict_base = {'num': 0,
                     'z_min': 999, 'z_max': -999,
                     'dec_min': 999, 'dec_max': -999,
                     'avg_mag': np.array([]), 'avg_mag_err': np.array([])}
        csp, atlas, ztf, combined = dict_base.copy(), dict_base.copy(), dict_base.copy(), dict_base.copy()
        for path in all_data_paths:
            name = path.split('/')[-1][:-4]
            if (name.replace('CSP', '').replace('ATLAS', '').replace('ZTF', '')
                    not in tb_uncut['objname']):
                continue

            if name[0] == "2":
                tb = tb_uncut[tb_uncut['objname'] == name]
                combined['num'] += 1
                atlas['num'] += 1
                ztf['num'] += 1

                # Set redshifts
                if tb['z'] < combined['z_min']: combined['z_min'] = tb['z'][0]
                if tb['z'] > combined['z_max']: combined['z_max'] = tb['z'][0]
                if tb['z'] < atlas['z_min']: atlas['z_min'] = tb['z'][0]
                if tb['z'] > atlas['z_max']: atlas['z_max'] = tb['z'][0]
                if tb['z'] < ztf['z_min']: ztf['z_min'] = tb['z'][0]
                if tb['z'] > ztf['z_max']: ztf['z_max'] = tb['z'][0]

                # Set declinations
                if tb['dec'] < combined['dec_min']: combined['dec_min'] = tb['dec'][0]
                if tb['dec'] > combined['dec_max']: combined['dec_max'] = tb['dec'][0]
                if tb['dec'] < atlas['dec_min']: atlas['dec_min'] = tb['dec'][0]
                if tb['dec'] > atlas['dec_max']: atlas['dec_max'] = tb['dec'][0]
                if tb['dec'] < ztf['dec_min']: ztf['dec_min'] = tb['dec'][0]
                if tb['dec'] > ztf['dec_max']: ztf['dec_max'] = tb['dec'][0]

                # Set mags
                data = np.genfromtxt(path, skip_header=2, dtype=str, delimiter=', ')
                tmax = tb['Tmax'][0]
                mjd = data[:, 5].astype(float)
                mag, dmag = data[:, 6].astype(float)[np.abs(mjd-tmax) < 5], data[:, 7].astype(float)[np.abs(mjd-tmax) < 5]
                combined['avg_mag'] = np.append(combined['avg_mag'], np.average(mag))
                combined['avg_mag_err'] = np.append(combined['avg_mag_err'], np.average(dmag))
                atlas['avg_mag'] = np.append(atlas['avg_mag'], np.average(mag))
                atlas['avg_mag_err'] = np.append(atlas['avg_mag_err'], np.average(dmag))
                ztf['avg_mag'] = np.append(ztf['avg_mag'], np.average(mag))
                ztf['avg_mag_err'] = np.append(ztf['avg_mag_err'], np.average(dmag))


            elif name[:3] == "CSP":
                tb = tb_uncut[tb_uncut['objname'] == name[3:]]
                csp['num'] += 1

                # Set redshifts
                if tb['z'] < csp['z_min']: csp['z_min'] = tb['z'][0]
                if tb['z'] > csp['z_max']: csp['z_max'] = tb['z'][0]

                # Set declinations
                if tb['dec'] < csp['dec_min']: csp['dec_min'] = tb['dec'][0]
                if tb['dec'] > csp['dec_max']: csp['dec_max'] = tb['dec'][0]

                # Set mags
                mjd, mag, dmag = np.array([]), np.array([]), np.array([])
                with open(path, 'r') as f:
                    for line in f.readlines():
                        line_arr = line.split(' ')
                        if len(line_arr) == 3:
                            mjd = np.append(mjd, float(line_arr[0])+53000)
                            mag = np.append(mag, float(line_arr[1]))
                            dmag = np.append(dmag, float(line_arr[2][:-1]))

                tmax = tb['Tmax'][0]
                mag, dmag = mag[np.abs(mjd-tmax) < 5], dmag[np.abs(mjd-tmax) < 5]
                csp['avg_mag'] = np.append(csp['avg_mag'], np.average(mag))
                csp['avg_mag_err'] = np.append(csp['avg_mag_err'], np.average(dmag))

            elif name[:5] == "ATLAS":
                tb = tb_uncut[tb_uncut['objname'] == name[5:]]
                atlas['num'] += 1

                # Set redshifts
                if tb['z'] < atlas['z_min']: atlas['z_min'] = tb['z'][0]
                if tb['z'] > atlas['z_max']: atlas['z_max'] = tb['z'][0]

                # Set declinations
                if tb['dec'] < atlas['dec_min']: atlas['dec_min'] = tb['dec'][0]
                if tb['dec'] > atlas['dec_max']: atlas['dec_max'] = tb['dec'][0]

                # Set mags
                data = np.genfromtxt(path, skip_header=1, dtype=str, delimiter=',')
                tmax = tb['Tmax'][0]
                mjd = data[:, 0].astype(float)
                mag, dmag = data[:, 1].astype(float)[np.abs(mjd-tmax) < 5], data[:, 2].astype(float)[np.abs(mjd-tmax) < 5]
                atlas['avg_mag'] = np.append(atlas['avg_mag'], np.average(mag))
                atlas['avg_mag_err'] = np.append(atlas['avg_mag_err'], np.average(dmag))

            elif name[:3] == "ZTF":
                tb = tb_uncut[tb_uncut['objname'] == name[3:]]
                ztf['num'] += 1

                # Set redshifts
                if tb['z'] < ztf['z_min']: ztf['z_min'] = tb['z'][0]
                if tb['z'] > ztf['z_max']: ztf['z_max'] = tb['z'][0]

                # Set declinations
                if tb['dec'] < ztf['dec_min']: ztf['dec_min'] = tb['dec'][0]
                if tb['dec'] > ztf['dec_max']: ztf['dec_max'] = tb['dec'][0]

                # Set mags
                data = np.genfromtxt(path, skip_header=57, dtype=str, delimiter=' ')
                tmax = tb['Tmax'][0]
                mjd = data[:, 22].astype(float) - 2400000.5
                flux, dflux, zpts = data[:, 24], data[:, 25], data[:, 20]
                mag = (-2.5 * np.log10(flux[flux != 'null'].astype(float))) + zpts[flux != 'null'].astype(float)
                dmag = np.abs(-1.08573620476 * (dflux[dflux != 'null'].astype(float) / flux[flux != 'null'].astype(float)))
                mag, dmag = mag[np.abs(mjd[flux != 'null']-tmax) < 5], dmag[np.abs(mjd[flux != 'null']-tmax) < 5]
                ztf['avg_mag'] = np.append(ztf['avg_mag'], np.average(mag))
                ztf['avg_mag_err'] = np.append(ztf['avg_mag_err'], np.average(dmag))

         # Simpify mag arrays
        csp['avg_mag'] = np.average(csp['avg_mag'][~np.isnan(csp['avg_mag'])])
        csp['avg_mag_err'] = np.average(csp['avg_mag_err'][~np.isnan(csp['avg_mag_err'])])
        atlas['avg_mag'] = np.average(atlas['avg_mag'][~np.isnan(atlas['avg_mag'])])
        atlas['avg_mag_err'] = np.average(atlas['avg_mag_err'][~np.isnan(atlas['avg_mag_err'])])
        ztf['avg_mag'] = np.average(ztf['avg_mag'][~np.isnan(ztf['avg_mag'])])
        ztf['avg_mag_err'] = np.average(ztf['avg_mag_err'][~np.isnan(ztf['avg_mag_err'])])
        combined['avg_mag'] = np.average(combined['avg_mag'][~np.isnan(combined['avg_mag'])])
        combined['avg_mag_err'] = np.average(combined['avg_mag_err'][~np.isnan(combined['avg_mag_err'])])

        # Get all data
        all_sources = dict_base.copy()
        all_sources['num'] = len(tb_uncut)

        all_sources['z_min'] = np.min([csp['z_min'], atlas['z_min'], ztf['z_min'], combined['z_min']])
        all_sources['dec_min'] = np.min([csp['dec_min'], atlas['dec_min'], ztf['dec_min'], combined['dec_min']])
        all_sources['avg_mag'] = np.average([csp['avg_mag'], atlas['avg_mag'], ztf['avg_mag'], combined['avg_mag']])
        all_sources['avg_mag_err'] = np.average([csp['avg_mag_err'], atlas['avg_mag_err'], ztf['avg_mag_err'], combined['avg_mag_err']])

        all_sources['z_max'] = np.max([csp['z_max'], atlas['z_max'], ztf['z_max'], combined['z_max']])
        all_sources['dec_max'] = np.max([csp['dec_max'], atlas['dec_max'], ztf['dec_max'], combined['dec_max']])
        all_sources['avg_mag'] = np.average([csp['avg_mag'], atlas['avg_mag'], ztf['avg_mag'], combined['avg_mag']])
        all_sources['avg_mag_err'] = np.average([csp['avg_mag_err'], atlas['avg_mag_err'], ztf['avg_mag_err'], combined['avg_mag_err']])

        print("\\begin{table*}\n"
              "\t\\centering\n"
              "\t\\caption{Summary of SN~1991bg-like SNe Samples}\n"
              "\t\\begin{tabular}{lcrrrrrr}\n"
              "\t\t\\hline\n"
              "\t\t& Number of SNe & \multicolumn{2}{r}{Redshift Range}& \multicolumn{2}{r}{Declination Range} & Average Mag. & Average Mag. Err\\\\\n"
              "\t\t&&min.&max.&min.&max.&$\pm5$~days from peak&$\pm5$~days from peak\\\\\n"
              "\t\t\\hline\n"
              f"\t\tCSP & ${csp['num']:.0f}$ & ${csp['z_min']:.4f}$ & ${csp['z_max']:.4f}$ & ${csp['dec_min']:.4f}$ & ${csp['dec_max']:.4f}$ & ${csp['avg_mag']:.4f}$ & ${csp['avg_mag_err']:.4f}$\\\\\n"
              f"\t\tATLAS & ${atlas['num']:.0f}$ & ${atlas['z_min']:.4f}$ & ${atlas['z_max']:.4f}$ & ${atlas['dec_min']:.4f}$ & ${atlas['dec_max']:.4f}$ & ${atlas['avg_mag']:.4f}$ & ${atlas['avg_mag_err']:.4f}$\\\\\n"
              f"\t\tZTF & ${ztf['num']:.0f}$ & ${ztf['z_min']:.4f}$ & ${ztf['z_max']:.4f}$ & ${ztf['dec_min']:.4f}$ & ${ztf['dec_max']:.4f}$ & ${ztf['avg_mag']:.4f}$ & ${ztf['avg_mag_err']:.4f}$\\\\\n"
              f"\t\tATLAS-ZTF & ${combined['num']:.0f}$ & ${combined['z_min']:.4f}$ & ${combined['z_max']:.4f}$ & ${combined['dec_min']:.4f}$ & ${combined['dec_max']:.4f}$ & ${combined['avg_mag']:.4f}$ & ${combined['avg_mag_err']:.4f}$\\\\\n"
              f"\t\tAll Surveys & ${all_sources['num']:.0f}$ & ${all_sources['z_min']:.4f}$ & ${all_sources['z_max']:.4f}$ & ${all_sources['dec_min']:.4f}$ & ${all_sources['dec_max']:.4f}$ & ${all_sources['avg_mag']:.4f}$ & ${all_sources['avg_mag_err']:.4f}$\\\\\n"
              "\t\t\\hline\n"
              "\t\\end{tabular}\n"
              "\t\\label{tab:sample_summary}\n"
              "\\end{table*}"
              )

        # =============================================================================================================
        # Table 2. Selection Criteria
        # Ex.
        # \begin{table*}
        #     \centering
        #     \begin{tabular}{lrr}
        #         \hline
        #         Parameter & SNooPy & SALT3\\
        #         \hline
        #         Redshift & $z > 0.015$ & $z > 0.015$\\
        #         Color & $-0.2 < E(B-V)_{\rm host} < 0.3$& $-0.6 < c < 0.6$\\
        #         Color Err.& $\sigma_{E(B-V)_{\rm host}} < 0.1$ & $\sigma_{c} < 0.1 $\\
        #         Stretch & $s_{BV} < 1.0$& $-3.2 < x_1 < 3.2$\\
        #         Stretch Err.&$\sigma_{s_{BV}} < 0.1$&$\sigma_{x_{1}} < 1$\\
        #         Dist. Mod. Err.& $\sigma_{\mu} < 0.2$& $\sigma_{\mu} < 0.2$\\
        #         Time Err.& $\sigma_{t_0} < 1.0$& $\sigma_{t0} < 1.0$\\
        #         \hline
        #     \end{tabular}
        #     \caption{Selection Criteria}
        #     \label{tab:selection_summary}
        # \end{table*}
        tb_uncut = utils.default_open('results/91bg_snpy-salt.txt', True)
        tb_snpy = tb_uncut[tb_uncut['algo'] == 'SNPY']
        tb_salt = tb_uncut[tb_uncut['algo'] == 'SALT']
        criteria = {
            'z_cmb': [0.015, 999],
            'Tmax_err': [-999, 1.0],
            'mu_err': [-999, 0.1],
            'EBVhost': [-0.3, 0.3],
            'EBVhost_err': [-999, 0.1],
            'st': [0.00, 1.0],
            'st_err': [-999, 0.1],
            'c': [-0.6, 0.6],
            'c_err': [-999, 0.1],
            'x1': [-4.0, 4.0],
            'x1_err': [-999, 0.2]
        }
        chuv_lim, chuv_lim_num = dataManipulation.selection_criteria(snpy_path='results/91bg_snpy.txt',
                                                                     salt_path='results/91bg_salt.txt',
                                                                     save_loc='results/91bg_snpy-salt_cut.txt',
                                                                     criteria=criteria)

        print("\\begin{table*}")
        print("\t\\centering")
        print("\t\\begin{tabular}{lrrr}")
        print("\t\t\\hline")
        print("\t\tParameter & SNooPy & SALT3 & Number of SNe Removed\\\\")
        print("\t\t\\hline")

        # Redshift Cut
        num_rm_z = (len(tb_snpy[tb_snpy['z_cmb'] < criteria['z_cmb'][0]]) + len(tb_snpy[tb_snpy['z_cmb'] > criteria['z_cmb'][1]]))
        print(f"\t\tRedshift & "
              f"$z > {criteria['z_cmb'][0]}$ & "
              f"$z > {criteria['z_cmb'][0]}$ & "
              f"{num_rm_z}\\\\")

        # Color Cut
        cl_lb = ["E(B-V)_{\\rm host}", "c"]
        num_rm_cl = (len(tb_snpy[tb_snpy['color'] < criteria['EBVhost'][0]]) + len(tb_snpy[tb_snpy['color'] > criteria['EBVhost'][1]]) +
                     len(tb_salt[tb_salt['color'] < criteria['c'][0]]) + len(tb_salt[tb_salt['color'] > criteria['c'][1]]))
        print(f"\t\tColor & "
              f"${criteria['EBVhost'][0]} < {cl_lb[0]} < {criteria['EBVhost'][1]}$ & "
              f"${criteria['c'][0]} < {cl_lb[1]} < {criteria['c'][1]}$ & "
              f"{num_rm_cl}\\\\")

        # Color Error Cut
        dcl_lb = ["\sigma_{E(B-V)_{\\rm host}}", "\sigma_{c}"]
        num_rm_dcl = (len(tb_snpy[tb_snpy['color_err'] > criteria['EBVhost_err'][1]]) +
                     len(tb_salt[tb_salt['color_err'] > criteria['c_err'][1]]))
        print(f"\t\tColor Err. & "
              f"${dcl_lb[0]} < {criteria['EBVhost_err'][1]}$ & "
              f"${dcl_lb[1]} < {criteria['c_err'][1]}$ & "
              f"{num_rm_dcl}\\\\")

        # Stretch Cut
        st_lb = ["s_{BV}", "x_1"]
        num_rm_st = (len(tb_snpy[tb_snpy['stretch'] > criteria['st'][1]]) +
                     len(tb_salt[tb_salt['stretch'] < criteria['x1'][0]]) + len(tb_salt[tb_salt['stretch'] > criteria['x1'][1]]))
        print(f"\t\tStretch & "
              f"${st_lb[0]} < {criteria['st'][1]}$ & "
              f"${criteria['x1'][0]} < {st_lb[1]} < {criteria['x1'][1]}$ & "
              f"{num_rm_st}\\\\")

        # Stretch Error Cut
        dst_lb = ["\sigma_{s_{BV}}", "\sigma_{x_{1}}"]
        num_rm_dst = (len(tb_snpy[tb_snpy['stretch_err'] > criteria['st_err'][1]]) +
                      len(tb_salt[tb_salt['stretch_err'] > criteria['x1_err'][1]]))
        print(f"\t\tStretch Err. & "
              f"${dst_lb[0]} < {criteria['st_err'][1]}$ & "
              f"${dst_lb[1]} < {criteria['x1_err'][1]}$ & "
              f"{num_rm_dst}\\\\")

        # Distance Modulus Error Cut
        dmu_lb = ["\sigma_{\mu}", "\sigma_{\mu}"]
        num_rm_dmu = (len(tb_snpy[tb_snpy['mu_err'] > criteria['mu_err'][1]]) +
                      len(tb_salt[tb_salt['mu_err'] > criteria['mu_err'][1]]))
        print(f"\t\tDist. Mod. Err. & "
              f"${dmu_lb[0]} < {criteria['mu_err'][1]}$ & "
              f"${dmu_lb[1]} < {criteria['mu_err'][1]}$ & "
              f"{num_rm_dmu}\\\\")

        # Max Time Error Cut
        dmu_lb = ["\sigma_{t_0}", "\sigma_{t0}"]
        num_rm_dmu = (len(tb_snpy[tb_snpy['Tmax_err'] > criteria['Tmax_err'][1]]) +
                      len(tb_salt[tb_salt['Tmax_err'] > criteria['Tmax_err'][1]]))
        print(f"\t\tTime Err. & "
              f"${dmu_lb[0]} < {criteria['Tmax_err'][1]}$ & "
              f"${dmu_lb[1]} < {criteria['Tmax_err'][1]}$ & "
              f"{num_rm_dmu}\\\\")

        print("\t\t\\hline")

        # Chauvenet Limit
        dmu_lb = ["\\delta\\mu", "\\delta\\mu"]
        print(f"\t\tChauvenet's Criterion & "
              f"${dmu_lb[0]} < {round(chuv_lim, 4)}$ & "
              f"${dmu_lb[1]} < {round(chuv_lim, 4)}$ & "
              f"{round(chuv_lim_num, 4)}\\\\")


        print("\t\t\\hline")
        print("\t\\end{tabular}")
        print("\t\\caption{Selection Criteria}")
        print("\t\\label{tab:selection_summary}")
        print("\\end{table*}")

        # =============================================================================================================
        # Table 3. Summary of SN 1991bg-like SNe Samples After Selection Criteria
        tb_cut = utils.default_open('results/91bg_snpy-salt_cut.txt', True)

        all_data_paths = glob.glob('data/combined_91bg/*.txt')
        dict_base = {'num': 0,
                     'z_min': 999, 'z_max': -999,
                     'dec_min': 999, 'dec_max': -999,
                     'avg_mag': np.array([]), 'avg_mag_err': np.array([])}
        csp, atlas, ztf, combined = dict_base.copy(), dict_base.copy(), dict_base.copy(), dict_base.copy()
        for path in all_data_paths:
            name = path.split('/')[-1][:-4]
            if (name.replace('CSP', '').replace('ATLAS', '').replace('ZTF', '')
                    not in tb_cut['objname']):
                continue

            if name[0] == "2":
                tb = tb_cut[tb_cut['objname'] == name]
                combined['num'] += 1
                atlas['num'] += 1
                ztf['num'] += 1

                # Set redshifts
                if tb['z'] < combined['z_min']: combined['z_min'] = tb['z'][0]
                if tb['z'] > combined['z_max']: combined['z_max'] = tb['z'][0]
                if tb['z'] < atlas['z_min']: atlas['z_min'] = tb['z'][0]
                if tb['z'] > atlas['z_max']: atlas['z_max'] = tb['z'][0]
                if tb['z'] < ztf['z_min']: ztf['z_min'] = tb['z'][0]
                if tb['z'] > ztf['z_max']: ztf['z_max'] = tb['z'][0]

                # Set declinations
                if tb['dec'] < combined['dec_min']: combined['dec_min'] = tb['dec'][0]
                if tb['dec'] > combined['dec_max']: combined['dec_max'] = tb['dec'][0]
                if tb['dec'] < atlas['dec_min']: atlas['dec_min'] = tb['dec'][0]
                if tb['dec'] > atlas['dec_max']: atlas['dec_max'] = tb['dec'][0]
                if tb['dec'] < ztf['dec_min']: ztf['dec_min'] = tb['dec'][0]
                if tb['dec'] > ztf['dec_max']: ztf['dec_max'] = tb['dec'][0]

                # Set mags
                data = np.genfromtxt(path, skip_header=2, dtype=str, delimiter=', ')
                tmax = tb['Tmax'][0]
                mjd = data[:, 5].astype(float)
                mag, dmag = data[:, 6].astype(float)[np.abs(mjd - tmax) < 5], data[:, 7].astype(float)[np.abs(mjd - tmax) < 5]
                combined['avg_mag'] = np.append(combined['avg_mag'], np.average(mag))
                combined['avg_mag_err'] = np.append(combined['avg_mag_err'], np.average(dmag))
                atlas['avg_mag'] = np.append(atlas['avg_mag'], np.average(mag))
                atlas['avg_mag_err'] = np.append(atlas['avg_mag_err'], np.average(dmag))
                ztf['avg_mag'] = np.append(ztf['avg_mag'], np.average(mag))
                ztf['avg_mag_err'] = np.append(ztf['avg_mag_err'], np.average(dmag))


            elif name[:3] == "CSP":
                tb = tb_cut[tb_cut['objname'] == name[3:]]
                csp['num'] += 1

                # Set redshifts
                if tb['z'] < csp['z_min']: csp['z_min'] = tb['z'][0]
                if tb['z'] > csp['z_max']: csp['z_max'] = tb['z'][0]

                # Set declinations
                if tb['dec'] < csp['dec_min']: csp['dec_min'] = tb['dec'][0]
                if tb['dec'] > csp['dec_max']: csp['dec_max'] = tb['dec'][0]

                # Set mags
                mjd, mag, dmag = np.array([]), np.array([]), np.array([])
                with open(path, 'r') as f:
                    for line in f.readlines():
                        line_arr = line.split(' ')
                        if len(line_arr) == 3:
                            mjd = np.append(mjd, float(line_arr[0]) + 53000)
                            mag = np.append(mag, float(line_arr[1]))
                            dmag = np.append(dmag, float(line_arr[2][:-1]))

                tmax = tb['Tmax'][0]
                mag, dmag = mag[np.abs(mjd - tmax) < 5], dmag[np.abs(mjd - tmax) < 5]
                csp['avg_mag'] = np.append(csp['avg_mag'], np.average(mag))
                csp['avg_mag_err'] = np.append(csp['avg_mag_err'], np.average(dmag))

            elif name[:5] == "ATLAS":
                tb = tb_cut[tb_cut['objname'] == name[5:]]
                atlas['num'] += 1

                # Set redshifts
                if tb['z'] < atlas['z_min']: atlas['z_min'] = tb['z'][0]
                if tb['z'] > atlas['z_max']: atlas['z_max'] = tb['z'][0]

                # Set declinations
                if tb['dec'] < atlas['dec_min']: atlas['dec_min'] = tb['dec'][0]
                if tb['dec'] > atlas['dec_max']: atlas['dec_max'] = tb['dec'][0]

                # Set mags
                data = np.genfromtxt(path, skip_header=1, dtype=str, delimiter=',')
                tmax = tb['Tmax'][0]
                mjd = data[:, 0].astype(float)
                mag, dmag = data[:, 1].astype(float)[np.abs(mjd - tmax) < 5], data[:, 2].astype(float)[
                    np.abs(mjd - tmax) < 5]
                atlas['avg_mag'] = np.append(atlas['avg_mag'], np.average(mag))
                atlas['avg_mag_err'] = np.append(atlas['avg_mag_err'], np.average(dmag))

            elif name[:3] == "ZTF":
                tb = tb_cut[tb_cut['objname'] == name[3:]]
                ztf['num'] += 1

                # Set redshifts
                if tb['z'] < ztf['z_min']: ztf['z_min'] = tb['z'][0]
                if tb['z'] > ztf['z_max']: ztf['z_max'] = tb['z'][0]

                # Set declinations
                if tb['dec'] < ztf['dec_min']: ztf['dec_min'] = tb['dec'][0]
                if tb['dec'] > ztf['dec_max']: ztf['dec_max'] = tb['dec'][0]

                # Set mags
                data = np.genfromtxt(path, skip_header=57, dtype=str, delimiter=' ')
                tmax = tb['Tmax'][0]
                mjd = data[:, 22].astype(float) - 2400000.5
                flux, dflux, zpts = data[:, 24], data[:, 25], data[:, 20]
                mag = (-2.5 * np.log10(flux[flux != 'null'].astype(float))) + zpts[flux != 'null'].astype(float)
                dmag = np.abs(
                    -1.08573620476 * (dflux[dflux != 'null'].astype(float) / flux[flux != 'null'].astype(float)))
                mag, dmag = mag[np.abs(mjd[flux != 'null'] - tmax) < 5], dmag[np.abs(mjd[flux != 'null'] - tmax) < 5]
                ztf['avg_mag'] = np.append(ztf['avg_mag'], np.average(mag))
                ztf['avg_mag_err'] = np.append(ztf['avg_mag_err'], np.average(dmag))

        # Simpify mag arrays
        csp['avg_mag'] = np.average(csp['avg_mag'][~np.isnan(csp['avg_mag'])])
        csp['avg_mag_err'] = np.average(csp['avg_mag_err'][~np.isnan(csp['avg_mag_err'])])
        atlas['avg_mag'] = np.average(atlas['avg_mag'][~np.isnan(atlas['avg_mag'])])
        atlas['avg_mag_err'] = np.average(atlas['avg_mag_err'][~np.isnan(atlas['avg_mag_err'])])
        ztf['avg_mag'] = np.average(ztf['avg_mag'][~np.isnan(ztf['avg_mag'])])
        ztf['avg_mag_err'] = np.average(ztf['avg_mag_err'][~np.isnan(ztf['avg_mag_err'])])
        combined['avg_mag'] = np.average(combined['avg_mag'][~np.isnan(combined['avg_mag'])])
        combined['avg_mag_err'] = np.average(combined['avg_mag_err'][~np.isnan(combined['avg_mag_err'])])

        # Get all data
        all_sources = dict_base.copy()
        all_sources['num'] = len(tb_cut)

        all_sources['z_min'] = np.min([csp['z_min'], atlas['z_min'], ztf['z_min'], combined['z_min']])
        all_sources['dec_min'] = np.min([csp['dec_min'], atlas['dec_min'], ztf['dec_min'], combined['dec_min']])
        all_sources['avg_mag'] = np.average([csp['avg_mag'], atlas['avg_mag'], ztf['avg_mag'], combined['avg_mag']])
        all_sources['avg_mag_err'] = np.average(
            [csp['avg_mag_err'], atlas['avg_mag_err'], ztf['avg_mag_err'], combined['avg_mag_err']])

        all_sources['z_max'] = np.max([csp['z_max'], atlas['z_max'], ztf['z_max'], combined['z_max']])
        all_sources['dec_max'] = np.max([csp['dec_max'], atlas['dec_max'], ztf['dec_max'], combined['dec_max']])
        all_sources['avg_mag'] = np.average([csp['avg_mag'], atlas['avg_mag'], ztf['avg_mag'], combined['avg_mag']])
        all_sources['avg_mag_err'] = np.average(
            [csp['avg_mag_err'], atlas['avg_mag_err'], ztf['avg_mag_err'], combined['avg_mag_err']])

        print("\\begin{table*}\n"
              "\t\\centering\n"
              "\t\\caption{Summary of Selected SN~1991bg-like SNe Samples}\n"
              "\t\\begin{tabular}{lcrrrrrr}\n"
              "\t\t\\hline\n"
              "\t\t& Number of SNe & \multicolumn{2}{r}{Redshift Range}& \multicolumn{2}{r}{Declination Range} & Average Mag. & Average Mag. Err\\\\\n"
              "\t\t&&min.&max.&min.&max.&$\pm5$~days from peak&$\pm5$~days from peak\\\\\n"
              "\t\t\\hline\n"
              f"\t\tCSP & ${csp['num']:.0f}$ & ${csp['z_min']:.4f}$ & ${csp['z_max']:.4f}$ & ${csp['dec_min']:.4f}$ & ${csp['dec_max']:.4f}$ & ${csp['avg_mag']:.4f}$ & ${csp['avg_mag_err']:.4f}$\\\\\n"
              f"\t\tATLAS-ZTF & ${combined['num']:.0f}$ & ${combined['z_min']:.4f}$ & ${combined['z_max']:.4f}$ & ${combined['dec_min']:.4f}$ & ${combined['dec_max']:.4f}$ & ${combined['avg_mag']:.4f}$ & ${combined['avg_mag_err']:.4f}$\\\\\n"
              f"\t\tAll Surveys & ${all_sources['num']:.0f}$ & ${all_sources['z_min']:.4f}$ & ${all_sources['z_max']:.4f}$ & ${all_sources['dec_min']:.4f}$ & ${all_sources['dec_max']:.4f}$ & ${all_sources['avg_mag']:.4f}$ & ${all_sources['avg_mag_err']:.4f}$\\\\\n"
              "\t\t\\hline\n"
              "\t\\end{tabular}\n"
              "\t\\label{tab:sample_summary_cut}\n"
              "\\end{table*}"
              )
        pass
    else:
        print(f"[!!!] '{style}' is not a valid style ['terminal'/'readme'/'latex']")
        return
    return
# ANALYZE =============================================================================================================
def compare_old_new_params(new_path: str = 'results/SN1991bglikeParameters_snpy-salt_cut.txt',
                           old_path: str = 'results/old/merged_params_cut.txt'):
    import numpy as np
    tb_new = utils.default_open(new_path, True)
    tb_old = utils.default_open(old_path, True)

    i = 0
    missing = []
    new_better, old_better = 0, 0
    for n in tb_old['objname']:
        if n in list(tb_new['objname']):
            i += 1
            print(i)
            print(f"Old --- "
                  f"{tb_old[tb_old['objname'] == n]['objname'][0]}: "
                  f"mu = {round(tb_old[tb_old['objname'] == n]['mu'][0], 3)} "
                  f"+/- = {round(tb_old[tb_old['objname'] == n]['mu_err'][0], 3)} "
                  f"({tb_old[tb_old['objname'] == n]['algo'][0]})")
            print(f"New --- "
                  f"{tb_new[tb_new['objname'] == n]['objname'][0]}: "
                  f"mu = {round(tb_new[tb_new['objname'] == n]['mu'][0], 3)} "
                  f"+/- = {round(tb_new[tb_new['objname'] == n]['mu_err'][0], 3)} "
                  f"({tb_new[tb_new['objname'] == n]['algo'][0]})")
            if tb_new[tb_new['objname'] == n]['mu_err'][0] > tb_old[tb_old['objname'] == n]['mu_err'][0]:
                print('Better Error: -Old-')
                old_better += 1
            else:
                print('Better Error: +New+')
                new_better += 1
            print('================================================================================================')
            # break
        else:
            missing.append(n)
    print(f"{len(tb_new)} / {len(tb_old)}")
    print(f"New v. Old: {new_better} v. {old_better}")
    print(f"Missing: {missing}")
    return
def brout_scholnic_expected_RMS():
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.table import Table

    # Get average SALT color
    tb_91bg = utils.default_open('results/91bg_snpy-salt_cut.txt', True, delimiter=', ')
    sn91bg_avg_c = np.average(tb_91bg[tb_91bg['algo'] == 'SALT']['color'])

    data = np.genfromtxt("txts/b21_fig5_data.csv", delimiter=',', dtype=str, skip_header=1)
    tb = Table(names=data[0, :], data=data[1:, :], dtype=[float, float])

    coefficients = np.polyfit(tb['c'], tb['sigma'], 2)
    a, b, c = coefficients
    quadratic_function = lambda x: a * x ** 2 + b * x + c

    x_fit = np.linspace(min(tb['c']), sn91bg_avg_c, 100)
    y_fit = quadratic_function(x_fit)

    plt.scatter(tb['c'], tb['sigma'], label='Original Data')
    plt.plot(x_fit, y_fit, label=f'Quadratic Fit: y = {a:.2f}x^2 + {b:.2f}x + {c:.2f}')
    plt.xlabel('c')
    plt.ylabel('sigma')
    plt.legend()
    plt.show()

    print(f"Expected RMS for c = {round(sn91bg_avg_c, 4)}: "
          f"sigma = {round(y_fit[np.where(x_fit == sn91bg_avg_c)[0][0]], 4)}")
    return
def check_SN_in_ATLAS_ZTF(objname: str):
    import requests

    # URL of the TNS object
    url = f"https://www.wis-tns.org/object/{objname}"

    # Set User-Agent to mimic a browser (TNS blocks bots)
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/119.0.0.0 Safari/537.36"}

    # Send GET request
    response = requests.get(url, headers=headers)

    # Check response status
    if response.status_code == 200:
        html_dump = str(response.text)
        if ("ATLAS" in html_dump) and ("ZTF" in html_dump):
            print(f"[+++] CANIDATE [{objname}]")
        elif "ATLAS" in html_dump:
            print(f"[~~~] ATLAS ONLY [{objname}]")
        elif "ZTF" in html_dump:
            print(f"[~~~] ZTF ONLY [{objname}]")
        else:
            print(f"[!!!] NEITHER [{objname}]")
    elif response.status_code == 429:
        raise ValueError("Exsiting...")
    else:
        print(f"[!!!] Failed to fetch data. Status Code: {response.status_code}")
    return
# RUNTIME =============================================================================================================
@run_time
def maunal_fitting():
    # Indivisually fititng desired sne
    ## CSP
    csp_names = "2008bt SN 2008bd SN 2007N SN 2007ba SN 2006gt SN 2006bd SN 2005bl".split(" SN ")
    sne_csp, failed_csp = [], []
    for i, n in enumerate(csp_names):
        print(f"[{i+1} / {len(csp_names)}] ========================================================================")
        sn_snpy = fitter.fit(data_loc=f'data/CSP-91bg/CSP{n}.txt', algo='snpy', rewrite=True)[0]
        sn_salt = fitter.fit(data_loc=f'data/CSP-91bg/CSP{n}.txt', algo='salt', rewrite=True)[0]
        for sn in [sn_snpy, sn_salt]:
            if ("mu" not in sn.params) or ("hostMass" not in sn.params):
                sn = None
            elif ("mu" in sn.params) and (sn.params['mu']['value'] < 0):
                sn = None
            elif ("hostMass" in sn.params) and (sn.params['hostMass']['value'] < 0):
                sn = None

        if (sn_snpy is None) and (sn_salt is None):
            failed_csp.append(n)
        elif (sn_salt is None):
            sne_csp.append(sn_snpy)
        else:
            sne_csp.append(sn_salt)
    print(f"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n"
          f"[+++] ({len(sne_csp) - len(failed_csp)} / {len(csp_names)})\n"
          f"      {failed_csp}\n"
          f"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n")

    ## ATLAS
    atlas_names = ("2024zaj SN 2023jah SN 2022ywc SN 2022yv SN 2022vse SN 2021mab SN 2024jhk").split(" SN ")
    sne_atlas, failed_atlas = [], []
    for i, n in enumerate(atlas_names):
        print(f"[{i+1} / {len(atlas_names)}] ========================================================================")
        sn_snpy = fitter.fit(data_loc=f'data/ATLAS-91bg/ATLAS{n}.txt', algo='snpy', rewrite=True)[0]
        sn_salt = fitter.fit(data_loc=f'data/ATLAS-91bg/ATLAS{n}.txt', algo='salt', rewrite=True)[0]
        for sn in [sn_snpy, sn_salt]:
            if ("mu" not in sn.params) or ("hostMass" not in sn.params):
                sn = None
            elif ("mu" in sn.params) and (sn.params['mu']['value'] < 0):
                sn = None
            elif ("hostMass" in sn.params) and (sn.params['hostMass']['value'] < 0):
                sn = None

        if (sn_snpy is None) and (sn_salt is None):
            failed_atlas.append(n)
        elif (sn_salt is None):
            sne_atlas.append(sn_snpy)
        else:
            sne_atlas.append(sn_salt)
    print(f"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n"
          f"[+++] ({len(sne_atlas) - len(failed_atlas)} / {len(atlas_names)})\n"
          f"      {failed_atlas}\n"
          f"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n")

    ## ZTF
    ztf_names = ["2021afur"]
    # 2024jhk file is in wrong format
    # 2019ecx failed SNPY & SALT
    sne_ztf, failed_ztf = [], []
    for i, n in enumerate(ztf_names):
        print(f"[{i+1} / {len(ztf_names)}] ========================================================================")
        sn_snpy = fitter.fit(data_loc=f'data/ZTF-91bg/ZTF{n}.txt', algo='snpy', rewrite=True)[0]
        if ("mu" not in sn_snpy.params) or ("hostMass" not in sn_snpy.params):
            sn_snpy = None
        elif ("mu" in sn_snpy.params) and (sn_snpy.params['mu']['value'] < 0):
            sn_snpy = None
        elif ("hostMass" in sn_snpy.params) and (sn_snpy.params['hostMass']['value'] < 0):
            sn_snpy = None

        sn_salt = fitter.fit(data_loc=f'data/ZTF-91bg/ZTF{n}.txt', algo='salt', rewrite=True)[0]
        if ("mu" not in sn_salt.params) or ("hostMass" not in sn_salt.params):
            sn_salt = None
        elif ("mu" in sn_salt.params) and (sn_salt.params['mu']['value'] < 0):
            sn_salt = None
        elif ("hostMass" in sn_salt.params) and (sn_salt.params['hostMass']['value'] < 0):
            sn_salt = None

        if (sn_snpy is None) and (sn_salt is None):
            failed_ztf.append(n)
        elif (sn_salt is None):
            sne_ztf.append(sn_snpy)
        elif (sn_snpy is None):
            sne_ztf.append(sn_salt)
        else:
            if sn_snpy.params['mu']['err'] > sn_salt.params['mu']['err']:
                sne_ztf.append(sn_salt)
            else:
                sne_ztf.append(sn_snpy)
        # print(f"\n{type(sn_snpy)} | {type(sn_salt)}\n")
    print(f"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n"
          f"[+++] ({len(sne_ztf) - len(failed_ztf)} / {len(ztf_names)})\n"
          f"      {failed_ztf}\n"
          f"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n")
    # print(len(sne_ztf))
    # print(len(failed_ztf))

    ## COMBINED
    combined_names = ("2024bjb SN 2023fwb SN 2023fot SN 2023cvq SN 2022zvu SN 2022ubt SN 2022rjs "
                      "SN 2022rey SN 2022dsu SN 2022aaok SN 2021wzb SN 2021uve SN 2021uef SN 2021qqr SN 2021jvp "
                      "SN 2021fnr SN 2021bls SN 2021abzd SN 2020vae SN 2020nta SN 2020acbx SN 2020abpe SN 2019op "
                      "SN 2019moq SN 2019luv SN 2019ecx SN 2019cp SN 2019cdc SN 2019be SN 2018lph SN 2018jag "
                      "SN 2018eyi").split(" SN ")
    sne_combined, failed_combined = [], []
    for i, n in enumerate(combined_names):
        print(f"[{i + 1} / {len(combined_names)}] ========================================================================")
        sn_snpy = fitter.fit(data_loc=f'data/combined_91bg/{n}.txt', algo='snpy', rewrite=True)[0]
        if ("mu" not in sn_snpy.params) or ("hostMass" not in sn_snpy.params):
            sn_snpy = None
        elif ("mu" in sn_snpy.params) and (sn_snpy.params['mu']['value'] < 0):
            sn_snpy = None
        elif ("hostMass" in sn_snpy.params) and (sn_snpy.params['hostMass']['value'] < 0):
            sn_snpy = None

        sn_salt = fitter.fit(data_loc=f'data/combined_91bg/{n}.txt', algo='salt', rewrite=True)[0]
        if ("mu" not in sn_salt.params) or ("hostMass" not in sn_salt.params):
            sn_salt = None
        elif ("mu" in sn_salt.params) and (sn_salt.params['mu']['value'] < 0):
            sn_salt = None
        elif ("hostMass" in sn_salt.params) and (sn_salt.params['hostMass']['value'] < 0):
            sn_salt = None

        if (sn_snpy is None) and (sn_salt is None):
            failed_combined.append(n)
        elif (sn_salt is None):
            sne_combined.append(sn_snpy)
        elif (sn_snpy is None):
            sne_combined.append(sn_salt)
        else:
            if sn_snpy.params['mu']['err'] > sn_salt.params['mu']['err']:
                sne_combined.append(sn_salt)
            else:
                sne_combined.append(sn_snpy)
        print(f"\n{type(sn_snpy)} | {type(sn_salt)}\n")
    print(f"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n"
          f"[+++] ({len(sne_combined) - len(failed_combined)} / {len(combined_names)})\n"
          f"      {failed_combined}\n"
          f"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n")

    # Combined resultant object lists and make a txt file of parameters
    dataManipulation.make_param_file(sne=sne_csp+sne_atlas+sne_ztf+sne_ztf+sne_combined,
                                     save_loc=f'manual_fit_test.txt')

    # Apply selection criteria to parameter file
    dataManipulation.selection_criteria(path=f'manual_fit_test.txt',
                                        subtype='91bg',
                                        save_loc=f'manual_fit_test_cut.txt')

    return
@run_time
def main():
    # File acquisition
    # run_queryCSP()  # Get CSP data from CSP site
    # run_queryATLAS()  # Query ATLAS using target list of known RA & DEC
    # run_queryZTF()  # Query ZTF using target list of known RA & DEC
    # queryZTF.download_email(tar_list="txts/target_files/sn1991bglike_tarlist.csv", save_loc = 'data/ZTF-91bg/')

    # File verification
    # utils.get_twomass()  # Checks for twomass velocities data (large file)
    # run_file_verification()  # Using target file, checks if SN exists in the proper directory

    # Data combination
    # run_data_combination()  # Combines data file if SN exsists in multiple surveys

    # Fit light curves
    # run_fitter(subtype = '91bg', algo = 'snpy', rewrite = True)
    # run_fitter(subtype = '91bg', algo = 'salt', rewrite = True)
    # run_fitter(subtype = 'norm', algo = 'snpy', rewrite = False)
    # run_fitter(subtype = 'norm', algo = 'salt', rewrite = False)

    # Selection Criteria
    # criteria_91bg = {
    #     'z_cmb':       [0.015, 999],
    #     'Tmax_err':    [-999, 1.0],
    #     'mu_err':      [-999, 0.1],
    #     'EBVhost':     [-1.0, 1.0],
    #     'EBVhost_err': [-999, 0.1],
    #     'st':          [-999, 999],
    #     'st_err':      [-999, 999],
    #     'c':           [-999, 999],
    #     'c_err':       [-999, 999],
    #     'x1':          [-999, 999],
    #     'x1_err':      [-999, 999]
    # }
    criteria_91bg = {
        'z_cmb':       [0.015, 999],
        'Tmax_err':    [-999, 1.0],
        'mu_err':      [-999, 0.1],
        'EBVhost':     [-0.3, 0.3],
        'EBVhost_err': [-999, 0.1],
        'st':          [0.00, 1.0],
        'st_err':      [-999, 0.1],
        'c':           [-0.6, 0.6],
        'c_err':       [-999, 0.1],
        'x1':          [-3.0, 3.0],
        'x1_err':      [-999, 0.2]
    }
    # dataManipulation.combine_snpy_salt(snpy_path='results/91bg_snpy.txt', salt_path='results/91bg_salt.txt',
    #                                    save_loc='results/91bg_snpy-salt.txt')
    dataManipulation.selection_criteria(snpy_path='results/91bg_snpy.txt', salt_path='results/91bg_salt.txt',
                                        save_loc='results/91bg_snpy-salt_cut.txt', criteria=criteria_91bg)
    # dataManipulation.selection_criteria(snpy_path='results/91bg_snpy.txt',
    #                                     save_loc='results/91bg_snpy-salt_cut.txt', criteria=criteria_91bg)
    # dataManipulation.selection_criteria(salt_path='results/91bg_salt.txt',
    #                                     save_loc='results/91bg_snpy-salt_cut.txt', criteria=criteria_91bg)

    criteria_norm = {
        'z_cmb':       [.015, 999], # hold
        'Tmax_err':    [-999, 2.0], # hold, 2.0
        'mu_err':      [-999, 0.1], # hold
        'EBVhost':     [-0.1, 0.2], #
        'EBVhost_err': [-999, 0.1], #
        'st':          [0.75, 1.18],# hold
        'st_err':      [-999, 0.2], # 0.2
        'c':           [-1.0, 1.0], #
        'c_err':       [-999, 999], # no cut
        'x1':          [-3.0, 3.0], #
        'x1_err':      [-999, 1.0]  # 1.0
    }
    # dataManipulation.combine_snpy_salt(snpy_path='results/norm_snpy.txt', salt_path='results/norm_salt.txt',
    #                                    save_loc='results/norm_snpy-salt.txt')
    # dataManipulation.selection_criteria(snpy_path='results/norm_snpy.txt', salt_path='results/norm_salt.txt',
    #                                     save_loc='results/norm_snpy-salt_cut.txt', criteria=criteria_norm)
    # dataManipulation.selection_criteria(snpy_path='results/norm_snpy.txt',
    #                                     save_loc='results/norm_snpy-salt_cut.txt', criteria=criteria_norm)
    # dataManipulation.selection_criteria(salt_path='results/norm_salt.txt',
    #                                     save_loc='results/norm_snpy-salt_cut.txt', criteria=criteria_norm)

    # Check lc of cut file
    # run_lc_checker('91bg')
    # run_lc_checker('norm')

    # Plot data
    # run_plotter()

    # Outout stats
    # run_stats('latex')
    return
@run_time
def experiment():

    # tb = utils.default_open('results/91bg_snpy-salt.txt', True)
    # resid = tb['mu'] - utils.current_cosmo().distmod(tb['z_cmb']).value
    # print(np.std(resid))


    return

if __name__ == '__main__':
    main()
    # experiment()
    # maunal_fitting()

