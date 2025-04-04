# M. D. Woods
# 03/19/2025
import glob
import os.path
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
import stats              # Outputs the stats of files

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
    # # # Get normal SNe
    # queryATLAS.download(tar_list="txts/target_files/normal_tarlist.csv", save_loc='data/ATLAS-norm/')
    return
def run_queryZTF():
    # Get SN1991bg-like SNe
    queryZTF.download(tar_list="txts/target_files/sn1991bglike_tarlist.csv", save_loc = 'data/ZTF-91bg/')
    # # Get normal SNe
    # queryZTF.download(tar_list="txts/target_files/normal_tarlist.csv", save_loc='data/ZTF-norm/')
    return
def run_file_verification():
    # SN1991bg-like SNe -- Verify where the SNe from the target file are in the directories
    utils.verify_downloads(combined_tarlist_path="txts/target_files/sn1991bglike_tarlist.csv", source='CSP', subtype='91bg')
    utils.verify_downloads(combined_tarlist_path="txts/target_files/sn1991bglike_tarlist.csv", source='ATLAS', subtype='91bg')
    utils.verify_downloads(combined_tarlist_path="txts/target_files/sn1991bglike_tarlist.csv", source='ZTF', subtype='91bg')

    # normal SNe -- Verify where the SNe from the target file are in the directories
    utils.verify_downloads(combined_tarlist_path="txts/target_files/normal_tarlist.csv", source='CSP', subtype='norm')
    utils.verify_downloads(combined_tarlist_path="txts/target_files/normal_tarlist.csv", source='ATLAS', subtype='norm')
    return
def run_data_combination(clear_old: bool = True):
    dataManipulation.combine_like_data("txts/target_files/sn1991bglike_tarlist.csv",
                                       'data/combined_91bg/',
                                       '91bg', clear_old=clear_old)
    dataManipulation.check_combined_stats('91bg')
    # dataManipulation.combine_like_data("txts/target_files/normal_tarlist.csv",
    #                                    'data/combined_norm/',
    #                                    'norm', clear_old=clear_old)
    # dataManipulation.check_combined_stats('norm')
    return
def run_algo_combination(subtype: str = '91bg'):
    dataManipulation.combine_snpy_salt(f'results/{subtype}_snpy.txt',
                                       f'results/{subtype}_salt.txt',
                                       f'results/{subtype}_snpy-salt.txt')
    dataManipulation.selection_criteria(f'results/{subtype}_snpy-salt.txt',
                                        subtype,
                                        f'results/{subtype}_snpy-salt_cut.txt')

    # dataManipulation.selection_criteria(f'results/{subtype}_snpy.txt',
    #                                     subtype,
    #                                     f'results/{subtype}_snpy_cut.txt')
    # dataManipulation.selection_criteria(f'results/{subtype}_salt.txt',
    #                                     subtype,
    #                                     f'results/{subtype}_salt_cut.txt')
    # dataManipulation.combine_snpy_salt(f'results/{subtype}_snpy.txt',
    #                                    f'results/{subtype}_salt.txt',
    #                                    f'results/{subtype}_snpy-salt.txt')
    # dataManipulation.combine_snpy_salt(f'results/{subtype}_snpy_cut.txt',
    #                                    f'results/{subtype}_salt_cut.txt',
    #                                    f'results/{subtype}_snpy-salt_cut.txt')
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

    # Apply selection criteria to parameter file
    dataManipulation.selection_criteria(path=f'results/{subtype}_{algo}.txt',
                                        subtype=subtype,
                                        save_loc=f'results/{subtype}_{algo}_cut.txt')

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
        n_o = f"fitting/{algo.lower()}-plots/{name}_lc.png"
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

    pm_redNorms = 'results/old/redNormSNe.txt'
    pm_dust = 'results/old/global_dust_params.txt'

    # WORKING =========================================================================================================
    plotter.resid_v_mass(path_91bg=pm_91bg_merged_cut,
                         path_norm=pm_norm_merged_cut,
                         save_loc=final_dir + 'resid_v_mass.png',
                         label = True)
    plotter.mu_v_z(path_91bg=pm_91bg_merged_cut,
                    path_norm=pm_norm_merged_cut,
                    save_loc=final_dir+'mu_v_z.png',
                    label = False)
    # Brout+Scolnic 2021 style Dust v. Scatter
    plotter.dust_v_scatter(path_91bg = pm_91bg_salt_cut, path_norm = pm_norm_salt_cut, path_dust = pm_dust,
                           bin_num = 50, bin_bounds = [0.1, 6.3], hist_bins = 60,
                           label = True, save_loc = final_dir+'dust_v_scatter.png')

    # SALT3 Plots
    plotter.alpha_beta(path_91bg=pm_91bg_salt_cut,
                        path_norm=pm_norm_salt_cut,
                        save_loc=final_dir+'alpha_beta.png')

    # Dust Plots
    plotter.abs_mag_v_color(path_91bg=pm_91bg_salt_cut, path_norm=pm_norm_salt_cut, path_dust=pm_dust,
                            save_loc=final_dir+'absMag_v_color.png')
    plotter.dust_hist(path_91bg=pm_91bg_salt_cut, path_red_norm=pm_redNorms, path_dust=pm_dust,
                       save_loc=final_dir+'dust_params.png')

    # Paramater Histograms
    plotter.param_hist(snpy_91bg_path=pm_91bg_snpy_cut, snpy_norm_path=pm_norm_snpy_cut,
                       salt_91bg_path=pm_91bg_salt_cut, salt_norm_path=pm_norm_salt_cut,
                       save_loc=final_dir + 'param_hist_cut.png', line_type='median')
    plotter.param_hist(snpy_91bg_path=pm_91bg_snpy_uncut, snpy_norm_path=pm_norm_snpy_uncut,
                       salt_91bg_path=pm_91bg_salt_uncut, salt_norm_path=pm_norm_salt_uncut,
                       save_loc=final_dir + 'param_hist_uncut.png', line_type='median')

    # Brout+Scolnic 2021 style Paramaters (Color & Stretch) v. Scatter
    plotter.params_v_scatter(path_snpy_91bg=pm_91bg_snpy_uncut, path_snpy_norm=pm_norm_snpy_uncut,
                             path_salt_91bg=pm_91bg_salt_uncut, path_salt_norm=pm_norm_salt_uncut,
                             bin_nums=[[10, 10], [10, 10], [10, 10], [10, 10]],
                             bin_bounds=[[0.17, 1.074], [-0.299, 0.4], [-3.291, 2.838], [-0.368, 0.67]],
                             label=True, save_loc=final_dir + 'params_v_scatter.png')

    # # REDUNDANT =======================================================================================================
    # plotter.color_v_scatter(path_snpy_91bg=pm_91bg_snpy_cut, path_snpy_norm=pm_norm_snpy_cut,
    #                         path_salt_91bg=pm_91bg_salt_cut, path_salt_norm=pm_norm_salt_cut,
    #                         bin_nums=[[10, 10], [10, 10]], bin_bounds=[[-0.3, 0.6], [-0.3, 0.6]], label=True,
    #                         save_loc=final_dir + 'color_v_scatter.png')
    # plotter.resid_v_mass_dust(path_91bg=pm_91bg_merged_cut, path_norm=pm_norm_merged_cut, path_dust=pm_dust,
    #                            save_loc=final_dir+'dust_resid_v_mass.png')

    # BROKEN ==========================================================================================================


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
        tb_91bg = utils.default_open('results/91bg_snpy-salt_cut.txt', True)
        tb_91bg['resid_mu'] = tb_91bg['mu'] - utils.current_cosmo().distmod(tb_91bg['z_cmb']).value
        tb_91bg['resid_mu_err'] = np.copy(tb_91bg['mu_err'])
        mass_step_10_91bg, resid_10_91bg = plotter.mass_step_calc(tb_91bg['mu'], tb_91bg['mu_err'],
                                                                  tb_91bg['resid_mu'], tb_91bg['hostMass'],
                                                                  tb_91bg['z_cmb'], cut = 10.0)
        mass_step_med_91bg, resid_med_91bg = plotter.mass_step_calc(tb_91bg['mu'], tb_91bg['mu_err'],
                                                                    tb_91bg['resid_mu'], tb_91bg['hostMass'],
                                                                    tb_91bg['z_cmb'], cut=np.median(tb_91bg['hostMass']))
        print("% Variable Results\n"
              "% SN~1991bg-like SNe Stats\n"
              "\def \msfOurNumSNe {"+f"{len(tb_91bg)}"+"}\n"
              "\def \msfOurMassStepLogTen {"+f"{round(mass_step_10_91bg['value'], 3)}"+"}\n"
              "\def \msfOurMassStepLogTenErr {"+f"{round(mass_step_10_91bg['err'], 3)}"+"}\n"
              "\def \msfOurMassStepMedian {"+f"{round(mass_step_med_91bg['value'], 3)}"+"}\n"
              "\def \msfOurMassStepMedianErr {"+f"{round(mass_step_med_91bg['err'], 3)}"+"}")
        tb_norm = utils.default_open('results/norm_snpy-salt_cut.txt', True)
        tb_norm['resid_mu'] = tb_norm['mu'] - utils.current_cosmo().distmod(tb_norm['z_cmb']).value
        tb_norm['resid_mu_err'] = np.copy(tb_norm['mu_err'])
        mass_step_10_norm, resid_10_norm = plotter.mass_step_calc(tb_norm['mu'], tb_norm['mu_err'],
                                                                  tb_norm['resid_mu'], tb_norm['hostMass'],
                                                                  tb_norm['z_cmb'], cut=10.0)
        mass_step_med_norm, resid_med_norm = plotter.mass_step_calc(tb_norm['mu'], tb_norm['mu_err'],
                                                                    tb_norm['resid_mu'], tb_norm['hostMass'],
                                                                    tb_norm['z_cmb'], cut=np.median(tb_91bg['hostMass']))
        print("% Normal SNe Stats\n"
              "\def \msfNormNumSNe {"+f"{len(tb_norm)}"+"}\n"
              "\def \msfNormMassStepLogTen {"+f"{round(mass_step_10_norm['value'], 3)}"+"}\n"
              "\def \msfNormMassStepLogTenErr {"+f"{round(mass_step_10_norm['err'], 3)}"+"}\n"
              "\def \msfNormMassStepMedian {"+f"{round(mass_step_med_norm['value'], 3)}"+"}\n"
              "\def \msfNormMassStepMedianErr {"+f"{round(mass_step_med_norm['err'], 3)}"+"}")

        # =============================================================================================================
        # Table 1. Summary of SN 1991bg-like SNe Samples
        # Ex.
        # \begin{table*}
        #     \centering
        #     \caption{Summary of SN~1991bg-like SNe Samples}
        #     \begin{tabular}{lcrrrrrr}
        #         \hline
        #          & Number of SNe & \multicolumn{2}{r}{Redshift Range}& \multicolumn{2}{r}{Declination Range} & Average Mag. & Average Mag. Err\\
        #          &&min.&max.&min.&max.&$\pm5$~days from peak&$\pm5$~days from peak\\
        #         \hline
        #         CSP & $11$ & $0.0039$ & $0.0424$ & $-24.9441$ & $20.5262$ & $17.1799$ & $0.0259$\\
        #         ATLAS & $10$ & $0.0149$& $0.0594$ & $-44.2625$ & $74.8300$& $17.3974$ & $0.0870$\\
        #         ZTF & $20$ & $0.0160$ & $0.1002$ & $-22.4477$ & $65.4675$ & $18.8785$ & $0.1230$\\
        #         ATLAS-ZTF$^{\ast}$ & $61$ & $0.0055$ & $0.0706$ & $-27.2149$ & $86.9325$ & $18.2714$ & $0.1098$\\
        #         All Surveys &  $102$ & $0.0039$ & $0.1002$ & $-44.2625$ & $86.9325$ & $18.2427$ & $0.1049$\\
        #         \hline
        #     \end{tabular}
        #
        #     *ATLAS and ZTF have 61 SNe that overlap, so these data contains the \textit{c}, \textit{o}, \textit{g}, \textit{r}, and \textit{i}-bands.
        #     \label{tab:sample_summary_uncut}
        # \end{table*}
        tb1 = utils.default_open('results/91bg_snpy-salt.txt', True)
        tb1_csp = tb1[tb1['origin'] == 'CSP']
        tb1_atlas = tb1[tb1['origin'] == 'ATLAS']
        tb1_ztf = tb1[tb1['origin'] == 'ZTF']
        tb1_combined = tb1[tb1['origin'] == 'COMBINED']

        # Average Mag. & Mag. Err ±5 days from peak
        avg_mag_dict = {}
        for tb, source, addon in zip([tb1_csp, tb1_atlas, tb1_ztf, tb1_combined],
                                     ['CSP-91bg', 'ATLAS-91bg', 'ZTF-91bg', 'combined_91bg'],
                                     ['CSP', 'ATLAS', 'ZTF', '']):
            arr, arr_err = np.array([]), np.array([])
            for i in range(len(tb['objname'])):
                name = tb['objname'][i]
                f = f"data/{source}/{addon}{name}.txt"
                sn = fitter.sneObj(source.lower(), 'snpy', f)
                avg_mag = np.average(sn.mag[np.abs(sn.mjd - tb['Tmax'][i]) < 5])
                avg_mag_err = np.average(sn.dmag[np.abs(sn.mjd - tb['Tmax'][i]) < 5])

                arr = np.append(arr, avg_mag)
                arr_err = np.append(arr_err, avg_mag_err)

            avg_mag_dict.update({source[:-5].lower(): {'value': np.average(arr[~np.isnan(arr)]),
                                               'err': np.average(arr_err[~np.isnan(arr_err)])}})

        avg_mag_dict.update({'all': {'value': np.average([avg_mag_dict['csp']['value'],
                                                          avg_mag_dict['atlas']['value'],
                                                          avg_mag_dict['ztf']['value'],
                                                          avg_mag_dict['combined']['value']]),
                                     'err': np.average([avg_mag_dict['csp']['err'],
                                                        avg_mag_dict['atlas']['err'],
                                                        avg_mag_dict['ztf']['err'],
                                                        avg_mag_dict['combined']['err']])}})

        combined_label = "ATLAS-ZTF$^{\\ast}$"
        print("\\begin{table*}\n"
              "\t\\centering\n"
              "\t\\caption{Summary of SN~1991bg-like SNe Samples}\n"
              "\t\\begin{tabular}{lcrrrrrr}\n"
              "\t\t\\hline\n"
              "\t\t& Number of SNe & \multicolumn{2}{r}{Redshift Range}& \multicolumn{2}{r}{Declination Range} & Average Mag. & Average Mag. Err\\\\\n"
              "\t\t&&min.&max.&min.&max.&$\pm5$~days from peak&$\pm5$~days from peak\\\\\n"
              "\t\t\\hline\n"
              f"\t\tCSP & ${len(tb1_csp)}$ & ${round(np.min(tb1_csp['z_cmb']), 4)}$ & ${round(np.max(tb1_csp['z_cmb']), 4)}$ & ${round(np.min(tb1_csp['dec']), 4)}$ & ${round(np.max(tb1_csp['dec']), 4)}$ & ${round(avg_mag_dict['csp']['value'], 3)}$ & ${round(avg_mag_dict['csp']['err'], 3)}$\\\\\n"
              f"\t\tATLAS & ${len(tb1_atlas)}$ & ${round(np.min(tb1_atlas['z_cmb']), 4)}$ & ${round(np.max(tb1_atlas['z_cmb']), 4)}$ & ${round(np.min(tb1_atlas['dec']), 4)}$ & ${round(np.max(tb1_atlas['dec']), 4)}$ & ${round(avg_mag_dict['atlas']['value'], 3)}$ & ${round(avg_mag_dict['atlas']['err'], 3)}$\\\\\n"
              f"\t\tZTF & ${len(tb1_ztf)}$ & ${round(np.min(tb1_ztf['z_cmb']), 4)}$ & ${round(np.max(tb1_ztf['z_cmb']), 4)}$ & ${round(np.min(tb1_ztf['dec']), 4)}$ & ${round(np.max(tb1_ztf['dec']), 4)}$ & ${round(avg_mag_dict['ztf']['value'], 3)}$ & ${round(avg_mag_dict['ztf']['err'], 3)}$\\\\\n"
              f"\t\t{combined_label} & ${len(tb1_combined)}$ & ${round(np.min(tb1_combined['z_cmb']), 4)}$ & ${round(np.max(tb1_combined['z_cmb']), 4)}$ & ${round(np.min(tb1_combined['dec']), 4)}$ & ${round(np.max(tb1_combined['dec']), 4)}$ & ${round(avg_mag_dict['combined']['value'], 3)}$ & ${round(avg_mag_dict['combined']['err'], 3)}$\\\\\n"
              f"\t\tAll Surveys & ${len(tb1)}$ & ${round(np.min(tb1['z_cmb']), 4)}$ & ${round(np.max(tb1['z_cmb']), 4)}$ & ${round(np.min(tb1['dec']), 4)}$ & ${round(np.max(tb1['dec']), 4)}$ & ${round(avg_mag_dict['all']['value'], 3)}$ & ${round(avg_mag_dict['all']['err'], 3)}$\\\\\n"
              "\t\t\\hline\n"
              "\t\\end{tabular}\n"
              "\n"
              "\t*ATLAS and ZTF have "+f"{len(tb1_combined)}"+" SNe that overlap, so these data contains the \\textit{c}, \\textit{o}, \\textit{g}, \\textit{r}, and \\textit{i}-bands.\n"
              "\t\\label{tab:sample_summary_uncut}\n"
              "\\end{table*}\n"
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
        # pass

        # =============================================================================================================
        # Table 3. Summary of SN 1991bg-like SNe Samples After Selection Criteria
        tb3 = utils.default_open('results/91bg_snpy-salt_cut.txt', True)
        tb3_csp = tb3[tb3['origin'] == 'CSP']
        tb3_atlas = tb3[tb3['origin'] == 'ATLAS']
        tb3_ztf = tb3[tb3['origin'] == 'ZTF']
        tb3_combined = tb3[tb3['origin'] == 'COMBINED']

        # Average Mag. & Mag. Err ±5 days from peak
        avg_mag_dict = {}
        for tb, source, addon in zip([tb3_csp, tb3_atlas, tb3_ztf, tb3_combined],
                                     ['CSP-91bg', 'ATLAS-91bg', 'ZTF-91bg', 'combined_91bg'],
                                     ['CSP', 'ATLAS', 'ZTF', '']):
            arr, arr_err = np.array([]), np.array([])
            for i in range(len(tb['objname'])):
                name = tb['objname'][i]
                f = f"data/{source}/{addon}{name}.txt"
                sn = fitter.sneObj(source.lower(), 'snpy', f)
                avg_mag = np.average(sn.mag[np.abs(sn.mjd - tb['Tmax'][i]) < 5])
                avg_mag_err = np.average(sn.dmag[np.abs(sn.mjd - tb['Tmax'][i]) < 5])

                arr = np.append(arr, avg_mag)
                arr_err = np.append(arr_err, avg_mag_err)

            avg_mag_dict.update({source[:-5].lower(): {'value': np.average(arr[~np.isnan(arr)]),
                                                       'err': np.average(arr_err[~np.isnan(arr_err)])}})

        avg_mag_dict.update({'ztf': {'value': 0.00, 'err': 0.00}}) # No ZTF-only in cut data
        avg_mag_dict.update({'all': {'value': np.average([avg_mag_dict['csp']['value'],
                                                          avg_mag_dict['atlas']['value'],
                                                          # avg_mag_dict['ztf']['value'],
                                                          avg_mag_dict['combined']['value']]),
                                     'err': np.average([avg_mag_dict['csp']['err'],
                                                        avg_mag_dict['atlas']['err'],
                                                        # avg_mag_dict['ztf']['err'],
                                                        avg_mag_dict['combined']['err']])}})
        print("\\begin{table*}\n"
              "\t\\centering\n"
              "\t\\caption{Summary of SN~1991bg-like SNe Samples}\n"
              "\t\\begin{tabular}{lcrrrrrr}\n"
              "\t\t\\hline\n"
              "\t\t& Number of SNe & \multicolumn{2}{r}{Redshift Range}& \multicolumn{2}{r}{Declination Range} & Average Mag. & Average Mag. Err\\\\\n"
              "\t\t&&min.&max.&min.&max.&$\pm5$~days from peak&$\pm5$~days from peak\\\\\n"
              "\t\t\\hline\n"
              f"\t\tCSP & ${len(tb3_csp)}$ & ${round(np.min(tb3_csp['z_cmb']), 4)}$ & ${round(np.max(tb3_csp['z_cmb']), 4)}$ & ${round(np.min(tb3_csp['dec']), 4)}$ & ${round(np.max(tb3_csp['dec']), 4)}$ & ${round(avg_mag_dict['csp']['value'], 3)}$ & ${round(avg_mag_dict['csp']['err'], 3)}$\\\\\n"
              f"\t\tATLAS & ${len(tb3_atlas)}$ & ${round(np.min(tb3_atlas['z_cmb']), 4)}$ & ${round(np.max(tb3_atlas['z_cmb']), 4)}$ & ${round(np.min(tb3_atlas['dec']), 4)}$ & ${round(np.max(tb3_atlas['dec']), 4)}$ & ${round(avg_mag_dict['atlas']['value'], 3)}$ & ${round(avg_mag_dict['atlas']['err'], 3)}$\\\\\n"
              # f"\t\tZTF & ${len(tb3_ztf)}$ & ${round(np.min(tb3_ztf['z_cmb']), 4)}$ & ${round(np.max(tb3_ztf['z_cmb']), 4)}$ & ${round(np.min(tb3_ztf['dec']), 4)}$ & ${round(np.max(tb3_ztf['dec']), 4)}$ & ${round(avg_mag_dict['ztf']['value'], 3)}$ & ${round(avg_mag_dict['ztf']['err'], 3)}$\\\\\n"
              f"\t\tZTF & $0$ & $0.0000$ & $0.0000$ & $0.0000$ & $0.0000$ & $0.0000$ & $0.0000$\\\\\n"
              f"\t\tATLAS-ZTF & ${len(tb3_combined)}$ & ${round(np.min(tb3_combined['z_cmb']), 4)}$ & ${round(np.max(tb3_combined['z_cmb']), 4)}$ & ${round(np.min(tb3_combined['dec']), 4)}$ & ${round(np.max(tb3_combined['dec']), 4)}$ & ${round(avg_mag_dict['combined']['value'], 3)}$ & ${round(avg_mag_dict['combined']['err'], 3)}$\\\\\n"
              f"\t\tAll Surveys & ${len(tb3)}$ & ${round(np.min(tb3['z_cmb']), 4)}$ & ${round(np.max(tb3['z_cmb']), 4)}$ & ${round(np.min(tb3['dec']), 4)}$ & ${round(np.max(tb3['dec']), 4)}$ & ${round(avg_mag_dict['all']['value'], 3)}$ & ${round(avg_mag_dict['all']['err'], 3)}$\\\\\n"
              "\t\t\\hline\n"
              "\t\\end{tabular}\n"
              "\t\\label{tab:sample_summary_cut}\n"
              "\\end{table*}"
              )
    else:
        print(f"[!!!] '{style}' is not a valid style ['terminal'/'readme'/'latex']")
        return
    return
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
    # maunal_fitting()

    # # File acquisition
    # run_queryCSP()  # Get CSP data from CSP site
    # run_queryATLAS()  # Query ATLAS using target list of known RA & DEC
    # run_queryZTF()  # Query ZTF using target list of known RA & DEC
    # queryZTF.download_email(tar_list="txts/target_files/sn1991bglike_tarlist.csv", save_loc = 'data/ZTF-91bg/')

    # # File verification
    # utils.get_twomass()  # Checks for twomass velocities data (large file)
    # run_file_verification()  # Using target file, checks if SN exists in the proper directory

    # # Data combination
    # run_data_combination()  # Combines data file if SN exsists in multiple surveys

    # Fit light curves
    # run_fitter(subtype = '91bg', algo = 'snpy', rewrite = True)
    # run_fitter(subtype = '91bg', algo = 'salt', rewrite = True)
    # run_fitter(subtype = 'norm', algo = 'snpy', rewrite = True)
    # run_fitter(subtype = 'norm', algo = 'salt', rewrite = True)

    # # # Combine algorithms
    # run_algo_combination('91bg')
    # run_algo_combination('norm')

    # # Check lc of cut file
    # run_lc_checker('91bg')
    # run_lc_checker('norm')

    # # Plot data
    # run_plotter()

    # # Outout stats
    # run_stats('latex')
    return
@run_time
def experiment():
    tb_snpy = utils.default_open('results/91bg_snpy.txt', True)
    tb_salt = utils.default_open('results/91bg_salt.txt', True)

    for n in tb_snpy['objname']:
        if n in list(tb_salt['objname']):
            print(list(tb_salt[tb_salt['objname'] == n]))
            print(list(tb_snpy[tb_snpy['objname'] == n]))
        break




    # overlap_files = glob.glob('data/combined_91bg/2*.txt')
    # all_names, all_mu_diff, all_mu_err_diff = [], [], []
    # combined_atlas_mu_diff, combined_atlas_mu_err_diff = [], []
    # combined_ztf_mu_diff, combined_ztf_mu_err_diff = [], []
    # for index, file in enumerate(overlap_files):
    #     name = file.split('/')[-1][:-4]
    #     hdr = f"[{index+1} / {len(overlap_files)}] {name} ============================================================"
    #     clr = "="*len(hdr)
    #     print(hdr)
    #
    #     algo = 'snpy'
    #     sn_atlas = fitter.fit(f'data/ATLAS-91bg/ATLAS{name}.txt', algo, False)[0]
    #     if sn_atlas is None:
    #         print(f"\n-=-=-=-=-=-=-=- [!!!] No match/fit for {name} in ATLAS-91bg data with {algo}! -=-=-=-=-=-=-=- \n")
    #         systime.sleep(1)
    #         continue
    #     sn_ztf = fitter.fit(f'data/ZTF-91bg/ZTF{name}.txt', algo, False)[0]
    #     if sn_ztf is None:
    #         print(f"\n-=-=-=-=-=-=-=- [!!!] No match/fit for {name} in ZTF-91bg data with {algo}! -=-=-=-=-=-=-=-\n")
    #         systime.sleep(1)
    #         continue
    #     sn_combined = fitter.fit(file, algo, False)[0]
    #     if sn_combined is None:
    #         print(f"\n-=-=-=-=-=-=-=- [!!!] No match/fit for {name} in COMBINED data with {algo}! -=-=-=-=-=-=-=-\n")
    #         systime.sleep(1)
    #         continue
    #
    #     # sn_atlas = fitter.fit(f'data/ATLAS-91bg/ATLAS{name}.txt', 'snpy', True)[0]
    #     # if (sn_atlas is None) or (sn_atlas.params['mu']['value'] < 0):
    #     #     print(f"\n-=-=-=-=-=-=-=-=-= [~~~] SNooPy failed, attemping SALT3 for {name} -=-=-=-=-=-=-=-=-=\n")
    #     #     sn_atlas = fitter.fit(f'data/ATLAS-91bg/ATLAS{name}.txt', 'salt', True)[0]
    #     #     if (sn_atlas is None) or (sn_atlas.params['mu']['value'] < 0):
    #     #         print(f"\n-=-=-=-=-=-=-=-=-= [!!!] No match for {name} in ATLAS-91bg data! -=-=-=-=-=-=-=-=-=\n")
    #     #         continue
    #
    #     # sn_ztf = fitter.fit(f'data/ZTF-91bg/ZTF{name}.txt', 'snpy', True)[0]
    #     # if (sn_ztf is None) or (sn_ztf.params['mu']['value'] < 0):
    #     #     print(f"\n-=-=-=-=-=-=-=-=-= [~~~] SNooPy failed, attemping SALT3 for {name} -=-=-=-=-=-=-=-=-=\n")
    #     #     sn_ztf = fitter.fit(f'data/ZTF-91bg/ZTF{name}.txt', 'salt', True)[0]
    #     #     if (sn_ztf is None) or (sn_ztf.params['mu']['value'] < 0):
    #     #         print(f"\n-=-=-=-=-=-=-=-=-= [!!!] No match for {name} in ZTF-91bg data! -=-=-=-=-=-=-=-=-=\n")
    #     #         continue
    #
    #     mu_dif = np.abs(sn_atlas.params['mu']['value'] - sn_ztf.params['mu']['value'])
    #     mu_err_diff = np.abs(sn_atlas.params['mu']['err'] - sn_ztf.params['mu']['err'])
    #     all_mu_diff.append(mu_dif)
    #     all_mu_err_diff.append(mu_err_diff)
    #     all_names.append(name)
    #
    #     combined_atlas_mu_diff.append(np.abs(sn_atlas.params['mu']['value'] - sn_combined.params['mu']['value']))
    #     combined_atlas_mu_err_diff.append(np.abs(sn_atlas.params['mu']['err'] - sn_combined.params['mu']['err']))
    #     combined_ztf_mu_diff.append(np.abs(sn_ztf.params['mu']['value'] - sn_combined.params['mu']['value']))
    #     combined_ztf_mu_err_diff.append(np.abs(sn_ztf.params['mu']['err'] - sn_combined.params['mu']['err']))
    #
    #     print(f"\n-=-=-=-=-=-=-= [{index+1} / {len(overlap_files)}] {name}: {mu_dif}+/-{mu_err_diff} =-=-=-=-=-=-=-\n")
    #     print(clr)
    #     # systime.sleep(2)
    #
    #     # if index > 10:
    #     #     break
    #
    # list11, list12, list21, list22, list31, list32 = [], [], [], [], [], []
    # for i, d_mu in enumerate(all_mu_diff):
    #     if all_mu_err_diff[i] > 100:
    #         print(f"{all_names[i]}: Error")
    #     elif d_mu > 1 or all_mu_err_diff[i] > 0.1:
    #         print(f"{all_names[i]}: {d_mu}+/-{all_mu_err_diff[i]}")
    #     else:
    #         list11.append(d_mu)
    #         list12.append(all_mu_err_diff[i])
    # print("===========================================================================================================")
    #
    # for i, d_mu in enumerate(combined_atlas_mu_diff):
    #     if combined_atlas_mu_err_diff[i] > 100:
    #         print(f"{all_names[i]}: Error")
    #     elif d_mu > 1 or combined_atlas_mu_err_diff[i] > 0.1:
    #         print(f"{all_names[i]}: {d_mu}+/-{combined_atlas_mu_err_diff[i]}")
    #     else:
    #         list21.append(d_mu)
    #         list22.append(combined_atlas_mu_err_diff[i])
    # print("===========================================================================================================")
    #
    # for i, d_mu in enumerate(combined_ztf_mu_diff):
    #     if combined_ztf_mu_err_diff[i] > 100:
    #         print(f"{all_names[i]}: Error")
    #     elif d_mu > 1 or combined_ztf_mu_err_diff[i] > 0.1:
    #         print(f"{all_names[i]}: {d_mu}+/-{combined_ztf_mu_err_diff[i]}")
    #     else:
    #         list31.append(d_mu)
    #         list32.append(combined_ztf_mu_err_diff[i])
    #
    # print("===========================================================================================================")
    # print(f"Average Dist. Mod. ATLAS-ZTF Difference: {np.average(list11)}")
    # print(f"Average Dist. Mod. Err. ATLAS-ZTF Difference: {np.average(list12)}")
    # print(len(list11), np.average(list11), np.average(list12))
    # print("===========================================================================================================")
    # print(f"Average Dist. Mod. (ATLAS+ZTF)-ATLAS Difference: {np.average(list21)}")
    # print(f"Average Dist. Mod. Err. (ATLAS+ZTF)-ATLAS Difference: {np.average(list22)}")
    # print(len(list21), np.average(list21), np.average(list22))
    # print("===========================================================================================================")
    # print(f"Average Dist. Mod. (ATLAS+ZTF)-ZTF Difference: {np.average(list31)}")
    # print(f"Average Dist. Mod. Err. (ATLAS+ZTF)-ZTF Difference: {np.average(list32)}")
    # print(len(list31), np.average(list31), np.average(list32))
    # print("===========================================================================================================")


    # # Sigma Clipping
    # from astropy.stats import sigma_clip, sigma_clipped_stats
    #
    # sigma_mask = [all_mu_diff == sigma_clip(data=all_mu_diff, masked=False, sigma=3)]
    # all_mu_diff = all_mu_diff[sigma_mask]
    # all_mu_err_diff = all_mu_err_diff[sigma_mask]
    #
    # sigma_mask = [combined_atlas_mu_diff == sigma_clip(data=combined_atlas_mu_diff, masked=False, sigma=3)]
    # combined_atlas_mu_diff = combined_atlas_mu_diff[sigma_mask]
    # combined_atlas_mu_err_diff = combined_atlas_mu_err_diff[sigma_mask]
    #
    # sigma_mask = [combined_ztf_mu_diff == sigma_clip(data=combined_ztf_mu_diff, masked=False, sigma=3)]
    # combined_ztf_mu_diff = combined_ztf_mu_diff[sigma_mask]
    # combined_ztf_mu_err_diff = combined_ztf_mu_err_diff[sigma_mask]

    # print("===========================================================================================================")
    # print(f"Average Dist. Mod. ATLAS-ZTF Difference: {np.average(all_mu_diff)}")
    # print(f"Average Dist. Mod. Err. ATLAS-ZTF Difference: {np.average(all_mu_err_diff)}")
    # print(len(all_mu_diff), np.average(all_mu_diff), np.average(all_mu_err_diff))
    # print("===========================================================================================================")
    # print(f"Average Dist. Mod. (ATLAS+ZTF)-ATLAS Difference: {np.average(combined_atlas_mu_diff)}")
    # print(f"Average Dist. Mod. Err. (ATLAS+ZTF)-ATLAS Difference: {np.average(combined_atlas_mu_err_diff)}")
    # print(len(combined_atlas_mu_diff), np.average(combined_atlas_mu_diff), np.average(combined_atlas_mu_err_diff))
    # print("===========================================================================================================")
    # print(f"Average Dist. Mod. (ATLAS+ZTF)-ZTF Difference: {np.average(combined_ztf_mu_diff)}")
    # print(f"Average Dist. Mod. Err. (ATLAS+ZTF)-ZTF Difference: {np.average(combined_ztf_mu_err_diff)}")
    # print(len(combined_ztf_mu_diff), np.average(combined_ztf_mu_diff), np.average(combined_ztf_mu_err_diff))
    # print("===========================================================================================================")


    # for name, mu_diff, mu_err_diff in zip(all_names, all_mu_diff, all_mu_err_diff):
    #     print(f"{name}: {mu_diff}+/-{mu_err_diff}")

    # sne_diffs = """
    # 2022omn: 0.43746447116605935+/-0.10508640074681364
    # 2024rkh: 0.1811057899193571+/-0.16183921497712436
    # 2022dsu: 1.282140844002349+/-0.05534293442307786
    # 2023abdv: 1.7267826002825046+/-2.9424837089284868
    # 2024xdl: 3.744972594249841+/-0.4728160316645956
    # 2023fwb: 0.3298418043455129+/-0.07125915267218907
    # 2021fnr: 0.39562462461530146+/-0.024385910984806698
    # 2019luv: 0.7743605895695538+/-0.08963125510210432
    # 2022kbc: 0.4920627391841137+/-0.15939638351438445
    # 2021zqs: 1.18306706523758+/-0.713787144789337
    # 2019cdc: 0.5044361406172229+/-0.0018202308293510383
    # 2022oux: 1.1263491947170081+/-0.02564679442170792
    # 2021zsz: 1.0041593491934435+/-0.12848467851006637
    # 2021abzd: 20.976178158231065+/-0.11606900291769076
    # 2019noh: 3.119025466329788+/-3.1056105601134347
    # 2018ame: 0.7905112192257491+/-0.16235914078155106
    # 2022aecb: 0.5394240782791329+/-1.3077178730571617
    # 2020ecn: 3.4590599609611488+/-2.5556582505465704
    # 2021xfm: 1.5087279632335324+/-0.4935881123291739
    # 2024luo: 1.393763812857685+/-0.01629514395398103
    # 2022uxl: 2.1792299755037305+/-0.281333638352215
    # 2021pom: 4.111762035322819+/-0.4517131718931864
    # 2020abpe: 0.1523178989100913+/-0.04613776688072789
    # 2024zls: 2.537549092277395+/-0.9181065700315219
    # 2022zsp: 3.7566796288328703+/-4.455983836547228
    # 2018ast: 3.322327568777496+/-0.16766175157696542
    # 2022xkq: 1.17835512971741+/-0.22624300277487716
    # 2024wdg: 0.1160199036689491+/-0.289387614325861
    # 2024bjb: 0.42663986122726527+/-0.21836054112342052
    # 2020aclx: 0.8510916535443371+/-0.7871856941306603
    # 2023cvq: 3.7936590804384522+/-0.6181935965503224
    # 2020mfd: 0.06553900042540306+/-0.28584533572277004
    # 2024fid: 121.65228205097502+/-1.7292152916950199
    # 2018ciw: 0.06275780466481962+/-0.4169477463671217
    # 2019moq: 1.9588193405262473+/-0.20847788060954706
    # 2018gro: 1.646541200338973+/-4.113524692707137
    # 2024jih: 28.300880203375847+/-0.11290546552259632
    # 2022ihz: 0.23193531305316384+/-1.0047190876353176
    # 2021cad: 2.2350115033078524+/-0.5011313547889173
    # 2019cp: 3.790027163320829+/-0.3357068246315346
    # 2020nta: 0.5918173165850149+/-0.011703888711282334
    # 2020acoo: 3.983524949356209+/-0.6273732291146854
    # 2024jgw: 7.008262183866272+/-0.5098458959865377
    # 2024zwb: 0.28358585193144137+/-0.13090574864679205
    # 2024vcj: 4.670173621348283+/-1.0769188897902942
    # 2024jhk: 2.7555246334027217+/-24.020932093856512
    # 2021qqr: 0.6784759073238504+/-0.18581526975466722
    # 2018baz: 5.4370102342053315+/-0.6249260231145987
    # 2021twa: 1.2893063858166158+/-0.3803775345197379
    # 2018efn: 0.056464113880686284+/-0.11521143097094837
    # 2023yrs: 22.061243727029627+/-0.9392550824433331
    # 2022rey: 0.2200958286350314+/-6.092244033355324
    # 2018jag: 0.5658984523600026+/-0.042599637683314415
    # 2020acbx: 1.3686940336514297+/-10.378631088557956
    # 2022abom: 0.6118704420858876+/-1.0124729929124263
    # 2019be: 0.9660875824609931+/-0.024936133563702396
    # 2021uef: 2.264527180085338+/-0.07667296580235243
    # 2023sps: 0.8317217494537701+/-0.5726733380010692
    # 2022zvu: 2.085160176404756+/-0.006593893793095176
    # 2022xnq: 1.3940951192716255+/-0.45590049739170607
    # 2018eyi: 0.5475073468061993+/-5.648428777786931
    # 2022skw: 20.479946890546795+/-11.377080457704805
    # 2020abmg: 1.3668110108165266+/-0.21497037137354508
    # 2023vjh: 2.250233017093585+/-0.15548518670877687
    # 2024xhs: 1.543916390949036+/-1.522749278877741
    # 2021wzb: 2.5176111202001437+/-0.24322145028177128
    # 2021jvp: 1.5017694038441576+/-0.27870211794182187
    # 2022aaok: 0.19273231830782578+/-0.02648165656876797
    # 2020vae: 3.8158049765461755+/-0.058799433622846
    # 2021agej: 1.0345037004930688+/-21.375446405119174
    # 2018lph: 1.464545085264099+/-0.08415039456938098
    # 2024abwg: 2.339943504613423+/-0.6807205293287651
    # 2021bls: 0.26311892955337157+/-0.0697471118345738
    # 2021uve: 0.37595558383133465+/-0.030078344440739528
    # 2022bsi: 2.658011145579138+/-3.338268787147059
    # 2022fjx: 1.787626903995708+/-0.13706721988866313
    # 2021qvv: 1.1680176720104285+/-3.5733505835008077
    # 2024sag: 0.9618996441052587+/-0.18151142841371615
    # """

    # all_mu_diff, all_mu_err_diff = np.array([]), np.array([])
    # split_sne_diffs = sne_diffs.split('\n')[1:]
    # for i in range(int(len(split_sne_diffs)/2)):
    #     name = split_sne_diffs[i].split(": ")[0][4:]
    #     mu_diff, mu_err_diff = split_sne_diffs[i].split(": ")[-1].split('+/-')
    #     all_mu_diff = np.append(all_mu_diff, float(mu_diff))
    #     all_mu_err_diff = np.append(all_mu_err_diff, float(mu_err_diff))
    # print(f"Average Dist. Mod. Difference: {np.average(all_mu_diff)}")
    # print(f"Average Dist. Mod. Err. Difference: {np.average(all_mu_err_diff)}")
    # print(len(all_mu_diff), np.average(all_mu_diff), np.average(all_mu_err_diff))

    # NaN fitting: ['ZTF2019cp']
    # TNS Error: ['2016iuh', '2018awi', '2016ije', '2016brx']
    # sne = fitter.fit('data/ZTF-91bg/*.txt', 'salt', True)

    # run_data_combination()


    # # Replacing Unknown discvery dates
    # import glob
    # from astropy.table import Table
    # import numpy as np
    # data = np.genfromtxt("txts/target_files/normal_tarlist.csv", delimiter=',', skip_header=1, dtype=str)
    # tar_tb = Table(names=data[0, :], data=data[1:, :])
    # with open("txts/target_files/normal_tarlist_new.csv", 'w') as f:
    #     print("# Object Type: Normal SNe Ia; Date: 2004 to 2025; Number SNe: 500; CSP: ~ | ATLAS: ~ | ZTF: ~ | ATLAS+ZTF: ~\n"
    #           "Name,RA,DEC,Redshift,Discovery Date (UT),Data Source(s)", file=f)
    #     comb_norms = glob.glob("data/combined_norm/*.txt")
    #     previous_done = ['2008go', '2006hx', '2005bo', '2009I', '2007bm', '2008hv', '2004gc', '2008gl', '2006D', '2009al', '2008O', '2005mc', '2007A', '2007le', '2006ej', '2008hu', '2004gs', '2009Y', '2007af', '2005lu', '2007ol', '2008fl', '2006py', '2004eo', '2006kf', '2008ia', '2007on', '2005ag', '2007as', '2006fw', '2004ey', '2008cd', '2007jd', '2005ku', '2005kc', '2008bc', '2007jg', '2007hx', '2007sr', '2008cf', '2008bf', '2006lu', '2008bq', '2007hj', '2006bh', '2008cc', '2007st', '2005el', '2005hj', '2009ds', '2008ar', '2007so', '2006ob', '2006br', '2005ki', '2005ir', '2006os', '2006ax', '2005iq', '2005hc', '2008bz', '2008fu', '2006ev', '2007nq', '2007ca', '2008R', '2005bg', '2007bd', '2009ag', '2006X', '2006dd', '2008gg', '2008gp', '2008hj', '2006gj', '2005A', '2009cz', '2005be', '2009ad', '2006et', '2006is', '2009P', '2007bc', '2008C', '2007cg', '2005am', '2005al', '2009aa', '2006ef', '2006eq', '2008fr', '2004ef', '2008fp', '2009D', '2005na']
    #     print(len(previous_done))
    #     for i, file in enumerate(comb_norms):
    #         n = file.split('/')[-1][:-4].split('CSP')[-1].split('ATLAS')[-1].split('ZTF')[-1]
    #         if n in previous_done: continue
    #
    #         print(f"[{i+1}/{len(comb_norms)}] {file}")
    #         r = utils.query_tns(objname=n)
    #
    #         line = ''
    #         for col in tar_tb[tar_tb['Name'] == f"SN {n}"].as_array()[0]:
    #             line += f'{col}, '
    #         print(line[:-2].replace("Unknown", f"{r['discoverydate']}"), file=f)
    #         previous_done.append(n)
    #         print(previous_done)
    #         # break




    # import glob
    # import os


    # old_names, new_names = [], []
    # for o_file in old_ztf: old_names.append(o_file.split('/')[-1][3:-4])
    # for n_file in new_ztf: new_names.append(n_file.split('/')[-1][3:-4])

    # unique_names = []
    # for n in new_names:
    #     if n not in old_names:
    #         unique_names.append(n)
    # print(unique_names)

    # print(len(unique_names)) # 112
    # for n in unique_names:
    #     shutil.copy(f"data/ZTF-91bg/ZTF{n}.txt", f"data/ZTF-91bg-v6/ZTF{n}.txt")





    # """
    # Update ZTF names to TNS names
    # """
    # import glob
    # import numpy as np
    # from astropy.time import Time as astrotime
    # new_files = glob.glob("/Users/mekhiwoods/Downloads/03252025_ztfdata/downloaded_data/*.txt")
    # for i, file in enumerate(new_files):
    #     if i+1 < (73):
    #         print(f"[~~~] ({i+1}/{len(new_files)}) File already checked, skipping...")
    #         continue
    #
    #     with open(file, 'r') as f:
    #         for n in range(3): f.readline()
    #         ra = f.readline().split(' ')[-2]
    #         dec = f.readline().split(' ')[-2]
    #
    #     data = np.genfromtxt(file, delimiter=' ', skip_header=57)
    #     mjd = data[:, 22].astype(float) - 2400000.5
    #     mjds, mjde = np.min(mjd), np.max(mjd)
    #     year_s = int(str(astrotime(mjds, format='mjd', scale='utc').to_datetime())[:4])
    #     year_e = int(str(astrotime(mjde, format='mjd', scale='utc').to_datetime())[:4])
    #
    #     r = utils.query_tns(coords=[ra, dec])
    #     r_year, r_name = int(r['objname'][:4]), r['objname']
    #
    #     if r_year <= year_e and r_year >= year_s:
    #         print(f"[+++] ({i+1}/{len(new_files)}) Copied {r_name} to... {file}")
    #         shutil.copy(file, f'data/ZTF-91bg/ZTF{r_name}.txt')
    #     else:
    #         print(f"[---] ({i+1}/{len(new_files)}) Date mismatch for {r_name}!")


    # from astropy.time import Time as astrotime
    # from astropy.table import Table
    # import matplotlib.pyplot as plt
    #
    # new_ztf = glob.glob('new_ztf/renamed/*.txt')
    # old_ztf = glob.glob('data/ZTF-91bg/*.txt')
    #
    # new_ztf_names, old_ztf_names = [], []
    # for f in new_ztf: new_ztf_names.append(f.split('/')[-1][3:-4])
    # for f in old_ztf: old_ztf_names.append(f.split('/')[-1][3:-4])
    #
    # for n in new_ztf_names:
    #     if n not in old_ztf_names:
    #         print(f"SN {n}")


    # new_ztf = glob.glob('new_ztf/*.txt')
    # for n, ztf_file in enumerate(new_ztf):
    #     if n < 23:
    #         continue
    #     print(f"[{n+1} / {len(ztf_file)}] {ztf_file}")
    #
    #     # Check if empty
    #     if len(np.genfromtxt(ztf_file)) == 0:
    #         print("[!!!] Empty file...")
    #         continue
    #
    #     # Read file
    #     with open(ztf_file, 'r') as f:
    #         for i in range(3): f.readline() # Skip first lines
    #         ra = f.readline().split(' ')[-2]
    #         dec = f.readline().split(' ')[-2]
    #         for i in range(50): f.readline() # Skip rest of header
    #         hdr = f.readline()[1:-1].split(', ')
    #     data = np.genfromtxt(ztf_file, delimiter=' ', skip_header=57, dtype=str)
    #     tb = Table(names=hdr, data=data)
    #
    #     # Get TNS data
    #     # utc = str(astrotime(np.median(jd), format='jd', scale='utc').to_datetime())[:4]
    #     r = utils.query_tns(coords=[ra, dec])
    #     objname, discdate = r['objname'], r['discoverydate']
    #     discdate = astrotime(discdate, format='iso', scale='utc').mjd
    #
    #     # Save copy to file
    #     print(f"[~~~] {ztf_file} -> new_ztf/renamed/ZTF{objname}.txt")
    #     shutil.copy(ztf_file, f"new_ztf/renamed/ZTF{objname}.txt")


        # # Clean table
        # tb = tb[tb['forcediffimflux'] != 'null']
        # tb['jd'] = tb['jd'].astype(float) - 2400000.5
        # tb.rename_column('jd', 'mjd')
        # tb = tb[(tb['mjd'].astype(float) < discdate+100)]
        # tb = tb[(tb['mjd'].astype(float) > discdate-100)]
        #
        # plt.title(f"{objname}")
        # colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
        # for i, filter in enumerate(np.unique(tb['filter'])):
        #     plt.errorbar(tb['mjd'].astype(float)[tb['filter'] == filter],
        #                  tb['forcediffimflux'].astype(float)[tb['filter'] == filter],
        #                  yerr=tb['forcediffimfluxunc'].astype(float)[tb['filter'] == filter],
        #                  fmt='o', color=colors[i], label=filter)
        # plt.legend()
        # plt.axvline(discdate)
        # plt.show()
        #
        #
        #
        #
        # input("Next?...")
        # break



    # ztf_all = ("2024jhk SN 2024bjb SN 2023fwb SN 2023fot SN 2023cvq SN 2022zvu SN 2022ubt SN 2022rjs SN 2022rey "
    #            "SN 2022dsu SN 2022aaok SN 2021wzb SN 2021uve SN 2021uef SN 2021qqr SN 2021jvp SN 2021fnr SN 2021bls "
    #            "SN 2021afur SN 2021abzd SN 2020vae SN 2020nta SN 2020acbx SN 2020abpe SN 2019op SN 2019moq SN 2019luv "
    #            "SN 2019ecx SN 2019cp SN 2019cdc SN 2019be SN 2018lph SN 2018jag SN 2018eyi SN 2018ame SN 2018baz "
    #            "SN 2018hkw SN 2019ahh SN 2020abmg SN 2020acoo SN 2020ais SN 2020ecn SN 2020mfd SN 2020yo SN 2021agej "
    #            "SN 2021agnf SN 2021bmu SN 2021pom SN 2021xfm SN 2021zsz SN 2022bsi SN 2022omn SN 2022xhh SN 2022xnq "
    #            "SN 2022zsp SN 2023abdv SN 2023sps SN 2023vjh SN 2024abwg SN 2024fid SN 2024jih SN 2024rkh SN 2024sag "
    #            "SN 2024vcj SN 2024wdg SN 2024xdl SN 2024xhs SN 2024yhg SN 2024zls SN 2024zwb SN 2018ast SN 2018ciw "
    #            "SN 2018efn SN 2018gro SN 2019noh SN 2021aare SN 2021cad SN 2021gel SN 2021qvv SN 2021twa SN 2021zqs "
    #            "SN 2022abom SN 2022fjx SN 2022ihz SN 2022kbc SN 2022oux SN 2022uxl SN 2022vxf SN 2022xkq SN 2023acdv "
    #            "SN 2023bhm SN 2023dk SN 2023omo SN 2023yrs SN 2024jgw SN 2024luo SN 2024nev SN 2024pbd "
    #            "SN 2025nn").split(" SN ")
    # ztf_cur = glob.glob("data/ZTF-91bg/*.txt")
    #
    # ztf_cur_names = []
    # for f in ztf_cur: ztf_cur_names.append(f.split("/")[-1][3:-4])
    #
    # for n in ztf_all:
    #     if n not in ztf_cur_names:
    #         print(n)

    # Plotting the problem children from the cross reference sheet
    # from astropy.table import Table
    # import matplotlib.pyplot as plt
    # h_l = {'mjd': 'mjd', 'mag': 'mag', 'dmag': 'dmag', 'filter': 'filter'}
    # problems = (("2018eyi SN 2018jag SN 2018lph SN 2019cp SN 2019luv SN 2019moq SN 2019op SN 2020abpe "
    #             "SN 2020acbx SN 2020vae SN 2021abzd SN 2021bls SN 2021fnr SN 2021uef SN 2021uve SN 2021wzb "
    #             "SN 2022dsu SN 2022rey SN 2022rjs SN 2022ubt SN 2022zvu SN 2023cvq SN 2023fwb SN 2024jhk").split(" SN "))
    # for i, n in enumerate(problems):
    #     # sn = fitter.fit(data_loc=f'data/combined_91bg/{n}.txt', algo='snpy', rewrite=True)
    #
    #     f = f'data/combined_91bg/{n}.txt'
    #     tb = Table.read(f, format="csv", delimiter=',', header_start=1)
    #
    #     # Magnitude Error cut, dmag < 1mag
    #     tb = tb[tb[h_l['dmag']] < 1]
    #
    #     print(f"[{i+1}/{len(problems)}] SN{n} | {len(tb)}pts @ {f}")
    #
    #     plt.figure(figsize=(12, 8))
    #     colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    #     for i, f in enumerate(np.unique(tb[h_l['filter']])):
    #         plt.errorbar(tb[h_l['mjd']][tb[h_l['filter']]==f], tb[h_l['mag']][tb[h_l['filter']]==f],
    #                      yerr=tb[h_l['dmag']][tb[h_l['filter']]==f],
    #                      fmt='o', color=colors[i], label=f)
    #     plt.title(f"SN{n} | {len(tb)}pts @ {f}")
    #     plt.legend()
    #     plt.gca().invert_yaxis()
    #     plt.show()
    #
    #
    #     # break
    #     input()



    # tb = utils.default_open('results/91bg_snpy-salt.txt', True)
    # for n in tb['objname']:
    #     print("SN", n)


    # from astropy.table import Table
    #
    # data = np.genfromtxt("txts/error_report_test.txt", skip_header=7, dtype=str, delimiter=';')
    # report_tbl = Table(names=data[0, :], data=data[1:, :], dtype=[str, list, list, list, list, list, list, list])
    #
    # for l in ["99999", "124", "107", "434", "222", "333", "404"]:
    #     err_col = []
    #     for j in report_tbl[l]:
    #         if len(j) > 2:
    #             for z in j[1:-1].replace("'", "").split(', '):
    #                 # print(z)
    #                 if z not in err_col:
    #                     err_col.append(z)
    #     if l == '404':
    #         for e in err_col:
    #             print("SN", e)
    # #     print(err_col)
    # #     # print(len(report_tbl[l]))
    # #     # print(report_tbl[l])



    # # # Run only ZTF data
    # # sne_ztf = fitter.fit(data_loc=f'data/ZTF-91bg/ZTF*.csv', algo='snpy', rewrite=True)
    #
    # import matplotlib.pyplot as plt
    # from astropy.table import Table
    # from astropy.stats import sigma_clip, sigma_clipped_stats
    #
    # ztf_files = glob.glob('data/ZTF-91bg/*.csv')
    # atlas_files = glob.glob('data/ATLAS-91bg/*.txt')
    #
    # h_l = {'mjd': 'mjd', 'mag': 'mag', 'dmag': 'magerr', 'filter': 'filtercode', 'limitmag': 'limitmag'}
    # # h_l = {'mjd': 'MJD', 'mag': 'm', 'dmag': 'dm', 'filter': 'F', 'limitmag': 'mag5sig'}
    # for f in ztf_files:
    # # for f in atlas_files:
    #     tb = Table.read(f, format="csv")
    #     # Magnitude Error cut, dmag < 1mag
    #     tb = tb[tb[h_l['dmag']] < 1]
    #
    #     # Magnitude Limit cut, +/- 3mag of magnitude limit
    #     limitmag, magtol = abs(tb[h_l['limitmag']] - tb[h_l['mag']]), 3
    #     tb = tb[limitmag < magtol]
    #
    #     print("SN", f.split('/')[-1][:-4], len(tb), f)
    #
    #     plt.figure(figsize=(12, 8))
    #     colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    #     for i, f in enumerate(np.unique(tb[h_l['filter']])):
    #         plt.errorbar(tb[h_l['mjd']][tb[h_l['filter']]==f], tb[h_l['mag']][tb[h_l['filter']]==f],
    #                      yerr=tb[h_l['dmag']][tb[h_l['filter']]==f],
    #                      fmt='o', color=colors[i], label=f)
    #     plt.title(f"SN{f.split('/')[-1][:-4]}")
    #     plt.legend()
    #     plt.gca().invert_yaxis()
    #     plt.show()
    #
    #     # break
    #     input()


    # Test if only has one filter and average error on points
    # ztf_files = glob.glob('data/ZTF-91bg/*.csv')

    # # Checking ATLAS light curves
    # atlas_files = glob.glob('data/ATLAS-91bg/*.txt')
    # bad_sne = 0
    # all_avg_dmag = []
    # for i, f in enumerate(atlas_files):
    #     filter_low_n = False
    #     data = np.genfromtxt(f, delimiter=',', dtype=str)
    #     hdr, data = list(data[0, :]), data[1:, :]
    #
    #     filters = data[:, hdr.index('F')]
    #     unq_filters = np.unique(filters)
    #     filter_line = "[~~~] Filters: "
    #     for unq_f in unq_filters:
    #         filter_line += f"{unq_f}({len(filters[filters == unq_f])}) "
    #         if len(filters[filters == unq_f]) < 10: filter_low_n = True
    #
    #     dmag = data[:, hdr.index('dm')].astype(float)
    #     avg_dmag = np.average(dmag[~np.isnan(dmag)])
    #     all_avg_dmag.append(avg_dmag)
    #
    #     if (len(unq_filters) == 1) or (filter_low_n) or (avg_dmag > 3.3):
    #         bad_sne += 1
    #         continue
    #     else:
    #         # shutil.copy(f, f"data/ATLAS-91bg_good/{f.split('/')[-1]}")
    #
    #         import matplotlib.pyplot as plt
    #         from astropy.stats import sigma_clip, sigma_clipped_stats
    #
    #         mjd = data[:, hdr.index('MJD')].astype(float)
    #         mag = data[:, hdr.index('m')].astype(float)
    #         dmag = data[:, hdr.index('dm')].astype(float)
    #         dmag = data[:, hdr.index('dm')].astype(float)
    #         mag5sig = data[:, hdr.index('mag5sig')].astype(float)
    #
    #         # Cut on error
    #         _, _, dmag_std = sigma_clipped_stats(dmag[~np.isnan(dmag)])
    #         dmag_limit = dmag_std*3
    #         mag5sig = mag5sig[dmag < dmag_limit]
    #         mjd = mjd[dmag < dmag_limit]
    #         mag = mag[dmag < dmag_limit]
    #         dmag = dmag[dmag < dmag_limit]
    #
    #         # Cut on within 1mag of 5-simga mag limit
    #         limitmag, magtol = abs(mag5sig - mag), 1
    #         mjd = mjd[limitmag < magtol]
    #         mag = mag[limitmag < magtol]
    #         dmag = dmag[limitmag < magtol]
    #
    #
    #         plt.figure(figsize=(12, 8))
    #         plt.errorbar(mjd, mag, yerr=dmag, fmt='o')
    #         plt.gca().invert_yaxis()
    #         plt.show()
    #     input("Next? ")




            # bad_sne += 1
            # print(f"[{i + 1} / {len(atlas_files)}] =========================================================================")
            # print(filter_line)
            # print(f"[~~~] Average dmag: {avg_dmag}")

    # print(f"Bad SNe: {bad_sne}/{len(atlas_files)}")

    # ATLAS
    # Bad SNe: 36/120
    # (len(unq_filters) == 1)   -> 11/120
    # (filter_low_n)            -> 11/120
    # (avg_dmag > 3.3)          -> 21/120



    # LIKELY CANIDATES: 2024zaj, 2024yhg, 2024ulb, 2024bjb, 2023omo, 2023abdv, 2022ihz, 2022fjx, 2022abom, 2021zsz, 2021aare, 2023sps, 2021zqs, 2018ast
    # ATLAS ONLY: 2023jah, 2022ywc, 2022vse, 2022an, 2020dq, 2020aejj, 2019so,
    #   2024stt, 2024pkh, 2024lhg, 2023mkp, 2022yv, 2022skw, 2022aecb, 2021mab, 2021jbp, 2020aclx, 2019exc, 2023ex,
    # ZTF ONLY: 2020yo, 2021afur, 2020ecn, 2018baz
    # NO ATLAS OR ZTF: 2017hle, 2018awi, 2017fzw, 2016iuh, 2016ije, 2016brx, 2017avj, 2016eqa, 2016brv,

    # missing_ztf = ("2023jah, 2022zsp, 2022ywc, 2022vse, 2022an, 2020dq, 2020aejj, 2019so, 2017hle, 2017avj, 2016eqa").split(', ')
    # for n in missing_ztf:
    #     check_SN_in_ATLAS_ZTF(n)

    # queryZTF.download(tar_list="txts/target_files/sn1991bglike_tarlist_atlastest.csv", save_loc='data/ZTF-91bg/')
    # queryATLAS.download(tar_list="txts/target_files/sn1991bglike_tarlist_atlastest.csv", save_loc='data/ATLAS-91bg/')
    # queryZTF.download(tar_list="txts/target_files/sn1991bglike_tarlist_missingZTF.csv", save_loc='data/ZTF-91bg/')
    # queryATLAS.download(tar_list="txts/target_files/sn1991bglike_tarlist_missingATLAS.csv", save_loc='data/ATLAS-91bg/')



    # # Target file
    # from astropy.table import Table
    # tar_data = np.genfromtxt("txts/target_files/sn1991bglike_tarlist.csv", delimiter=',', dtype=str)
    # tar_tb = Table(names=tar_data[0, :], data=tar_data[1:, :])
    # tar_names = []
    # for n in tar_tb['Name']: tar_names.append(n[3:])
    #
    # # Old names - goal
    # old_params = utils.default_open('results/old/merged_params_cut.txt', True)
    # old_params.remove_row(np.where(old_params['objname'] == '2022ydr')[0][0])  # Remove the SNIa (2022ydr)
    # old_names = list(old_params['objname'])
    #
    # # Current data names
    # csp_files = glob.glob('data/CSP-91bg/*.txt')
    # atlas_files = glob.glob('data/ATLAS-91bg/*.txt')
    # ztf_files = glob.glob('data/ZTF-91bg/*.csv')
    # csp_f_names, atlas_f_names, ztf_f_names = [], [], []
    # for names, files in zip([csp_f_names, atlas_f_names, ztf_f_names], [csp_files, atlas_files, ztf_files]):
    #     for f in files:
    #         n = (str(f.split('/')[-1][:-4]).
    #              replace("CSP", "").replace("ATLAS", "").replace("ZTF", ""))
    #         names.append(n)
    # new_names = csp_f_names+atlas_f_names+ztf_f_names
    #
    # # Check
    # csp_c, atlas_c, ztf_c, atlas_ztf_c = 0, 0, 0, 0
    # for n in tar_names:
    #     if (n in csp_f_names):
    #         csp_c += 1
    #     if (n in atlas_f_names) and (n in ztf_f_names):
    #         atlas_ztf_c += 1
    #     if (n in atlas_f_names) and (n not in ztf_f_names):
    #         atlas_c += 1
    #     if (n in ztf_f_names) and (n not in atlas_f_names):
    #         ztf_c += 1
    # print(csp_c, atlas_c, ztf_c, atlas_ztf_c)

        # line = f"[+++] {n} in"
        # if n in csp_f_names:
        #     line += " CSP"
        # if n in atlas_f_names:
        #     line += " ATLAS"
        # if n in ztf_f_names:
        #     line += " ZTF"
        #
        # if line == f"[+++] {n} in":
        #     print(f"[!!!] No matches {n}")
        # else:
        #     print(line)

    # # Check
    # for n in tar_names:
    #     line = f"[+++] {n} in"
    #     if n in csp_f_names:
    #         line += " CSP"
    #     if n in atlas_f_names:
    #         line += " ATLAS"
    #     if n in ztf_f_names:
    #         line += " ZTF"
    #
    #     if line == f"[+++] {n} in":
    #         print(f"[!!!] No matches {n}")
    #     else:
    #         print(line)





    # sn = fitter.fit('')
    #

    # import numpy as np
    # failed = ("2023jah, 2022zsp, 2022ywc, 2022vse, 2022an, 2020dq, 2020aejj, 2019so, 2017hle, 2017avj, 2016eqa, "
    #           "2016brv, 2024zaj, 2024yhg, 2024ulb, 2024stt, 2024rkh, 2024pkh, 2024lhg, 2024bjb, 2023omo, 2023mkp, "
    #           "2023fot, 2023cvq, 2023abdv, 2022yv, 2022skw, 2022ihz, 2022fjx, 2022bsi, 2022aecb, 2022abom, 2021zsz, "
    #           "2021mab, 2021jvp, 2021jbp, 2021gel, 2021bls, 2021aare, 2020acoo, 2020aclx, 2019op, 2019exc, 2019cp, "
    #           "2019be, 2018jag, 2018hkw, 2018awi, 2009F, 2008bt, 2008bi, 2008bd, 2007N, 2007ba, 2007ax, 2007al, "
    #           "2006mr, 2006gt, 2006bd, 2005ke, 2005bl, 2023ex, 2018lph, 2017fzw, 2016iuh, 2016ije, 2016brx").split(', ')
    # tarlist = np.genfromtxt('txts/target_files/sn1991bglike_tarlist.csv', str, delimiter=',')
    # for n in failed:
    #     print(tarlist[np.where(tarlist[:, 0] == f"SN {n}")[0][0]][-1])
        # r = utils.query_tns(objname=n)
        # print(r)
        # break




    # # 2023jah, 2024yhg, 2007ba
    # ebv, st = [], []
    # ebv_e, st_e = [], []
    # mu, mu_e = [], []
    # Tmax, Tmax_e = [], []
    # for n in ["ATLAS2023jah", "ATLAS2024yhg", "CSP2007ba"]:
    #     sn = fitter.fit(f'data/combined_91bg/{n}.txt', 'snpy', rewrite=False)
    #     # sn[0].plot()
    #     # ebv.append(sn[0].params['EBVhost']['value'])
    #     # ebv_e.append(sn[0].params['EBVhost']['err'])
    #     # st.append(sn[0].params['st']['value'])
    #     # st_e.append(sn[0].params['st']['err'])
    #     # mu.append(sn[0].params['mu']['value'])
    #     # mu_e.append(sn[0].params['mu']['err'])
    #     Tmax.append(sn[0].params['Tmax']['value'])
    #     Tmax_e.append(sn[0].params['Tmax']['err'])
    #     # input('\nNext?\n')
    # # print(ebv)
    # # print(ebv_e)
    # # print(st)
    # # print(st_e)
    # print(Tmax)
    # print(Tmax_e)




    # import glob
    # import numpy as np
    # ztf_files = glob.glob('data/ZTF-91bg/*.csv')
    # og_tarfile = np.genfromtxt('txts/target_files/sn1991bglike_tarlist.csv', delimiter=',', dtype=str, skip_header=2)
    # for file in ztf_files:
    #     data = np.genfromtxt(file, delimiter=',', dtype=str)
    #     if len(data.shape) > 1:
    #         print(data[:, 0])
    #     # if len(data) < 10:
    #     #     name = f"SN {file.split('/')[-1][3:-4]}"
    #     #     name_pos = np.where(og_tarfile[:, 0] == name)[0][0]
    #     #     with open('temp_target_file.txt', 'a') as f:
    #     #         line = ''
    #     #         for c in og_tarfile[name_pos, :]:
    #     #             line += f"{c},"
    #     #         print(line[:-1], file=f)
    #
    # # queryZTF.download('temp_target_file.txt', 'data/ZTF-91bg/')




    # 2023fot
    # queryATLAS.download(tar_list="txts/target_files/sn1991bglike_tarlist_2023fot.csv", save_loc='data/ATLAS-91bg/')
    # dataManipulation.combine_like_data('txts/target_files/sn1991bglike_tarlist_2023fot.csv', 'data/combined_91bg/', '91bg', clear_old = False)
    # sn = fitter.fit('data/combined_91bg/2023fot.txt', 'snpy', rewrite=True)
    # print(sn[0].params['mu']['value'], sn[0].params['mu']['err'])

    # 2019moq
    # queryATLAS.download(tar_list="txts/target_files/sn1991bglike_tarlist_2019moq.csv", save_loc='data/ATLAS-91bg/')
    # dataManipulation.combine_like_data('txts/target_files/sn1991bglike_tarlist_2019moq.csv', 'data/combined_91bg/', '91bg', clear_old = False)
    # sn = fitter.fit('data/combined_91bg/2019moq.txt', 'snpy', rewrite=True)
    # sn[0].plot('flux')

    # 2021uef
    # queryATLAS.download(tar_list="txts/target_files/sn1991bglike_tarlist_2021uef.csv", save_loc='data/ATLAS-91bg/')
    # dataManipulation.combine_like_data('txts/target_files/sn1991bglike_tarlist_2021uef.csv', 'data/combined_91bg/', '91bg', clear_old = False)
    # sn = fitter.fit('data/combined_91bg/2021uef.txt', 'snpy', rewrite=True)
    # sn[0].plot('flux')




    # dataManipulation.selection_criteria('results/91bg_snpy-salt.txt')
    # compare_old_new_params(new_path = 'results/91bg_snpy-salt_cut.txt',
    #                        old_path = 'results/old/merged_params_cut.txt')
    # compare_old_new_params(new_path = 'results/91bg_snpy-salt.txt',
    #                        old_path = 'results/old/merged_params_cut.txt')
    # dataManipulation.selection_criteria('temp.txt', '91bg', 'temp_cut.txt')
    # compare_old_new_params(new_path='temp.txt',
    #                        old_path='results/old/merged_params_cut.txt')

    # Missing: ['2019be', '2020acbx', '2021bls', '2021mab', '2024jhk', '2018lph']
    # Just failed - snpy:                   2021bls, 2021mab, 2024jhk
    #             - salt: 2019be, 2020acbx, 2021bls, 2021mab, 2024jhk, 2018lph
    # sn = fitter.sneObj('combined_91bg', 'snpy', 'data/combined_91bg/2021mab.txt')
    # sn.plot('mag')
    # sn = fitter.fit('data/combined_91bg/2021bls.txt', 'snpy', rewrite=True)
    # sn = fitter.fit('data/combined_91bg/2021bls.txt', 'salt', rewrite=True)

    # dataManipulation.make_param_file(sn, 'n_temp.txt')

    # import glob
    # for file in glob.glob(f'data/combined_91bg/2*.txt'):
    #     sn = fitter.sneObj('combined_91bg', 'salt', file)
    #     sn.plot()
    #     input("\nNext?...\n")
    #     # break


    # sn = fitter.fit('data/combined_91bg/2020abmg.txt', 'salt', rewrite=True)
    # sn[0].plot()

    # # Create TNS Header
    # headers = {'Authorization': f'Token {utils.get_apikeys()["atlas_key"]}', 'Accept': 'application/json'}
    #
    # # # Get known names
    # # known_names = []
    # # for k in glob.glob(save_loc + "*.txt"):
    # #     known_names.append(k.split('/')[-1].split('ATLAS')[-1].split('.txt')[0])
    #
    # # Download light curves for each set of coordinates
    # failed_names = []
    # objnames, ra, dec, disc_date = ['2020abmg'], ['123.852270875'], ['-6.957704925'], ['2020-12-03 12:47:31.200']
    # for i, n in enumerate(zip(objnames, ra, dec, disc_date)):
    #     # hdr = f"[{i + 1} / {len(objnames)}] Downloading {n[0]} ({round(n[1], 3)}, {round(n[2], 3)})... ==================="
    #     # csr = f"{'=' * len(hdr)}"
    #     #######################
    #     # print(hdr)
    #
    #     # Check if already downloaded
    #     # if n[0] in known_names:
    #     #     print(f'[---] Already downloaded! Skipping...\n{csr}')
    #     #     continue
    #
    #     # Attempt download
    #     task_url = queryATLAS.initate_download(n[1], n[2], n[3], headers)
    #     result_url = queryATLAS.check_download(task_url, headers)
    #     if result_url is None:
    #         print(f"[!!!] Result unavailable! Skipping...\n{csr}")
    #         failed_names.append(n[0])
    #         continue
    #
    #     # Successful Download
    #     queryATLAS.save_download(result_url, headers, f'{save_loc}ATLAS{n[0]}.txt')
    #
    # # Report failed downloads
    # print(f"[~~~] The following resulting in errors...\n[", end='')
    # for n in failed_names:
    #     print(f"{n}", end=', ')
    # print(']\n')




    # tb_new = utils.default_open('results/SN1991bglikeParameters_snpy-salt_cut.txt', True)
    # tb_old = utils.default_open('results/old/merged_params_cut.txt', True)
    # tb_old_fix = utils.default_open('old_params_fix.txt', True)
    #
    # resid_fix = tb_old_fix['mu'].astype(float) - utils.current_cosmo().distmod(tb_old_fix['z_cmb'].astype(float)).value
    # # print(list(tb_old['objname']))
    #
    # dataManipulation.selection_criteria('results/SN1991bglikeParameters_snpy-salt.txt',
    #                                     '',
    #                                     'results/SN1991bglikeParameters_snpy-salt_cut.txt')
    # plotter.resid_v_mass('results/SN1991bglikeParameters_snpy-salt_cut.txt',
    #                      'results/old/aaronDo_salt2_params_cut.txt', label=True)



    # names = list(tb_old['objname'])
    # for i, n in enumerate(names):
    #     if i < 27: continue
    #     print(f"[{i+1} / {len(names)}] {n}: ", end='')
    #     r = utils.query_tns(objname=str(n))
    #     print(f"{r['object_type']['name']}")




    # dataManipulation.selection_criteria('results/91bg_snpy-salt.txt', 'norm', 'sn91bg_test.txt')
    # dataManipulation.selection_criteria('results/norm_snpy-salt.txt', 'norm', 'norm_test.txt')
    # plotter.resid_v_mass('sn91bg_test.txt', 'norm_test.txt', label = True)

    # brout_scholnic_expected_RMS()

    # print("Old")
    # tb = utils.default_open('results/old/merged_params.txt', True)
    # resid = tb['mu'].astype(float) - utils.current_cosmo().distmod(tb['z_cmb'].astype(float)).value
    # print(np.std(resid[abs(resid) < 0.4]))

    # print("New")
    # tb = utils.default_open('results/91bg_snpy-salt.txt', True)
    # resid = tb['mu'].astype(float) - utils.current_cosmo().distmod(tb['z_cmb'].astype(float)).value
    # # print(np.std(resid[abs(resid) < 0.4]))
    #
    # # print(len(np.array(abs(resid) < 0.4)))
    # tb_bad = tb[~(abs(resid) < 0.4)]
    # tb_best = tb[(abs(resid) < 0.4)]
    # tol = 0.05
    #
    # print(f"'z_cmb': [{min(tb_best['z_cmb'])-tol}, "
    #       f"{max(tb_best['z_cmb'])+tol}],")
    # print(f"'Tmax_err': [{min(tb_best['Tmax_err'])-tol}, "
    #       f"{max(tb_best['Tmax_err'])+tol}],")
    # print(f"'mu_err': [{min(tb_best['mu_err'])-tol}, "
    #       f"{max(tb_best['mu_err'])+tol}],")
    # print(f"'EBVhost': [{min(tb_best[tb_best['algo'] == 'SNPY']['color'])-tol}, "
    #       f"{max(tb_best[tb_best['algo'] == 'SNPY']['color'])+tol}],")
    # print(f"'EBVhost_err': [{min(tb_best[tb_best['algo'] == 'SNPY']['color_err'])-tol}, "
    #       f"{max(tb_best[tb_best['algo'] == 'SNPY']['color_err'])+tol}],")
    # print(f"'st': [{min(tb_best[tb_best['algo'] == 'SNPY']['stretch'])-tol}, "
    #       f"{max(tb_best[tb_best['algo'] == 'SNPY']['stretch'])+tol}],")
    # print(f"'st_err': [{min(tb_best[tb_best['algo'] == 'SNPY']['stretch_err'])-tol}, "
    #       f"{max(tb_best[tb_best['algo'] == 'SNPY']['stretch_err'])+tol}],")
    # print(f"'c': [{min(tb_best[tb_best['algo'] == 'SALT']['color'])-tol}, "
    #       f"{max(tb_best[tb_best['algo'] == 'SALT']['color'])+tol}],")
    # print(f"'c_err': [{min(tb_best[tb_best['algo'] == 'SALT']['color_err'])-tol}, "
    #       f"{max(tb_best[tb_best['algo'] == 'SALT']['color_err'])+tol}],")
    # print(f"'x1': [{min(tb_best[tb_best['algo'] == 'SALT']['stretch'])-tol}, "
    #       f"{max(tb_best[tb_best['algo'] == 'SALT']['stretch'])+tol}],")
    # print(f"'x1_err': [{min(tb_best[tb_best['algo'] == 'SALT']['stretch_err'])-tol}, "
    #       f"{max(tb_best[tb_best['algo'] == 'SALT']['stretch_err'])+tol}]")


    # # New cutting func ==========================================================================================
    # from astropy.table import Table, vstack
    # from astropy.stats import sigma_clip, sigma_clipped_stats
    # # criteria = {
    # #     'z_cmb': [-0.03506500785421299, 0.1298494505246801],
    # #     'Tmax_err': [-0.049855317069886954, 4.35280159722459],
    # #     'mu_err': [-0.04201354847310502, 0.7405862611109835],
    # #     'EBVhost': [-0.3212796775028634, 0.7554503028153857],
    # #     'EBVhost_err': [-0.033601535432013656, 0.18699488981537687],
    # #     'st': [0.21629154916757753, 1.1562290779572633],
    # #     'st_err': [-0.03930799632829331, 0.343594415861289],
    # #     'c': [0.04277177483489296, 1.2634417737451562],
    # #     'c_err': [-0.04763315383815581, 0.25434594036056063],
    # #     'x1': [-4.247329661115002, 0.1270202119216745],
    # #     'x1_err': [-0.03582847954485664, 1.7780040047944603]
    # # }
    # # criteria = {
    # #     'z_cmb':       [-0.0350, 0.1298],
    # #     'Tmax_err':    [-0.0498, 4.3528],
    # #     'mu_err':      [-0.0420, 0.7405],
    # #     'EBVhost':     [-0.3212, 0.7554],
    # #     'EBVhost_err': [-0.0336, 0.1869],
    # #     'st':          [0.2162, 1.1562],
    # #     'st_err':      [-0.0393, 0.3435],
    # #     'c':           [0.0427, 1.2634],
    # #     'c_err':       [-0.0476, 0.2543],
    # #     'x1':          [-4.2473, 0.1270],
    # #     'x1_err':      [-0.0358, 1.7780]
    # # }
    # criteria = {
    #     'z_cmb': [0.015, 999],
    #     'Tmax_err': [0, 1.0],
    #     'mu_err': [0, 0.2],
    #     'EBVhost': [-1, 1],
    #     'EBVhost_err': [0, 0.1],
    #     'st': [0, 1.5],
    #     'st_err': [0, 0.1],
    #     'c': [0, 1],
    #     'c_err': [0, 0.1],
    #     'x1': [-4.0, 0.0],
    #     'x1_err': [0, 0.1]
    # }
    # criteria = {
    #     'z_cmb':        [0.015, 999],
    #     'EBVhost':      [-0.2, 0.3],
    #     'EBVhost_err':  [-999, 0.1],
    #     'st':           [-999, 1.0],
    #     'st_err':       [-999, 0.1],
    #     'c':            [-0.6, 0.6],
    #     'c_err':        [-999, 0.1],
    #     'x1':           [-4.2, 3.2],
    #     'x1_err':       [-999, 1.0],
    #     'Tmax_err':     [-999, 1.0],
    #     't0_err':       [-999, 1.0],
    #     'mu_err':       [-999, 0.1],
    # }
    # # Load data
    # tb = utils.default_open('results/91bg_snpy-salt.txt', True)
    #
    # # Sigma Clipping
    # resid = (tb['mu'].astype(float) - utils.current_cosmo().distmod(tb['z_cmb'].astype(float)).value)
    # clipped_data = sigma_clip(resid, sigma=3, maxiters=5)
    # tb_clipped = Table(names=tb.colnames, dtype=[object] * len(tb.colnames))
    # for i in range(len(clipped_data)):
    #     if ~clipped_data.mask[i]:
    #         tb_clipped.add_row(tb[i])
    # resid_clipped = (tb_clipped['mu'].astype(float) - utils.current_cosmo().distmod(tb_clipped['z_cmb'].astype(float)).value)
    # print(f"[~~~] Sigma Clipping on Hubble Residual: σ = 3... SNe = {len(resid)} ---> {len(resid_clipped)} "
    #       f"(Hubble Residual Scatter = {round(np.std(resid), 3)} --> {round(np.std(resid_clipped), 3)})")
    # tb = tb_clipped.copy()
    #
    # # Basic cuts
    # for p in ['z_cmb', 'Tmax_err', 'mu_err']:
    #     print(f"{p}: {len(tb)} -> ", end='')
    #     tb = tb[(tb[p] > criteria[p][0]) & (tb[p] <= criteria[p][1])]
    #     print(f"{len(tb)} ", end='')
    #     print(f"({np.std(tb['mu'].astype(float) - utils.current_cosmo().distmod(tb['z_cmb'].astype(float)).value)})")
    #
    # # Seperate SNPY & SALT
    # tb_snpy = tb[tb['algo'] == 'SNPY']
    # tb_salt = tb[tb['algo'] == 'SALT']
    #
    # # Cut on SNPY params
    # for f_name, d_name in zip(['EBVhost', 'st'], ['color', 'stretch']):
    #     for n_end in ['', '_err']:
    #         print(f"{f_name+n_end}: {len(tb_snpy)+len(tb_salt)} -> ", end='')
    #         tb_snpy = tb_snpy[(tb_snpy[d_name+n_end] >= criteria[f_name+n_end][0]) &
    #                           (tb_snpy[d_name+n_end] <= criteria[f_name+n_end][1])]
    #         print(f"{len(tb_snpy)+len(tb_salt)} ", end='')
    #
    #         tb_temp = vstack([tb_snpy, tb_salt])
    #         print(f"({np.std(tb_temp['mu'].astype(float) - utils.current_cosmo().distmod(tb_temp['z_cmb'].astype(float)).value)})")
    #
    # # Cut on SALT params
    # for f_name, d_name in zip(['c', 'x1'], ['color', 'stretch']):
    #     for n_end in ['', '_err']:
    #         print(f"{f_name+n_end}: {len(tb_snpy)+len(tb_salt)} -> ", end='')
    #         tb_salt = tb_salt[(tb_salt[d_name+n_end] >= criteria[f_name+n_end][0]) &
    #                           (tb_salt[d_name+n_end] <= criteria[f_name+n_end][1])]
    #         print(f"{len(tb_snpy)+len(tb_salt)} ", end='')
    #
    #         tb_temp = vstack([tb_snpy, tb_salt])
    #         print(f"({np.std(tb_temp['mu'].astype(float) - utils.current_cosmo().distmod(tb_temp['z_cmb'].astype(float)).value)})")
    #
    # # Recombine SNPY & SALT
    # tb_combined = vstack([tb_snpy, tb_salt])


    # dataManipulation.selection_criteria('results/91bg_snpy-salt.txt')


    # print(f"{p}: {len(tb)} -> ", end='')
    # tb = tb[(tb['color'] > criteria['EBVhost'][0]) & (tb['color'] < criteria['EBVhost'][1])]
    # print(f"{len(tb)} ", end='')
    # print(f"({np.std(tb['mu'].astype(float) - utils.current_cosmo().distmod(tb['z_cmb'].astype(float)).value)})")



    # norm_path = 'results/norm_snpy-salt.txt'
    # # norm_path = 'results/old/aaronDo_salt2_params_cut.txt'
    # sn91bg_path = 'results/91bg_snpy-salt.txt'
    # # dataManipulation.selection_criteria(norm_path, 'norm_test.txt')
    # dataManipulation.selection_criteria(sn91bg_path, 'sn91bg_test.txt')
    # plotter.resid_v_mass(path_91bg='sn91bg_test.txt',
    #                      path_norm='norm_test.txt',
    #                      label=True)

    # criteria = {
    #     'z_cmb':        [0.015, 999],
    #     'EBVhost':      [-0.2, 0.3],
    #     'EBVhost_err':  [-999, 0.1],
    #     'st':           [-999, 1.0],
    #     'st_err':       [-999, 0.1],
    #     'c':            [-0.6, 0.6],
    #     'c_err':        [-999, 0.1],
    #     'x1':           [-3.2, 3.2],
    #     'x1_err':       [-999, 1.0],
    #     'Tmax':         [-999, 999],
    #     'Tmax_err':     [-999, 1.0],
    #     't0_err':       [-999, 1.0],
    #     'mu':           [-999, 999],
    #     'mu_err':       [-999, 0.2],
    #     'hostMass':     [-999, 999],
    #     'hostMass_err': [-999, 999],
    # }
    # sn_snpy = fitter.fit(data_loc=f'data/combined_91bg/2022vxf.txt', algo='snpy', rewrite=True)
    # print(sn_snpy[0].objname)
    # for p in sn_snpy[0].params:
    #     if p == "chisquare": continue
    #     print(f"{p}: {sn_snpy[0].params[p]['value']} +/- {sn_snpy[0].params[p]['err']} "
    #           f"({criteria[p][0]} < {p} < {criteria[p][1]}) "
    #           f"({criteria[p + '_err'][0]} < {p + '_err'} < {criteria[p + '_err'][1]})")


    # compare_old_new_params(new_path='results/manual_params.csv',)

    # dataManipulation.selection_criteria('results/manual_params.csv', 'results/manual_params_cut.csv')

    # sne_ztf_91bg = fitter.fit(data_loc='data/combined_91bg/ZTF*.txt', algo='snpy', rewrite=False)

    # r = utils.query_tns(coords=["174.61858", "20.52622"])
    # r = utils.query_tns(objname="2006bd")
    # print(r)

    # a = 'snpy'
    # for n in ['2005ke.txt', 'ATLAS2017fof.txt', 'ZTF2016brx.txt']:
    #     sn = fitter.fit(data_loc=f'data/combined_91bg/{n}', algo=a, rewrite=True)
    #     input("Next? ...")
    # sn = fitter.fit(data_loc=f'data/combined_91bg/2005ke.txt', algo=a, rewrite=True)
    # sn = fitter.fit(data_loc=f'data/combined_91bg/ATLAS2017fof.txt', algo=a, rewrite=True)
    # sn = fitter.fit(data_loc=f'data/combined_91bg/ZTF2016brx.txt', algo=a, rewrite=True)

    # # for n in ['2019be', '2020acbx', '2021mab', '2006bd']:
    # for n in ['2006bd']:
    #     sn_snpy = fitter.fit(data_loc=f'data/combined_91bg/{n}.txt', algo='snpy', rewrite=True)
    #     sn_salt = fitter.fit(data_loc=f'data/combined_91bg/{n}.txt', algo='salt', rewrite=True)
    #
    #     print("=======================================================")
    #     print(f"SNPY: {sn_snpy[0].params['mu']['value']} +/- {sn_snpy[0].params['mu']['err']}")
    #     print(f"SALT: {sn_salt[0].params['mu']['value']} +/- {sn_salt[0].params['mu']['err']}")
    #     print("=======================================================")
    #
    #     break

    # for n, a in zip(['2019be', '2020acbx', '2006bd'],
    #                 ['salt', 'salt', 'salt']):
    #     sn = fitter.fit(data_loc=f'data/combined_91bg/{n}.txt', algo=a, rewrite=True)
    #     dataManipulation.make_param_file(sn, f'{n}.txt')


    # print(f"{sn_2018jag[0].objname}, COMBINED, 91BG, SNPY, {sn_2018jag[0].coords[0]}, {sn_2018jag[0].coords[1]}, "
    #       f"{sn_2018jag[0].z}, {sn_2018jag[0].z_cmb}, {min(sn_2018jag[0].mjd)}, {max(sn_2018jag[0].mjd)}, "
    #       f"{sn_2018jag[0].params['mu']['value']}, {sn_2018jag[0].params['mu']['err']}, "
    #       f"{sn_2018jag[0].params['hostMass']['value']}, {sn_2018jag[0].params['hostMass']['err']},"
    #       f"{sn_2018jag[0].params['Tmax']['value']}, {sn_2018jag[0].params['Tmax']['value']}")

    # sn_2019be = fitter.fit(data_loc=f'data/combined_91bg/2018jag.txt', algo='salt', rewrite=True)
    # sn_2022aaok = fitter.fit(data_loc=f'data/combined_91bg/2022aaok.txt', algo='salt', rewrite=True)
    # sn_ATLAS2024jhk = fitter.fit(data_loc=f'data/combined_91bg/ATLAS2024jhk.txt', algo='snpy', rewrite=True)

    # objname	origin	subtype	algo	ra	dec	z	z_cmb	MJDs	MJDe	mu	mu_err	hostMass	hostMass_err	Tmax	Tmax_err	stretch	stretch_err	color	color_err	amplitude	amplitude_err	peak_mag	peak_mag_err

    # sne_ztf_91bg = fitter.fit(data_loc='data/combined_91bg/ZTF*.txt', algo='snpy', rewrite=False)
    # sn = fitter.sneObj('atlas-91bg', 'salt', f'data/combined_91bg/ATLAS2024rkh.txt')
    # sn = fitter.sneObj('ztf-91bg', 'salt', f'data/combined_91bg/ZTF2016brx.txt')
    # sn = fitter.sneObj('combined_91bg', 'salt', f'data/combined_91bg/2023acdv.txt')
    # sn = fitter.sneObj('atlas-91bg', 'salt', f'data/combined_91bg/ATLAS2024luo.txt')
    # sn.plot(y_type='mag')

    # import glob
    # algo = 'snpy'
    # # # ATLAS
    # # for p in glob.glob('data/combined_91bg/ATLAS*.txt'):
    # #     sn = fitter.sneObj('atlas-91bg', algo, p)
    # #     sn.plot('flux', fmt='-o')
    # #     input("Next? ")
    #
    # # # ZTF
    # # for p in glob.glob('data/combined_91bg/ZTF*.txt'):
    # #     sn = fitter.sneObj('ztf-91bg', algo, p)
    # #     sn.plot('flux', fmt='-o')
    # #     input("Next? ")
    #
    # # COMBINED
    # for p in glob.glob('data/combined_91bg/2*.txt'):
    #     sn = fitter.sneObj('combined_91bg', algo, p)
    #     sn.plot('flux', fmt='-o')
    #     input("Next? ")

    # import numpy as np
    # from astropy.table import Table
    # import matplotlib.pyplot as plt
    #
    # ztf_data = np.genfromtxt('data/ZTF-91bg/ZTF2022fjx.txt', skip_header=56, delimiter=' ', dtype=str)
    # with open('data/ZTF-91bg/ZTF2022fjx.txt') as f:
    #     for i in range(55): f.readline()
    #     ztf_hdr = f.readline()[1:].rstrip('\n').split(', ')
    # ztf_tbl = Table(data=ztf_data, names=ztf_hdr)
    #
    # atlas_data = np.genfromtxt('data/ATLAS-91bg/ATLAS2022fjx.txt', delimiter=',', dtype=str)
    # atlas_hdr, atlas_data = atlas_data[0, :], atlas_data[1:, :]
    # atlas_tbl = Table(data=atlas_data, names=atlas_hdr)
    #
    # atlas_mjd = atlas_tbl['MJD'].astype(float)
    # atlas_flux = atlas_tbl['uJy'].astype(float)
    # atlas_dflux = atlas_tbl['duJy'].astype(float)
    #
    # ztf_mjd = ztf_tbl[ztf_tbl['forcediffimflux'] != 'null']['jd'].astype(float) - 2400000.5
    # ztf_flux = ztf_tbl[ztf_tbl['forcediffimflux'] != 'null']['forcediffimflux'].astype(float)
    # ztf_dflux = ztf_tbl[ztf_tbl['forcediffimfluxunc'] != 'null']['forcediffimfluxunc'].astype(float)
    #
    # plt.errorbar(ztf_mjd, ztf_flux, yerr=ztf_dflux, color='green', fmt='o')
    # # plt.errorbar(atlas_mjd, atlas_flux, yerr=atlas_dflux, color='blue', fmt='o')
    # plt.show()

    # print(np.max(tbl[tbl['source'] == 'ATLAS']['mjd'].astype(float)))
    # print(np.max(tbl[tbl['source'] == 'ZTF']['mjd'].astype(float)))
    # print(
    #     np.max(tbl[tbl['source'] == 'ZTF']['mjd'].astype(float)) -
    #     np.max(tbl[tbl['source'] == 'ATLAS']['mjd'].astype(float))
    # )

    # print(sn.coords)

    # plotter.resid_v_mass('results/SN1991bglikeParameters_snpy-salt_cut.txt',
    #                      'results/old/aaronDo_salt2_params.txt')

    # import numpy as np
    # from astropy.table import Table
    # import glob
    #
    # tarlist = np.genfromtxt("txts/target_files/all_possible_91bgs.csv", delimiter=',', dtype=str, skip_header=1)
    # tb = Table(names=tarlist[0], data=tarlist[1:])
    #
    # downloaded = []
    # for p in glob.glob('data/ZTF-91bg/*.txt'): downloaded.append("SN " + p.split('/')[-1][3:-4])
    #
    # print(len(downloaded), len(tb['Name']))

    # queryATLAS.download(tar_list="txts/target_files/manual_ATLAS.csv", save_loc='data/ATLAS-91bg/')

    # result_url = 'https://fallingstar-data.com/forcedphot/static/results/job2293296.txt'
    # name = '2023jah'
    # headers = {'Authorization': f'Token {utils.get_apikeys()["atlas_key"]}', 'Accept': 'application/json'}
    # queryATLAS.save_download('result_url', headers, f'data/ATLAS-91bg/ATLAS{name}.txt')

    # plotter.resid_v_mass('results/SN1991bglikeParameters_snpy-salt_cut.txt',
    #                      'results/old/aaronDo_salt2_params.txt')

    # fitter.fit(f"data/combined_91bg/ATLAS2016ajf.txt", 'salt', True)
    # sn = fitter.sneObj('atlas-91bg', 'salt', "data/ATLAS-91bg/ATLAS2006bd.txt")
    # sn = fitter.sneObj('atlas-91bg', 'salt', f"data/combined_91bg/ATLAS{p[2:]}.txt")

    # import numpy as np
    # # paths = ("SN2024rkh, SN2020dq, SN2019so, SN2024luo, SN2016brv,<br/> SN2017dzs, SN2017fof, SN2017hle, SN2022uxl, "
    # #          "SN2015bo,<br/> SN2024iyx, SN2024zls, SN2022zsp, SN2022jou, SN2024wdg,<br/> SN2016gkd, SN2024eny, "
    # #          "SN2020aejj, SN2024jgw, SN2024zaj,<br/> SN2021oyx, SN2024yhg, SN2017avj, SN2018efn, SN2024vcj,<br/> "
    # #          "SN2016ajm, SN2009al, SN2022ywc, SN2022vse, "
    # #          "SN2024xhs,<br/>").replace('<br/>', '').replace(',', '').split(' ')
    # paths = ("SN2024rkh, SN2020dq, SN2019so, SN2016brv, SN2017fof,<br/> SN2018chl, SN2022uxl, SN2015bo, SN2022jou, "
    #          "SN2024wdg,<br/> SN2016gkd, SN2024eny, SN2020aejj, SN2018efn, SN2024vcj,<br/> SN2016ajm, SN2009al, "
    #          "SN2022ywc, SN2022vse").replace('<br/>', '').replace(',', '').split(' ')
    # for p in paths:
    #     print(f"====================\ndata/combined_91bg/ATLAS{p[2:]}.txt\n====================")
    #     # sn = fitter.sneObj('atlas-91bg', 'salt', f"data/combined_91bg/ATLAS{p[2:]}.txt")
    #     # try:
    #     #     sn.plot(y_type='flux')
    #     # except Exception as e:
    #     #     print(e)
    #     fitter.fit(f"data/combined_91bg/ATLAS{p[2:]}.txt", 'salt', True)
    #     # break
    #     input('Next? ')

    # fitter.fit(f"data/combined_91bg/ATLAS2016ajf.txt", 'salt', True)

    # no_cut = ['2018ame', '2018baz', '2018ciw', '2018hkw', '2019be', '2019bwi', '2019cp', '2019op', '2020acbx', '2020acoo', '2020fhs', '2020vae', '2021aare', '2021diu', '2021gel', '2021mab', '2021oyx', '2021pom', '2021twa', '2021xfm', '2021zsz', '2022aaok', '2022aecb', '2022bsi', '2022ihz', '2022oux', '2022rjs', '2022ydr', '2022yv', '2023abdv', '2023bhm', '2023dk', '2023ex', '2023fwb', '2023jah', '2023mkp', '2024fid', '2024jhk', '2017fzw', '2017hle', '2018gro', '2020abmg', '2020aclx', '2020ecn', '2020yo', '2022skw', '2022vxf', '2024luo']
    # cut = ['2008bd', '2018eyi', '2018jag', '2019be', '2019cp', '2019moq', '2019op', '2020acbx', '2020vae', '2021jvp', '2021mab', '2021uef', '2022aaok', '2022rjs', '2022ydr', '2023fwb', '2024jhk', '2018lph']
    #
    # all_names = []
    # for n in no_cut+cut:
    #     if n not in all_names:
    #         all_names.append(n)
    # all_names = sorted(all_names)
    #
    # import glob
    # downloaded_data = glob.glob('data/combined_91bg/*.txt')
    # downloaded_names = []
    # for d in downloaded_data:
    #     downloaded_names.append(
    #         d.split('/')[-1][:-4].replace('CSP', '').replace('ATLAS', '').replace('ZTF', '')
    #     )
    #
    # missing_from_downloaded = []
    # for n in all_names:
    #     if n not in downloaded_names:
    #         missing_from_downloaded.append(n)
    #         # break
    # print(missing_from_downloaded)
    # # ['2017fzw', '2019bwi', '2020fhs', '2021diu', '2021twa', '2022ihz', '2022ydr', '2023jah']
    return

if __name__ == '__main__':
    # main()
    experiment()
