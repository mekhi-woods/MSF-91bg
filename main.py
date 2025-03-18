# M. D. Woods
# 02/19/2025
import glob
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
    # Get SN1991bg-like SNe
    queryATLAS.download(tar_list="txts/target_files/sn1991bglike_tarlist.csv", save_loc='data/ATLAS-91bg/')
    # # Get normal SNe
    queryATLAS.download(tar_list="txts/target_files/normal_tarlist.csv", save_loc='data/ATLAS-norm/')
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
    # Combine UNCUT paramater file
    dataManipulation.combine_snpy_salt(f'results/{subtype}_snpy.txt',
                                       f'results/{subtype}_salt.txt',
                                       f'results/{subtype}_snpy-salt.txt')
    # Combine CUT paramater file
    dataManipulation.combine_snpy_salt(f'results/{subtype}_snpy_cut.txt',
                                       f'results/{subtype}_salt_cut.txt',
                                       f'results/{subtype}_snpy-salt_cut.txt')
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
    dataManipulation.make_param_file(sne=sne_csp+sne_atlas+sne_ztf+sne_combined,
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
def run_lc_checker():
    tb = utils.default_open("results/91bg_snpy-salt_cut.txt", True)
    for name, algo, source, subtype in zip(tb['objname'], tb['algo'], tb['origin'], tb['subtype']):
        n_o = f"fitting/{algo.lower()}-plots/{name}_lc.png"
        n_f = f"results/final_lc/{name}_lc.png"
        shutil.copyfile(n_o, n_f)
        print(f"[~~~] Copied {name}: '{n_o}' -> '{n_f}'...")
    return
def run_plotter(final_dir: str = 'plots/'):

    # OLDEST PARAMS
    pm_norm_salt_cut = 'results/old/aaronDo_salt2_params_cut.txt'
    pm_norm_salt_uncut = 'results/old/aaronDo_salt2_params.txt'
    pm_norm_snpy_cut = 'results/old/dr3_params.txt'  # Needs to be cut?, no mu
    pm_norm_snpy_uncut = 'results/old/dr3_params.txt'
    pm_norm_merged_cut = 'results/old/aaronDo_salt2_params_cut.txt'  # only contains salt fitted
    pm_norm_merged_uncut = 'results/old/aaronDo_salt2_params.txt'  # only contains salt fitted
    #
    # pm_91bg_salt_cut = 'results/old/combiend__salt_params_cut.txt'
    # pm_91bg_salt_uncut = 'results/old/combiend__salt_params.txt'
    # pm_91bg_snpy_cut = 'results/old/combiend__snpy_params_cut.txt'
    # pm_91bg_snpy_uncut = 'results/old/combiend__snpy_params.txt'
    # pm_91bg_merged_cut = 'results/old/merged_params_cut.txt'
    # pm_91bg_merged_uncut = 'results/old/merged_params.txt'

    # OLD PARAMS
    pm_norm_salt_cut = 'results/NormalParameters_salt_cut.txt'
    pm_norm_salt_uncut = 'results/NormalParameters_salt.txt'
    pm_norm_snpy_cut = 'results/NormalParameters_snpy_cut.txt'
    pm_norm_snpy_uncut = 'results/NormalParameters_snpy.txt'
    pm_norm_merged_cut = 'results/NormalParameters_snpy-salt_cut.txt'
    pm_norm_merged_uncut = 'results/NormalParameters_snpy-salt.txt'

    pm_91bg_salt_cut = 'results/SN1991bglikeParameters_salt_cut.txt'
    pm_91bg_salt_uncut = 'results/SN1991bglikeParameters_salt.txt'
    pm_91bg_snpy_cut = 'results/SN1991bglikeParameters_snpy_cut.txt'
    pm_91bg_snpy_uncut = 'results/SN1991bglikeParameters_snpy.txt'
    pm_91bg_merged_cut = 'results/SN1991bglikeParameters_snpy-salt_cut.txt'
    pm_91bg_merged_uncut = 'results/SN1991bglikeParameters_snpy-salt.txt'

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
    # plotter.resid_v_mass(path_91bg=pm_91bg_merged_cut,
    #                      path_norm=pm_norm_merged_cut,
    #                      save_loc=final_dir + 'resid_v_mass.png',
    #                      label = True)
    # plotter.mu_v_z(path_91bg=pm_91bg_merged_cut,
    #                 path_norm=pm_norm_merged_cut,
    #                 save_loc=final_dir+'mu_v_z.png',
    #                 label = False)

    # # Brout+Scolnic 2021 style Dust v. Scatter
    # plotter.dust_v_scatter(path_91bg = pm_91bg_salt_cut, path_norm = pm_norm_salt_cut, path_dust = pm_dust,
    #                        bin_num = 50, bin_bounds = [0.1, 6.3], hist_bins = 60,
    #                        label = True, save_loc = final_dir+'dust_v_scatter.png')
    #
    # # SALT3 Plots
    # plotter.alpha_beta(path_91bg=pm_91bg_salt_cut,
    #                     path_norm=pm_norm_salt_cut,
    #                     save_loc=final_dir+'alpha_beta.png')
    #
    # # Dust Plots
    # # plotter.dust_hist(path_91bg=pm_91bg_salt_cut, path_red_norm=pm_redNorms, path_dust=pm_dust,
    # #                    save_loc=final_dir+'dust_params.png')
    # plotter.resid_v_mass_dust(path_91bg=pm_91bg_merged_cut, path_norm=pm_norm_merged_cut, path_dust=pm_dust,
    #                            save_loc=final_dir+'dust_resid_v_mass.png')
    # # plotter.abs_mag_v_color(path_91bg=pm_91bg_salt_cut, path_norm=pm_norm_salt_cut, path_dust=pm_dust,
    # #                         save_loc=final_dir+'absMag_v_color.png')
    #
    # # Paramater Histograms
    # plotter.param_hist(snpy_91bg_path=pm_91bg_snpy_cut, snpy_norm_path=pm_norm_snpy_cut,
    #                    salt_91bg_path=pm_91bg_salt_cut, salt_norm_path=pm_norm_salt_cut,
    #                    save_loc=final_dir + 'param_hist_cut.png', line_type='median')
    # plotter.param_hist(snpy_91bg_path=pm_91bg_snpy_uncut, snpy_norm_path=pm_norm_snpy_uncut,
    #                    salt_91bg_path=pm_91bg_salt_uncut, salt_norm_path=pm_norm_salt_uncut,
    #                    save_loc=final_dir + 'param_hist_uncut.png', line_type='median')
    #
    # # Brout+Scolnic 2021 style Paramaters (Color & Stretch) v. Scatter -- just need to run for norms first
    # plotter.params_v_scatter(path_snpy_91bg=pm_91bg_snpy_uncut, path_snpy_norm=pm_norm_snpy_uncut,
    #                          path_salt_91bg=pm_91bg_salt_uncut, path_salt_norm=pm_norm_salt_uncut,
    #                          bin_nums=[[10, 10], [10, 10], [10, 10], [10, 10]],
    #                          bin_bounds=[[0.17, 1.074], [-0.299, 0.4], [-3.291, 2.838], [-0.368, 0.67]],
    #                          label=True, save_loc=final_dir + 'params_v_scatter.png')
    # plotter.color_v_scatter(path_snpy_91bg=pm_91bg_snpy_cut, path_snpy_norm=pm_norm_snpy_cut,
    #                         path_salt_91bg=pm_91bg_salt_cut, path_salt_norm=pm_norm_salt_cut,
    #                         bin_nums=[[10, 10], [10, 10]], bin_bounds=[[-0.3, 0.6], [-0.3, 0.6]], label=True,
    #                         save_loc=final_dir + 'color_v_scatter.png')

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
        pass
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
@run_time
def main():
    # # File acquisition
    # run_queryCSP()  # Get CSP data from CSP site
    # run_queryATLAS()  # Query ATLAS using target list of known RA & DEC
    # run_queryZTF()  # Query ZTF using target list of known RA & DEC

    # # File verification
    # utils.get_twomass()  # Checks for twomass velocities data (large file)
    # run_file_verification()  # Using target file, checks if SN exists in the proper directory
    #
    # # Data combination
    # run_data_combination()  # Combines data file if SN exsists in multiple surveys
    #
    # # Fit light curves
    # run_fitter(subtype = '91bg', algo = 'snpy', rewrite = True)
    # run_fitter(subtype = '91bg', algo = 'salt', rewrite = True)
    # run_fitter(subtype = 'norm', algo = 'snpy', rewrite = False)
    # run_fitter(subtype = 'norm', algo = 'salt', rewrite = False)
    #
    # # Combine algorithms
    # run_algo_combination('91bg')
    # run_algo_combination('norm')
    #
    # # Check lc of cut file
    # run_lc_checker()
    #
    # # Plot data
    # run_plotter()

    # # Outout stats
    # run_stats('readme')
    return
def experiment():

    import numpy as np
    failed = ("2023jah, 2022zsp, 2022ywc, 2022vse, 2022an, 2020dq, 2020aejj, 2019so, 2017hle, 2017avj, 2016eqa, "
              "2016brv, 2024zaj, 2024yhg, 2024ulb, 2024stt, 2024rkh, 2024pkh, 2024lhg, 2024bjb, 2023omo, 2023mkp, "
              "2023fot, 2023cvq, 2023abdv, 2022yv, 2022skw, 2022ihz, 2022fjx, 2022bsi, 2022aecb, 2022abom, 2021zsz, "
              "2021mab, 2021jvp, 2021jbp, 2021gel, 2021bls, 2021aare, 2020acoo, 2020aclx, 2019op, 2019exc, 2019cp, "
              "2019be, 2018jag, 2018hkw, 2018awi, 2009F, 2008bt, 2008bi, 2008bd, 2007N, 2007ba, 2007ax, 2007al, "
              "2006mr, 2006gt, 2006bd, 2005ke, 2005bl, 2023ex, 2018lph, 2017fzw, 2016iuh, 2016ije, 2016brx").split(', ')
    tarlist = np.genfromtxt('txts/target_files/sn1991bglike_tarlist.csv', str, delimiter=',')
    for n in failed:
        print(tarlist[np.where(tarlist[:, 0] == f"SN {n}")[0][0]][-1])
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
    main()
    # experiment()
