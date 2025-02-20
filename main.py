# M. D. Woods
# 02/19/2025
import time as systime

import numpy as np

import utils              # Utility functions for various common actions accross project; Here to be initiated for later.
import queryCSP           # Retrieves "Photometric Data release 3" tarball from Krisciunas et al. 2017 hosted on
                          # https://csp.obs.carnegiescience.edu/data, then seperates data using target file.
import queryATLAS         # Retrieves data from ATLAS server and stores it in data\ -- requires 'atlas_key' in 'txts/api_keys.txt'.
import queryZTF           # Retrieves data from ZTF server sends out request to be send to an email.
import fitter             # Uses data in data\ to create SNooPy & SALT fits.
import plotter            # Plots data found in files in results\.
import dataManipulation

def run_queryCSP():
    utils.verify_downloads(combined_tarlist_path="txts/target_files/tarlist_CSP-ATLAS_norm.csv",
                           source='CSP', subtype='norm')
    utils.verify_downloads(combined_tarlist_path="txts/target_files/tarlist_CSP-ATLAS-ZTF_91bg.csv",
                           source='CSP', subtype='91bg')
    queryCSP.download(save_loc='data/')  # location of CSP-norm & CSP-91bg
    queryCSP.seperate_data(tar_list='txts/target_files/tarlist_CSP_91bg.csv',
                           data_loc='data/',
                           save_loc='data/CSP-91bg')  # Seperate 1991bg-like out of all CSP-I data
    queryCSP.seperate_data(tar_list='txts/target_files/tarlist_CSP_norm.csv',
                           data_loc='data/',
                           save_loc='data/CSP-norm')  # Seperate Normals out of all CSP-I data
    return
def run_queryATLAS():
    utils.verify_downloads(combined_tarlist_path="txts/target_files/tarlist_CSP-ATLAS-ZTF_91bg.csv",
                           source='ATLAS', subtype='91bg')
    utils.verify_downloads(combined_tarlist_path="txts/target_files/tarlist_CSP-ATLAS_norm.csv",
                           source='ATLAS', subtype='norm')
    queryATLAS.download(tar_list="txts/target_files/tarlist_CSP-ATLAS-ZTF_91bg.csv", save_loc='data/ATLAS-91bg/',
                        combined_tarlist=True)
    queryATLAS.download(tar_list="txts/target_files/tarlist_CSP-ATLAS_norm.csv", save_loc='data/ATLAS-norm/',
                        combined_tarlist=True)
    return
def run_queryZTF():
    # utils.verify_downloads(combined_tarlist_path="txts/target_files/tarlist_CSP-ATLAS-ZTF_91bg.csv",
    #                        source='ZTF', subtype='91bg')
    # queryZTF.download(tns_targets_path = "txts/target_files/tarlist_ZTF_91bg.csv", save_loc = 'data/ZTF-91bg/')
    # queryZTF.download(tns_targets_path = "txts/target_files/tarlist_ZTF_norm.csv", save_loc = 'data/ZTF-norm/')
    return
def run_dataManipulation():
    # dataManipulation.selection_criteria(path='results/old/merged_params.txt', save_loc='results/old/merged_params_cut.txt')
    # dataManipulation.selection_criteria(path='results/91bg_csp+atlas+ztf_snpy.txt')
    # dataManipulation.selection_criteria(path='results/old/merged_params.txt', save_loc='results/old/merged_params_cut.txt')
    # dataManipulation.selection_criteria(path='results/old/combiend__snpy_params.txt')
    # dataManipulation.selection_criteria(path='results/old/combiend__salt_params.txt')
    return
def run_fitter(algo: str = 'snpy', rewrite: bool = True):
    """
    :param algo: Algorithm to use for fitting
    :param rewrite: Whether to rewrite the fit data or not
    :return: None
    """
    # All 91bg-like Batch Run Prompts
    sne_csp_91bg = fitter.fit(data_loc = 'data/CSP-91bg/*.txt', algo = algo, rewrite = rewrite)
    sne_atlas_91bg = fitter.fit(data_loc = 'data/ATLAS-91bg/*.txt', algo = algo, rewrite = rewrite)
    sne_ztf_91bg = fitter.fit(data_loc = 'data/ZTF-91bg/*.txt', algo = algo, rewrite = rewrite)
    dataManipulation.make_param_file(sne=sne_csp_91bg+sne_atlas_91bg+sne_ztf_91bg,
                                     save_loc=f'results/SN1991bglikeParameters_{algo}.txt')
    dataManipulation.selection_criteria(f'results/SN1991bglikeParameters_{algo}.txt',
                                        f'results/SN1991bglikeParameters_{algo}_cut.txt')

    # All Normal Batch Run Prompts
    sne_csp_norm = fitter.fit(data_loc = 'data/CSP-norm/*.txt', algo = algo, rewrite = rewrite)
    sne_atlas_norm = fitter.fit(data_loc = 'data/ATLAS-norm/*.txt', algo = algo, rewrite = rewrite)
    dataManipulation.make_param_file(sne=sne_csp_norm+sne_atlas_norm,
                                     save_loc=f'results/NormalParameters_{algo}.txt')
    dataManipulation.selection_criteria(f'results/NormalParameters_{algo}.txt',
                                        f'results/NormalParameters_{algo}_cut.txt')
    return
def run_plotter():
    final_dir = 'plots/'

    pm_norm_salt_cut = 'results/old/aaronDo_salt2_params_cut.txt'
    pm_norm_salt_uncut = 'results/old/aaronDo_salt2_params.txt'
    pm_norm_snpy_cut = 'results/old/dr3_params.txt'  # Needs to be cut?, no mu
    pm_norm_snpy_uncut = 'results/old/dr3_params.txt'
    pm_norm_merged_cut = 'results/old/aaronDo_salt2_params_cut.txt'  # only contains salt fitted
    pm_norm_merged_uncut = 'results/old/aaronDo_salt2_params.txt'  # only contains salt fitted

    pm_91bg_salt_cut = 'results/old/combiend__salt_params_cut.txt'
    pm_91bg_salt_uncut = 'results/old/combiend__salt_params.txt'
    pm_91bg_snpy_cut = 'results/old/combiend__snpy_params_cut.txt'
    pm_91bg_snpy_uncut = 'results/old/combiend__snpy_params.txt'
    pm_91bg_merged_cut = 'results/old/merged_params_cut.txt'
    pm_91bg_merged_uncut = 'results/old/merged_params.txt'

    pm_redNorms = 'results/old/redNormSNe.txt'
    pm_dust = 'results/old/global_dust_params.txt'

    # plotter.resid_v_mass(path_91bg=pm_91bg_merged_cut,
    #                       path_norm=pm_norm_merged_cut,
    #                       save_loc=final_dir+'resid_v_mass.png',
    #                       label = False)
    # plotter.mu_v_z(path_91bg=pm_91bg_merged_cut,
    #                 path_norm=pm_norm_merged_cut,
    #                 save_loc=final_dir+'mu_v_z.png',
    #                 label = False)
    #
    # ## SALT3 Plots
    # plotter.alpha_beta(path_91bg=pm_91bg_salt_cut,
    #                     path_norm=pm_norm_salt_cut,
    #                     save_loc=final_dir+'alpha_beta.png')
    # plotter.abs_mag_v_color(path_91bg=pm_91bg_salt_cut, path_red_norm=pm_redNorms, path_dust=pm_dust,
    #                         save_loc=final_dir+'absMag_v_color.png')

    ## Dust Plots
    # plotter.dust_hist(path_91bg=pm_91bg_salt_cut, path_red_norm=pm_redNorms, path_dust=pm_dust,
    #                    save_loc=final_dir+'dust_params.png')
    # plotter.resid_v_mass_dust(path_91bg=pm_91bg_merged_cut, path_norm=pm_norm_merged_cut, path_dust=pm_dust,
    #                            save_loc=final_dir+'dust_resid_v_mass.png')
    ## Brout+Scolnic 2021 style Dust v. Scatter
    plotter.dust_v_scatter(path_91bg = pm_91bg_merged_cut, path_norm = pm_redNorms, path_dust = pm_dust,
                           bin_num = 50, bin_bounds = [0.1, 6.3], hist_bins = 60,
                           label = True, save_loc = final_dir+'dust_v_scatter.png')

    # ## Paramater Histograms
    # plotter.param_hist(snpy_91bg_path=pm_91bg_snpy_uncut, snpy_norm_path=pm_norm_snpy_uncut,
    #                    salt_91bg_path=pm_91bg_salt_uncut, salt_norm_path=pm_norm_salt_uncut,
    #                    save_loc=final_dir + 'param_hist_uncut.png', line_type='median')
    # plotter.param_hist(snpy_91bg_path=pm_91bg_snpy_cut, snpy_norm_path=pm_norm_snpy_cut,
    #                    salt_91bg_path=pm_91bg_salt_cut, salt_norm_path=pm_norm_salt_cut,
    #                    save_loc=final_dir + 'param_hist_cut.png', line_type='median')
    #
    # ## Brout+Scolnic 2021 style Paramaters (Color & Stretch) v. Scatter
    # plotter.params_v_scatter(path_snpy_91bg=pm_91bg_snpy_cut, path_snpy_norm='results/old/norm_snpy_params_cut.txt',
    #                          path_salt_91bg=pm_91bg_salt_cut, path_salt_norm=pm_norm_salt_cut,
    #                          bin_nums=[[10, 10], [10, 10], [10, 10], [10, 10]],
    #                          bin_bounds=[[0.17, 1.074], [-0.299, 0.4], [-3.291, 2.838], [-0.368, 0.67]],
    #                          label=True, save_loc=final_dir + 'params_v_scatter.png')
    # plotter.color_v_scatter(path_snpy_91bg=pm_91bg_snpy_cut, path_snpy_norm='results/old/norm_snpy_params_cut.txt',
    #                         path_salt_91bg=pm_91bg_salt_cut, path_salt_norm=pm_norm_salt_cut,
    #                         bin_nums=[[10, 10], [10, 10]], bin_bounds=[[-0.3, 0.6], [-0.3, 0.6]], label=True,
    #                         save_loc=final_dir + 'color_v_scatter.png')
    #

    return

if __name__ == '__main__':
    start = systime.time()  # Runtime tracker
    # utils.get_twomass()  # Checks for twomass velocities data (large file)
    # run_queryCSP()
    # run_queryATLAS()
    # run_queryZTF()
    # run_dataManipulation()

    # run_fitter(algo = 'snpy', rewrite = False)
    # run_fitter(algo = 'salt', rewrite = False)
    # dataManipulation.combine_snpy_salt('results/SN1991bglikeParameters_snpy.txt',
    #                                    'results/SN1991bglikeParameters_salt.txt',
    #                                    'results/SN1991bglikeParameters_snpy-salt.txt')
    # dataManipulation.combine_snpy_salt('results/SN1991bglikeParameters_snpy.txt',
    #                                    'results/SN1991bglikeParameters_salt.txt',
    #                                    'results/SN1991bglikeParameters_snpy-salt.txt')
    # run_plotter()

    # import os
    #
    # valid = []
    # for sn_name in ("SN2023mkp, SN2021fnr, SN2021zsz, SN2023fot, SN2021aare, SN2021agej, SN2024pbd, SN2020abmg, "
    #                 "SN2021abzd, SN2022aecb, SN2023dk, SN2023abdv, SN2022ywc, SN2021agnf, SN2024ulb, SN2024xhs, "
    #                 "SN2021wzb, SN2021bls, SN2021jbp, SN2024sag, SN2021gel, SN2021bmu").split(', '):
    #     dataset = 'atlas-91bg'
    #     path = f'data/ATLAS-91bg/ATLAS{sn_name[2:]}.txt'
    #     algo = 'snpy'
    #     class_loc = f"classes/{dataset}/{dataset}_{path.split('/')[-1].split('.txt')[0]}_{algo}_class.txt"
    #     sn_obj = fitter.sneObj('class', algo, class_loc)
    #     if len(np.unique(sn_obj.filter)) == 1:
    #         pass
    #     elif len(sn_obj.mag[sn_obj.filter == np.unique(sn_obj.filter)[0]]) > 5 and len(sn_obj.mag[sn_obj.filter == np.unique(sn_obj.filter)[1]]) > 5:
    #         valid.append(f"{sn_name}+{len(sn_obj.mag[sn_obj.filter == np.unique(sn_obj.filter)[0]])}+{len(sn_obj.mag[sn_obj.filter == np.unique(sn_obj.filter)[1]])}")
    #
    #     print(f"{sn_name} filters: {np.unique(sn_obj.filter)} ====================")
    #     for f in np.unique(sn_obj.filter):
    #         print(f"| {f}: {len(sn_obj.mag[sn_obj.filter == f])}", end=' ')
    #     print('|\n ================================================')
    #
    #
    # print(valid)
    #
    #

    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
