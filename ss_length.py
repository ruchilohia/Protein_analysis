from itertools import groupby
import numpy as np
from scipy.spatial.distance import pdist
import scipy.cluster.hierarchy as hac
import numpy as np
from sklearn.decomposition import PCA
from scipy.spatial.distance import squareform
import fabfile
import subprocess
import math
import random
from sklearn.cluster import KMeans
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import seaborn as sns
matplotlib.use('Agg')
from matplotlib.colors import LinearSegmentedColormap

import matplotlib.colors as colors


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    function taken from
    https://stackoverflow.com/questions/7404116/...
        ...defining-the-midpoint-of-a-colormap-in-matplotlib
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {  'red': [],  'green': [], 'blue': [],  'alpha': []  }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


def get_my_base(base):

    class DPCA(base):
        """
        Does PCA and makes clusters and finds contacts
        """

        @fabfile._decorator_calc_loop
        def calc_dihedral_angle(self, func_dir_name):
                """calculates the phi psi and omega dihedral angles at every residue
                # this needs to be run only once at all required temperatures"""
                self.population_calc()
                self._get_angle_index()

                return

        def get_a(self, func_dir_name, psi, phi):
                """reads phi and psi values and stores them as integers depending on the SS structure"""

                inp_cov_1 = psi
                inp_cov_1[inp_cov_1 > 50] = 500
                inp_cov_1[inp_cov_1 < -120] = 500
                inp_cov_1[inp_cov_1 != 500] = 1  # helix
                inp_cov_1[inp_cov_1 == 500] = -1  # beta

                inp_cov_phi = phi
                inp_cov_phi[inp_cov_phi >= 0] = 500  # L-helix
                inp_cov_phi[inp_cov_phi < -90] = 300  # beta
                # the remaining region is pp2 region

                my_matrix = np.empty(inp_cov_1.shape, dtype=int)
                for (x, y), value in np.ndenumerate(my_matrix):
                    if inp_cov_1[x, y] == 1:
                        if inp_cov_phi[x, y] != 500:
                            my_matrix[x, y] = 1  # helix
                        else:
                            my_matrix[x, y] = 3  # lh
                    else:
                        if inp_cov_phi[x, y] == 500:
                            my_matrix[x, y] = 3  # lh
                        elif inp_cov_phi[x, y] == 300:
                            my_matrix[x, y] = -1  # beta
                        else:
                            my_matrix[x, y] = -2  # pp2 , instead of -2, I am writing pp2 as beta

                np.savetxt('%s/my_matrix.out' %
                           self.repl_dir_parent, my_matrix, fmt='%d')

                return ('%s/my_matrix.out' % self.repl_dir_parent)

        def ss_stretch(self, func_dir_name):
            """writes the length of ss_strch formed at each residue id"""
            for form in self.form_list:
                for i in range(0, self.calc_replica):
                    self.i = i
                    for s in range(self.parent_cluster_no):
                        self.parent_cluster_id = s
                        self.form = form
                        self.repl_dir_parent = self.repl_dir_parent_method(func_dir_name)
                        self.population_calc()
                        self.load_angle_values()
                        _reslist = range(self.resid_st, self.resid_end)
                        inp_cov_1 = pd.read_csv((self.get_a(func_dir_name, self.psi, self.phi)), sep=" ", names=[str(i) for i in _reslist])  # in these dataframes, the time frame is the index
                        ex_beta = open("%s/ex_beta.out" % self.repl_dir_parent, "w")
                        ex_lh = open("%s/ex_lh.out" % self.repl_dir_parent, "w")
                        ex_pp2 = open("%s/ex_pp2.out" % self.repl_dir_parent, "w")
                        ex_helix = open("%s/ex_helix.out" % self.repl_dir_parent, "w")

                        for my_file in [ex_beta, ex_lh, ex_pp2, ex_helix]:
                            my_file.write(" ".join(map(str, _reslist)))
                            my_file.write("\n")

                        for index, time_frame in inp_cov_1.iterrows():

                            for bit, group in groupby(time_frame):

                                y = list(group)  # breaks down the time frame into lists of consecutive identical numbers
                                for i in range(len(y)):
                                    if bit == 1:  # helix
                                        ex_helix.write(str(len(y)) + " ")
                                        ex_beta.write(str(0) + " ")
                                        ex_lh.write(str(0) + " ")
                                        ex_pp2.write(str(0) + " ")

                                    elif bit == -1:  # beta
                                        ex_beta.write(str(len(y)) + " ")
                                        ex_helix.write(str(0) + " ")
                                        ex_lh.write(str(0) + " ")
                                        ex_pp2.write(str(0) + " ")

                                    elif bit == -2:  # ppII
                                        ex_pp2.write(str(len(y)) + " ")
                                        ex_helix.write(str(0) + " ")
                                        ex_lh.write(str(0) + " ")
                                        ex_beta.write(str(0) + " ")

                                    else:  # lh helix
                                        ex_beta.write(str(0) + " ")
                                        ex_helix.write(str(0) + " ")
                                        ex_lh.write(str(len(y)) + " ")
                                        ex_pp2.write(str(0) + " ")

                            ex_helix.write("\n")
                            ex_beta.write("\n")
                            ex_lh.write("\n")
                            ex_pp2.write("\n")

                        ex_helix.close()
                        ex_beta.close()
                        ex_lh.close()
                        ex_pp2.close()

            return self

        def generate_data_frame(self, func_dir_name, types, clusname):
            """creates combined python dataframe"""

            _reslist = range(self.resid_st, self.resid_end)
            for i in ['form', 'replica_id', 'cluster_index', 'time_frame', 'Temperature', 'rgyr']:
                _reslist.append(i)

            helix = pd.DataFrame(columns=map(str, _reslist))

            for form in self.form_list:
                # for i in range(0, self.calc_replica):
                for i in range(0, self.calc_replica):
                    self.i = i
                    for s in range(self.parent_cluster_no):
                        self.parent_cluster_id = s
                        self.form = form
                        self.population_calc(clusname)
                        self.repl_dir_parent = self.repl_dir_parent_method(func_dir_name)

                        _helix = pd.read_csv("%s/ex_%s.out" % (self.repl_dir_parent, types), sep=" ", index_col=False)
                        rgyr_path = ("/".join(map(str, [self.output, self.form, "rgyr", str(self.i), "gyrate.xvg"])))
                        _helix['rgyr'] = np.loadtxt(rgyr_path)[:, 1]
                        _helix['form'] = self.form
                        _helix['replica_id'] = self.demux_file_per_temperature()
                        _helix['cluster_index'] = self.clus_list_index
                        _helix['time_frame'] = self.time_frame
                        _helix['Temperature'] = self.curr_temp_list
                        helix = helix.append(_helix)
            helix.to_csv("%s/%s/%s/pandas_format_%s" % (self.output, self.form, func_dir_name, types), sep='\t')
            return helix

        def generate_data_frame_66_cluster(self, func_dir_name, clusname="None"):
            """creates combined python dataframe"""

            for form in self.form_list:
                _reslist = ['helix', 'beta', 'lh', 'pp2']
                for i in ['form', 'replica_id', 'cluster_index', 'time_frame', 'Temperature', 'rgyr']:
                    _reslist.append(i)

                helix = pd.DataFrame(columns=map(str, _reslist))
                # for i in range(0, self.calc_replica):
                for i in range(0, self.calc_replica):
                    self.i = i
                    for s in range(self.parent_cluster_no):
                        self.parent_cluster_id = s
                        self.form = form
                        self.population_calc(clusname)
                        self.repl_dir_parent = self.repl_dir_parent_method(func_dir_name)
                        rgyr_path = ("/".join(map(str, [self.output, self.form, "rgyr", str(self.i), "gyrate.xvg"])))
                        _helix = pd.read_csv("%s" % (rgyr_path), sep="\t\t|\t|\s+|", index_col=False, names=['rgyr'], engine='python', usecols=[1])
                        _helix['helix'] = pd.read_csv("%s/ex_helix.out" % (self.repl_dir_parent), sep=" ", index_col=False)['66']
                        _helix['beta'] = pd.read_csv("%s/ex_beta.out" % (self.repl_dir_parent), sep=" ", index_col=False)['66']
                        _helix['lh'] = pd.read_csv("%s/ex_lh.out" % (self.repl_dir_parent), sep=" ", index_col=False)['66']
                        _helix['pp2'] = pd.read_csv("%s/ex_pp2.out" % (self.repl_dir_parent), sep=" ", index_col=False)['66']
                        _helix['form'] = self.form
                        _helix['replica_id'] = self.demux_file_per_temperature()
                        _helix['cluster_index'] = self.clus_list_index
                        _helix['time_frame'] = self.time_frame
                        _helix['Temperature'] = self.curr_temp_list
                        helix = helix.append(_helix)
                helix.to_csv("%s/%s/%s/pandas_format_cluster_66" % (self.output,  form, func_dir_name), sep='\t')
            return helix

        def hist_length_plot(self, func_dir_name, _ncol=1, _nrow=1, types="helix", c_max=5, st_time=0, clusname="None"):
                """heat map of residue id vs ss structure length at every temperature"""
                input_file = self.generate_data_frame(func_dir_name, types, clusname)
                st_frame = st_time / self.xtc_freq
                helix = input_file.loc[input_file['time_frame'] >= st_frame]
                reslist = range(self.resid_st, self.resid_end)
                reslist_min = np.amin(reslist)
                reslist_max = np.amax(reslist)
                bin_length_y = np.amax(reslist) - np.amin(reslist)
                helix_length_range = [helix[map(str, reslist)].values.min() + 2, helix[map(str, reslist)].values.max() + 1]
                bin_length_x = max(helix_length_range) - min(helix_length_range)

                # Particular form for all the given temeprature together
                grouped_form = helix.groupby(['form'])

                ncol = 2
                nrow = int(np.ceil(grouped_form.ngroups / ncol))

                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)

                for (key, ax) in zip(grouped_form.groups.keys(), axes.flatten()):
                            form_input = grouped_form.get_group(key)

                            Ht, xedges, yedges = np.histogram2d((form_input[map(str, reslist)]).astype(float).values.flatten(), (reslist * len(form_input)), bins=[bin_length_x, bin_length_y], range=[helix_length_range, [reslist_min, reslist_max]], normed=False, weights=None)
                            #np.savetxt("tmp.csv", Ht)
                            #tmp = pd.read_csv("tmp", sep="\t\t|\t|\s+|", index_col=False,, engine='python', usecols=[1])
                            H = Ht * 100.0 / len(form_input.index)
                            H[H == 0] = 'Nan'
                            figpath = self.fig_path_method(func_dir_name)
                            figname = "%s/%s.pdf" % (figpath, types)
                            imgp = ax.imshow(H.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', vmax=3, vmin=0, cmap=plt.get_cmap('Blues'))
                            cbar = ax.figure.colorbar(imgp, ax=ax)
                            # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
                            # self.ax.set_yticks(np.arange(reslist_min, reslist_max, 10))
                            # self.ax.set_xticks(np.arange(helix_length_range[0], helix_length_range[1], 2))
                            fig.set_size_inches(ncol * 8, nrow * 8)
                            ax.set_title(key)

                plt.savefig(figname, bbox_inches='tight', dpi=100)

                # Particular form for all the given temeprature together
                grouped_form = helix.groupby(['form'])

                ncol = 2
                nrow = int(np.ceil(grouped_form.ngroups / ncol))

                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)

                for (key, ax) in zip(grouped_form.groups.keys(), axes.flatten()):
                            form_input = grouped_form.get_group(key)
                            print form_input.head()
                            #print form_input.groupby(['24']).get_group(3)
                            #print form_input.groupby(['24'])['replica_id'].nunique()
                            #print form_input.groupby(['24'])['replica_id'].nunique()
                            #print form_input.groupby(['24', ])['replica_id'].nunique()
                            #pd_xyz = form_input.groupby(['24']).agg({"24": 'count' , "replica_id": lambda x: x.nunique()}).reset_index()
                            #print pd_xyz

                            #pd_xyz['count'] = pd_xyz['count'].transform(lambda x: x * 100.0 / y['time_frame'].nunique())

                            Ht, xedges, yedges = np.histogram2d((form_input[map(str, reslist)]).astype(float).values.flatten(), (reslist * len(form_input)), bins=[bin_length_x, bin_length_y], range=[helix_length_range, [reslist_min, reslist_max]], normed=False, weights=None)
                            H = Ht * 100.0 / len(form_input.index)
                            H[H == 0] = 'Nan'
                            figpath = self.fig_path_method(func_dir_name)
                            figname = "%s/%s.pdf" % (figpath, types)
                            imgp = ax.imshow(H.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', vmax=3, vmin=0, cmap=plt.get_cmap('Blues'))
                            cbar = ax.figure.colorbar(imgp, ax=ax)
                            # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
                            # self.ax.set_yticks(np.arange(reslist_min, reslist_max, 10))
                            # self.ax.set_xticks(np.arange(helix_length_range[0], helix_length_range[1], 2))
                            fig.set_size_inches(ncol * 8, nrow * 8)
                            ax.set_title(key)

                plt.savefig(figname, bbox_inches='tight', dpi=100)


                fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 4), sharey=True)
                form_input = grouped_form.get_group(self.form_list[0])
                print grouped_form.groups.keys()[0]
                if types == 'helix':
                    len_min = 8
                    max_helix = 19
                if types == 'beta':
                    len_min = 5
                    max_helix = 12
                Ht, xedges, yedges = np.histogram2d((reslist * len(form_input)), (form_input[map(str, reslist)]).astype(float).values.flatten(), bins=[bin_length_y, [1, len_min, max_helix]], normed=False, weights=None)
                H_0 = Ht * 100.0 / len(form_input.index)

                form_input = grouped_form.get_group(self.form_list[1])
                Ht, xedges, yedges = np.histogram2d((reslist * len(form_input)), (form_input[map(str, reslist)]).astype(float).values.flatten(), bins=[bin_length_y, [1, len_min, max_helix]], normed=False, weights=None)
                H_1 = Ht * 100.0 / len(form_input.index)

                H = H_0 - H_1
                H[H == 0] = 'Nan'

                figpath = self.fig_path_method(func_dir_name)
                figname = "%s/%s_diff_cum.pdf" % (figpath, types)
                vmax= 4.0
                #cmap = LinearSegmentedColormap.from_list('mycmap', [(0 / vmax, 'darkorange'), (2 / vmax, 'white'), (4 / vmax, 'steelblue')])
                #cmap = LinearSegmentedColormap.from_list('mycmap', [(0 / vmax, 'g'), (2 / vmax, 'white'), (4 / vmax, 'darkorange')])
                cmap = LinearSegmentedColormap.from_list('mycmap', [(0 / vmax, 'red'), (2 / vmax, 'white'), (4 / vmax, 'steelblue')])
                #cmap = LinearSegmentedColormap.from_list('mycmap', [(0 / vmax, 'g'), (2 / vmax, 'white'), (4 / vmax, 'red')])
                imgp = ax.imshow(H.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', vmax=4, vmin=-4, cmap=cmap)
                cbar = ax.figure.colorbar(imgp, ax=ax)
                # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
                # self.ax.set_yticks(np.arange(reslist_min, reslist_max, 10))
                # self.ax.set_xticks(np.arange(helix_length_range[0], helix_length_range[1], 2))
                ax.set_yticks([1, max_helix])
                fig.set_size_inches(16, 4)
                # ax.set_title(key)
                plt.savefig(figname, bbox_inches='tight', dpi=100)

                fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 4), sharey=True)
                form_input = grouped_form.get_group(grouped_form.groups.keys()[0])
                print grouped_form.groups.keys()[0]
                if types == 'helix':
                    len_min = 4
                if types == 'beta':
                    len_min = 5
                Ht, xedges, yedges = np.histogram2d((reslist * len(form_input)), (form_input[map(str, reslist)]).astype(float).values.flatten(), bins=[bin_length_y, [1, len_min, max(helix_length_range)]], normed=False, weights=None)
                H_0 = Ht * 100.0 / len(form_input.index)

                form_input = grouped_form.get_group(grouped_form.groups.keys()[1])
                Ht, xedges, yedges = np.histogram2d((reslist * len(form_input)), (form_input[map(str, reslist)]).astype(float).values.flatten(), bins=[bin_length_y, [1, len_min, max(helix_length_range)]], normed=False, weights=None)
                H_1 = Ht * 100.0 / len(form_input.index)

                H = H_0 - H_1
                H[H == 0] = 'Nan'

                figpath = self.fig_path_method(func_dir_name)
                figname = "%s/replica_id_%s_diff_cum.pdf" % (figpath, types)
                imgp = ax.imshow(H.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', vmax=12, vmin=-12, cmap=plt.get_cmap('PiYG'))
                cbar = ax.figure.colorbar(imgp, ax=ax)
                # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
                # self.ax.set_yticks(np.arange(reslist_min, reslist_max, 10))
                # self.ax.set_xticks(np.arange(helix_length_range[0], helix_length_range[1], 2))
                ax.set_yticks([1, max(helix_length_range)])
                fig.set_size_inches(16, 4)
                # ax.set_title(key)
                plt.savefig(figname, bbox_inches='tight', dpi=100)


                fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 4), sharey=True)
                form_input = grouped_form.get_group(grouped_form.groups.keys()[0])
                print grouped_form.groups.keys()[0]
                if types == 'helix':
                    len_min = 5
                if types == 'beta':
                    len_min = 5
                Ht, xedges, yedges = np.histogram2d((form_input[map(str, reslist)]).astype(float).values.flatten(), (reslist * len(form_input)), bins=[bin_length_x, bin_length_y], normed=False, weights=None, range=[helix_length_range, [reslist_min, reslist_max]])
                H_0 = Ht * 100.0 / len(form_input.index)

                form_input = grouped_form.get_group(grouped_form.groups.keys()[1])
                Ht, xedges, yedges = np.histogram2d((form_input[map(str, reslist)]).astype(float).values.flatten(), (reslist * len(form_input)), bins=[bin_length_x, bin_length_y], normed=False, weights=None, range=[helix_length_range, [reslist_min, reslist_max]])
                H_1 = Ht * 100.0 / len(form_input.index)

                H = H_0 - H_1
                H[H == 0] = 'Nan'

                figpath = self.fig_path_method(func_dir_name)
                figname = "%s/%s_diff_per_box.pdf" % (figpath, types)
                imgp = ax.imshow(H.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', vmax=5, vmin=-5, cmap=plt.get_cmap('PiYG'))
                cbar = ax.figure.colorbar(imgp, ax=ax)
                # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
                # self.ax.set_yticks(np.arange(reslist_min, reslist_max, 10))
                # self.ax.set_xticks(np.arange(helix_length_range[0], helix_length_range[1], 2))
                #ax.set_xticks([1, max(helix_length_range)])
                fig.set_size_inches(8, 8)
                # ax.set_title(key)
                plt.savefig(figname, bbox_inches='tight', dpi=100)

                ncol = 2
                nrow = int(np.ceil(grouped_form.ngroups / ncol))

                grouped_form = helix.groupby(['form'])
                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)
                if types == 'helix':
                    len_min = 5
                if types == 'beta':
                    len_min = 5

                for (key, ax) in zip(grouped_form.groups.keys(), axes.flatten()):
                            form_input = grouped_form.get_group(key)

                            Ht, xedges, yedges = np.histogram2d((form_input[map(str, reslist)]).astype(float).values.flatten(), (reslist * len(form_input)), bins=[bin_length_x, bin_length_y], normed=False, weights=None, range=[[len_min, 18], [reslist_min, reslist_max]])
                            H = Ht * 100.0 / len(form_input.index)
                            H[H == 0] = 'Nan'
                            figpath = self.fig_path_method(func_dir_name)
                            figname = "%s/%s_minimum_r.pdf" % (figpath, types)
                            imgp = ax.imshow(H.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', vmin=0, vmax=10, cmap=plt.get_cmap('Blues'))
                            cbar = ax.figure.colorbar(imgp, ax=ax)
                            # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
                            # self.ax.set_yticks(np.arange(reslist_min, reslist_max, 10))
                            # self.ax.set_xticks(np.arange(helix_length_range[0], helix_length_range[1], 2))
                            fig.set_size_inches(ncol * 8, nrow * 8)
                            ax.set_title(key)

                plt.savefig(figname, bbox_inches='tight', dpi=100)
                #  Plots the kde and histogram distribution of Rg for the given Rg distribution
                # helix  = helix.loc[helix['form'] == "M_hid_11"]

                # length of SS at reidue 66 vs Rg bins
                bin_length_y_rg = 20
                bin_length_x_sum = 20
                helix_length_range_rg = [helix[map(str, reslist)].values.min(), helix[map(str, reslist)].values.max()]
                # Particular form for all the given temeprature together
                # helix['sum'] = (helix[map(str, reslist)]).astype(float).sum(axis=1)
                helix['sum'] = (helix[['66']]).astype(float).sum(axis=1)
                grouped_form = helix.groupby(['form'])

                ncol = 2
                nrow = int(np.ceil(grouped_form.ngroups / ncol))

                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)

                for (key, ax) in zip(grouped_form.groups.keys(), axes.flatten()):
                            form_input = grouped_form.get_group(key)

                            Ht, xedges, yedges = np.histogram2d(form_input['sum'], (form_input['rgyr']).astype(float), bins=[bin_length_x_sum, bin_length_y_rg], range=[[helix['sum'].min(), 20], [helix['rgyr'].min(), helix['rgyr'].max()]], normed=False, weights=None)
                            H = Ht * 100.0 / len(form_input.index)
                            H[H == 0] = 'Nan'
                            figpath = self.fig_path_method(func_dir_name)
                            figname = "%s/66_rgyr_%s.pdf" % (figpath, types)
                            imgp = ax.imshow(H.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', vmax=5, vmin=0, cmap=plt.get_cmap('Blues'))
                            cbar = ax.figure.colorbar(imgp, ax=ax)
                            # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
                            # self.ax.set_yticks(np.arange(reslist_min, reslist_max, 10))
                            # self.ax.set_xticks(np.arange(helix_length_range[0], helix_length_range[1], 2))
                            fig.set_size_inches(ncol * 8, nrow * 8)
                            ax.set_title(key)

                plt.savefig(figname, bbox_inches='tight', dpi=100)

                ncol = 2
                nrow = int(np.ceil(grouped_form.ngroups / ncol))

                grouped_form_list_o = helix
                #grouped_form_list = helix.loc[(helix['93'] > 5) | (helix['91'] > 5)]  
                grouped_form_list = helix.loc[(helix['91'] > 7)]              
                grouped_form_o = grouped_form_list_o.groupby(['form'])
                grouped_form = grouped_form_list.groupby(['form'])
                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)

                for (key, ax) in zip(grouped_form.groups.keys(), axes.flatten()):
                            form_input = grouped_form.get_group(key)

                            Ht, xedges, yedges = np.histogram2d((form_input[map(str, reslist)]).astype(float).values.flatten(), (reslist * len(form_input)), bins=[bin_length_x, bin_length_y], range=[helix_length_range, [reslist_min, reslist_max]], normed=False, weights=None)
                            H = Ht * 100.0 / len((grouped_form_o.get_group(key)).index)
                            H[H == 0] = 'Nan'
                            figpath = self.fig_path_method(func_dir_name)
                            figname = "%s/%s_93_r_inverted.pdf" % (figpath, types)
                            imgp = ax.imshow(H, interpolation='nearest', origin='low', extent=[yedges[0], yedges[-1], xedges[0], xedges[-1]], aspect='auto', vmin=0, vmax=3, cmap=plt.get_cmap('Blues'))
                            cbar = ax.figure.colorbar(imgp, ax=ax)
                            # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
                            # self.ax.set_yticks(np.arange(reslist_min, reslist_max, 10))
                            # self.ax.set_xticks(np.arange(helix_length_range[0], helix_length_range[1], 2))
                            fig.set_size_inches(ncol * 16, nrow * 4)
                            #ax.set_title(key)

                plt.savefig(figname, bbox_inches='tight', dpi=100)

                ncol = 2
                nrow = int(np.ceil(grouped_form.ngroups / ncol))

                grouped_form_list = helix.loc[helix['66'] > 7]
                grouped_form_list_o = helix
                grouped_form_o = grouped_form_list_o.groupby(['form'])
                grouped_form = grouped_form_list.groupby(['form'])
                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)

                for (key, ax) in zip(grouped_form.groups.keys(), axes.flatten()):
                            form_input = grouped_form.get_group(key)

                            Ht, xedges, yedges = np.histogram2d((form_input[map(str, reslist)]).astype(float).values.flatten(), (reslist * len(form_input)), bins=[bin_length_x, bin_length_y], range=[helix_length_range, [reslist_min, reslist_max]], normed=False, weights=None)
                            H = Ht * 100.0 / len((grouped_form_o.get_group(key)).index)
                            H[H == 0] = 'Nan'
                            figpath = self.fig_path_method(func_dir_name)
                            figname = "%s/%s_66_r.pdf" % (figpath, types)
                            imgp = ax.imshow(H, interpolation='nearest', origin='low', extent=[yedges[0], yedges[-1], xedges[0], xedges[-1]], aspect='auto', vmin=0, vmax=3, cmap=plt.get_cmap('Blues'))
                            cbar = ax.figure.colorbar(imgp, ax=ax)
                            # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
                            # self.ax.set_yticks(np.arange(reslist_min, reslist_max, 10))
                            # self.ax.set_xticks(np.arange(helix_length_range[0], helix_length_range[1], 2))
                            fig.set_size_inches(ncol * 16, nrow * 4)
                            ax.set_title(key)

                plt.savefig(figname, bbox_inches='tight', dpi=100)


                #  Plots the kde and histogram distribution of Rg
                #grouped_form_list = helix.loc[helix['66'] > 5]
                grouped_form = helix.groupby(['form'])

                ncol = 1
                nrow = 1
                fig, ax = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)
                for key, mc, my_label in zip(grouped_form.groups.keys(), ['teal', 'saddlebrown', 'steelblue', 'darkorange'], ['M(p)', 'V(p)', 'M', 'V' ]):
                        form_input = grouped_form.get_group(key)
                        form_input_oro = form_input[form_input['replica_id'] != 100]
                        form_input = form_input[form_input['66'] < 8]
                        # form_input = form_input[form_input['replica_id'] != 15]
                        #form_input = form_input[form_input['replica_id'] != 9]
                        #form_input = form_input[form_input['replica_id'] != 36]
                        #form_input = form_input[form_input['replica_id'] != 27]
                        #form_input = form_input[form_input['replica_id'] != 2]
                        #form_input = form_input[form_input['replica_id'] != 2]
                        #form_input = form_input[form_input['replica_id'] != 49]
                        hist, bins = np.histogram(form_input['rgyr'], bins=60, range=(1.2, 4))
                        #ax.bar(bins[:-1], hist.astype(np.float32) / hist.sum(), width=(bins[1]-bins[0]), alpha=0.5)
                        center = ((bins[:-1] + bins[1:]) / 2)
                        hist = [item * 1.0 / form_input_oro['rgyr'].shape[0] * 1.0 for item in hist]
                        ax.plot(center, hist, lw=2, label=key, color=mc)   

                        #ax.bar(bins[:-1], hist, width=(bins[1] - bins[0]), alpha=0.5)
                        #form_input['rgyr'].hist(ax=ax, alpha=0.5, label=key, bins=50, normed = 1)
                        # form_input['rgyr'].plot.kde(ax=ax, label=key, lw=2)   # this plots the probablity distribution
                ax.legend(loc='best')
                ax.set_ylabel("Frequency")
                ax.set_xlabel("Rg (nm)")
                figpath = self.fig_path_method(func_dir_name)
                figname = "%s/%s.pdf" % (figpath, "fig_gen_histogram")
                fig.set_size_inches(ncol * 8, nrow * 8)
                ax.set_title('rgyr')
                plt.savefig(figname, bbox_inches='tight', dpi=100)


                # overlapping plot of V and M at each replcia id Rg- this should show if the Rg of every replica is changing with time
                grouped_form_list = helix.loc[helix['66'] > 8]

                grouped_form = grouped_form_list.groupby(['replica_id'])
                ncol = 2
                nrow = int(np.ceil(grouped_form.ngroups / ncol))    

                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)  

                for (key, ax) in zip(grouped_form.groups.keys(), axes.flatten()):
                        ax.set_color_cycle(['steelblue', 'teal', 'darkorange', 'saddlebrown'])
                        for key2, form_input in grouped_form.get_group(key).groupby('form'):
                            form_input_oro = helix[helix['replica_id'] != 100]
                            hist, bins = np.histogram(form_input['rgyr'], bins=30, range=(1, 5))
                            ax.bar(bins[:-1], hist.astype(np.float32) / form_input_oro['rgyr'].shape[0], width=(bins[1] - bins[0]), alpha=0.5)
                            center = ((bins[:-1] + bins[1:]) / 2)
                            hist = [item * 1.0 / form_input_oro['rgyr'].shape[0] * 1.0 for item in hist]
                            ax.plot(center, hist, lw=2, label=key2) 

                            figpath = self.fig_path_method(func_dir_name)
                            figname = "%s/%s.pdf" % (figpath, "fig_gen_replica_hist")
                            ax.set_title(key)
                ax.legend(loc='best')
                fig.set_size_inches(ncol * 8, nrow * 8)
                plt.savefig(figname, bbox_inches='tight', dpi=100)


                # overlapping plot of V and M at each replcia id Rg- this should show if the Rg of every replica is changing with time
                grouped_form_list = helix.loc[helix['66'] > 8]

                grouped_form = grouped_form_list.groupby(['replica_id'])
                ncol = 2
                nrow = int(np.ceil(grouped_form.ngroups / ncol))    

                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)  

                for (key, ax) in zip(grouped_form.groups.keys(), axes.flatten()):
                        ax.set_color_cycle(['steelblue', 'teal', 'darkorange', 'saddlebrown'])
                        for key2, form_input in grouped_form.get_group(key).groupby('form'):    

                            figpath = self.fig_path_method(func_dir_name)
                            figname = "%s/%s.pdf" % (figpath, "fig_gen_replica")
                            form_input.sort_values(by=['time_frame'])['rgyr'].plot(ax=ax, lw=2, use_index=False, label=key2)
                            ax.set_title(key)
                ax.legend(loc='best')
                fig.set_size_inches(ncol * 8, nrow * 8)
                plt.savefig(figname, bbox_inches='tight', dpi=100)

                for my_resid in ['65, 66']:

                    if types == 'helix':
                        grouped_form_list = [helix.loc[helix['91'] == 0], helix.loc[helix['91'] == 1], helix.loc[helix['91'] > 2]]
                        grouped_form_list = [helix.loc[helix['63'] == 0], helix.loc[helix['63'] == 1], helix.loc[helix['63'] == 2], helix.loc[helix['63'] == 3], helix.loc[helix['63'] == 4], helix.loc[helix['63'] == 5], helix.loc[helix['63'] > 6] ]
                        grouped_form_list = [helix.loc[helix['63'] == 1]]

                    else:
                        grouped_form_list = [helix.loc[helix[my_resid] != 100], helix.loc[helix['resid'] == 1], helix.loc[helix['66'] == 2], helix.loc[helix['66'] == 3], helix.loc[helix['66'] > 3]]

                    ncol = 7
                    nrow = 1
                    fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)
                    for (key, ax) in zip(grouped_form_list, axes.flatten()):
                            form_input_pre = key.groupby(['form'])
                            for key in form_input_pre.groups.keys():
                                form_input = form_input_pre.get_group(key)
                                form_input_oro = form_input[form_input['replica_id'] != 100]
                                print key
                                y = form_input['time_frame']
                                y.to_csv("%s/%s/%s/time_frame_63_1h" % (self.output, key, func_dir_name), sep='\t')
                                hist, bins = np.histogram(form_input['rgyr'], range=(1, 5), bins=40)
                                #ax.bar(bins[:-1], hist.astype(np.float32) / hist.sum(), width=(bins[1]-bins[0]), alpha=0.5)
                                center = ((bins[:-1] + bins[1:]) / 2)
                                hist = [item * 1.0 / helix.loc[helix['form'] == key].shape[0] * 1.0 for item in hist]
                                ax.bar(bins[:-1], hist, width=(bins[1] - bins[0]), alpha=0.5)
                                #hist = [item * 1.0 / form_input['rgyr'].shape[0] * 1.0 for item in hist]
                                popu = form_input.shape[0] * 100.0 / helix.loc[helix['form'] == key].shape[0]
                                ax.plot(center, hist, lw=2, label="key_%s" % popu )
                                ax.legend(loc='best')
                    ax.set_ylabel("Frequency")
                    ax.set_xlabel("Rg (nm)")
                    figpath = self.fig_path_method(func_dir_name)
                    figname = "%s/fig_gen_histogram_%s_jji.pdf" % (figpath, types)
                    fig.set_size_inches(ncol * 8, nrow * 8)
                    ax.set_title('rgyr')
                    plt.savefig(figname, bbox_inches='tight', dpi=100)

                #  Plots the kde and histogram distribution of Rg on the basis of presence or absence of specific contacts
                # helix  = helix.loc[helix['form'] == "M_hid_11"]

                contact_file_list = ["time_frame_95_63", "time_frame_63_34", "time_frame_66_63", "time_frame_70_66", "time_frame_96_70", "time_frame_70_33"]
                contact_file_list = ["time_frame_93_66", "time_frame_66_32", "time_frame_95_66", "time_frame_66_34", "time_frame_66_29"]
                contact_file_list = ["time_frame_78_88", "time_frame_38_103", "time_frame_93_88", "time_frame_93_68", "time_frame_93_64", "combined_long_short" ]
                helix['time_frame'] = helix['time_frame'].astype(float)

                ncol = 3
                nrow = 6
                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4))

                for (contact_file, ax_row) in zip(contact_file_list, axes):
                        form_input_pre = helix.groupby(['form'])
                        for (key, ax) in zip(form_input_pre.groups.keys(), ax_row):
                            form_input = form_input_pre.get_group(key)
                            ax.set_color_cycle(['dimgray', 'lime'])
                            contact_frames = np.loadtxt("%s/%s/sb_6.5/%s" % (self.output, key, contact_file))[:, 1]
                            for form_input_contacts in [form_input[~form_input.time_frame.isin(contact_frames)], form_input[form_input.time_frame.isin(contact_frames)]]:
                                    hist, bins = np.histogram(form_input_contacts['rgyr'], bins=60, range=(1.2, 4))
                                    #ax.bar(bins[:-1], hist.astype(np.float32) / hist.sum(), width=(bins[1]-bins[0]), alpha=0.5)
                                    center = ((bins[:-1] + bins[1:]) / 2)
                                    #hist = [item * 1.0 / form_input_contacts['rgyr'].shape[0] * 1.0 for item in hist]
                                    hist = [item * 1.0 / form_input['rgyr'].shape[0] * 1.0 for item in hist]

                                    #ax.legend(['MD', 'NMR'])
                                    ax.bar(bins[:-1], hist, width=(bins[1] - bins[0]), alpha=0.6)
                                    # form_input_contacts['rgyr'].plot.kde(ax=ax, label=key, lw=2)   # this plots the probablity distribution
                                    #form_input_contacts['rgyr'].plot.hist(ax=ax, label=key, lw=2, bins=50, normed = 1)
                                    popu = form_input_contacts.shape[0] * 100.0 / form_input.shape[0]
                                    #ax.plot(center, hist, lw=4, label="%s %s" % (key[0], int(popu)))
                                    ax.legend(loc='best')
                # ax.set_ylabel("Frequency")
                #ax.set_xlabel("Rg (nm)")
                figpath = self.fig_path_method(func_dir_name)
                figname = "%s/fig_gen_histogram_%s_contacts_interchange_colors.pdf" % (figpath, types)
                fig.set_size_inches(ncol * 10, nrow * 10)
                # ax.set_title('rgyr')
                plt.savefig(figname, bbox_inches='tight', dpi=100)

                #  Plots the kde and histogram distribution of Rg for the given Rg distribution
                # helix  = helix.loc[helix['form'] == "M_hid_11"]

                if types == 'helix':
                    grouped_form_list = [helix.loc[helix['91'] == 0], helix.loc[helix['91'] == 1], helix.loc[helix['91'] > 2]]
                    grouped_form_list = [helix.loc[helix['63'] == 0], helix.loc[helix['63'] == 1], helix.loc[helix['63'] == 2], helix.loc[helix['63'] == 3], helix.loc[helix['63'] == 4], helix.loc[helix['63'] == 5], helix.loc[helix['63'] > 6] ]
                    grouped_form_list = [helix.loc[helix['63'] == 1]]

                else:
                    grouped_form_list = [helix.loc[helix['66'] != 100], helix.loc[helix['66'] == 1], helix.loc[helix['66'] == 2], helix.loc[helix['66'] == 3], helix.loc[helix['66'] > 3]]

                ncol = 7
                nrow = 1
                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)
                for (key, ax) in zip(grouped_form_list, axes.flatten()):
                        form_input_pre = key.groupby(['form'])
                        for key in form_input_pre.groups.keys():
                            form_input = form_input_pre.get_group(key)
                            form_input_oro = form_input[form_input['replica_id'] != 100]
                            print key
                            if key == "V_hid_11":

                                form_input = form_input[form_input['replica_id'] != 2]
                                form_input = form_input[form_input['replica_id'] != 31]
                                form_input = form_input[form_input['replica_id'] != 100]
                            else:
                                form_input = form_input[form_input['replica_id'] != 14]
                                form_input = form_input[form_input['replica_id'] != 54]
                                form_input = form_input[form_input['replica_id'] != 15]

                            y = form_input['time_frame']
                            y.to_csv("%s/%s/%s/time_frame_63_1h" % (self.output, key, func_dir_name), sep='\t')
                            hist, bins = np.histogram(form_input['rgyr'], range=(1, 5), bins=40)
                            #ax.bar(bins[:-1], hist.astype(np.float32) / hist.sum(), width=(bins[1]-bins[0]), alpha=0.5)
                            center = ((bins[:-1] + bins[1:]) / 2)
                            hist = [item * 1.0 / helix.loc[helix['form'] == key].shape[0] * 1.0 for item in hist]
                            ax.bar(bins[:-1], hist, width=(bins[1] - bins[0]), alpha=0.5)
                            #hist = [item * 1.0 / form_input['rgyr'].shape[0] * 1.0 for item in hist]
                            popu = form_input.shape[0] * 100.0 / helix.loc[helix['form'] == key].shape[0]
                            ax.plot(center, hist, lw=2, label="key_%s" % popu )
                            ax.legend(loc='best')
                ax.set_ylabel("Frequency")
                ax.set_xlabel("Rg (nm)")
                figpath = self.fig_path_method(func_dir_name)
                figname = "%s/fig_gen_histogram_%s_jji.pdf" % (figpath, types)
                fig.set_size_inches(ncol * 8, nrow * 8)
                ax.set_title('rgyr')
                plt.savefig(figname, bbox_inches='tight', dpi=100)



                # length of SS at reidue 66 vs Rg bins
                bin_length_y_rg = 20
                bin_length_x_sum = 20
                helix_length_range_rg = [helix[map(str, reslist)].values.min(), helix[map(str, reslist)].values.max()]
                # Particular form for all the given temeprature together
                # helix['sum'] = (helix[map(str, reslist)]).astype(float).sum(axis=1)
                helix['sum'] = (helix[['66']]).astype(float).sum(axis=1)
                grouped_form = helix.groupby(['form'])

                ncol = 2
                nrow = int(np.ceil(grouped_form.ngroups / ncol))

                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)

                for (key, ax) in zip(grouped_form.groups.keys(), axes.flatten()):
                        form_input = grouped_form.get_group(key)
                        # form_input_contacts = form_input.loc[form_input['63'] >= 4]
                        form_input_contacts = form_input.loc[form_input['rgyr'] < 1.6]
                        hist, bins = np.histogram(form_input_contacts['rgyr'], range=(1, 5), bins=30)
                        # ax.bar(bins[:-1], hist.astype(np.float32) / hist.sum(), width=(bins[1]-bins[0]), alpha=0.5)
                        center = ((bins[:-1] + bins[1:]) / 2)
                        hist = [item * 1.0 / form_input['rgyr'].shape[0] * 1.0 for item in hist]
                        popu = form_input_contacts.shape[0] * 100.0 / form_input.shape[0]
                        ax.plot(center, hist, lw=1, label="%s %s" % (key[0], int(popu)))
                        ax.legend(loc='best')

                figpath = self.fig_path_method(func_dir_name)
                figname = "%s/fig_gen_histogram_%s_66.pdf" % (figpath, types)
                fig.set_size_inches(ncol * 10, nrow * 10)

                plt.savefig(figname, bbox_inches='tight', dpi=100)

                #  Plots the kde and histogram distribution of Rg on the basis of presence or absence of specific contacts
                # helix  = helix.loc[helix['form'] == "M_hid_11"]

                contact_file_list = ["time_frame_96_63", "time_frame_63_34"]
                helix['time_frame'] = helix['time_frame'].astype(float)

                ncol = 3
                nrow = 2
                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4))
                for (contact_file, ax_row) in zip(contact_file_list, axes):
                        form_input_pre = helix.groupby(['form'])
                        for key in form_input_pre.groups.keys():
                            form_input = form_input_pre.get_group(key)
                            contact_frames = np.loadtxt("%s/%s/all_pairs/%s" % (self.output, key, contact_file))[:, 1]
                            for form_input_contacts, ax in zip([form_input[~form_input.time_frame.isin([0])], form_input[~form_input.time_frame.isin(contact_frames)], form_input[form_input.time_frame.isin(contact_frames)]], ax_row):
                                    hist, bins = np.histogram(form_input_contacts['rgyr'], range=(1, 5), bins=30)
                                    # ax.bar(bins[:-1], hist.astype(np.float32) / hist.sum(), width=(bins[1]-bins[0]), alpha=0.5)
                                    center = ((bins[:-1] + bins[1:]) / 2)
                                    hist = [item * 1.0 / form_input_contacts['rgyr'].shape[0] * 1.0 for item in hist]
                                    popu = form_input_contacts.shape[0] * 100.0 / form_input.shape[0]
                                    ax.plot(center, hist, lw=1, label="%s %s" % (key[0], int(popu)))
                                    ax.legend(loc='best')
                # ax.set_ylabel("Frequency")
                #ax.set_xlabel("Rg (nm)")
                figpath = self.fig_path_method(func_dir_name)
                figname = "%s/fig_gen_histogram_%s_contacts.pdf" % (figpath, types)
                fig.set_size_inches(ncol * 10, nrow * 10)
                # ax.set_title('rgyr')
                plt.savefig(figname, bbox_inches='tight', dpi=100)

                # combined frames

                # contact_file_list_pre = ["time_frame_95_63", "time_frame_63_34", "time_frame_66_63", "time_frame_70_66", "time_frame_96_70", "time_frame_70_33"]
                contact_file_list = [["time_frame_95_63", "time_frame_63_34", "time_frame_96_70", "time_frame_70_33"], ["time_frame_96_70", "time_frame_70_33"]]
                helix['time_frame'] = helix['time_frame'].astype(float)
                ncol = 3
                nrow = 6
                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4))
                for (contact_file, ax_row) in zip(contact_file_list, axes):
                        form_input_pre = helix.groupby(['form'])
                        for (key, ax) in zip(form_input_pre.groups.keys(), ax_row):
                            form_input = form_input_pre.get_group(key)
                            contact_frames_pre = [np.loadtxt("%s/%s/all_pairs/%s" % (self.output, key, i))[:, 1] for i in contact_file]
                            contact_frames = np.concatenate(contact_frames_pre)
                            for form_input_contacts in [form_input[~form_input.time_frame.isin([0])], form_input[~form_input.time_frame.isin(contact_frames)], form_input[form_input.time_frame.isin(contact_frames)]]:
                                    hist, bins = np.histogram(form_input_contacts['rgyr'], range=(1, 5), bins=30)
                                    #ax.bar(bins[:-1], hist.astype(np.float32) / hist.sum(), width=(bins[1]-bins[0]), alpha=0.5)
                                    center = ((bins[:-1] + bins[1:]) / 2)
                                    #hist = [item * 1.0 / form_input_contacts['rgyr'].shape[0] * 1.0 for item in hist]
                                    hist = [item * 1.0 / form_input['rgyr'].shape[0] * 1.0 for item in hist]
                                    # form_input_contacts['rgyr'].plot.kde(ax=ax, label=key, lw=2)   # this plots the probablity distribution
                                    #form_input_contacts['rgyr'].plot.hist(ax=ax, label=key, lw=2, bins=50, normed = 1)
                                    popu = form_input_contacts.shape[0] * 100.0 / form_input.shape[0]
                                    ax.plot(center, hist, lw=2, label="%s %s" % (key[0], int(popu)))
                                    ax.legend(loc='best')
                # ax.set_ylabel("Frequency")
                #ax.set_xlabel("Rg (nm)")
                figpath = self.fig_path_method(func_dir_name)
                figname = "%s/fig_gen_histogram_%s_contacts_interchange_colors_combined.pdf" % (figpath, types)
                fig.set_size_inches(ncol * 10, nrow * 10)
                # ax.set_title('rgyr')
                plt.savefig(figname, bbox_inches='tight', dpi=100)

                contact_file_list = ["time_frame_96_63", "time_frame_63_34"]
                helix['time_frame'] = helix['time_frame'].astype(float)

                ncol = 2
                nrow = 2
                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)
                for (contact_file, ax_row) in zip(contact_file_list, axes):
                        form_input_pre = helix.groupby(['form'])
                        for (key, ax) in zip(form_input_pre.groups.keys(), ax_row):
                            form_input = form_input_pre.get_group(key)
                            contact_frames = np.loadtxt("%s/%s/all_pairs/%s" % (self.output, key, contact_file))[:, 1]
                            # for form_input_contacts, ax in zip([form_input[~form_input.time_frame.isin(contact_frames)], form_input[form_input.time_frame.isin(contact_frames)]], ax_row):
                            form_input_contacts = form_input[form_input.time_frame.isin(contact_frames)]
                            Ht, xedges, yedges = np.histogram2d((form_input_contacts[map(str, reslist)]).astype(float).values.flatten(), (reslist * len(form_input_contacts)), bins=[bin_length_x, bin_length_y], range=[helix_length_range, [reslist_min, reslist_max]], normed=False, weights=None)
                            H = Ht * 100.0 / len(form_input_contacts.index)
                            H[H == 0] = 'Nan'
                            figpath = self.fig_path_method(func_dir_name)
                            imgp = ax.imshow(H.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', vmax=c_max, vmin=0, cmap=plt.get_cmap('Blues'))
                            cbar = ax.figure.colorbar(imgp, ax=ax)
                fig.set_size_inches(ncol * 24, nrow * 24)
                figname = "%s/%s_hist_contacts.pdf" % (figpath, types)
                plt.savefig(figname, bbox_inches='tight', dpi=100)

                ncol = 2
                nrow = 3
                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)

                for (key, ax_row) in zip(grouped_form_list, axes):
                        form_input_pre = key.groupby(['form'])
                        for (key, ax) in zip(form_input_pre.groups.keys(), ax_row):
                            form_input = form_input_pre.get_group(key)
                            Ht, xedges, yedges = np.histogram2d((form_input[map(str, reslist)]).astype(float).values.flatten(), (reslist * len(form_input)), bins=[bin_length_x, bin_length_y], range=[helix_length_range, [reslist_min, reslist_max]], normed=False, weights=None)
                            H = Ht * 100.0 / len(form_input.index)
                            H[H == 0] = 'Nan'
                            figpath = self.fig_path_method(func_dir_name)
                            imgp = ax.imshow(H.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', vmax=c_max, vmin=0, cmap=plt.get_cmap('Blues'))
                            cbar = ax.figure.colorbar(imgp, ax=ax)
                            # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
                            # self.ax.set_yticks(np.arange(reslist_min, reslist_max, 10))
                            # self.ax.set_xticks(np.arange(helix_length_range[0], helix_length_range[1], 2))
                            ax.set_title(key)

                fig.set_size_inches(ncol * 24, nrow * 24)
                figname = "%s/%s_hist.pdf" % (figpath, types)
                plt.savefig(figname, bbox_inches='tight', dpi=100)

                fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 4), sharey=True)
                form_input = grouped_form.get_group(grouped_form.groups.keys()[0])
                print grouped_form.groups.keys()[0]
                if types == 'helix':
                    len_min = 6
                if types == 'beta':
                    len_min = 4
                Ht, xedges, yedges = np.histogram2d((form_input[map(str, reslist)]).astype(float).values.flatten(), (reslist * len(form_input)), bins=[[1, len_min, max(helix_length_range)], bin_length_y], range=[helix_length_range, [reslist_min, reslist_max]], normed=False, weights=None)
                H_0 = Ht * 100.0 / len(form_input.index)

                form_input = grouped_form.get_group(grouped_form.groups.keys()[1])
                Ht, xedges, yedges = np.histogram2d((form_input[map(str, reslist)]).astype(float).values.flatten(), (reslist * len(form_input)), bins=[[1, len_min, max(helix_length_range)], bin_length_y], range=[helix_length_range, [reslist_min, reslist_max]], normed=False, weights=None)
                H_1 = Ht * 100.0 / len(form_input.index)

                H = H_0 - H_1
                H[H == 0] = 'Nan'

                figpath = self.fig_path_method(func_dir_name)
                figname = "%s/%s_diff.pdf" % (figpath, types)
                imgp = ax.imshow(H.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', vmax=4, vmin=-4, cmap=plt.get_cmap('PiYG'))
                cbar = ax.figure.colorbar(imgp, ax=ax)
                # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
                # self.ax.set_yticks(np.arange(reslist_min, reslist_max, 10))
                # self.ax.set_xticks(np.arange(helix_length_range[0], helix_length_range[1], 2))
                fig.set_size_inches(8, 8)
                # ax.set_title(key)
                plt.savefig(figname, bbox_inches='tight', dpi=100)

                # Particular form for a particular replica_id for all temerature combined
                grouped_form = helix.groupby(['form', 'replica_id'])
                ncol = 2
                nrow = int(np.ceil(grouped_form.ngroups / ncol))

                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)

                for (key, ax) in zip(grouped_form.groups.keys(), axes.flatten()):
                            form_input = grouped_form.get_group(key)

                            Ht, xedges, yedges = np.histogram2d((form_input[map(str, reslist)]).values.flatten(), (reslist * len(form_input)), bins=[bin_length_x, bin_length_y], range=[helix_length_range, [reslist_min, reslist_max]], normed=False, weights=None)
                            H = Ht * 100.0 / len(form_input.index)

                            H[H == 0] = 'Nan'

                            figpath = self.fig_path_method(func_dir_name)
                            figname = "%s/%s_rep_id.pdf" % (figpath, types)

                            imgp = ax.imshow(H.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', vmax=c_max, vmin=0.0)
                            fig.set_size_inches(ncol * 8, nrow * 8)
                            ax.set_title(key)

                plt.savefig(figname, bbox_inches='tight', dpi=100)

                # Particular form for a particular temperature
                grouped_form = helix.groupby(['form', 'Temperature'])
                ncol = 2
                nrow = int(np.ceil(grouped_form.ngroups / ncol))

                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)

                for (key, ax) in zip(grouped_form.groups.keys(), axes.flatten()):
                            form_input = grouped_form.get_group(key)

                            Ht, xedges, yedges = np.histogram2d((form_input[map(str, reslist)]).values.flatten(), (reslist * len(form_input)), bins=[bin_length_x, bin_length_y], range=[helix_length_range, [reslist_min, reslist_max]], normed=False, weights=None)
                            H = Ht * 100.0 / len(form_input.index)

                            H[H == 0] = 'Nan'

                            figpath = self.fig_path_method(func_dir_name)
                            figname = "%s/%s_temperature.pdf" % (figpath, types)

                            imgp = ax.imshow(H.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', vmax=c_max, vmin=0.0)

                            fig.set_size_inches(ncol * 8, nrow * 8)
                            ax.set_title(key)
                plt.savefig(figname, bbox_inches='tight', dpi=100)

                return types