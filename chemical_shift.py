import subprocess
import fabfile
import numpy as np
import math
import paramiko
import colorsys
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.patches as mpatches
from os.path import exists as path_exists
import pylab as plot



def get_my_base(base):

    class sparta_calc(base):
        @fabfile._decorator_calc_loop
        # calculating shift every 500ps, amount of ps in each block
        def calc_chemical_shifts(self, func_dir_name, my_xtc_freq=100, custom_st_time=0):
                    """calculates the sparta_calc shifts from MD trajectories
                    custom start time is to save time from calculating the chemical shift for the same trajectory"""

                    self.population_calc()
                    end_time = self.tot_frames * self.xtc_freq
                    """for every temprature"""
                    for counter in range(custom_st_time, end_time, my_xtc_freq):

                        subprocess.call("echo 1 > inp ; gmx trjconv -f %s -s %s -dump %s -o input_3.pdb < inp ; sed -i s/HID/HIS/ input_3.pdb && sed -i s/HIE/HIS/ input_3.pdb && sed -i s/HIP/HIS/ input_3.pdb && grep -v 'NME\|ACE' input_3.pdb > inp_4.pdb; mkdir struct_tab ; mkdir pred_tab" % (
                            self.xtc, self.tpr, counter), shell=True, cwd=self.repl_dir)
                        subprocess.call("/u1/home/lohia/data/Programs/SPARTA+/sparta+ -in inp_4.pdb -out pred_tab/pred_%s.tab -outS struct_tab/struct_%s.tab ; rm -rf \#*" % (
                            counter, counter), shell=True, cwd=self.repl_dir)

        def generate_data_frame(self, func_dir_name, st_time, my_xtc_freq):
            """concatenated the calculated MD shifts from every frame into one frame"""
            headlist = ['RESID', 'RESNAME', 'ATOMNAME', 'SS_SHIFT_SPARTA', 'SS_SHIFT', 'SHIFT', 'RC_SHIFT_SPARTA', 'HM_SHIFT', 'EF_SHIFT', 'SIGMA', 'RC_SHIFT_POLSUN']
            for i in ['form', 'replica_id', 'cluster_index', 'time_frame', 'Temperature']:
                headlist.append(i)
            helix = pd.DataFrame(columns=map(str, headlist))
            for form in self.form_list:
                inp_dir_form_rc = ("/u1/home/lohia/Dropbox/lab_work/output_02/%s/nmr/rc.tab" % (form[0]))
                if not path_exists(inp_dir_form_rc):
                    inp_dir_form_rc = ("/u1/home/lohia/Dropbox/lab_work/output_02/%s/nmr/rc.tab" % ("V"))  # by default V nmr shifts are taken, this should be only take for ff input
                    print "taking ff"
                for i in range(0, self.calc_replica):
                    self.i = i
                    for s in range(self.parent_cluster_no):
                        self.parent_cluster_id = s
                        self.form = form
                        self.population_calc()
                        end_time = self.tot_frames * self.xtc_freq
                        end_time = 2000000 # if i want to use a different end time
                        #end_time = 900000 # if i want to use a different end time
                        self.repl_dir_parent = self.repl_dir_parent_method(func_dir_name)
                        for counter in range(st_time, end_time, my_xtc_freq):

                            inp_dir_form = self.repl_dir_parent + "/pred_tab/pred_%s.tab" % counter
                            _helix = pd.read_csv(inp_dir_form, sep="\t\t|\t|\s+|", index_col=False, skiprows=27, names=['RESID', 'RESNAME', 'ATOMNAME', 'SS_SHIFT_SPARTA', 'SHIFT', 'RC_SHIFT_SPARTA', 'HM_SHIFT', 'EF_SHIFT', 'SIGMA'], engine='python')
                            time_frame = counter * 1.0 / self.xtc_freq
                            _helix['form'] = self.form
                            _helix['Temperature'] = self.temp_ticks[self.replica_curr]
                            # _helix['Time'] = counter
                            _helix['time_frame'] = time_frame
                            _helix['cluster_index'] = self.clus_curr
                            _helix['replica_id'] = self.demux_file_per_temperature()[int(time_frame)]

                            # calculates SS_shift with respect to polsun random coil chemical shifts for disorederd proteins
                            _helix_rc = pd.read_csv(inp_dir_form_rc, sep="\t\t|\t|\s+|", index_col=False, skiprows=2, names=['RESID', 'CA', 'CB', 'C', 'N', 'HN', 'HA'], engine='python')
                            _helix_rc['RESID'] = _helix_rc['RESID'].apply(lambda x: x + 21)
                            _helix_rc_melted = pd.melt(_helix_rc, id_vars=['RESID'], value_vars=['CA', 'CB', 'C', 'N', 'HN', 'HA'], var_name='ATOMNAME', value_name='RC_SHIFT_POLSUN')
                            _merged_data = pd.merge(_helix, _helix_rc_melted, on=['RESID', 'ATOMNAME'])
                            _merged_data['RC_SHIFT_POLSUN'] = _merged_data['RC_SHIFT_POLSUN'].apply(lambda x: float(x))
                            _merged_data['SS_SHIFT'] = _merged_data['SHIFT'] - _merged_data['RC_SHIFT_POLSUN']

                            # concatenates the data from every frame
                            helix = helix.append(_merged_data)

            return helix

        def fig_SS_chemical_shifts_plot(self, func_dir_name, my_xtc_freq=100, st_time=0):
            """plots secondary chemical shifts where random coil shifts are used from polsun prediction"""
            input_file = self.generate_data_frame(func_dir_name, st_time, my_xtc_freq)
            for form in ["V", "M"]:
                inp_dir_form = ("/u1/home/lohia/Dropbox/lab_work/output_02/%s/nmr/%s.tab" % (form, form))
                _helix = pd.read_csv(inp_dir_form, sep="\t\t|\t|\s+|", index_col=False, skiprows=6, names=['RESID', 'RESNAME', 'ATOMNAME', 'SHIFT'], engine='python')
                inp_dir_form_rc = ("/u1/home/lohia/Dropbox/lab_work/output_02/%s/nmr/rc.tab" % (form))
                _helix_rc = pd.read_csv(inp_dir_form_rc, sep="\t\t|\t|\s+|", index_col=False, skiprows=2, names=['RESID', 'CA', 'CB', 'C', 'N', 'HN', 'HA'], engine='python')

                # calculates SS_shift with respect to polsun random coil chemical shifts for disorederd proteins
                _helix['RESID'] = _helix['RESID'].apply(lambda x: x + 21)
                _helix_rc['RESID'] = _helix_rc['RESID'].apply(lambda x: x + 21)
                _helix_rc_melted = pd.melt(_helix_rc, id_vars=['RESID'], value_vars=['CA', 'CB', 'C', 'N', 'HN', 'HA'], var_name='ATOMNAME', value_name='RC_SHIFT_POLSUN')
                _merged_data = pd.merge(_helix, _helix_rc_melted, on=['RESID', 'ATOMNAME'])
                _merged_data['RC_SHIFT_POLSUN'] = _merged_data['RC_SHIFT_POLSUN'].apply(lambda x: float(x))
                _merged_data['SS_SHIFT'] = _merged_data['SHIFT'] - _merged_data['RC_SHIFT_POLSUN']

                _merged_data['Temperature'] = 280
                _merged_data['form'] = "NMR_%s" % form
                _merged_data['time_frame'] = 3000000
                input_file = input_file.append(_merged_data)

            st_frame = st_time / self.xtc_freq
            helix = input_file[input_file['time_frame'] >= st_frame]  # this input file has chemical shifts for all the possible forms
            grouped = helix.groupby(['ATOMNAME'])
            ncol = 1
            nrow = int(np.ceil(grouped.ngroups / ncol))

            for form in self.form_list:
                curr_form = "NMR_%s" % form[0]  # for every form in the form list, I substact the given form NMR chemical shift

                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)

                for (key, ax) in zip(grouped.groups.keys(), axes.flatten()):
                    if key in ['CA', 'CB', 'C']:
                        print key
                        print form

                        grouped_input = grouped.get_group(key)
                        ss_shift = grouped_input.groupby(['RESID', 'form'])['SS_SHIFT'].mean().unstack()

                        # calculate rmsd with respect to nmr
                        #ss_shift = ss_shift.dropna()

                        shift_diff = ss_shift.subtract(ss_shift[curr_form], axis='index')

                        shift_diff_squared = shift_diff.apply(lambda x: x ** 2)
                        rmsd = shift_diff_squared.mean() ** .5
                        print rmsd
                        ax.set_color_cycle(['red', 'black'])
                        if form == "V_hid_11":
                            ax.set_color_cycle(['steelblue', 'black'])
                            ax.legend(['MD', 'NMR'])
                        else:
                            ax.set_color_cycle(['steelblue', 'black'])
                            ax.legend(['MD', 'NMR'])
                        ss_shift[form].plot(ax=ax, lw=2)

                        #ss_shift[form].plot.scatter(ax=ax, lw=2)
                        #ss_shift[form].plot.scatter(ax=ax, c='black')
                        ss_shift[curr_form].plot(ax=ax, lw=2, c='black' )
                        #print ss_shift.reset_index()['RESID']
                        ax.fill_between(ss_shift.reset_index()['RESID'], ss_shift.reset_index()[curr_form], ss_shift.reset_index()[form], where=(ss_shift.reset_index()[curr_form] - ss_shift.reset_index()[form]).abs() > 0, alpha=0.3)
                        ax.fill_between(ss_shift.reset_index()['RESID'], ss_shift.reset_index()[curr_form] + 0.5, 3, facecolor='gainsboro', alpha=0.8)
                        ax.fill_between(ss_shift.reset_index()['RESID'], ss_shift.reset_index()[curr_form] - 0.5, -3, facecolor='gainsboro', alpha=0.8)
                        ax.fill_between(ss_shift.reset_index()['RESID'], ss_shift.reset_index()[curr_form] + 0.5, ss_shift.reset_index()[form], where=(ss_shift.reset_index()[curr_form] + 0.5 - ss_shift.reset_index()[form]) <= 0, alpha=1, facecolor='steelblue', interpolate=True)
                        ax.fill_between(ss_shift.reset_index()['RESID'], ss_shift.reset_index()[curr_form] - 0.5, ss_shift.reset_index()[form], where=(ss_shift.reset_index()[curr_form] - 0.5 - ss_shift.reset_index()[form]) >= 0, alpha=1, facecolor='steelblue', interpolate=True)
                        ax.set_color_cycle(['red', 'black'])
                        if form == "V_hid_11":
                            ax.set_color_cycle(['darkorange', 'black'])
                            ax.legend(['MD', 'NMR'])
                        else:
                            ax.set_color_cycle(['steelblue', 'black'])
                            ax.legend(['MD', 'NMR'])
                        figpath = self.fig_path_method(func_dir_name)
                        figname = "%s/%s_raw.pdf" % (figpath, form)
                        # ax.legend(loc='best')
                        #fig.set_size_inches(ncol * 10, nrow * 4)
                        # ax.set_title(key)
                        ax.set_ylim([-2, 2])
                        ax.set_xlim([self.resid_st, self.resid_end])
                        ax.set_xticks(np.arange(24, 113, 10))
                        fig.set_size_inches(ncol * 16, nrow * 4)
                        ax.set_xlabel("Residue")
                        #ax.set_ylabel("(ppm)" %key)
                        if key == 'CA':
                            ax.set_ylabel(r'$\Delta \delta C_\alpha (ppm)$')
                        elif key == 'CB':
                            ax.set_ylabel(r'$\Delta \delta C_\beta (ppm)$')
                        else:
                            ax.set_ylabel(r'$\Delta \delta C\' (ppm)$')
                plt.tight_layout()
                plt.savefig(figname, bbox_inches='tight', dpi=100)

                for nmr_form, m_c, shift in zip([form, curr_form], ['steelblue', 'darkorange'], [0, 0.2]):
                    ss_shift_per_atom = helix[helix['form'] == nmr_form].groupby(['RESID', 'ATOMNAME'])['SS_SHIFT'].mean().unstack()
                    ss_shift_per_atom['ca-cb'] = ss_shift_per_atom['CA'] - ss_shift_per_atom['CB']
                    axes.flatten()[-1].bar(ss_shift_per_atom.index + shift, ss_shift_per_atom['ca-cb'], label=nmr_form, width=1, color=m_c, alpha=0.5)
                    axes.flatten()[-1].set_title('CA-CB')
                    # axes.flatten()[-1].legend(loc='best')
                plt.tight_layout()
                plt.savefig(figname, bbox_inches='tight')

            # plots the V vs M ca and cb for both the forms, # V vs M for replica id

            nrow = 8

            fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12, 4), sharey=True)

            for form in self.form_list:

                for (key, ax) in zip(grouped.groups.keys(), axes.flatten()):
                    if key in ['CA', 'CB', 'C']:
                        grouped_input = grouped.get_group(key)
                        ss_shift = grouped_input.groupby(['RESID', 'form'])['SS_SHIFT'].mean().unstack()
                        # calculate rmsd with respect to nmr
                        ss_shift[form].plot(ax=ax, lw=2)
                        ax.legend(loc='best')
                        fig.set_size_inches(ncol * 10, nrow * 4)
                        ax.set_title(key)
                        ax.set_ylim([-2, 2])
                        ax.set_xlim([self.resid_st, self.resid_end])

            for form in ["NMR_V", "NMR_M"]:

                for (key, ax) in zip(['CA', 'CB', 'C'], [axes.flatten()[3], axes.flatten()[4], axes.flatten()[5]]):
                    if key in ['CA', 'CB', 'C']:
                        grouped_input = grouped.get_group(key)
                        ss_shift = grouped_input.groupby(['RESID', 'form'])['SS_SHIFT'].mean().unstack()
                        ss_shift[form].plot(ax=ax, lw=2)
                        ax.legend(loc='best')
                        fig.set_size_inches(ncol * 10, nrow * 4)
                        ax.set_title(key)
                        ax.set_ylim([-2, 2])
                        ax.set_xlim([self.resid_st, self.resid_end])

            for nmr_form, m_c, shift in zip(["NMR_V", "NMR_M"], ['darkorange', 'steelblue'], [0, 0.2]):
                ss_shift_per_atom = helix[helix['form'] == nmr_form].groupby(['RESID', 'ATOMNAME'])['SS_SHIFT'].mean().unstack()
                ss_shift_per_atom['ca-cb'] = ss_shift_per_atom['CA'] - ss_shift_per_atom['CB']
                axes.flatten()[6].bar(ss_shift_per_atom.index + shift, ss_shift_per_atom['ca-cb'], label=nmr_form[-1], width=1, color=m_c, alpha=0.7)
                ax.set_xticks(np.arange(24, 113, 10))
                axes.flatten()[6].set_title(r'$\Delta \delta C_\alpha - \Delta \delta C_\beta $')
                axes.flatten()[6].legend(loc='best')
                axes.flatten()[6].set_xlabel("Residue")
                axes.flatten()[6].set_ylabel(r'$ \delta (ppm)$')
                # axes.flatten()[6].set_title('CA-CB')

            for nmr_form, m_c, shift in zip(self.form_list, ['steelblue', 'darkorange'], [0, 0.2]):
                ss_shift_per_atom = helix[helix['form'] == nmr_form].groupby(['RESID', 'ATOMNAME'])['SS_SHIFT'].mean().unstack()
                ss_shift_per_atom['ca-cb'] = ss_shift_per_atom['CA'] - ss_shift_per_atom['CB']
                axes.flatten()[7].bar(ss_shift_per_atom.index + shift, ss_shift_per_atom['ca-cb'], label=nmr_form[0], width=1, color=m_c, alpha=0.7)
                ax.set_xticks(np.arange(24, 113, 10))
                axes.flatten()[7].set_ylabel(r'$\Delta \delta C_\alpha - \Delta \delta C_\beta (ppm)$')
                axes.flatten()[7].legend(loc='best')
                # axes.flatten()[7].set_title('CA-CB')

            figname = "%s/VvsM.pdf" % (figpath)
            for ax in axes.flatten():
                ax.legend(loc='best')
                fig.set_size_inches(ncol * 10, nrow * 4)
                ax.set_ylim([-2, 2])
                ax.set_xlim([self.resid_st, self.resid_end])

            plt.tight_layout()
            plt.savefig(figname, bbox_inches='tight')
            return

    return sparta_calc
