import fabfile
from os.path import expanduser
import mindist
import cischeck
import demux
import Rgyr
import secondary_structure
import clustering_cutoff
import contacts
import h_bond_1
import free_energy_resid
import ss_length
import vmd_rep
import rmsd
import no_h_bonds
import sasa
import vmd_rep
import mdmat
import chemical_shift
import chemical_shift_ff
import hydropro_calc
import end_dist
import contacts_rg
#import polystat

"""This is the place from which all the analysis functions are defined and run"""

def all_initialization():
    """initialization of user specific paths"""
    home = expanduser("~")
    fabfile.Calc.dropbox_path = "%s/Dropbox/lab_work/output_30qd_long/" % home
    fabfile.Calc.dropbox_path_mac = "%s/Dropbox/lab_work/output_30qd_long/" % home
    fabfile.Calc.git_dir = "/u1/home/lohia/Dropbox/replex_md/analysis_01/"
    fabfile.Calc.git_dir_mac = "/Users/lohia/Dropbox/replex_md/analysis_01/"
    fabfile.Calc.output = "/g3/home/lohia/T-rex-NMR-region/trex-04/output/"
    fabfile.Calc.input_path = "/g3/home/lohia/T-rex-NMR-region/trex-04/"
    fabfile.Calc.gro_file_name = "resid_23-113-capped.gro"
    fabfile.Calc.xtc_freq = 100
    fabfile.Calc.st_time = 0
    fabfile.Calc.end_time = 3000000  # 3us
    fabfile.Calc.resid_st = 23
    fabfile.Calc.resid_end = 113
    fabfile.Calc.num_r = 64
    fabfile.Calc.calc_replica = fabfile.Calc.num_r
    fabfile.Calc.form_list = ["V_hid_11", "M_hid_11"]
    fabfile.Calc.end_append_list = [50, 53]
    fabfile.Calc.nmr_list = ["19358.str.txt", "19357.str.txt"]
    fabfile.Calc.xtc_file_name = "trajout.xtc"
    # xtc file name before ramining minidst frames
    fabfile.Calc.xtc_file_name_before_mindist = "trajout.xtc"
    fabfile.Calc.tpr_file_name = "ex.tpr"
    fabfile.Calc.cont_file_folder = "/"
    fabfile.Calc.cont_file_output = "/"
    fabfile.Calc.temp_file_name = "rex_replica_temp.xvg"

def common_steps(analysis_name, base_class):
    myfunc = analysis_name
    return myfunc.get_my_base(base_class)

def _chemical_shift(Temperature=1, base_class=fabfile.Calc):  # This function file is also made public
    x = common_steps(chemical_shift, base_class)
    x(1).calc_chemical_shifts("sparta_calc", my_xtc_freq=100, custom_st_time = 0)
    x(1).fig_SS_chemical_shifts_plot_diff_NMR("sparta_calc", st_time=0)
    x(1).fig_SS_chemical_shifts_plot("sparta_calc", my_xtc_freq=100, st_time=800000)

def _mindist(Temperature=1, base_class=fabfile.Calc):

    x = common_steps(mindist, base_class)
    x(Temperature).calc("mindist")
    x(Temperature).hist_length_pandas("mindist")
    x(Temperature).fig_gen("mindist", st_time=800000)


def _vmd_rep(Temperature=1, base_class=fabfile.Calc):

    x = common_steps(vmd_rep, base_class)
    x(Temperature).temp("loading")


def _rmsd(Temperature=1, base_class=fabfile.Calc):

    x = common_steps(rmsd, base_class)
    x(Temperature).calc("rmsd")
    x(Temperature).fig_gen("rmsd", "rmsd.xvg", _ncol=1, _nrow=5)


def _sasa(Temperature=1, base_class=fabfile.Calc):

    x = common_steps(sasa, base_class)
    x(Temperature).calc("sasa")
    x(Temperature).generate_data_frame("sasa")
    x(Temperature).fig_gen("sasa", st_time=800000)


def _cischeck(Temperature=1, base_class=fabfile.Calc):

    x = common_steps(cischeck, base_class)
    x(Temperature).calc("cischeck")
    x(Temperature).fig_gen_all("cischeck", "cischeck.xvg", _nrow=1, _ncol=1)


def _demux(Temperature=1, base_class=fabfile.Calc):
    x = common_steps(demux, base_class)
    x(1).calc("demux")
    x(1).hist_length_pandas("demux")
    x(1).fig_gen("demux", st_time=800000)
    x(num_r).rep_mix("demux", "rex_replica_temp.xvg", _ncol=4, _nrow=64 / 4)
    x(1).round_trip("demux", st_time=800000)


def _Rgyr(Temperature=1, base_class=fabfile.Calc, my_clus="67_only", clus_id=0, cont_class=False):
    x = common_steps(Rgyr, base_class)
    x(Temperature).calc("rgyr")
    x(Temperature).hist_length_pandas("rgyr")
    x(Temperature).fig_gen("rgyr", st_time=800000)

def _hydropro_calc(Temperature=1, base_class=fabfile.Calc):
    x = common_steps(hydropro_calc, base_class)
    x(1).calc_hydrorpo_radius("hydropro_calc", my_xtc_freq=100, custom_st_time=0)
    x(1).fig_SS_chemical_shifts_plot_diff_NMR("sparta_calc", st_time=0)
    x(1).fig_hydropro_plot("hydropro_calc", my_xtc_freq=500, st_time=0)


def _secondary_structure(Temperature=1, base_class=fabfile.Calc, dpca_name="dpca_67", resid_1=67, resid_3=67):
    x = common_steps(secondary_structure, base_class)

    if base_class == fabfile.Clusters:
        x(Temperature).cluster_prop_calc(
            "ss_analysis", "ss.xvg", fabfile.Clusters.clustername)
        x(Temperature).cluster_percent_population("ss_analysis",
                                                  "ss.xvg", _nrow=1, _ncol=fabfile.Clusters.cutoff)
        x(Temperature).ss_vmd("ss_analysis", "ss.xvg", _nrow=base_class.cutoff,
                              _ncol=2, extra="all_helix", mfigname="all_helix", ylim_u=40)
    else:
        x(Temperature).calc_vmd("ss_analysis")
        x(Temperature).ss_vmd("ss_analysis", st_time=800000)


def _contacts(Temperature=1, base_class=fabfile.Calc):

    # this is always calculated for the frames with with unprocessed mindist
    x = common_steps(contacts, base_class)
    x(Temperature).all_contacts_ca("all_pairs_8.5_ca", 8.5)


def _free_energy_resid(Temperature=1, base_class=fabfile.Calc, resid=66, cont_class=False):
    x = common_steps(free_energy_resid, base_class)
    x(Temperature, cont_class=cont_class).calc("free_energy_%s" % resid, resid)
    x(Temperature, cont_class=cont_class).gibbs_single("free_energy_%s" %resid, "gibbs_%s.dat" % resid, _nrow=base_class.cutoff, _ncol=1)


def _mdmat(Temperature=1, base_class=fabfile.Calc):

    x = common_steps(mdmat, base_class)
    x(Temperature).calc("mdmat")
    x(Temperature).gibbs_single("mdmat", "diff.dat", _nrow=Temperature, _ncol=1)


def _ss_length(Temperature=1, base_class=fabfile.Calc):

    x = common_steps(test, base_class)

    x(Temperature).calc_dihedral_angle("ss_stretch") # this has to be run only once
    x(Temperature).ss_stretch("ss_stretch")
    x(Temperature).generate_data_frame_66_cluster("ss_stretch")
    x(Temperature).hist_length_plot("ss_stretch", _ncol=1, _nrow=1, types="helix", c_max=4,  st_time=800000)


def _h_bond_1(Temperature=1, base_class=fabfile.Calc, rep_list=range(0, fabfile.Calc.num_r)):
    """This step prepares the file for drawing in cytoscape"""
    x = common_steps(h_bond_1, base_class)
    x(Temperature).hist_length_pandas("all_pairs_8.5_ca", "contact_all.dat")
    x(Temperature).contacts_correlation_jaccard_domain("all_pairs_8.5_ca", st_frame=8000,given_contact=contact_list,diff=diff)
    x(Temperature).hist_length_plot("sample", st_time=800000)

def _cytoscape_draw(Temperature=1, base_class=fabfile.Calc):
    import cyto_draw
    x = common_steps(cyto_draw, base_class)
    rep_list = range(0, fabfile.Calc.num_r)
    x(Temperature).cyto_plotting_pandas("all_pairs_8.5_ca")


def _polystat(Temperature=1, base_class=fabfile.Calc, my_clus="67_only", clus_id=0, cont_class=False):
    x = common_steps(polystat, base_class)
    x(Temperature).calc_rij("polystat", st_time=800000)
    x(Temperature).fig_gen("polystat")




# code execution

for form_list, end_append_list, num_r in zip([["M_hid_11", "V_hid_11"]], [['117', '123']], [64]):

    all_initialization()
    fabfile.Calc.num_r = num_r
    fabfile.Calc.end_append_list = end_append_list
    fabfile.Calc.form_list = form_list
    fabfile.Calc.st_time = 0


    home = expanduser("~")
    contact_list = ['95_66', '93_90']
    _h_bond_1(Temperature=1, base_class=fabfile.Calc, rep_list=range(0, num_r))
