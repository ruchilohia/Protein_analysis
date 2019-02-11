import subprocess
from subprocess import Popen, PIPE  # noqa
from functools import wraps
import pandas as pd
import paramiko
import matplotlib
import inspect
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import colorsys
import os
import errno
from os.path import expanduser
from os.path import exists as path_exists
from os.path import join as path_join

"""This is the base file defining generic properties for all other files"""

font = {
    'weight': 'bold',
    'size': 20.0}

axes = {'labelweight': 'bold',
        'grid': True,
        'titlesize': 30.0,
        'labelsize': 20.0,
        'titleweight': 'bold'}

lines = {'linewidth': 10,
         }

figure = {'facecolor': '0.75',
          'titlesize': 25.0,
          'titleweight': 'bold'}

legend = {'fontsize': 12,
          'labelspacing': 0.25}

grid = {'linewidth': 0.5
        # ,'alpha' : 0.2
}

plt.rc('font', **font)
plt.rc('axes', **axes)
plt.rc('lines', **lines)
#plt.rc ('figure' , **figure)
plt.rc('grid', **grid)
plt.rc('legend', **legend)

#params = {'legend.fontsize': 20,
#          'legend.handlelength': 1.2}
#plot.rcParams.update(params)

def put_file(path):
    """opens a remote file and makes a local copy on my mac"""
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    #ssh.connect("taranis", username="lohia")
    ssh.connect("kestrel.ccib.rutgers.edu", username="rl487")
    sftp = ssh.open_sftp()
    try:
        sftp.open(path, 'r')
    except IOError:
        pass
    f = sftp.get(path, "/Users/lohia/remote_del/tmp")
    ssh.close()
    return "/Users/lohia/remote_del/tmp"


def get_file(sour, dest):
    """copies the local file to the remote server"""
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    # ssh.connect("taranis", username="lohia")
    ssh.connect("kestrel.ccib.rutgers.edu", username="rl487")
    sftp = ssh.open_sftp()
    try:
        sftp.open(sour, 'r')
    except IOError:
        pass
    f = sftp.put(sour, dest)
    ssh.close()


def temprature(n, t_low, t_high):

    return [int((t_low * 1.0) * ((t_high * 1.0) / t_low) ** (x / float(n - 1)) for x in range(0, n))]


def node_color(gro):
        """color of each node for scatter plot depending on the residue type
        :param gro: .gro file path
        :returns: color and resid list
        """
        #gro = "/g3/home/lohia/T-rex-short-region/trex-01/output/V/0/resid_23-113-capped.gro"
        f1 = open("tmp", 'w')
        f1.write("#!/usr/bin/tclsh\n")
        f1.write("mol new %s  waitfor -1\n" % gro)
        f1.write("set f [open hp "'"w"'"]\n")
        f1.write(
            "puts $f "'"[[atomselect top  "'"hydrophobic and type CA"'" ] get resid]"'"\n")
        f1.write("close $f\n")
        f1.write("set f [open aci "'"w"'"]\n")
        f1.write(
            "puts $f "'"[[atomselect top  "'"acidic and type CA"'" ] get resid]"'"\n")
        f1.write("close $f\n")
        f1.write("set f [open base "'"w"'"]\n")
        f1.write(
            "puts $f "'"[[atomselect top  "'"basic and type CA"'" ] get resid]"'"\n")
        f1.write("close $f\n")
        f1.write("set f [open xaxis "'"w"'"]\n")
        f1.write(
            "puts $f "'"[[atomselect top  "'"type CA"'" ] get resid]"'"\n")
        f1.write("close $f\n")
        f1.write("exit\n")
        f1.close()
        subprocess.call("vmd -e tmp", shell=True)
        color = []
        hp = np.loadtxt("hp")
        aci = np.loadtxt("aci")
        base = np.loadtxt("base")
        xaxis = np.loadtxt("xaxis")
        for i in [int(x) for x in xaxis]:
            if i in aci:
                color.append("red")
            elif i in base:
                color.append("blue")
            elif i in hp:
                color.append("white")
            else:
                color.append("green")
        return color, [int(x) for x in xaxis]


def node_color(res_name_list, reslist):
        """color of each node for scatter plot depending on the residue type"""
        color = []
        hp = ['MET', 'VAL', 'TRP', 'PHE', 'ILE', 'LEU', 'ALA', 'TYR', 'HID']
        base = ['ARG', 'LYS', 'HIP']
        aci = ['GLU', 'ASP']
        for i in res_name_list:
            if i in aci:
                color.append("red")
            elif i in base:
                color.append("blue")
            elif i in hp:
                color.append("white")
            else:
                color.append("green")
        return color, [int(x) for x in reslist]


def node_prop(gro):
        """color of each node for scatter plot depending on the residue type
        :param gro: .gro file path
        :returns: color and resid list

        It also outputs helicity,hdrophobicity and amino acid name
        """
        f1 = open("tmp", 'w')
        f1.write("#!/usr/bin/tclsh\n")
        f1.write("mol new %s  waitfor -1\n" % gro)
        f1.write("set f [open xaxis "'"w"'"]\n")
        f1.write(
            "puts $f "'"[[atomselect top  "'"type CA"'" ] get resid]"'"\n")
        f1.write("close $f\n")
        f1.write("set f [open resname "'"w"'"]\n")
        f1.write(
            "puts $f "'"[[atomselect top  "'"type CA"'" ] get resname]"'"\n")
        f1.write("close $f\n")
        f1.write("set f [open hp "'"w"'"]\n")
        f1.write(
            "puts $f "'"[[atomselect top  "'"hydrophobic and type CA"'" ] get resid]"'"\n")
        f1.write("close $f\n")
        f1.write("set f [open aci "'"w"'"]\n")
        f1.write(
            "puts $f "'"[[atomselect top  "'"acidic and type CA"'" ] get resid]"'"\n")
        f1.write("close $f\n")
        f1.write("set f [open base "'"w"'"]\n")
        f1.write(
            "puts $f "'"[[atomselect top  "'"basic and type CA"'" ] get resid]"'"\n")
        f1.write("close $f\n")
        f1.write("exit\n")
        f1.close()
        subprocess.call("vmd -e tmp", shell=True)
        color = []
        xaxis = np.loadtxt("xaxis")
        hp = np.loadtxt("hp")
        aci = np.loadtxt("aci")
        base = np.loadtxt("base")
        my_resname = np.loadtxt("resname", dtype='str')

        helicity = {'ALA': 0, 'LEU': 0.21, 'ARG': 0.21, 'MET': 0.24, 'LYS': 0.26, 'GLN': 0.39, 'GLU': 0.40, 'ILE': 0.41, 'TRP': 0.49,
                    'SER': 0.50, 'TYR': 0.53, 'PHE': 0.54, 'VAL': 0.61, 'HID': 0.61, 'ASN': 0.65, 'THR': 0.66, 'CYS': 0.68, 'ASP': 0.69, 'GLY': 1, 'PRO': 1}

        hydrophobicity = {'ALA': 1.800, 'LEU': 3.8, 'MET': 1.900, 'ILE': 4.5, 'TYR': -1.3, 'PHE': 2.8, 'VAL': 4.2, 'PRO': -1.6, 'GLY': -0.4, 'THR': -0.7, 'TRP':-0.9, 'SER':-0.8, 'HIS': -3.2, 'GLU': -3.5, 'GLN': -3.5, 'ASP':-3.5, 'ASN':-3.5 , 'LYS': -3.9, 'ARG': -4.5}

        amino = {'ALA': 'A', 'LEU': 'L', 'ARG': 'R', 'MET': 'M', 'LYS': 'K', 'GLN': 'Q', 'GLU': 'E', 'ILE': 'I', 'TRP': 'W', 'SER': 'S',
                 'TYR': 'Y', 'PHE': 'F', 'VAL': 'V', 'HID': 'H', 'ASN': 'N', 'THR': 'T', 'CYS': 'C', 'ASP': 'D', 'GLY': 'G', 'PRO': 'P'}

        for i, res in zip([int(x) for x in xaxis], my_resname):
            if i in aci:
                color.append("red")
                print i, helicity[res], 100.0, amino[res]
            elif i in base:

                print i, helicity[res], -100.0, amino[res]
            elif i in hp:
                print i, helicity[res], hydrophobicity[res], amino[res]
            else:
                print i, helicity[res], -50.0, amino[res]

        return color, [int(x) for x in xaxis]


def line_color(parent_cluster_no, parent_cluster_id):
    """every parent replica is given color on the basis of hot to cold
    :param parent_cluster_no: total no of temperatures
    :param parent_cluster_id: present temperature id
    :returns : the color for a given parent_cluser_id

    """

    if parent_cluster_no > 1:
        n = parent_cluster_no * 1.0
        HSV_tuples = [(x * 0.7 / n, 1, 1) for x in range(parent_cluster_no)]
        RGB_tuples_temp = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
        RGB_tuples = []
        for i in reversed(RGB_tuples_temp):
            RGB_tuples.append(i)
        return RGB_tuples[parent_cluster_id]
    else:
        return "red"


def line_width(parent_cluster_no, parent_cluster_id):
    """every parent replica is given color on the basis of hot to cold"""
    print "******"
    if parent_cluster_no > 1:
        n = parent_cluster_no
        m_lw = [5.0 * (5.0 / 5.0) ** (x / float(n - 1)) for x in range(0, n)]
        RGB_tuples = []
        for i in reversed(m_lw):
            RGB_tuples.append(i)
        return RGB_tuples[parent_cluster_id]
        # return m_lw[parent_cluster_id]
    else:
        return 10


def line_alpha(parent_cluster_no, parent_cluster_id):
    """every parent replica is given transparency on the basis of hot to cold"""
    print "******"
    if parent_cluster_no > 1:
        n = parent_cluster_no
        m_lw = [0.1 * (1.0 / 0.1) ** (x / float(n - 1)) for x in range(0, n)]
        RGB_tuples = []
        for i in reversed(m_lw):
            RGB_tuples.append(i)
        return RGB_tuples[parent_cluster_id]
    else:
        return 10

"""remove the archive once the genralized script works"""

def _decorator_calc_loop(func):
        """
        creates a loop for every form and every parent_cluster no and every replica id
        generates path for inp/out files and directories
        *** I still don't understand the wraps function
        """
        @wraps(func)
        def common_func(self, func_dir_name, *args, **kwargs):
            for form, end_append in zip(self.form_list, self.end_append_list):
                self.form = form
                self.end_append = end_append
                for s in range(self.parent_cluster_no):
                    self.parent_cluster_id = s
                    func_dir_form = self.output_func_dir(func_dir_name)
                    self.func_dir_form = func_dir_form
                    if not path_exists(func_dir_form):
                        os.makedirs(func_dir_form)
                    subprocess.call(
                        "rm -rf  \#*", cwd=func_dir_form, shell=True)
                    #for i in range(self.calc_replica):
                    #for i in range(46, 47):
                    #for i in [31,63]:
                    for i in [0]:
                        self.i = i
                        self.xtc = self.xtc_file(self.xtc_file_name)
                        self.tpr = (
                            "/".join(map(str, [self.output, self.form, self.i, self.tpr_file_name])))
                        self.gro = (
                            "/".join(map(str, [self.output, self.form, self.i, self.gro_file_name])))
                        self.repl_dir = (
                            "/".join(map(str, [self.func_dir_form, self.i])))
                        self.repl_dir_parent = self.repl_dir_parent = self.repl_dir_parent_method(func_dir_name)
                        if not path_exists(self.repl_dir):
                            os.makedirs(self.repl_dir)
                        subprocess.call("mkdir %s" % self.repl_dir, shell=True)
                        func(self, func_dir_name, *args, **kwargs)
        return common_func

class Calc(object):

    """with this class we assingn all the inpt variables
    param: num_r : total number of replica,this will always remain fixed
           calc_replica : total number of replcas, this is always changing in every loop
           xtc_freq : xtc frequecny in picoseconds
           end_append_list: id of the last part which is being appended
           xtc_file_name: name of the xtc file on  which all further calulcations will be performed
            """
    home = expanduser("~")
    dropbox_path = "%s/Dropbox/lab_work/output_02/" % home
    dropbox_path_mac = "%s/Dropbox/lab_work/output_02/" % home
    git_dir = "/u1/home/lohia/Dropbox/replex_md/analysis_01/"
    git_dir_mac = "/Users/lohia/Dropbox/replex_md/analysis_01/"
    output = "/g3/home/lohia/T-rex-NMR-region/trex-02/output/"
    input_path = "/g3/home/lohia/T-rex-NMR-region/trex-02/"
    gro_file_name = "resid_23-113-capped.gro"
    num_r = 60
    xtc_freq = 50
    st_time = 100000
    end_time = 250000
    resid_st = 23
    resid_end = 113
    form_list = ["V", "M"]
    end_append_list = [50, 50]
    nmr_list = ["19358.str.txt", "19357.str.txt"]
    xtc_file_name = "trajout.xtc"
    # xtc file name before ramining minidst frames
    tpr_file_name = "ex.tpr"
    cont_file_folder = "cont_traj"
    cont_file_output = "rep_id"
    temp_file_name = "rex_replica_temp.xvg"
    calc_replica = 60
    color_list = ["black", "green", "blue", "red",
                  "cyan", "magenta", "peru", "slategray", "olive"]
    # this is arbitary, since there are no parent cluster for this class
    parent_cluster_no = 1
    class_name = "Calc"
    temp_ticks = [float(300.0 * (420 / 300.0) ** (x / float(60 - 1))) for x in range(0, 60)]
    cutoff = 1
    clustername = "none"
    cont_file_folder = "cont_traj"
    cont_file_output = "rep_id"
    angle_values = "dpca_phi_psi_omega"

    def __init__(self, calc_replica, cont_class=False):  # cont_class is only a place holder for now. it makes sense for cluster class
        self.fs = self.st_time / self.xtc_freq
        self.ft = self.end_time / self.xtc_freq
        self.tot_frames = (self.ft - self.fs) + 1
        self.calc_replica = calc_replica
        self.mylist_formlist = [
            [0 for i in range(self.calc_replica)] for j in self.form_list]
        self.tot_rep = calc_replica
        self.mylist_empty = [0 for i in range(self.calc_replica)]
        self.mylist_empty_2 = [0 for i in range(self.calc_replica)]
        self.form_dict = {}
        self.pandas_csv_files_dir = self.output + "/" + "pandas_csv_files/"

    def __call__(self):
        print "dropbox_path is %s" % self.dropbox_path
        print "dropbox_path_mac is %s " % self.dropbox_path_mac
        print "git_dir is %s " % self.git_dir
        print "git_dir mac is %s" % self.git_dir_mac
        print "output is %s " % self.output
        print "input_path is %s " % self.input_path
        print "gro_file is %s " % self.gro_file_name
        print "num_r is %s " % self.num_r
        print "xtc_freq is %s " % self.xtc_freq
        print "st_time is %s " % self.st_time
        print "end_time is %s " % self.end_time
        print "resid_st is %s " % self.resid_st
        print "resid_end is %s " % self.resid_end
        print "form_list is %s " % self.form_list
        print "calc_replica is %s " % self.calc_replica

    def pre_analysis(self):
        """concatenates all the tajectories and makes them as whole
        update the part number for the outputs
        the trajectories start at 0 and end at the given end time"""
        for form, end_append in zip(self.form_list, self.end_append_list):
            inp = "%s/%s" % (self.input_path, form)

            subprocess.call("echo 1 > pre_analysis", cwd="%s" %
                            inp, shell=True)
            for i in range(self.calc_replica):
                inp2 = "%s/%s/%s/" % (self.output, form, i)
                # subprocess.call("gmx trjcat -f ex.part0{001..022}.xtc ex.part0{024..%s}.xtc -e %s -o cat_out.xtc  ; rm -rf \#*" % (
                #    end_append, self.end_time), shell=True, cwd="%s/%s" % (inp, i))  # this is just for M_hip_11 pre_analysis
                subprocess.call("gmx trjcat -f ex.part0{001..%s}.xtc -e %s -o cat_out.xtc  ; rm -rf \#*" % (
                    end_append, self.end_time), shell=True, cwd="%s/%s" % (inp, i))
                subprocess.call(
                    "gmx trjconv -s ex.tpr -f cat_out.xtc -pbc whole -o cat_out_whole.xtc < ../pre_analysis ; rm -rf \#* ; rm -rf cat_out.xtc ", shell=True, cwd="%s/%s" % (inp, i))
                subprocess.call("mkdir %s ; mkdir %s/%s ; mkdir %s ; cp cat_out_whole.xtc %s/%s ; cp %s %s ; cp ex.tpr  %s " % (self.output, self.output, form, inp2, inp2, self.xtc_file_name_before_mindist, self.gro_file_name, inp2, inp2), shell=True, cwd="%s/%s" % (inp, i))
                subprocess.call("rm -rf cat_out_whole.xtc", shell=True, cwd="%s/%s" % (inp, i))

    def energy_cat(self):
        """concatenates energy file"""
        for form, end_append in zip(self.form_list, self.end_append_list):
            inp = "%s/%s" % (self.input_path, form)

            subprocess.call("echo 1 > pre_analysis", cwd="%s" %
                            inp, shell=True)
            for i in range(self.calc_replica):
                inp2 = "%s/%s/%s/" % (self.output, form, i)
                subprocess.call("gmx eneconv -f ex.part00{01..%s}.edr -dt 50 -b %s -e %s -o all_energy.edr " % (
                    end_append, self.st_time, self.end_time), shell=True, cwd="%s/%s" % (inp, i))
                subprocess.call("cp all_energy.edr %s/all_energy.edr " %
                                (inp2), shell=True, cwd="%s/%s" % (inp, i))

    def initial_struct_energy(self):
        """finds initial energy for each structure"""
        for form in self.form_list:
            inp = "%s/%s" % (self.input_path, form)

            subprocess.call("echo 1 > pre_analysis", cwd="%s" %
                            inp, shell=True)
            for i in range(self.calc_replica):
                inp2 = "%s/%s/%s/" % (self.output, form, i)
                subprocess.call("cp equil.part0001.edr %s/equil.part0001.edr" %
                                (inp2), shell=True, cwd="%s/%s" % (inp, i))

    def output_form(self, form, filetype='/'):
        # gives the path of output and form, we use form instead of self.form because it is being called from another class beofre self.form is declared
        # it will be a good idea to get rid of form and replace it with self.from in all the methods
        _file = (
            "/".join(map(str, [self.output, form, self.cont_file_output, filetype])))
        return _file

    def output_form_temperature(self, form, i, filetype='/'):
        # gives the path of output and form and temperature
        _file = (
            "/".join(map(str, [self.output, form, i, self.cont_file_output, filetype])))
        return _file

    def xtc_file(self, filetype):
        """filetype = name of the xtc file
        gives the xtc ,tpr , cont_traj , gro_file for given clus_id"""
        _file = (
            "/".join(map(str, [self.output, self.form, self.i, filetype])))
        return _file

    def output_func_dir(self, filetype):
        """creates and returns the location for the directory of the function folders"""
        output_func = ("/".join(map(str, [self.output, self.form, self.cont_file_output, filetype])))

        return output_func

    def fig_path_method(self, func_dir_name):
        """creates and returns the location for the directory of the figure folders"""
        self.fig_path = (
            "/".join(map(str, [self.dropbox_path_mac, "figures", self.cont_file_output, func_dir_name])) + "/")
        if not path_exists(self.fig_path):
            os.makedirs(self.fig_path)
        #self.figname = "%s/replica_%s_clus_%s_form_%s.pdf" % (self.fig_path, self.replica_curr, self.clus_curr, self.form)
        return self.fig_path

    def vmd_file(self, f1):
        """loading gro and xtc file in vmd without equilibration"""
        #f1 = open("tmp", 'w')
        f1.write("#!/usr/bin/tclsh\n")
        f1.write("mol new %s  waitfor -1\n" % self.gro)
        f1.write("mol addfile %s type {xtc}  first %s waitfor -1\n" % (self.xtc, str(self.fs)))
        return f1

    def angle_values_directory(self):
        """creates a directory and stores dihedral angles"""
        if not path_exists(self.output_func_dir(self.angle_values)):
            os.mkdir(self.output_func_dir(self.angle_values))
        self.angle_values_dir = self.output_func_dir(self.angle_values) + "/" + str(self.replica_curr)
        if not path_exists(self.angle_values_dir):
            os.mkdir(self.angle_values_dir)
        return self.angle_values_dir

    def _get_angle_index(self):
        """calculates phi psi and omega values for every frame, like pre_analysis, it is immune to start time ( self.st_time) """
        tclfile = "index"
        repl_dir = self.angle_values_directory()
        f1 = open("tmp", 'w')
        f1 = self.vmd_file(f1)
        f1.write("source %s/dPCA/%s.tcl\n" % (self.git_dir, tclfile))
        f1.write("%s %s %s %s \n " % (tclfile, repl_dir, self.resid_st, (self.resid_end - 1)))
        f1.write("exit\n")
        f1.close()
        subprocess.call("vmd -e tmp", shell=True)
        subprocess.call("cat diangle_phi.ndx > diangle.ndx", shell=True, cwd=repl_dir)
        subprocess.call("echo '[phi]' | cat - diangle.ndx > temp && mv temp diangle.ndx", shell=True, cwd=repl_dir)
        subprocess.call("gmx angle -f %s -n diangle.ndx -ov dihedral -type dihedral -all" % self.xtc, shell=True, cwd=repl_dir)
        subprocess.call("tail -n +17 dihedral.xvg > dihedral_phi ; rm -rf \#*", shell=True, cwd=repl_dir)
        subprocess.call("cat diangle_psi.ndx  > diangle.ndx", shell=True, cwd=repl_dir)
        subprocess.call("echo '[psi]' | cat - diangle.ndx > temp && mv temp diangle.ndx", shell=True, cwd=repl_dir)
        subprocess.call("gmx angle -f %s -n diangle.ndx -ov dihedral -type dihedral -all" % self.xtc, shell=True, cwd=repl_dir)
        subprocess.call("tail -n +17 dihedral.xvg  > dihedral_psi ; rm -rf \#*", shell=True, cwd=repl_dir)
        subprocess.call("cat diangle_omega.ndx  > diangle.ndx", shell=True, cwd=repl_dir)
        subprocess.call("echo '[omega]' | cat - diangle.ndx > temp && mv temp diangle.ndx", shell=True, cwd=repl_dir)
        subprocess.call("gmx angle -f %s -n diangle.ndx -ov dihedral -type dihedral -all" % self.xtc, shell=True, cwd=repl_dir)
        subprocess.call("tail -n +17 dihedral.xvg  > dihedral_omega ; rm -rf \#*", shell=True, cwd=repl_dir)

        return repl_dir

    def load_angle_values(self):
        repl_dir = self.angle_values_directory()
        self.psi = (np.loadtxt("%s/dihedral_psi" % repl_dir)[:, 2::])
        self.phi = (np.loadtxt("%s/dihedral_phi" % repl_dir)[:, 2::])
        self.omega = (np.loadtxt("%s/dihedral_omega" % repl_dir)[:, 2::])
        return

#    def get_angle_index_loadfile(self):

    def repl_dir_parent_method(self, func_dir_name):
        """creates and returns the location for the directory of the function folders with temperature"""
        repl_dir_parent = ("/".join(map(str, [self.output, self.form, self.cont_file_output, func_dir_name, str(self.i)])))
        return repl_dir_parent

    def population_calc(self, clusname="None"):
        "loads properties which are differnt in clusters"

        self.replica_curr = self.i
        self.clus_curr = 0
        self.tot_frames = self.demux_file_method(self.form).shape[0] + 1  # demux file has the actual 0th frame missing, it has 1 less frame
        self.population = self.tot_frames
        self.frame_list = [i + 1  for i in range(0, self.tot_frames)]
        self.frame_list_index = [i for i in range(0, self.tot_frames)]
        self.curr_temp_list = [self.temp_ticks[self.replica_curr]] * self.tot_frames
        self.time_frame = range(0, self.tot_frames)
        Calc.pt_file_method(self, clusname)

    def demux_file_method(self, form):
        """gives the output for demux file"""
        path_tmp = path_join(self.output, form, "demux/0/replica_index_50_s.xvg")
        if not path_exists(path_tmp):
            path_tmp = put_file(path_join(self.output, form, "demux/0/replica_index_50_s.xvg"))
        self.demux_file = np.loadtxt(path_tmp)
        return self.demux_file

    def demux_frames_index_dict(self, curr_temperature, replica_index):

        #  outputs all the frames INDEX belonging to a particular replica at a given temeperature.
        #  this is valid when the output frame is written from time 0 from the trajout.xtc file, if the output is written from vmd file be careful - the gro file should not be included.
        # use this when the indexes have to be matched
        #  param: demux_f: replica index file

        return [i + 1 for i, j in enumerate(self.demux_file[:, curr_temperature + 1]) if j == replica_index]  # the frame index is shifted, frame order at time 0 in index file is the index for frame 1 in

    def demux_frames_dict(self, curr_temperature, replica_index):

        #  outputs all the frames belonging to a particular replica at a given temeperature.
        #  this is valid when the gromacs analysis tool is used and the output frmae is written from time 0
        #  this is used when the gromacs frames at time 0 are written starting from 1
        # use this when the frames have to be matched
        #  param: demux_f: replica index file

        return [i + 2 for i, j in enumerate(self.demux_file[:, curr_temperature + 1]) if j == replica_index]  # +1 is added two time, 1: naturally shifted, we add it all the time, 2: the frame numbers counting in the contact file is starting from 1, hence we have to add 1 as well.

    def replica_index_dict(self, curr_temperature, time_frame):

        # this will return a replica_id for a given time frame
        # the input time frame should start from 0 without any clusters
        if time_frame == 0:
            return curr_temperature  # the 0 time frame is the actual replica we start with
        else:

            return self.demux_file[time_frame - 1, curr_temperature + 1]  # 1 is substracted to get the correct replica

    def demux_file_per_temperature(self):

        # this will return a replica_id for a given time frame
        # the input time frame should start from 0 without any clusters
        my_path = self.output_form_temperature(self.form, self.replica_curr)
        replica_id = open("%s/replica_id.out" % my_path, "w")
        _list = self.demux_file[:, 0 + 1]
        _list2 = np.insert(_list, 0, 0)
        replica_id.write("\n".join(map(str, _list2)))
        replica_id.close()
        return _list2

    def pt_file_method(self, clusname="None"):
        cluster_dir = Calc.output_form(self, self.form, clusname)
        pt_file = ("/".join(map(str, [cluster_dir, self.i, "pt.dat"])))
        if not path_exists(pt_file):
            try:
                _pt = put_file(pt_file)
            except IOError:
                self.clus_list_index = [0] * self.tot_frames
        else:
            self.clus_list_index = np.loadtxt(_pt, skiprows=1)
        return

class Clusters(Calc):
    """
    this is for making clusters
    we don't need to change the start time for clusters , since the gromacs calculation with -b option already doesnt find them ( because the trajcout find we are using doesn't have them initial time step probably). however for vmd calculations they have to be updated or we can use a new method
    param: calc_replica: is the number of clusters for which each property is calculated this is equilvaent to total no of clusters in the parent file
    parent_cluster_no: Total number of temperature
    It always inherits from Calc class, the constructor takes the variable cont_class, if yoou want to make clusters on contnious replcias
    """

    st_time = 0
    parent_cluster_id = 0
    class_name = "Clusters"
    cutoff = 4
    clustername = "67"

    def __init__(self, parent_cluster_no, cont_class=False):
        self.cont_class = cont_class
        if cont_class:   # if you want the clusters to be calculated for continiuos trajectories
                self.cont_file_folder = "cont_traj"
                self.cont_file_output = "rep_id"
        else:
                self.cont_file_folder = "/"
                self.cont_file_output = "/"
        self.parent_cluster_no = parent_cluster_no
        self.calc_replica = self.cutoff
        self.fs = self.st_time / self.xtc_freq
        self.ft = self.end_time / self.xtc_freq
        self.mylist_formlist = [
            [0 for i in range(self.parent_cluster_no)] for j in self.form_list]
        self.tot_rep = parent_cluster_no
        self.mylist_empty = [0 for j in range(self.parent_cluster_no)]
        self.mylist_empty2 = [
            [0 for i in range(self.cutoff)] for j in self.form_list]
        self.mylist_empty_2 = [0 for i in range(self.parent_cluster_no)]
        self.mylist_empty3 = [
            [0 for i in range(self.cutoff)] for j in self.form_list]
        # self.mylist_formlist_tuple = [([0 for j in range(self.parent_cluster_no)],[0 for j in range(self.parent_cluster_no)],[0 for j in range(self.parent_cluster_no)],[0 for j in range(self.parent_cluster_no)],[0 for j in range(self.parent_cluster_no)]),([0 for j in range(self.parent_cluster_no)],[0 for j in range(self.parent_cluster_no)],[0 for j in range(self.parent_cluster_no)],[0 for j in range(self.parent_cluster_no)],[0 for j in range(self.parent_cluster_no)])]

    def pre_pre_analysis(self):
        """creates the cluster_dir/parent_id directory"""
        for form in self.form_list:
            for self.parent_cluster_id in range(self.parent_cluster_no):

                cluster_dir = Calc.output_form(self, form, filetype=self.clustername)
                subprocess.call(" if [ ! -d %s ] ; then mkdir %s ; fi " %
                                (cluster_dir, cluster_dir), shell=True)

                # making cluster directory with parent id

                cluster_dir_output = cluster_dir + "/" + str(self.parent_cluster_id)
                subprocess.call(" if [ ! -d %s ] ; then mkdir %s ; fi " %
                                (cluster_dir_output, cluster_dir_output), shell=True)
                print "making"
                print cluster_dir_output

    def cluster_prop_calc(self, func_dir_name, load_file_name, clusname):
        # func_dir_name,load_file_name are the name of the parent cluster properties files
        # this will only work if the property is calculated for every frame
        # saves the property for each cluster from the parent cluster
        # this is customized for calculating secondary struture properties
        for form in self.form_list:
            for self.parent_cluster_id in range(self.parent_cluster_no):
                for i in range(self.calc_replica):
                    cluster_dir = Calc.output_form(self, form, filetype=self.clustername)
                    cluster_dir_output = cluster_dir + "/" + str(self.parent_cluster_id)
                    subprocess.call(
                        "mkdir %s/%s" % (cluster_dir_output, func_dir_name), shell=True)
                    subprocess.call(
                        "mkdir %s/%s/%s" % (cluster_dir_output, func_dir_name, i), shell=True)
                    f = np.loadtxt("%s/pt.dat" % cluster_dir_output)
                    # +1 for i bcoz the clus_id is stored from 1 , #the resid no. at the begning of the file thus +1 is added to the frame no
                    frame_no = [k + 1 + self.fs for k,
                                j in enumerate(f) if j == i + 1]
                    property_dir = Calc.output_form(self, form, filetype=func_dir_name)
                    property_path = "/%s/%s" % (
                        property_dir, self.parent_cluster_id, load_file_name)
                    subprocess.call("mkdir %s/%s" %
                                    (cluster_dir_output, i), shell=True)
                    property_array = np.loadtxt(property_path, dtype='str')

                    a = property_array[:, frame_no]
                    np.savetxt("%s/%s/%s/%s" % (cluster_dir_output,
                                                func_dir_name, i, load_file_name), a, fmt='%s')

    def cluster_prop_calc_rg(self, func_dir_name, load_file_name, clusname):
        # func_dir_name,load_file_name are the name of the parent cluster properties files
        # this will only work if the property is calculated for every frame
        # saves the property for each cluster from the parent cluster
        # this is customized for rg calculation for clusters
        for form in self.form_list:
            for self.parent_cluster_id in range(self.parent_cluster_no):
                for i in range(self.calc_replica):
                    cluster_dir = Calc.output_form(self, form, filetype=self.clustername)
                    cluster_dir_output = cluster_dir + "/" + str(self.parent_cluster_id)
                    subprocess.call(
                        "mkdir %s/%s" % (cluster_dir_output, func_dir_name), shell=True)
                    subprocess.call(
                        "mkdir %s/%s/%s" % (cluster_dir_output, func_dir_name, i), shell=True)
                    f = np.loadtxt("%s/pt.dat" % cluster_dir_output)
                    # +1 for i bcoz the clus_id is stored from 1
                    frame_no = [k + self.fs for k,
                                j in enumerate(f) if j == i + 1]
                    property_dir = Calc.output_form(self, form, filetype=func_dir_name)
                    property_path = "/%s/%s/%s" % (
                        property_dir, self.parent_cluster_id, load_file_name)
                    subprocess.call("mkdir %s/%s" %
                                    (cluster_dir_output, i), shell=True)
                    property_array = np.loadtxt(property_path)
                    print i
                    a = property_array[frame_no, :]
                    np.savetxt("%s/%s/%s/%s" % (cluster_dir_output,
                                                func_dir_name, i, load_file_name), a, fmt='%s')

    def pre_analysis(self, gro_dir, cluster_gro_index):
        """creates all the cluster xtc files"""
        for form in self.form_list:
            for self.parent_cluster_id in range(self.parent_cluster_no):
                cluster_dir = Calc.output_form(self, form, filetype=self.clustername)
                # making cluster directory with parent id
                cluster_dir_output = cluster_dir + "/" + str(self.parent_cluster_id)
                subprocess.call("mkdir %s" % (cluster_dir), shell=True)
                subprocess.call("mkdir %s" % (cluster_dir_output), shell=True)
                parent_path = (
                    "/".join([self.output, form, str(self.parent_cluster_id)]))
                subprocess.call(
                    "cp %s/%s ." % (parent_path, self.tpr_file_name), shell=True, cwd=cluster_dir_output)
                subprocess.call(
                    "cp %s/%s ." % (parent_path, self.gro_file_name), shell=True, cwd=cluster_dir_output)
                if self.cont_class:
                    xtc_file = ("/".join(map(str, [self.output, form, self.cont_file_folder])) + "/%s_%s" % (str(self.parent_cluster_id), self.xtc_file_name))
                else:
                    xtc_file = ("/".join(map(str, [self.output, form, str(self.parent_cluster_id), self.xtc_file_name])))
                _out_path = Calc.output_form(self, form)
                subprocess.call("echo 1 > pro ; gmx trjconv -f %s -s %s -sub %s < pro ; rm -rf \#*" % (xtc_file, self.tpr_file_name, ("/".join(map(str, [_out_path, gro_dir, self.parent_cluster_id, cluster_gro_index])))), shell=True, cwd=cluster_dir_output)

    def xtc_file(self, filetype):
        """gives the xtc file for given clus_id"""
        cluster_dir = Calc.output_form(self, self.form, filetype=self.clustername)
        _file = (("/".join(map(str, [cluster_dir, self.parent_cluster_id, self.i + 1]))) + ".xtc")
        return _file

    def output_func_dir(self, filetype):
        """creates and returns the location for the directory of the function folders"""
        cluster_dir = Calc.output_form(self, self.form, filetype=self.clustername)
        output_func = (
            "/".join(map(str, [cluster_dir, self.parent_cluster_id, filetype])))

        return output_func

    def fig_path_method(self, func_dir_name):
        "creates and returns the location for the directory of the function folders"
        tmp = "/".join(map(str, [self.dropbox_path_mac, "figures", self.cont_file_output, self.clustername]))
        if not path_exists(tmp):
            os.makedirs(tmp)
        tmp = "/".join(map(str, [self.dropbox_path_mac, "figures", self.cont_file_output,
                                 self.clustername, self.parent_cluster_id]))
        if not path_exists(tmp):
            os.makedirs(tmp)
        self.fig_path = (
            "/".join(map(str, [self.dropbox_path_mac, "figures", self.cont_file_output, self.clustername, self.parent_cluster_id, func_dir_name])) + "/")
        self.figname = "%s/%s_replica%s_clus%s_%s.pdf" % (self.fig_path, func_dir_name, self.replica_curr, self.clus_curr, self.form)

        return

    def vmd_file(self, f1):
        """loading gro and xtc file in vmd without equilibration"""
        # f1 = open("tmp", 'w')
        f1.write("#!/usr/bin/tclsh\n")
        f1.write("mol new %s  waitfor -1\n" % self.gro)
        f1.write("mol addfile %s type {xtc} waitfor -1\n" % self.xtc)
        return f1

    def population_calc(self):
        "calculates the perecentage population of each cluster for reweighting"
        self.replica_curr = self.parent_cluster_id
        self.clus_curr = self.i
        self.demux_file_method(self.form)
        _pt = self.output_func_dir("pt.dat")
        if not path_exists(_pt):
            _pt = put_file(self.output_func_dir("pt.dat"))
        pt = np.loadtxt(_pt, skiprows=1)
        self.tot_frames = len(pt.tolist())
        self.population = ((pt == (self.clus_curr + 1)).sum())
        self.frame_list = [i + self.fs + 1 for i, j in enumerate(pt) if int(j) == self.clus_curr + 1]
        self.frame_list_index = [i + self.fs for i, j in enumerate(pt) if int(j) == self.clus_curr + 1]
        self.clus_list_index = pt
        self.curr_temp_list = self.temp_ticks[self.replica_curr] * self.tot_frames
        self.time_frame = range(0, self.tot_frames)

    def repl_dir_parent_method(self, func_dir_name):
        """creates and returns the location for the directory of the function folders with temeprature"""
        repl_dir_parent = ("/".join(map(str, [self.output, self.form, self.cont_file_output, func_dir_name, str(self.parent_cluster_id)])))
        return repl_dir_parent


class Sub_clusters(Clusters):
    """
    this is for making subclusters from clusters
    we don't need to change the start time for clusters , since the gromacs calculation with -b option already doesnt find them ( because the trajcout find we are using doesn't have them initial time step probably). however for vmd calculations they have to be updated or we can use a new method
    """

    st_time = 0
    end_time = 250000

    def __init__(self, calc_replica, clustername, top_cluster_dir, parent_cluster_no, parent_cluster_id):
        "calc_replica is the number of clusters for which each property is calculated"
        self.calc_replica = calc_replica
        self.clustername = clustername
        self.parent_cluster_id = parent_cluster_id
        # self.clustername/self.parent_cluster_id from the above cluster
        self.top_cluster_dir = top_cluster_dir
        # total no of parents clusters for each of them the subcluter properties are calculated
        self.parent_cluster_no = parent_cluster_no
        self.fs = self.st_time / self.xtc_freq
        self.ft = self.end_time / self.xtc_freq

    def pre_analysis(self, gro_dir, cluster_gro_index):
        for form in self.form_list:
            cluster_dir = "%s/%s/%s/%s/" % (self.output,
                                            form, self.top_cluster_dir, self.clustername)
            subprocess.call("mkdir %s" % (cluster_dir), shell=True)
            # making cluster directory with parent id
            cluster_dir_output = cluster_dir + str(self.parent_cluster_id)
            subprocess.call("mkdir %s" % (cluster_dir_output), shell=True)

            parent_path = ("/".join([self.output, form, self.top_cluster_dir]))

            subprocess.call("cp %s/%s ." % (parent_path,
                                            self.tpr_file_name), shell=True, cwd=cluster_dir_output)
            subprocess.call("cp %s/%s ." % (parent_path,
                                            self.gro_file_name), shell=True, cwd=cluster_dir_output)

            subprocess.call("echo 1 > pro ; g_trjconv -f %s/%s.xtc -s %s -sub %s < pro" % (parent_path, str(self.parent_cluster_id + 1), self.tpr_file_name, ("/".join(
                map(str, [self.output, form, self.top_cluster_dir, gro_dir, self.parent_cluster_id, cluster_gro_index])))), shell=True, cwd=cluster_dir_output)

    def xtc_file(self, filetype):
        """gives the xtc file for given clus_id"""
        _file = (("/".join(map(str, [self.output, self.form, self.top_cluster_dir,
                                     self.clustername, self.parent_cluster_id, self.i + 1]))) + ".xtc")
        return _file

    def output_func_dir(self, filetype):
        "creates and returns the location for the directory of the function folders"
        output_func = ("/".join(map(str, [self.output, self.form, self.top_cluster_dir,
                                          self.clustername, self.parent_cluster_id, filetype])))

        return output_func

    def fig_path_method(self, func_dir_name, home_dir):
        "creates and returns the location for the directory of the function folders"
        subprocess.call("mkdir %s" % (
            "/".join(map(str, [home_dir, "figures", self.top_cluster_dir, self.clustername]))), shell=True)
        subprocess.call("mkdir %s" % (
            "/".join(map(str, [home_dir, "figures", self.top_cluster_dir, self.clustername, self.parent_cluster_id]))), shell=True)
        output_func = ("/".join(map(str, [home_dir, "figures", self.top_cluster_dir,
                                          self.clustername, self.parent_cluster_id, func_dir_name])) + "/")
        self.figname = "%s/%s_%s_%s_%s.pdf" % (self.fig_path, func_dir_name, self.replica_curr, self.clus_curr, self.form)

        return output_func
if __name__ == "__main__":
    # node_prop()
    #print Calc(1).tot_frames
    Calc(60).energy_cat()
