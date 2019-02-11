import networkx as nx
import subprocess
import fabfile
import time
import numpy as np
from matplotlib.transforms import Bbox
import pandas as pd
from cyto_api import*

new_defaults = {

    # Network defaults
    'NETWORK_BACKGROUND_PAINT': 'white'
}


def get_my_base(base):
    class Hb_cyto(base):
        """makes the conatct networks in cytoscape"""

        def cyto_plotting_pandas(self, func_dir_name, cg_len=1, cut=4, rep_n=2, ss_type='sbsb'):
            # this combines all the temperature replicas togther
            """params: cg_len: coarse grained residues length
            cut: cutoff values for the maps to be loaded in cytoscape
            rep_n: cutoff value for number of replica a particular contact can exist in
            ss_type: ss_type file in the input file
            """

            for form in self.form_list:

                self.form = form
                for xyz in ['helix_6_93']:

                #for xyz in ['frames_66_32_sidechain', 'frames_66_92_sidechain', 'frames_65_39_sidechain', 'frames_31_70_sidechain', 'frames_62_66_sidechain', 'frames_77_93_sidechain']:

                    #input_file = "%s/%s/%s/pd_xyz_%s" % (self.output, self.form, func_dir_name, xyz) # saves the pandas format output
                    #input_file = "%s/%s/%s/pd_xyz_0" % (self.output, self.form, func_dir_name) # saves the pandas format output
                    input_file = "%s/%s/%s/pd_xyz_3a_3b_0" % (self.output, self.form, func_dir_name) # saves the pandas format output
                    #input_file = "%s/%s/%s/pd_xyz_0_domain_8_cytoscape_not_divided" % (self.output, self.form, func_dir_name) # saves the pandas format output
                    fname = fabfile.put_file(input_file)
                    print "loading following form"
                    print self.form

                    cy = CyRestClient(ip='127.0.0.1', port=1234)

                    f = fname

                    G = nx.read_edgelist(f, nodetype=int, delimiter="\t", create_using=nx.MultiDiGraph(), data=[
                                         ("weight", float), ("no_replica", float), ("diff", int)])
                    g_cy = loadcy(G, cy)
                    g_cy = nodecolortable(g_cy)
                    style1 = vis_sty(cy, name='A_contact_pairs')
                    g_cy = sapply(style1, g_cy, cy)
                    # style1 = update(new_defaults, style1)
                    response = cytoscape("layout", "attributes-layout", {"NodeAttribute": "key"})
                    response = cytoscape("view", "fit content")
                    response = cytoscape("network", "rename", {"name": "%s%s" % ('norm', form)})
                    #response = cytoscape("view", "export", method="HELP") 

                    time.sleep(2)

            return

        def cyto_plotting_pandas_diff(self, func_dir_name, cg_len=1, cut=4, rep_n=2, ss_type='sbsb'):
            # this combines all the temperature replicas togther
            """params: cg_len: coarse grained residues length
            cut: cutoff values for the maps to be loaded in cytoscape
            rep_n: cutoff value for number of replica a particular contact can exist in
            ss_type: ss_type file in the input file
            """

            form = self.form_list[0]
            self.form = form
            input_file = "%s/%s/%s/pd_xyz_diff_3a_3b_%s_%s_0" % (self.output, self.form, func_dir_name, self.form_list[0], self.form_list[1]) # saves the pandas format output
            #input_file = "%s/%s/%s/pd_xyz_diff_%s_%s_0" % (self.output, self.form, func_dir_name, self.form_list[0], self.form_list[1]) # saves the pandas format output
            print input_file
            fname = fabfile.put_file(input_file)

            cy = CyRestClient(ip='127.0.0.1', port=1234)

            f = fname

            G = nx.read_edgelist(f, nodetype=int, delimiter="\t", create_using=nx.MultiDiGraph(), data=[
                                 ("weight", float), ("diff", float), ("domain-domain", str)])
            g_cy = loadcy(G, cy)
            g_cy = nodecolortable(g_cy)
            style1 = vis_sty(cy, name='A_contact_pairs')
            g_cy = sapply(style1, g_cy, cy)
            # style1 = update(new_defaults, style1)
            response = cytoscape("layout", "attributes-layout", {"NodeAttribute": "key"})
            response = cytoscape("view", "fit content")
            response = cytoscape("network", "rename", {"name": "%s%s" % ('diff', form)})

            time.sleep(2)

            return