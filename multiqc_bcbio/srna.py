#!/usr/bin/env python
"""
small rna-seq custom reports
"""

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import re
from collections import defaultdict
import yaml

from multiqc import config, BaseMultiqcModule
from multiqc import plots


class parse(BaseMultiqcModule):

    def __init__(self):
        self.mirs = None
        self.iso = None
        mirbase = self.bcbio_mirna_stats()
        self.general = self.bcbio_general_stats()

    def bcbio_general_stats(self):
        bcbio_data = list()
        fns = self.find_log_files(config.sp['bcbio']['seqbuster'])
        data = defaultdict(dict)
        for f in fns:
            s_name = self.clean_s_name(f['fn'], root=None)
            with open(os.path.join(f['root'], f['fn'])) as in_handle:
                for line in in_handle:
                    cols = line.strip().split()
                    if cols[0] == "mirs":
                        data[s_name]["miRNAs"] = int(cols[1])
                    if cols[0] == "isomirs":
                        data[s_name]["isomiRs"] = int(cols[1])
        return data

    def bcbio_mirna_stats(self):
        bcbio_data = list()
        fns = self.find_log_files(config.sp['bcbio']['seqbuster'])
        mirs_data = defaultdict(dict)
        mirs_key = OrderedDict()
        iso_data = defaultdict(dict)
        iso_key = OrderedDict()
        for f in fns:
            s_name = self.clean_s_name(f['fn'], root=None)
            with open(os.path.join(f['root'], f['fn'])) as in_handle:
                for line in in_handle:
                    cols = line.strip().split()
                    if line.startswith("mirs_"):
                        mirs_key[cols[0]] = {'name': cols[0].replace("_", " ")}
                        mirs_data[s_name][cols[0]] = int(cols[1])
                    if line.startswith("iso_"):
                        iso_key[cols[0]] = {'name': cols[0].replace("_", " ")}
                        iso_data[s_name][cols[0]] = int(cols[1])
        self.write_data_file(mirs_data, "seqbuster_mirs")
        self.write_data_file(iso_data, "seqbuster_isomirs")
        if mirs_data:
            cnfg = {'ylab': '# of miRNAs'}
            cnfg['title'] = "Number of miRNAs with changes"
            self.mirs = plots.bargraph.plot(mirs_data, mirs_key, cnfg)
        if iso_data:
            cnfg = {'ylab': '# of isomiRs'}
            cnfg['title'] = "Number of isomiRs with changes"
            self.iso = plots.bargraph.plot(iso_data, iso_key, cnfg)


def add_srna_headers(dt):
    h = OrderedDict()
    for s in dt:
        for k in dt[s]:
            if k.find("miRNAs") > -1 or k.find("isomiRs") > -1:
                h[k] = {'title': k, 'min': 0, 'format': '{d}'}
    return h
