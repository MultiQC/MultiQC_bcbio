#!/usr/bin/env python

""" MultiQC module to parse output from bcbio"""

from __future__ import print_function
from collections import OrderedDict
import logging
import re
from collections import defaultdict

from multiqc import config, BaseMultiqcModule
try:
    from multiqc import plots
except ImportError:
    plots = None

# Initialise the logger
log = logging.getLogger(__name__)

INTRO_COVERAGE = """
                    <p>This section shows the percentage of regions with a certain
                    level of nucleotides covered by a given number of reads.</p>
                 """
INTRO_COVERAGE_AVG = """
                    <p>This section shows the percentage of nucleotides (y-axis)
                    covered by a given number of reads (x-axis).</p>
                 """
INTRO_VARIANT = """
                    <p>This section shows the percentage of variatns (y-axes) that a)
                    are covered by a givin number of reads or more (x-axis), or b)
                    %of GC content (x-axis).</p>
                """

def linegraph(self, data, config):
    if not plots:
        log.info("Version < 0.6")
        return self.plot_xy_data(data, config)
    else:
        return plots.linegraph.plot(data, config)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='bcbio',
        anchor='bcbio', target='bcbio',
        href='http://github.com/bchapman/bcbio-nextgen/',
        info="bcbio-nextgen makes calcultation of"\
        "coverage over regions and variants"\
        "if data is available during the analysis.")

        # Find and load any bcbio reports
        self.bcbio_data = dict()
        for f in self.find_log_files(config.sp['bcbio']['metrics']):
            parsed_data = self.parse_bcbio_report(f['f'])
            if parsed_data is not None:
                s_name = self.clean_s_name(f['fn'], root=None)
                if s_name in self.bcbio_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                self.add_data_source(f)
                self.bcbio_data[s_name] = parsed_data

        if len(self.bcbio_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.bcbio_data)))

        # Write parsed report data to a file
        self.write_data_file(self.bcbio_data, 'multiqc_bcbio_metrics')

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.metrics_stats_table()

        # Coverage plot
        # Only one section, so add to the intro

        self.sections = list()
        coverage_plot = self.bcbio_coverage_chart(config.sp['bcbio']['coverage'])
        coverage_avg_plot = self.bcbio_coverage_avg_chart(config.sp['bcbio']['coverage_avg'])
        variant_plot = self.bcbio_variants_chart(config.sp['bcbio']['variants'])
        if coverage_avg_plot:
            self.sections.append({
                'name': 'Coverage Profile',
                'anchor': 'bcbio-fraction-coverage-all',
                'content': INTRO_COVERAGE_AVG + coverage_avg_plot})
        if coverage_plot:
            self.sections.append({
                'name': 'Coverage Profile Along Regions',
                'anchor': 'bcbio-fraction-coverage',
                'content': INTRO_COVERAGE + coverage_plot})
        if variant_plot:
            self.sections.append({
                'name': 'Variant Profile',
                'anchor': 'bcbio-fraction-variants',
                'content': INTRO_VARIANT + variant_plot})


    def parse_bcbio_report (self, raw_data):
        """ Parse the bcbio log file. """
        parsed_data = {}
        for line in raw_data.split("\n"):
            fields = line.split("\t")
            if len(line) == 0:
                continue
            try:
                parsed_data[fields[0].strip()] = float(fields[1].strip())
            except:
                pass
        if len(parsed_data) == 0: return None
        return parsed_data


    def metrics_stats_table(self):
        """ Take the parsed stats from the bcbio report and add them to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers['offtarget_pct'] = {
            'title': '% Off-targets',
            'description': '% Off-targets reads',
            'max': 100,
            'min': 0,
            'modify': lambda x: x * 100,
            'scale': 'RdYlGn',
            'format': '{:.1f}%'
        }
        headers['Duplicates_pct'] = {
            'title': '% Dups',
            'description': '% of duplicated reads from samtools stats',
            'max': 100,
            'min': 0,
            'modify': lambda x: x * 100,
            'scale': 'RdYlGn',
            'format': '{:.1f}%'
        }
        headers['avg_coverage_per_region'] = {
            'title': 'Average coverage',
            'description': 'Average coverage per target',
            'scale': 'RdYlGn',
            'format': '{:.2f}%'
        }
        headers['Variations_total'] = {
            'title': 'Total Variations',
            'description': 'Numbers of Total Variations',
        }
        headers['Variations_in_dbSNP_pct'] = {
            'title': '% in dbSNP',
            'description': 'Numbers of Variations in dbSNP',
            'format': '{:.1f}%'
        }
        headers['Transition/Transversion'] = {
            'title': 'Ts/Tv',
            'description': 'Transition/Transversion Ratio',
            'scale': 'RdYlGn',
            'format': '{:.2f}%'
        }
        headers['Variations_heterozygous'] = {
            'title': 'Variations heterozygous',
            'description': 'Numbers of Heterozygous Variations',
        }
        headers['Variations_homozygous'] = {
            'title': 'Variations homozygous',
            'description': 'Numbers of Homozygous Variations',
        }
        if any(['avg_coverage_per_region' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['avg_coverage_per_region'] = {
                    'title': 'Avg Depth',
                    'description': 'Average read coverage on region'
                    }
        self.general_stats_addcols(self.bcbio_data, headers)


    def bcbio_coverage_chart (self, names) :
        """ Make the bcbio assignment rates plot """

        bcbio_data = []
        parsed_data = defaultdict(dict)
        seen = set()
        for f in self.find_log_files(names):
            s_name = self.clean_s_name(f['fn'], root=None)
            for line in f['f'].split("\n"):
                if not line.startswith("percent"):
                    continue
                fields = line.split("\t")
                x = 100 - float(fields[1])
                y = float(fields[2])
                if s_name not in parsed_data[fields[0]]:
                    parsed_data[fields[0]][s_name] = []
                parsed_data[fields[0]][s_name].append((x, y))
                seen.add(s_name)
            if s_name in seen:
                self.add_data_source(f)
        for ipct in [1, 5, 10, 20, 40, 50, 60, 70, 80, 100]:
            pct = "percentage%s" % ipct
            data_obj = {}
            for s in parsed_data[pct]:
                data_obj.update({s: dict(parsed_data[pct][s])})
            bcbio_data.append(data_obj)

        # Config for the plot
        config = {
                'data_labels': [
                    {'name': 'cutoff 1 nts'},
                    {'name': 'cutoff 5 nts'},
                    {'name': 'cutoff 10 nts'},
                    {'name': 'cutoff 20 nts'},
                    {'name': 'cutoff 40 nts'},
                    {'name': 'cutoff 50 nts'},
                    {'name': 'cutoff 60 nts'},
                    {'name': 'cutoff 70 nts'},
                    {'name': 'cutoff 80 nts'},
                    {'name': 'cutoff 100 nts'},
                    ],
            'id': 'bcbio_coverage_plot',
            'title': 'completeness',
            'ylab': '% nts in the regions covered',
            'xlab': '% regions',
            'ymin': 0,
            'ymax': 100,
        }

        if bcbio_data:
            return linegraph(self, bcbio_data, config)

    def bcbio_coverage_avg_chart (self, names) :
        """ Make the bcbio assignment rates plot """

        bcbio_data = list()
        parsed_data_depth = defaultdict(dict)
        for f in self.find_log_files(names):
            s_name = self.clean_s_name(f['fn'], root=None)
            for line in f['f'].split("\n"):
                if not line.startswith("percentage"):
                    continue
                fields = line.split("\t")
                y = float(fields[1])
                x = float(fields[0].replace("percentage", ""))
                parsed_data_depth[s_name][x] = y

            if s_name in parsed_data_depth:
                self.add_data_source(f)

        bcbio_data.append(parsed_data_depth)

        # Config for the plot
        config = {
                'xlab': "number of reads",
                'ylab': "pct of region covevered"
        }

        if bcbio_data:
            return linegraph(self, bcbio_data, config)

    def bcbio_variants_chart (self, names) :
        """ Make the bcbio assignment rates plot """

        bcbio_data = list()
        parsed_data_depth = defaultdict(dict)
        parsed_data_cg = defaultdict(dict)
        for f in self.find_log_files(names):
            s_name = self.clean_s_name(f['fn'], root=None)
            data_cg = defaultdict(float)
            for line in f['f'].split("\n"):
                if line.startswith("pct"):
                    continue
                fields = line.split("\t")
                if len(fields) < 3:
                    continue
                y = 100 - float(fields[0])
                x_depth = float(fields[1])
                x_cg = float(fields[2])
                parsed_data_depth[s_name][x_depth] = y
                data_cg[x_cg] = float(fields[0])
            init_cg = 0
            for cg in sorted(data_cg, key=data_cg.get):
                corrected = data_cg[cg] - init_cg
                init_cg = data_cg[cg]
                parsed_data_cg[s_name][cg] = corrected

            if s_name in parsed_data_depth:
                self.add_data_source(f)

        bcbio_data.append(parsed_data_depth)
        bcbio_data.append(parsed_data_cg)

        # Config for the plot
        config = {
                'ylab': 'pct of variants',
                'xlab': 'number of reads | CG content',
                'data_labels': [
                    {'name': 'Reads depth', 'ymax':100},
                    {'name': 'CG content', 'ymax':100, 'xmin':0, 'xmax':100}
                    ]
        }

        if bcbio_data:
            return linegraph(self, bcbio_data, config)
