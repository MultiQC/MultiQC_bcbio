#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" MultiQC module to parse output from bcbio"""

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import re
from collections import defaultdict
import yaml

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, table, heatmap
from multiqc_bcbio import srna

# Initialise the logger
log = logging.getLogger('multiqc.multiqc_bcbio')

INTRO_COVERAGE = """
                    <p>This section shows the percentage of regions with a certain
                    level of nucleotides covered by a given number of reads:</p>
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

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='bcbio',
        anchor='bcbio', target='bcbio',
        href='http://github.com/bchapman/bcbio-nextgen/',
        info="bcbio-nextgen calculates "\
        "coverage over target regions and variants "\
        "if data is available during the analysis.")

        # Find and load any bcbio reports
        self.bcbio_data = dict()
        for f in self.find_log_files('bcbio/metrics'):
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

        # detect mirna
        mirna_stats = srna.parse()
        self.bcbio_data.update(mirna_stats.general)
        # Write parsed report data to a file
        self.write_data_file(self.bcbio_data, 'multiqc_bcbio_metrics')

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.metrics_stats_table()

        # General target stats
        target_infos = list(self.find_log_files('bcbio/target'))
        if len(target_infos) == 1:
            add_project_info(yaml.load(target_infos[0]['f']))

        # Coverage plot
        # Only one section, so add to the intro

        coverage_plot = self.bcbio_coverage_chart('bcbio/coverage')
        coverage_avg_plot = self.bcbio_coverage_avg_chart('bcbio/coverage_avg')
        qsignature_plot = None # disable plotting for now
        #qsignature_plot = self.bcbio_qsignature_chart('bcbio/qsignature')
        for umi_section in self.bcbio_umi_stats("bcbio/umi"):
            self.add_section(**umi_section)

        viral_stats = self.get_viral_stats("bcbio/viral")
        if viral_stats:
            self.add_section(**viral_stats)

        damage_stats = self.get_damage_stats("bcbio/damage")
        if damage_stats:
            self.add_section(**damage_stats)

        if mirna_stats.mirs:
            self.add_section(
                name='miRNAs stats',
                anchor='bcbio-mirs',
                plot=mirna_stats.mirs)
        if mirna_stats.iso:
            self.add_section(
                name='Isomirs stats',
                anchor='bcbio-isomirs',
                plot=mirna_stats.iso)
        if coverage_avg_plot:
            self.add_section(
                name='Coverage Profile',
                anchor='bcbio-fraction-coverage-all',
                description=INTRO_COVERAGE_AVG,
                plot=coverage_avg_plot)
        if coverage_plot:
            self.add_section(
                name='Coverage Profile Along Regions',
                anchor='bcbio-fraction-coverage',
                description=INTRO_COVERAGE,
                plot=coverage_plot)
        if qsignature_plot:
            self.add_section(
                name='qSignature Profile',
                anchor='bcbio-fraction-qsignature',
                plot=qsignature_plot)

    def parse_bcbio_report(self, raw_data):
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

        try:
            from multiqc_bcbio import srna
        except ImportError:
            log.error("Cannot import MultiQC_bcbio sRNA.")
        else:
            headers.update(srna.add_srna_headers(self.bcbio_data))

        read_format = '{:,.1f}&nbsp;' + config.read_count_prefix
        if config.read_count_multiplier == 1:
            read_format = '{:,.0f}'
            
        if any(['Total_reads' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Total_reads'] = {
                'title': 'Reads',
                'description': 'Total read count in the output BAM file (including unmapped, etc.)',
                'min': 0,
                'modify': lambda x: x * config.read_count_multiplier,
                'shared_key': 'read_count',
                'format': read_format,
            }
        if any(['Mapped_reads' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Mapped_reads'] = {
                'title': 'Aln',
                'description': 'Total number of read alignments',
                'min': 0,
                'modify': lambda x: x * config.read_count_multiplier,
                'shared_key': 'read_count',
                'format': read_format,
                'hidden': True,
            }
        if any(['Mapped_reads_pct' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Mapped_reads_pct'] = {
                'title': 'Aln',
                'description': '% Mapped reads',
                'min': 0,
                'max': 100,
                'suffix': '%',
                'scale': 'RdYlGn',
                'format': '{:,.1f}',
            }
        if any(['Duplicates_pct' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Duplicates_pct'] = {
                'title': 'Dup',
                'description': '% Duplicated reads',
                'min': 0, 'max': 100, 'suffix': '%',
                'scale': 'RdYlGn',
                'format': '{:,.1f}'
            }
        if any(['Ontarget_pct' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Ontarget_pct'] = {
                'title': 'On-trg',
                'description': '% On-target (both mates, primary) mapped not-duplicate reads',
                'min': 0, 'max': 100, 'suffix': '%',
                'scale': 'RdYlGn',
                'format': '{:,.1f}'
            }
        if any(['Ontarget_padded_pct' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Ontarget_padded_pct'] = {
                'title': '±200bp',
                'description': '% Reads that overlap target regions extended by 200 bp. Expected to be 1-2% higher.',
                'min': 0, 'max': 100, 'suffix': '%',
                'scale': 'RdYlGn',
                'format': '{:,.1f}',
                'hidden': True,
            }
        if any(['Usable_pct' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Usable_pct'] = {
                'title': 'Usable',
                'description': '% Unique reads mapped on target in the total number of original reads.',
                'min': 0, 'max': 100, 'suffix': '%',
                'scale': 'RdYlGn',
                'format': '{:,.1f}'
            }
        if any(['Avg_coverage' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Avg_coverage'] = {
                'title': 'Cov',
                'description': 'Average target read coverage',
                'format': '{:,.1f}',
            }
        if any(['Average_insert_size' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Average_insert_size'] = {
                'title': 'Mean IS',
                'description': 'Mean insert size (samtools stats)',
                'format': '{:.0f}',
            }
        if any(['Disambiguated_ambiguous_reads' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers.update(_get_disambiguited(self.bcbio_data))

        if any(['rRNA_rate' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['rRNA_rate'] = {
                'title': 'rRNA pct',
                'description': '% alignments to rRNA. Depending on the library preparation methods used, the proportion '
                               'of rRNA sequences should be quite low. If large proportions of rRNA sequences are seen, '
                               'it is wise to consider if the depth of the remaining sequences is sufficient for further analyses.',
                'max': 100,
                'min': 0,
                'modify': lambda x: x * 100,
                'scale': 'RdYlGn',
                'format': '{:,.1f}'
            }
        if len(headers.keys()):
            self.general_stats_addcols(self.bcbio_data, headers)

    def bcbio_coverage_chart(self, names):
        """ Make the bcbio assignment rates plot """

        parsed_data = defaultdict(dict)
        seen = set()
        for f in self.find_log_files(names):
            s_name = self.clean_s_name(f['fn'], root=None)
            for line in f['f'].split("\n"):
                if not line.startswith("percent"):
                    continue
                cutoff_reads, region_pct, bases_pct, sample = line.split("\t")
                x = 100 - float(region_pct)
                y = float(bases_pct)
                if s_name not in parsed_data[cutoff_reads]:
                    parsed_data[cutoff_reads][s_name] = []
                parsed_data[cutoff_reads][s_name].append((x, y))
                seen.add(s_name)
            if s_name in seen:
                self.add_data_source(f)

        bcbio_data = []
        cutoffs = []
        for pct_key in sorted(parsed_data.keys(), key=lambda k: int(k.split("percentage")[1])):
            if any(any(v > 0 for v in dict(d).values()) for d in parsed_data[pct_key].values()):
                data_obj = {}
                for s in parsed_data[pct_key]:
                    data_obj.update({s: dict(parsed_data[pct_key][s])})
                bcbio_data.append(data_obj)
                cutoffs.append(int(pct_key.split("percentage")[1]))

        if bcbio_data and bcbio_data[0] and cutoffs:
            return linegraph.plot(bcbio_data, {
                'data_labels': [
                    {'name': str(c) + 'x'} for c in cutoffs
                ],
                'id': 'bcbio_coverage_plot',
                'title': 'Completeness',
                'xlab': '% regions',
                'ylab': '% bases in the regions covered',
                'ymin': 0,
                'ymax': 100,
            })

    def bcbio_coverage_avg_chart(self, names):
        """ Make the bcbio assignment rates plot """

        x_threshold = 0
        data = defaultdict(dict)
        for f in self.find_log_files(names):
            s_name = self.clean_s_name(f['fn'], root=None)
            for line in f['f'].split("\n"):
                if not line.startswith("percentage"):
                    continue
                cutoff_reads, bases_pct, sample = line.split("\t")
                y = float(bases_pct)
                x = int(cutoff_reads.replace("percentage", ""))
                data[s_name][x] = y
                if y > 1.0:
                    x_threshold = max(x_threshold, x)

            if s_name in data:
                self.add_data_source(f)

        if data:
            return linegraph.plot(data, {
                "xlab": "number of reads",
                "ylab": '% bases in the regions covered',
                "xmax": x_threshold,
            })

    def bcbio_umi_stats(self, names):
        """Prepare table of statistics on UMIs in the file.
        """
        parsed_data = defaultdict(dict)
        for f in self.find_log_files(names):
            with open(os.path.join(f['root'], f['fn'])) as in_handle:
                parsed_data.update(yaml.safe_load(in_handle))
        if len(parsed_data) > 0:
            umi_table = self._bcbio_umi_table(parsed_data)
            umi_count_plot = self._bcbio_umi_count_plot(parsed_data)
            return [umi_table, umi_count_plot]
        return []

    def _bcbio_umi_count_plot(self, parsed_data):
        plot_data = {}
        for s, info in parsed_data.items():
            plot_data[s] = info["umi_counts"]
        config = {'xlab': "Reads per UMI", 'ylab': "Count",
                  "xDecimals": False}
        return {'name': 'UMI count distribution',
                'anchor': 'umi-stats-counts',
                'plot': linegraph.plot([plot_data], config)}

    def _bcbio_umi_table(self, parsed_data):
        keys = OrderedDict()
        keys['umi_consensus_mapped'] = {'title': 'Consensus mapped',
                                        'description': 'Count of UMI consensus reads mapped',
                                        'format': '{:n}'}
        keys['umi_consensus_pct'] = {'title': 'Consensus reduction',
                                     'description': 'Percent of original reads removed by consensus',
                                     'suffix': '%',
                                     'format': '{:,.1f}'}
        keys['umi_baseline_all'] = {'title': 'Orig. total',
                                    'description': 'Total reads in the original BAM',
                                    'format': '{:n}'}
        keys['umi_baseline_mapped'] = {'title': "Orig. mapped",
                                       'description': 'Count of original mapped reads',
                                       'format': '{:n}'}
        keys['umi_baseline_duplicate_pct'] = {'title': 'Orig. dup',
                                              'description': 'Percentage original duplicates',
                                              'suffix': '%',
                                              'format': '{:,.1f}'}
        keys['umi_reduction_median'] = {'title': 'Dup. reduction (median)',
                                        'description': 'Reduction in duplicates per position by UMIs (median)',
                                        'suffix': 'x',
                                        'format': '{:n}'}
        keys['umi_reduction_max'] = {'title': 'Dup. reduction (max)',
                                     'description': 'Reduction in duplicates per position by UMIs (maximum)',
                                     'suffix': 'x',
                                     'format': '{:n}'}
        return {'name': 'UMI barcode statistics',
                'anchor': 'umi-stats',
                'plot': table.plot(parsed_data, keys)}

    def get_viral_stats(self, fnames):
        """Provide counts of top viral hits for samples.
        """
        to_show = 5
        data = {}
        for f in self.find_log_files(fnames):
            with open(os.path.join(f['root'], f['fn'])) as in_handle:
                sample_name = in_handle.readline().strip().split()[-1]
                counts = []
                for line in in_handle:
                    contig, count = line.strip().split("\t")
                    counts.append((count, contig))
                counts.sort(reverse=True)
                if counts:
                    data[sample_name] = {"counts": ", ".join(["%s (%s)" % (v, c) for (c, v) in counts[:to_show]])}
        keys = OrderedDict()
        keys["counts"] = {"title": "Virus (count)",
                          "description": "Top %s viral sequences, with counts, found in unmapped reads" % to_show}
        if data:
            return {"name": "Viral mapping read counts",
                    "anchor": "viral-counts",
                    "plot": table.plot(data, keys)}

    def get_damage_stats(self, fnames):
        """Summarize statistics on samples with DNA damage.
        """
        data = {}
        keys = set([])
        for f in self.find_log_files(fnames):
            with open(os.path.join(f['root'], f['fn'])) as in_handle:
                cur = yaml.safe_load(in_handle)
                keys = keys | set(cur["changes"].keys())
                data[cur["sample"]] = cur["changes"]
        if data:
            cols = OrderedDict()
            for k in sorted(list(keys), reverse=True):
                cols[k] = {"title": k, "format": "{:n}"}
            return {"name": "DNA damage and bias filtering",
                    "anchor": "damage-stats",
                    "plot": table.plot(data, cols)}

    def bcbio_qsignature_chart(self, names) :
        """ Make the bcbio assignment rates plot """

        hmdata = list()
        data = defaultdict(dict)
        for f in self.find_log_files(names):
            s_name = self.clean_s_name(f['fn'], root=None)
            for l in f['f'].splitlines():
                cols = l.strip().split()
                data[cols[0]][cols[1]] = float(cols[2])
                data[cols[1]][cols[0]] = float(cols[2])
                data[cols[0]][cols[0]] = 0
                data[cols[1]][cols[1]] = 0

        names = data.keys()
        for name in names:
            row = list()
            for name2 in names:
                row.append(data[name][name2])
            hmdata.append(row)

        return heatmap.plot(hmdata, names)


def _get_disambiguited(dt):
    """Get headers for disamb."""
    h = OrderedDict()
    for s in dt:
        for k in dt[s]:
            if k.startswith("Disambiguated"):
                h[k] = {
                    'title': k.replace("_", " ").replace("Disambiguated", "Disamb."),
                    'description': 'When samples are at risk of cross-species contamination (e.g. those '
                                   'derived from PDXs), an attempt is performed to remove reads that are '
                                   'assigned to the different species involved. This metric shows the number '
                                   'of removed reads.',
                    'min': 0,
                    'format': '{:.0f}',
                    'shared_key': 'read_count',
                    'scale': 'RdYlGn',
                }
    return h


def add_project_info(data):
    genome_info = data.get('genome_info')  # {name, size}
    coverage_bed_info = data.get('coverage_bed_info')
    variants_regions_info = data.get('variants_regions_info')

    if genome_info:
        config.report_header_info = []
        config.report_header_info.append({"Genome:": genome_info['name'] + ', ' + _format_decimal(genome_info['size'], 'bp')})
        if coverage_bed_info and coverage_bed_info['bed'] == variants_regions_info['bed']:
            config.report_header_info.append({"Target:": _format_bed_info(variants_regions_info, genome_info)})
        else:
            if coverage_bed_info:
                config.report_header_info.append({"Target for coverage:": _format_bed_info(coverage_bed_info, genome_info)})
            if variants_regions_info and variants_regions_info['bed']:
                config.report_header_info.append({"Target for var. calling:": _format_bed_info(variants_regions_info, genome_info)})


def _format_bed_info(d, genome_info):
    bed, size, regions, genes = d['bed'], \
                                _format_decimal(d.get('size')), \
                                _format_decimal(d.get('regions')),\
                                _format_decimal(d.get('genes'))
    bed_name = os.path.basename(bed)
    html = '<a href={bed}>{bed_name}</a>'
    if size is not None:
        percent = (100.0 * d['size'] / genome_info['size']) if genome_info.get('size') else 0
        html += ' ({size} bp'
        if percent >= 0.01:
            html += ' or {percent:.2f}% of the genome'
        if regions is not None:
            html += ', {regions} region' + ('s' if d.get('regions') != 1 else '')
        if genes is not None:
            html += ', {genes} gene' + ('s' if d.get('genes') != 1 else '')
        html += ')'
    return html.format(**locals())


def _format_decimal(value, unit=None):
    if value is None:
        return None
    unit_str = ('<span style="font-size: 50%; line-height: 1;">&nbsp;</span>' + unit) if unit else ''
    if value <= 9999:
        return '{value}{unit_str}'.format(**locals())
    else:
        v = '{value:,}{unit_str}'.format(**locals())
        return v.replace(',', '<span style="margin-left: 0.2em;"></span>')
