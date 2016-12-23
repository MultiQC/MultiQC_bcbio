#!/usr/bin/env python

""" MultiQC module to parse output from bcbio"""

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import re
from collections import defaultdict
import yaml

from multiqc import config, BaseMultiqcModule
from multiqc import plots
from multiqc_bcbio import srna

# Initialise the logger
log = logging.getLogger(__name__)

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

def linegraph(self, data, config):
    return plots.linegraph.plot(data, config)

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

        # detect mirna
        mirna_stats = srna.parse()
        self.bcbio_data.update(mirna_stats.general)
        # Write parsed report data to a file
        self.write_data_file(self.bcbio_data, 'multiqc_bcbio_metrics')

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.metrics_stats_table()

        # General target stats
        target_infos = list(self.find_log_files(config.sp['bcbio']['target']))
        if len(target_infos) == 1:
            add_project_info(yaml.load(target_infos[0]['f']))

        # Coverage plot
        # Only one section, so add to the intro

        self.sections = list()
        coverage_plot = self.bcbio_coverage_chart(config.sp['bcbio']['coverage'])
        coverage_avg_plot = self.bcbio_coverage_avg_chart(config.sp['bcbio']['coverage_avg'])
        variant_plot = self.bcbio_variants_chart(config.sp['bcbio']['variants'])
        qsignature_plot = None # disable plotting for now
        #qsignature_plot = self.bcbio_qsignature_chart(config.sp['bcbio']['qsignature'])
        variant_stats = self.bcbio_variants_stats(config.sp['bcbio']['vcfstats'])
        umi_stats = self.bcbio_umi_stats(config.sp["bcbio"]["umi"])

        if umi_stats:
            self.sections.extend(umi_stats)
        if mirna_stats.mirs:
            self.sections.append({'name': 'miRNAs stats',
                                  'anchor': 'bcbio-mirs',
                                  'content': mirna_stats.mirs})
        if mirna_stats.iso:
            self.sections.append({'name': 'Isomirs stats',
                                  'anchor': 'bcbio-isomirs',
                                  'content': mirna_stats.iso})
        if variant_stats:
            self.sections.append(variant_stats)
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
        if qsignature_plot:
            self.sections.append({
                'name': 'qSignature Profile',
                'anchor': 'bcbio-fraction-qsignature',
                'content': qsignature_plot})

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
        headers.update(srna.add_srna_headers(self.bcbio_data))

        if any(['Total_reads' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Total_reads'] = {
                'title': 'Reads',
                'description': 'Total sequences in the bam file',
                'min': 0,
                'modify': lambda x: x / 1000000,
                'shared_key': 'read_count',
                'format': '{:.2f} M',
            }
        if any(['Mapped_reads' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Mapped_reads'] = {
                'title': 'Mapped',
                'description': 'Mapped reads number',
                'min': 0,
                'modify': lambda x: x / 1000000,
                'shared_key': 'read_count',
                'format': '{:.2f} M',
                'hidden': True,
            }
        if any(['Mapped_reads_pct' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Mapped_reads_pct'] = {
                'title': '% Aln',
                'description': '% Mapped reads',
                'min': 0, 'max': 100, 'suffix': '%',
                'scale': 'RdYlGn',
                'format': '{:.1f}%',
            }
        if any(['Duplicates_pct' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Duplicates_pct'] = {
                'title': '% Dup',
                'description': '% Duplicated mapped reads',
                'min': 0, 'max': 100, 'suffix': '%',
                'scale': 'RdYlGn',
                'format': '{:.1f}%'
            }
        if any(['Ontarget_pct' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Ontarget_pct'] = {
                'title': '% On-trg',
                'description': '% On-target mapped not-duplicate reads',
                'min': 0, 'max': 100, 'suffix': '%',
                'scale': 'RdYlGn',
                'format': '{:.1f}%'
            }
        if any(['Ontarget_padded_pct' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Ontarget_padded_pct'] = {
                'title': '% On-trg+200bp',
                'description': '% Reads that overlap target regions extended by 200 bp. Expected to be 1-2% higher.',
                'min': 0, 'max': 100, 'suffix': '%',
                'scale': 'RdYlGn',
                'format': '{:.1f}%',
                'hidden': True,
            }
        if any(['Usable_pct' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Usable_pct'] = {
                'title': '% Usable',
                'description': '% Unique reads mapped on target in the total number of original reads.',
                'min': 0, 'max': 100, 'suffix': '%',
                'scale': 'RdYlGn',
                'format': '{:.1f}%'
            }
        if any(['Avg_coverage' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Avg_coverage'] = {
                'title': 'Depth',
                'description': 'Average target read coverage',
                'format': '{:.2f}',
            }

        if any(['Disambiguated_ambiguous_reads' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers.update(_get_disambiguited(self.bcbio_data))

        if any(['Variations_total' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Variations_total'] = {
                'title': 'Total Variations',
                'description': 'Numbers of Total Variations',
                'min': 0, 'format': '{:.0f}', 'share_key': 'variants'
            }
        if any(['Variations_in_dbSNP_pct' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Variations_in_dbSNP_pct'] = {
                'title': '% in dbSNP',
                'description': 'Numbers of Variations in dbSNP',
                'format': '{:.1f}%',
                'min': 0,
                'max': 100
            }
        if any(['Variations_heterozygous' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Variations_heterozygous'] = {
                'title': 'Variations heterozygous',
                'description': 'Numbers of Heterozygous Variations',
                'min': 0, 'format': '{d}', 'share_key': 'variants'
            }
        if any(['Variations_alt_homozygous' in self.bcbio_data[s] for s in self.bcbio_data]):
            headers['Variations_alt_homozygous'] = {
                'title': 'Variations homozygous',
                'description': 'Numbers of Homozygous Variations',
                'min': 0, 'format': '{d}', 'share_key': 'variants'
            }
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
                'format': '{:.1f}%'
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
            return linegraph(self, bcbio_data, {
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
            return linegraph(self, data, {
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

    def _bcbio_umi_count_plot(self, parsed_data):
        plot_data = {}
        for s, info in parsed_data.items():
            plot_data[s] = info["umi_counts"]
        config = {'xlab': "Reads per UMI", 'ylab': "Count",
                  "xDecimals": False}
        return {'name': 'UMI count distribution',
                'anchor': 'umi-stats-counts',
                'content': linegraph(self, [plot_data], config)}

    def _bcbio_umi_table(self, parsed_data):
        keys = OrderedDict()
        keys['umi_consensus_mapped'] = {'title': 'Consensus mapped',
                                        'description': 'Count of UMI consensus reads mapped',
                                        'format': '{:n}'}
        keys['umi_consensus_pct'] = {'title': 'Consensus reduction',
                                        'description': 'Percent of original reads removed by consensus',
                                        'format': '{:.1f}%'}
        keys['umi_baseline_mapped'] = {'title': "Original mapped",
                                        'description': 'Count of original mapped reads',
                                        'format': '{:n}'}
        keys['umi_baseline_duplicate_pct'] = {'title': 'Original duplicates',
                                                'description': 'Percentage original duplicates',
                                                'format': '{:.1f}%'}
        keys['umi_baseline_all'] = {'title': 'Original total',
                                    'description': 'Total reads in the original BAM',
                                    'format': '{:n}'}
        keys['umi_reduction_median'] = {'title': 'Duplicate reduction (median)',
                                        'description': 'Reduction in duplicates per position by UMIs (median)',
                                        'format': '{:n}x'}
        keys['umi_reduction_max'] = {'title': 'Duplicate reduction (max)',
                                        'description': 'Reduction in duplicates per position by UMIs (maximum)',
                                        'format': '{:n}x'}
        return {'name': 'UMI barcode statistics',
                'anchor': 'umi-stats',
                'content': plots.table.plot(parsed_data, keys)}

    def bcbio_variants_stats(self, names):
        """ Parsing stats from VCF files """

        parsed_data = defaultdict(dict)
        for f in self.find_log_files(names):
            with open(os.path.join(f['root'], f['fn'])) as in_handle:
                parsed_data.update(yaml.safe_load(in_handle))
        keys = OrderedDict()
        defaults = {
                'min': 0,
                'format': '{d}',
                'shared_key': 'variants'
                }
        keys['Variations (total)'] = dict(defaults, **{'description': 'Total variants detected'})
        keys['Variations (homozygous)'] = dict(defaults, **{'description': 'Total homozygous variants detected'})
        keys['Variations (alt homozygous)'] = dict(defaults, **{'description': 'Total alternative homozygous variants detected'})
        keys['Variations (heterozygous)'] = dict(defaults, **{'description': 'Total heterozygous variants detected'})
        keys['Variations (SNPs)'] = dict(defaults, **{'description': 'Total SNPs detected'})
        keys['Variations (indels)'] = dict(defaults, **{'description': 'Total indels detected'})
        keys['Variations (ts/tv)'] = dict(defaults, **{'description': 'TS/TV ratio'})

        if parsed_data:
            return {'name': 'Variant Summary Table (bcftools)',
                    'anchor': 'bcftools-stats',
                    'content': plots.table.plot(parsed_data, keys)}

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
                    {'name': 'Variants depth', 'ylab': 'pct of variants', 'xlab': 'number of reads', 'ymax':100},
                    {'name': 'CG content', 'ymax':100, 'ylab': 'pct of  variants', 'xlab': 'CG conetnt', 'xmin':0, 'xmax':100}
                    ]
        }

        if bcbio_data[0]:
            return linegraph(self, bcbio_data, config)

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

        return plots.heatmap.plot(hmdata, names)


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
