from __future__ import absolute_import
from distutils.version import StrictVersion

from .bcbio import MultiqcModule
import multiqc
from multiqc import config

config.sp['bcbio'] = {'metrics': {'fn': '*_bcbio.txt'},
                      'coverage': {'fn': '*_bcbio_coverage.txt'},
                      'coverage_avg': {'fn': '*_bcbio_coverage_avg.txt'},
                      'variants': {'fn': '*_bcbio_variants.txt'},
                      'target': {'fn': 'target_info.yaml'},
                      'qsignature': {'fn': '*bcbio_qsignature.ma'},
                      'vcfstats': {'fn': '*_bcbio_variants_stats.txt'},
                      'seqbuster': {'contents': 'seqbuster'},
                      'umi': {'fn': '*_umi_stats.yaml'},
                      }

config.fn_clean_exts.append({'type': 'regex', 'pattern': '_bcbio.*'})



for module, value_dict in {
    'FastQC': {
        'percent_duplicates': False,
        'total_sequences': False,
    },
    'QualiMap': {
        'percentage_aligned': False,
    },
    'Samtools Stats': {
        'non-primary_alignments': False,
        'reads_mapped': False,
        'reads_mapped_percent': False,
        'raw_total_sequences': False,
    }}.items():

    if module not in config.table_columns_visible:
        config.table_columns_visible[module] = dict()
    config.table_columns_visible[module].update(value_dict)

