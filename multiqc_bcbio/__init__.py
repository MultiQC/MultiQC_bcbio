from __future__ import absolute_import
from distutils.version import StrictVersion

from .bcbio import MultiqcModule
import multiqc
from multiqc import config
from multiqc.utils.config import update_dict
from multiqc.utils.report import get_filelist

update_dict(
    config.sp, {'bcbio/metrics': {'fn': '*_bcbio.txt'},
                'bcbio/coverage': {'fn': '*_bcbio_coverage.txt'},
                'bcbio/coverage_avg': {'fn': '*_bcbio_coverage_avg.txt'},
                'bcbio/variants': {'fn': '*_bcbio_variants.txt'},
                'bcbio/target': {'fn': 'target_info.yaml'},
                'bcbio/qsignature': {'fn': '*bcbio_qsignature.ma'},
                'bcbio/vcfstats': {'fn': '*_bcbio_variants_stats.txt'},
                'bcbio/seqbuster': {'contents': 'seqbuster'},
                'bcbio/umi': {'fn': '*_umi_stats.yaml'},
                'bcbio/viral': {'fn': '*viral*-counts.txt'},
                'bcbio/damage': {'fn': '*damage.yaml'},
                })


config.fn_clean_exts.append({'type': 'regex', 'pattern': '_bcbio.*'})
get_filelist()


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

