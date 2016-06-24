from __future__ import absolute_import
from distutils.version import StrictVersion

from .bcbio import MultiqcModule
import multiqc
from multiqc import config

config.sp['bcbio'] = {'metrics': {'fn': '*_bcbio.txt'},
                      'coverage': {'fn': '*_bcbio_coverage.txt'},
                      'coverage_avg': {'fn': '*_bcbio_coverage_avg.txt'},
                      'variants': {'fn': '*_bcbio_variants.txt'},
                      'qsignature': {'fn': '*bcbio_qsignature.ma'},
                      'vcfstats': {'fn': '*_bcbio_variants_stats.txt'},
                      }

config.fn_clean_exts.append({'type': 'regex', 'pattern': '_bcbio.*'})
