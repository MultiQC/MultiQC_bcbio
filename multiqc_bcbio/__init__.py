from __future__ import absolute_import

from .bcbio import MultiqcModule

from multiqc import config

config.sp['bcbio'] = {'metrics': {'fn': '*_bcbio.txt'},
                      'coverage': {'fn': '*_bcbio_coverage.txt'},
                      'coverage_avg': {'fn': '*_bcbio_coverage_avg.txt'},
                      'variants': {'fn': '*_bcbio_variants.txt'},
                      }
