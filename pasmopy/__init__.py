"""Patient-Specific Modeling in Python"""

__author__ = __maintainer__ = "Hiroaki Imoto"
__email__ = "himoto@protein.osaka-u.ac.jp"
'''
try:
    from setuptools_scm import get_version

    __version__ = get_version(root="..", relative_to=__file__)
    del get_version
except (LookupError, ImportError):
    import sys
    
    if sys.version_info[:2] < (3, 8):
        from importlib_metadata import version
    else:
        from importlib.metadata import version
    __version__ = version(__name__)
    del version
'''
__version__ = "0.0.1"

from .construction import Text2Model
from .individualization import Individualization
from .patient_model import PatientModelAnalyses, PatientModelSimulations
