"""dustline: measured extinction law (R_V) and 3D dust run A_X(D) along any sightline.

Public API:

    import dustline
    res = dustline.Sightline(ra=..., dec=...).run()
    res.rv                    # measured R_V (median, MAD, N)
    res.extinction("I")       # A_I(D) table (D_kpc, A_med, A_16, A_84, n_stars)
    res.law                   # band ratios A_X/A_I under the measured law
"""

from .api import ExtinctionResult, Sightline
from .userphot import UserPhotometry

__version__ = "0.7.0"
__all__ = ["Sightline", "ExtinctionResult", "UserPhotometry", "__version__"]
