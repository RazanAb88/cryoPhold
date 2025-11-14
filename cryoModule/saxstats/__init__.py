# Lightweight shim: import everything from saxstats.saxstats if present,
# otherwise provide safe stubs and clear error messages for missing functionality.
try:
    from .saxstats import *
except Exception as e:
    # provide informative attributes so imports like "from saxstats import *" don't crash at import time
    def _not_implemented(*args, **kwargs):
        raise ImportError("saxstats functionality is not installed. Replace this shim with the real saxstats package.")
    # example stub names; you can add more stubs as validation.py requires
    parse_saxs = _not_implemented
    compute_intensity = _not_implemented
    __all__ = ['parse_saxs', 'compute_intensity']
