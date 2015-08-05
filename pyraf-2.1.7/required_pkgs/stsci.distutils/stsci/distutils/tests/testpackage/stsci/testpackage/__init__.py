try:
    from .version import (__version__, __svn_revision__, __svn_full_info__,
                          __setup_datetime__)
except ImportError:
    __version__ = ''
    __svn_revision__ = ''
    __svn_full_info__ = ''
    __setup_datetime__ = None
