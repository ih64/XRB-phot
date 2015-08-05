#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distribute_setup import use_setuptools
    use_setuptools()
    from setuptools import setup

import sys

from pkg_resources import get_distribution, Requirement


# If stsci.distutils is being used to install another package in the stsci
# namespace package, we may need to first re-import the stsci package so that
# all the entries (including the current path) are added to stsci.__path__
# Deleting 'stsci' from sys.modules will force such a re-import.
if 'stsci' in sys.modules:
    del sys.modules['stsci']


# This is a workaround for http://bugs.python.org/setuptools/issue20; most
# packages that have this package as a setup-requirement also have d2to1 as
# a setup-requirement, which can lead to bugginess.
# See also http://mail.python.org/pipermail/distutils-sig/2011-May/017812.html
# for a description of the problem (in my example, package_A is d2to1 and
# package_B is stsci.distutils).
# This issue was fixed in distribute 0.6.17 and in setuptools 0.6c10, but
# leaving in support for older versions for now.
requirements = [Requirement.parse('setuptools<0.6c10'),
                Requirement.parse('distribute<0.6.19')]
# Distribution will actually convert a requirement for any setuptools version
# to a requirement for distribute, so if distribute is in use the first
# requirement is useless.  setuptools does something similar: yes, setuptools
# and distribute are actually antagonistic toward each other--ridiculous.
if requirements[0].key == requirements[1].key:
    del requirements[0]
try:
    # Note: If distribute is installed get_distribution('setuptools') returns
    # the installed distribute distribution
    has_issue205 = any([get_distribution('setuptools') in req
                        for req in requirements])
except:
    has_issue205 = False

if has_issue205:
    import sys
    from pkg_resources import working_set
    from setuptools import sandbox
    from setuptools.command import easy_install

    # Monkey patch setuptools so that subsequent calls to run_setup also
    # have this patch:
    _old_run_setup = sandbox.run_setup
    def run_setup(setup_script, args):
        save_entries = working_set.entries[:]
        save_entry_keys = working_set.entry_keys.copy()
        save_by_key = working_set.by_key.copy()
        save_modules = sys.modules.copy()
        try:
            _old_run_setup(setup_script, args)
        finally:
            working_set.entries = save_entries
            working_set.entry_keys = save_entry_keys
            working_set.by_key = save_by_key
            sys.modules.update(save_modules)
            for key in list(sys.modules):
                if key not in save_modules:
                    del sys.modules[key]
    sandbox.run_setup = run_setup
    easy_install.run_setup = run_setup

    # Patch the current call to run_setup
    save_entries = working_set.entries[:]
    save_entry_keys = working_set.entry_keys.copy()
    save_by_key = working_set.by_key.copy()
    save_modules = sys.modules.copy()
try:
    setup(
        setup_requires=['d2to1>=0.2.9'],
        namespace_packages=['stsci'], packages=['stsci'],
        d2to1=True,
    )
finally:
    if has_issue205:
        working_set.entries = save_entries
        working_set.entry_keys = save_entry_keys
        working_set.by_key = save_by_key
        sys.modules.update(save_modules)
        for key in list(sys.modules):
            if key not in save_modules:
                del sys.modules[key]
