"""module iraf.py -- home for all the IRAF tasks and basic access functions

$Id: iraf.py 1463 2011-06-24 22:58:30Z stsci_embray $

R. White, 1999 Jan 25
"""
from __future__ import division # confidence high

from iraffunctions import *

# a few CL tasks have modified names (because they start with '_')
import iraffunctions

_curpack = iraffunctions.curpack
_allocate = iraffunctions.clAllocate
_deallocate = iraffunctions.clDeallocate
_devstatus = iraffunctions.clDevstatus

del iraffunctions
