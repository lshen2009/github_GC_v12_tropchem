#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 12:28:15 2019

@author: lulushen
"""

import filecmp
import os
import glob
os.chdir('/Users/lulushen/Documents/GC_speedup/GEOS-Chem/github_GC_v12_tropchem/Code.GC12/GeosCore/')
files=glob.glob('*.F90')+glob.glob('*.F')

os.chdir('/Users/lulushen/Documents/GC_speedup/GEOS-Chem/github_GC_v12_tropchem')
for file in files:    
    if not filecmp.cmp('Code.GC12/GeosCore/'+file,'Code.GC12_Jan12/GeosCore/'+file):
        print(file)