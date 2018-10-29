#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 15:52:09 2018

@author: huyn
"""

a= None
def retrieve(n):
    global a
    a = n
retrieve(2)
print (a)