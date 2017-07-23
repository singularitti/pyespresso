#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# created at Jul 23, 2017 4:16 PM by Nil-Zil


import os
import re
import shelve
import shutil
import subprocess
import time
from math import exp
from random import randint

import atommass as at
import flow
import func_cij
import functions
import numpy as np
from func_cij import Ekl, apply_strain
from scipy.optimize import curve_fit
