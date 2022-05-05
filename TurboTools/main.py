# -*- coding: utf-8 -*-
"""
Created on Sun May  1 15:33:30 2022

@author: Karel Van den Borre

Main file to load, print and save a map in a different format
"""

import matplotlib.pyplot as plt
import TurboMap as tm
import pandas as pd
import os

# close any plots which are still open
plt.close('all')

# import the experimental data. Pandas allows to access all fields by name
# this test genertes the operating line of the machines, after an initial transient
testResults = pd.read_csv("experiments/Full007.rpt", delimiter="\t", skiprows=4, encoding = "ISO-8859-1")
# drop first rows to remove initial conditions experiment
testResults=testResults.iloc[2:, :].reset_index(drop=True)

# create folder to save figures
if not os.path.exists("figures"):
    os.mkdir("figures") 

# load the high pressure compressor, scale and plot it
HPC = tm.TurboMap("maps/HPC01_SI.Map", 20, 10, True, name="HPC")
HPC.scaleMap(1.0000321000000001, 0.78706299099999999, 0.98847553799999999)
HPC.writeMap("HPC01_SI.xml")
HPC.printMap()
HPC.plotDesingPoint(0.5, 1)
HPC.plotExperiment(testResults, "HPC", "rainbow")
plt.savefig("figures/HPC", dpi=300)

# load the high pressure compressor, scale and plot it
HPT = tm.TurboMap("maps/HPT01_SI.Map", 20, 5, False, name="HPT")
HPT.scaleMap(0.90707883899999997, 0.83650503799999998, 0.97804945799999998)
HPT.writeMap("HPT01_SI.xml")
HPT.printMap()
HPT.plotDesingPoint(0.5, 1)
HPT.plotExperiment(testResults, "HPT", "rainbow")
plt.savefig("figures/HPT", dpi=300)

# show plots, also when running from cmd
plt.show()
