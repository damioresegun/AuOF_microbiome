#!/usr/bin/env python3
###########################################################################################################################################################
#
#                               Script to carry out workflow control of microbiome analysis for the Antibiotics Under our Feet project
# 
# Requires: Via conda(fastqc, trim-galore, bmtagger, kraken2, bracken, krona, krakentools)
# Considerations: Remove need for kraken and bracken path arguments. Force installation via conda
# Author: Damilola Oresegun	                                                                                                                                 		              #
###########################################################################################################################################################
import enum
import os
import shutil
import sys
import argparse
import subprocess
from pathlib import Path
##########################################################################
# set the needed arguments
def get_args():
    parser