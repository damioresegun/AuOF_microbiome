#!/usr/bin/env python3
'''Script to hold various small tools and recurring commands
Author: Damilola Oresegun
'''

import os


def makeDirectory(directory):
    '''Function checks if a directory exists. If so, nothing is done,
    if the folder does not exist, the folder and any parent folder,
    is also made
    Input: Path to directory to be checked/made
    Output: None
    Usage: makeDirectory(path/to/directory)
    '''
    if os.path.exists(directory):
        pass
    else:
        os.makedirs(directory)


def prechecks(input, output):
    if os.path.exists(input):
        inp = "Good"
        pass
    else:
        inp = "Failed"
    if os.path.exists(output):
        out = "Good"
    else:
        out = "Make"
    return inp, out
        
