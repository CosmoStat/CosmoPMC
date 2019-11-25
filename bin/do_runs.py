#!/usr/bin/env python

"""do_runs.py

:Author: Martin Kilbinger

*Date*: 2018

:Package: CosmoPMC

"""


# Compability with python2.x for x>6
from __future__ import print_function


import sys
import os
import copy
import re

import numpy as np
import pylab as plt

from astropy.io import ascii
from astropy.table import Table, Column

from optparse import OptionParser
from optparse import OptionGroup

import mkstuff as stuff



def params_default():
    """Set default parameter values.

    Parameters
    ----------
    None

    Returns
    -------
    p_def: class stuff.param
        parameter values
    """

    p_def = stuff.param(
        nrun = 10,
        mode = 'rt'
    )

    return p_def



def parse_options(p_def):
    """Parse command line options.

    Parameters
    ----------
    p_def: class stuff.param
        parameter values

    Returns
    -------
    options: tuple
        Command line options
    args: string
        Command line string
    """

    usage  = "%prog [OPTIONS]"
    parser = OptionParser(usage=usage)

    parser.add_option('-N', '--Nrun', dest='nrun', default=p_def.nrun,
        help='number of runs, default = {}'.format(p_def.nrun))
    parser.add_option('-m', '--mode', dest='mode', default=p_def.mode,
        help='run mode, combination of chars \'r\' (run) \'t\' (tar) = {}'.format(p_def.mode))

    parser.add_option('-n', '--dry_run', dest='dry_run', action='store_true', help='do not run, only print commands')
    parser.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbose output')

    options, args = parser.parse_args()

    return options, args



def check_options(options):
    """Check command line options.

    Parameters
    ----------
    options: tuple
        Command line options

    Returns
    -------
    erg: bool
        Result of option check. False if invalid option value.
    """

    see_help = 'See option \'-h\' for help.'

    return True



def update_param(p_def, options):
    """Return default parameter, updated and complemented according to options.
    
    Parameters
    ----------
    p_def:  class stuff.param
        parameter values
    optiosn: tuple
        command line options
    
    Returns
    -------
    param: class stuff.param
        updated paramter values
    """

    param = copy.copy(p_def)

    # Update keys in param according to options values
    for key in vars(param):
        if key in vars(options):
            setattr(param, key, getattr(options, key))

    # Add remaining keys from options to param
    for key in vars(options):
        if not key in vars(param):
            setattr(param, key, getattr(options, key))

    # Do extra stuff 
    param.nrun = int(param.nrun)

    return param



def main(argv=None):
    """Main program.
    """

    # Set default parameters
    p_def = params_default()

    # Command line options
    options, args = parse_options(p_def)
    # Without option parsing, this would be: args = argv[1:]

    if check_options(options) is False:
        return 1

    param = update_param(p_def, options)

    niter = 8


    # Save calling command
    stuff.log_command(argv)
    if param.verbose:
        stuff.log_command(argv, name='sys.stderr')


    if param.verbose is True:
        print('Start program {}'.format(os.path.basename(argv[0])))


    ### Start main program ###

    # Run PMC
    if re.search('r', param.mode) is not None:
        for i in range(param.nrun):
            dir_name = 'run_{:02d}'.format(i)

            if param.verbose:
                print(dir_name)

            stuff.mkdir_p(dir_name)

            os.chdir(dir_name)

            stuff.ln_s('../config_pmc', 'config_pmc', verbose=param.verbose, force=False)

            stuff.ln_s('../prior_mvdens', 'prior_mvdens', verbose=param.verbose, force=False)

            test_file = 'evidence'
            if os.path.isfile(test_file) == False or os.path.getsize(test_file) == 0:
                cmd = 'cosmo_pmc.pl -p R -n 1'
                stuff.run_cmd(cmd, verbose=param.verbose, run=(not param.dry_run))
            else:
                if param.verbose:
                    print('File \'evidence\' exists, skipping...')

            os.chdir('..')


    # Create tar file with results
    if re.search('t', param.mode) is not None:
        to_tar = ''
        for i in range(param.nrun):
            to_tar = '{} run_{:02d}/iter_{}/all_cont2d.pdf'.format(to_tar, i, niter-1)

        tar_ball = 'contours.tgz'
        cmd = 'tar czf {} {}'.format(tar_ball, to_tar)
        stuff.run_cmd(cmd, verbose=param.verbose, run=(not param.dry_run))
        stuff.run_cmd('ls-for-scp.sh {}'.format(tar_ball), verbose=param.verbose)


    ### End main program

    if param.verbose is True:
        print('Finish program {}'.format(os.path.basename(argv[0])))

    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))

