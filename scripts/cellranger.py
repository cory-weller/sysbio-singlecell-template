#!/usr/bin/env python3
# coding: utf-8
import argparse
import sys
import os
print(os.environ)
sys.tracebacklimit = 0

try:
    from sysbio_sc import *
except ModuleNotFoundError:
    sys.path.append('src')
    from sysbio_sc import *


import argparse
parser = argparse.ArgumentParser(
    description="Run cellranger within python subprocess."
)

# Positional argument
#parser.add_argument("input", help="Input file path")

# Optional arguments
parser.add_argument("--library", default=None, help="library ID, i.e. key in cellranger-libraries.yaml")


args = parser.parse_args()

if args.library is None:
    raise ValueError('cellranger library argument missing!')


#===================================================================================================
# 01 Load config
#===================================================================================================

project_dir = get_project_dir()                                 # Returns PosixPath object
configfile = project_dir / 'config.yaml'                 # Default config name within project_dir
config = import_config(configfile)                              # Loads as DotDict; exits if cannot load
libraries = import_config('cellranger-libraries.yaml')
library = args.library
sample_prefixes = libraries[library]
if isinstance(sample_prefixes, str):
    sample_prefixes = [sample_prefixes]

# Parse create-bam argument
falsy = ['None', 'none', 'NA', 'N/A', 'na', 'n/a', False, 'False', 'false', 'F', 'f', 'N', 'N']
truthy = ['Yes', 'yes', 'Y', 'y', 'True', True, 'true', 'T']

# lower-case arg required by cellranger
if config.cellranger.create_bam in truthy:
    config.cellranger.create_bam = 'true'      
elif config.cellranger.create_bam in falsy:
    config.cellranger.create_bam = 'false'
else:
    raise ValueError('''Config for cellranger create_bam, should be "true" or "false"''')

#===================================================================================================
# 02 Pre-run Checks
#===================================================================================================
data_dir = require_path(project_dir / config.data_dir, label='data_dir', kind='dir', create=True)
#slurm_id = get_slurm_id()                                       # Exits if no SLURM_JOB_ID

temp_dir = require_path(os.environ['TMPDIR'], label='temporary workign directory', kind='dir', create=True)
output_dir = require_path(data_dir / 'CELLRANGER' / library, label='Cellranger output dir', kind='dir', create=True)
transcriptome_dir = require_path(config.cellranger.transcriptome, label='Transcriptome', kind='dir', create=False)  # Exits if does not exist
check_write_access(temp_dir)                                    # Exits if not writable
require_command('cellranger')                                   # Exits if command not in PATH


#===================================================================================================
#  03 Build --libraries csv for running cellranger
#===================================================================================================

os.chdir(temp_dir)

with open('libraries.csv', 'w') as outfile:
    outfile.write('fastqs,sample,library_type,\n')
    for i in sample_prefixes:
        outfile.write(f"{data_dir},{i},{config.cellranger.assay},\n")


#===================================================================================================
#  03 Run cellranger
#===================================================================================================

if config.cellranger.tool == 'count':
    cmd = ['cellranger','count',
        '--id', library,
        '--create-bam', config.cellranger.create_bam,
        '--libraries', 'libraries.csv',
        '--output-dir', library,
        '--transcriptome', transcriptome_dir,
        '--chemistry', config.cellranger.chemistry,
        '--disable-cell-annotation',
        '--nosecondary']
    print(cmd)
    try:
        cmd = [str(x) for x in cmd]
        print('Running command:')
        print(' '.join(cmd))
        subprocess.run(cmd)
        print('Cellranger complete, copying files...')
        # Copy outputs back to permanent dir after running
        shutil.copytree(src=f"{library}/outs",
                        dst=output_dir, 
                        dirs_exist_ok=True)
    except:
        raise
else:
    raise RuntimeError(f"cellranger tool set to {config.cellranger.tool} but code undefined")

