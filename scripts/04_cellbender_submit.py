#!/usr/bin/env python3
# coding: utf-8

import sys
sys.tracebacklimit = 0

try:
    from sysbio_sc import *
except ModuleNotFoundError:
    sys.path.append('src')
    from sysbio_sc import *


#===================================================================================================
# 01 Load config
#===================================================================================================

project_dir = get_project_dir()                                 # Returns PosixPath object
configfile = project_dir / 'config.yaml'                 # Default config name within project_dir
config = import_config(configfile)                              # Loads as DotDict; exits if cannot load
data_dir = project_dir / config.data_dir

# make generalizable format where config.step refers to this step in processing
step = 'cellbender'
config.step = config[step]

#===================================================================================================
# Pre-run checks
#===================================================================================================

# Fix importing 'None' or 'NA' as literal string instead of None object
falsy = ['None', 'none', 'NA', 'N/A', 'na', 'n/a', False, 'False', 'false']
if config.step.rerun_debug in falsy:
    config.step.rerun_debug = None

# for legibility when creating sbatch command
slurm = config.step.slurm

sbatch_file = require_path(path=project_dir / 'src' / 'sbatch-wrapper.sh', label='sbatch wrapper script', kind='file', create=False)



#===================================================================================================
#  RUN
#===================================================================================================

# Load library ID file
library_id_file = require_path(path = data_dir / 'libraryIDs.txt', label='LibraryID:N mapping file', kind='file', create=False)
df = pd.read_csv(library_id_file, sep='\t')


# Collect list of samples that haven't been done
ids_to_run = []
ids_finished = []


if config.step.rerun_debug is None:
    for libraryID in df['libraryID']:
        output = Path(data_dir / 'CELLBENDER' / libraryID /'output_filtered.h5')
        if output.exists():
            ids_finished.append(libraryID)
        elif not output.exists() or config.step.clobber is True:
            # Ensure input exists
            input = require_path(path=data_dir / 'CELLRANGER' / libraryID /'raw_feature_bc_matrix.h5', label='raw cellranger h5', kind='file', create=False)
            sample_N = int(df.loc[df['libraryID'] == libraryID, 'N'].iloc[0])
            ids_to_run.append(sample_N)
    print(f"INFO: {len(ids_finished)} samples already finished")
else:
    try:
        sample_N = int(df.loc[df['libraryID'] == config.step.rerun_debug, 'N'].iloc[0])
        ids_to_run.append(sample_N)
    except IndexError as e:
        raise RuntimeError(f"Tried to run explicitly with a single sample, {config.step.rerun_debug}, which does not exist in {config.data_dir}/libraryIDs.txt") from e


print(f"INFO: {len(ids_to_run)} sample(s) to run")

if len(ids_to_run) == 0:
    raise RuntimeError("All outputs exist, nothing to do unless cellbender clobber is set to True!")

# Convert list of sample numbers to ranges for job array
ranges_to_run = to_ranges(ids_to_run)
array_string = ','.join([f"{rng[0]}-{rng[1]}" for rng in ranges_to_run])
gres_string = ','.join([f"{x}:{slurm.gres[x]}" for x in slurm.gres])
slurm.modules = ' '.join(slurm.modules)         # string separate by spaces, in case of multiple args

python_script = require_path(project_dir / 'src' / 'cellbender-run.py', label = 'cellbender python script', kind='file', create=False)


custom_env = os.environ.copy()
custom_env["CONDA_SRC"] = config.conda_source
custom_env["CONDA_ENV"] = slurm.conda_env
custom_env["SCRIPT_PATH"] = python_script


sbatch_cmd = ['sbatch', 
                '--export=ALL',
                f"--array={array_string}%{slurm.array_limit}",
                f"--partition={slurm.partition}",
                f"--cpus-per-task={slurm.cpus}",
                f"--mem={slurm.mem}",
                f"--time={slurm.walltime}",
                f"--ntasks={slurm.ntasks}",
                f"--gres={gres_string}",
                sbatch_file,
                slurm.modules
            ]

sbatch_cmd = [str(x) for x in sbatch_cmd]

print("Running command:")
print(' '.join(sbatch_cmd))

subprocess.run(sbatch_cmd, env=custom_env)