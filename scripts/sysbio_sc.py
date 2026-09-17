#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import sys
import os
import yaml
import json
import subprocess
#import regex as re
import shutil

class DotDict(dict):
    '''
    define DotDict class for simpler config access
    instead of: mydict['key1']['key2']
    access by: mydict.key1.key2
    '''
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for key, value in self.items():
            self[key] = self._wrap(value)
    
    @classmethod
    def _wrap(cls, value):
        if isinstance(value, dict) and not isinstance(value, cls):
            return cls(value)
        if isinstance(value, (list, tuple)):
            return type(value)(cls._wrap(v) for v in value)
        return value
    
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(key)
    
    def __setattr__(self, key, value):
        self[key] = self._wrap(value)


def import_dotdict(yamlfile):
    '''
    Imports config yaml file and returns dot dict object. Exits if load is unsuccessful.
    '''
    try:
        with open(yamlfile, 'r') as f:
            contents = DotDict(yaml.safe_load(f))
        return contents
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Couldn't load file {yamlfile}, check paths?") from e

def get_project_dir():
    try:
        project_dir = Path(__file__).resolve().parent.parent.absolute()
        return project_dir
    except NameError as e:
        raise RuntimeError ("Couldn't automatically determine script or project location!") from e

def get_slurm_id():
    # Check running as slurm job, exits if not
    try:
        slurm_id = int(os.environ['SLURM_JOB_ID'])
        return slurm_id
    except KeyError as e:
        raise KeyError(f"No SLURM job ID detected within environment, check this is running via sbatch") from e

def get_slurm_array_id():
    # Check running as slurm job ARRAY
    try:
        slurm_array_id = int(os.environ['SLURM_ARRAY_TASK_ID'])
        return(slurm_array_id)
    except KeyError as e:
        raise KeyError(f"No SLURM job array task ID detected within environment, check this is running via sbatch --array") from e

def require_path(path, label, kind, create=False):
    """kind must be 'file' or 'dir'."""
    if kind not in ("file", "dir"):
        raise ValueError("kind must be 'file' or 'dir'")
    path = Path(path).resolve()
    if path.exists():
        # exists — verify it's the right type
        if kind == "file" and not path.is_file():
            raise FileNotFoundError(f"{label} {path} exists but is not a file")
        if kind == "dir" and not path.is_dir():
            raise NotADirectoryError(f"{label} {path} exists but is not a directory")
        return path.absolute()
    # doesn't exist
    if not create:
        raise FileNotFoundError(f"required {label} {kind} {path} does not exist")
    elif create:
        try:
            if kind == "file":
                path.parent.mkdir(parents=True, exist_ok=True)  # ensure parent dir
                path.touch()
            elif kind == 'dir':
                path.mkdir(parents=True, exist_ok=True)
            return path.absolute()
        except OSError as e:
            raise OSError(f"ERROR: could not create {label} {kind} {path}") from e


def require_command(name):
    path = shutil.which(name)
    if path is None:
        raise OSError(f"required bash command not found: {name}")

def check_write_access(path):
    '''
    Consider a directory (or if a file is provided, consider its parent directory).
    Create directory (and its parents) if necessary.
    Return the absolute path if the directory is writable, otherwise raise error.
    '''
    path = Path(path).resolve().absolute()
    if path.is_file():
        path = path.parent
    if not os.access(path, os.W_OK):
        try:
            os.makedirs(path, exist_ok = True)
            return(path)
        except:
            raise PermissionError(f"No write permissions for {path}")


def to_ranges(nums):
    '''
    Converts a set of integers to the minimum ranges required to represent the same integers.
    Range start and end points are both inclusive.
    e.g. to_ranges([1,2,3,4,5,6,8,9,10,11,12,13]) yields [(1, 6), (8, 13)]
    '''
    nums = sorted(set(nums))
    ranges = []
    start = prev = nums[0]
    for n in nums[1:]:
        if n == prev + 1:
            prev = n
        else:
            ranges.append((start, prev))
            start = prev = n
    ranges.append((start, prev))
    return ranges

def find_file(pattern):
    try:
        result = subprocess.run(['find', '.', '-type', 'f', '-name', pattern], capture_output=True, text=True).stdout.strip()
        if len(result.split()) > 1:
            raise RuntimeError("Multiple files found for pattern {pattern}, this pattern should be unique")
        return(Path(result).resolve().absolute())
    except:
        raise

