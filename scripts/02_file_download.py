#!/usr/bin/env python3
# coding: utf-8

import sys
import time
sys.tracebacklimit = 0

try:
    from sysbio_sc import *
except ModuleNotFoundError:
    try:
        sys.path.append('src')
        from sysbio_sc import *
    except:
        raise


#===================================================================================================
# 01 Load config
#===================================================================================================

project_dir = get_project_dir()                                 # Returns PosixPath object
configfile = project_dir / 'config.yaml'                 # Default config name within project_dir
config = import_config(configfile)                              # Loads as DotDict; exits if cannot load
data_dir = project_dir / config.data_dir


import synapseclient
import asyncio
import hashlib

from synapseclient.models import (
    File, Folder, Project, Table, EntityView, Dataset,
    DatasetCollection, MaterializedView, SubmissionView, VirtualTable
)


def get_md5(path:PosixPath):
    path = str(path.resolve().absolute())
    hash_md5 = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


async def main():
    # Authenticate to synapse
    if config.synapse.token.startswith('~'):
        token_file = require_path(Path(config.synapse.token).expanduser(), label='synapse auth token', kind='file', create=False)
    else:
        token_file = require_path(config.synapse.token, label='synapse auth token', kind='file', create=False)
    with open(token_file, 'r') as f:
        token = f.read().strip()

    syn = synapseclient.login(authToken=token)


    # Import reads metadata file
    metadata_filename = data_dir / config.synapse.metadata_summary
    df = pd.read_csv(metadata_filename, sep='\t')

    if config.synapse.clobber:
        print(f"WARNING: synapse.clobber is set to True. Files will be removed and overwritten! Waiting 10 seconds before continuing")
        time.sleep(10)
                  
    bad_files = []
    good_files = []
    #retries = []        # Tuple (synid, filepath, version, meta_md5)

    # FIRST check of files
    # iterate over rows of metadata:
    for i in range(df.shape[0]):
        x = df.iloc[i]
        filepath = data_dir / x['file_name']
        md5_filepath = Path(str(filepath) + '.md5')
        synid = x['synid']
        version = int(x['version_label'])
        meta_md5 = x['content_md5']
        
        # Delete files clobber=True and files exist
        if config.synapse.clobber:
            for file in [filepath, md5_filepath]:
                if file.exists():
                    file.unlink()
        
        # Download file if it does not exist
        if not filepath.exists():
            await syn.get_async(synid, download_file=True,
                            downloadLocation=filepath.parent,
                            if_collision="overwrite.local",
                            synapse_client=syn)
        
        # Calculate md5 if it does not exist
        if not md5_filepath.exists():
            file_md5 = get_md5(filepath)
        else:
            with open(md5_filepath, 'r') as infile:
                file_md5 = infile.read().strip()
        
        # Compare md5
        if file_md5 != meta_md5:
            print(f"{synid} Local file md5 ({file_md5}) does not match metadata md5 ({meta_md5}) for file {filepath} version {version}!")
            print(f"Removing file {filepath}")
            os.remove(filepath)
            os.remove(f"{filepath}.md5")
            bad_files.append((synid, filepath, version, meta_md5))
        
        # Only written if file md5 matches metadata md5
        elif file_md5 == meta_md5:
            good_files.append((synid, filepath, version, meta_md5))
            with open(f"{filepath}.md5", 'w') as outfile:
                outfile.write(file_md5 + '\n')
    print(f"Total files in {metadata_filename}: {df.shape[0]}")
    print(f"Good files: {len(good_files)}")
    print(f"Bad files: {len(bad_files)}")

    
    # Test if every file in metadata passed md5 check
    if df.shape[0] != len(good_files):
        raise RuntimeError(f"Not all files pass md5 check. {len(bad_files)} bad files:\n{'\n'.join(bad_files)}")


if __name__ == "__main__":
    asyncio.run(main())
    print("Completed successfully")
    sys.exit(0)