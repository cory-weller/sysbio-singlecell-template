#!/usr/bin/env python3
# coding: utf-8

import sys
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

from synapseclient.models import (
    File, Folder, Project, Table, EntityView, Dataset,
    DatasetCollection, MaterializedView, SubmissionView, VirtualTable
)


def write_mdata(mdata, fn):
    header = [x[0] for x in mdata[0]]
    values = []
    for x in mdata:
        values.append([i[1] for i in x])
    os.makedirs(Path(fn).parent, exist_ok=True)
    with open(fn, 'w') as outfile:
        outfile.write('\t'.join(header) + '\n')     # Write first half of tuple as header
        for v in values:
            outfile.write('\t'.join(v) + '\n')     # Write second half of tuple as values

def format_property(p):
    o = {}
    for key,value in p.items():
        if not isinstance(value, list):
            o[key] = str(value)
        if isinstance(value, list):
            o[key] = ';'.join([str(x) for x in value])
    return(o)

def update_files(synfolder, filetree='.'):
    # Update each file in the current node
    for synfile in synfolder.files:
        synfile.file_name = f"{filetree}/{synfile.file_handle.file_name}"
    # Recurse into nested folders
    for folder in synfolder.folders:
        update_files(folder, f"{filetree}/{folder.id}")

# def get_files(synfolder):
#     for synfile in synfolder.files:
#         yield synfile
#     for folder in synfolder.folders:
#         yield from get_files(folder)

def get_items(_dict, _keys):
    items = []
    for key in _keys:
        if key in _dict:
            items.append((key, _dict[key]))
        else:
            items.append((key, 'N/A'))
    return(items)

def extract_mdata(SYNFILE):
    file_handle_order = ['id',
                    'etag',
                    'created_by',
                    'created_on',
                    'modified_on',
                    'concrete_type',
                    'content_type',
                    'content_md5',
                    'storage_location_id',
                    'content_size',
                    'status',
                    'bucket_name',
                    'key',
                    'preview_id',
                    'is_preview',
                    'external_url']
    annotations_order = ['sex',
                        'assay',
                        'grant',
                        'organ',
                        'study',
                        'tissue',
                        'runType',
                        'species',
                        'cellType',
                        'dataType',
                        'platform',
                        'consortium',
                        'fileFormat',
                        'readLength',
                        'specimenID',
                        'dataSubtype',
                        'libraryPrep',
                        'individualID',
                        'resourceType',
                        'isModelSystem',
                        'isMultiSpecimen',
                        'nucleicAcidSource',
                        'individualIdSource']
    metadata = []
    synid = SYNFILE.id
    metadata =  [('synid',synid)]
    metadata += [('file_name', SYNFILE.file_name)]
    metadata += [('version_label', SYNFILE.version_label)]
    metadata += get_items(format_property(SYNFILE.file_handle.__dict__), file_handle_order)
    metadata += get_items(format_property(SYNFILE.annotations), annotations_order)
    return(metadata)


async def get_files(synid, syn):
    empty_entity = await syn.get_async(synid, download_file=False)
    if empty_entity.concreteType == 'org.sagebionetworks.repo.model.FileEntity':  # is File
        entity = await File(id=synid, download_file=False).get_async()    # Single file object with metadata
        entity.file_name = entity.id + '/' + entity.file_handle.file_name
        return([entity])
    elif empty_entity.concreteType == 'org.sagebionetworks.repo.model.Folder':     # is Folder
        entity = Folder(id=synid)
        entity = await entity.sync_from_synapse_async(download_file=False)       # Folder with list of file objects
        update_files(entity, entity.id)
        return(flatten_files(entity))

def flatten_files(node):
    files = []
    if node.files is not None:
        files += (node.files)
    for folder in node.folders:
        files += flatten_files(folder)
    return files


async def main():

    # Check if final output exists:
    metadata_output = data_dir / config.synapse.metadata_summary
    if metadata_output.exists():
        if config.synapse.clobber:
            print(f"Final output {metadata_output} exists, but re-running anyway because synapse.clobber=True")
        else:
            print(f"Final output {metadata_output} exists, stopping because synapse.clobber=False")
            return


    # Authenticate to synapse
    if config.synapse.token.startswith('~'):
        token_file = require_path(Path(config.synapse.token).expanduser(), label='synapse auth token', kind='file', create=False)
    else:
        token_file = require_path(config.synapse.token, label='synapse auth token', kind='file', create=False)
    with open(token_file, 'r') as f:
        token = f.read().strip()

    syn = synapseclient.login(authToken=token)


    # Recursively syncfolder data

    config_datasets = config.synapse.datasets
    # Ensure it is a list, even if length of 1
    if isinstance(config_datasets, str):
        datasets_ids = [config_datasets]
    else:
        datasets_ids = config_datasets

    print(f"Retrieving metadata for datasets: {' '.join(datasets_ids)}")

    combined_metadata = []
    for synapse_id in datasets_ids:
        mdata = await get_files(synapse_id, syn)
        combined_metadata += [extract_mdata(x) for x in mdata]


    # with open("combined_metadata.pkl", "wb") as file:
    #     pickle.dump(combined_metadata, file)

    # Write to final metadata file for sequencing reads to output
    write_mdata(combined_metadata, metadata_output)




if __name__ == "__main__":
    asyncio.run(main())
    print("Completed successfully")
    sys.exit(0)