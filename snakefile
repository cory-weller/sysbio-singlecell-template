#!/usr/bin/env python

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

project_dir = Path(os.getcwd())                                 # Returns PosixPath object
# config_file = project_dir / 'config.yaml'                 # Default config name within project_dir
# configfile: 'settings/config.yaml'                              # Loads as DotDict; exits if cannot load
config = DotDict(config)
#print(config)

LIBRARIES = import_dotdict('cellranger-libraries.yaml')
data_dir = project_dir / config.data_dir
synapse_metadata_summary = data_dir / config.synapse.metadata_summary

# Convert paths to posix strings
data_dir = data_dir.as_posix()
synapse_metadata_summary = synapse_metadata_summary.as_posix()

localrules: get_metadata, build_library_mapping, somalier_relate,
lscratch_tmpdir = "/lscratch/$SLURM_JOB_ID"

libraries = list(LIBRARIES.keys())
# VR297, not enough cells
bad_libraries = ['VR297']
for x in bad_libraries:
    libraries.remove(x)

all_cellbender_output = expand(f"{data_dir}/CELLBENDER/{{libraryID}}/output.h5", libraryID=libraries)
all_pileup_output = expand(f"{data_dir}/FREEMUXLET/{{libraryID}}/ALL_AF.vcf.gz", libraryID=libraries)
all_intersect_output = expand(f"{data_dir}/FREEMUXLET/{{libraryID}}/AF.filtered.bam", libraryID=libraries)
all_dsc_pileup_output = expand(f"{data_dir}/FREEMUXLET/{{libraryID}}/allsites_pileup.plp.gz", libraryID=libraries)
all_freemuxlet_output = expand(f"{data_dir}/FREEMUXLET/{{libraryID}}/demux.clust1.samples.gz", libraryID=libraries)
all_singlet_output=expand(f"{data_dir}/FREEMUXLET/{{libraryID}}/{{N}}-singlets.bam", libraryID=libraries, N=[0,1,2])
somalier_relate_output = f"{data_dir}/FREEMUXLET/somalier.pairs.tsv"

filter_bc_list = expand(f"{data_dir}/CELLBENDER/{{libraryID}}/qc-barcodes.txt", libraryID=libraries)
filtered_bam_list = expand(f"{data_dir}/CELLRANGER/{{libraryID}}/qc_bam.bam", libraryID=libraries)
qc_somalier_relate_output = f"{data_dir}/FREEMUXLET/qc-somalier.pairs.tsv"


def get_overlapping_libraries(EXCLUDE_LIBS):
    with open(f"{data_dir}/overlapping-library-pairs.txt", 'r') as infile:
        good_pairs = []
        for line in infile:
            pair = line.strip().split()
            if set(pair).isdisjoint(EXCLUDE_LIBS):
                good_pairs.append(pair)
        outputs = [f"{data_dir}/FREEMUXLET/FINGERPRINTS/{x[0]}_{x[1]}.fingerprints.txt" for x in good_pairs]
    return(outputs)



# ── Rules ──────────────────────────────────────────────────────
rule all: 
    input:         f'{data_dir}/FREEMUXLET/FINGERPRINTS/deconcolved.tsv', filtered_bam_list, all_cellbender_output, all_pileup_output, filter_bc_list

#get_overlapping_libraries(EXCLUDE_LIBS=bad_libraries)

rule get_sample_ids:
    input: f"{data_dir}/syn51753257/AMP-AD_DiverseCohorts_assay_multiome_metadata.csv"
    output: f"{data_dir}/donor_ids.txt"
    shell:
        '''
        awk -F',' 'NR > 1 {{print $1}}' {input} | cut -d '_' -f 1 | sort -u  > | awk '{{print $1,$1}}' {data_dir}/donor_ids.txt
        '''

rule get_synapse_metadata:
    """
    Iterates over synapse IDs defined in config.yaml, recursively,
    retrieving file metadata and writing to separate metadata.tsv files.
    Outputs metadata in nested directories, reflecting file structure on synapse.
    Requires a synapse authentication token, its location defined in config.yaml
    """
    output: synapse_metadata_summary
    conda: 'envs/sysbio_singlecell'
    shell: 
        """
        python3 src/01_get_synapse_metadata.py
        """


rule download_data:
    """
    Iterates over the combined metadata folder. For every file, in the list,
    if the file does not exist: it is downloaded using synapse get.
    If the file exists,
        If file.md5 has not been generated, it is generated.
        file.md5 is compared to the metadata table.
        if the md5 values do not match,
            the file and file.md5 are deleted.
            The bad file is logged as an error.
    Only if no errors occur does .download-check get generated.
    """
    input: synapse_metadata_summary
    output: f"{data_dir}/download_data.done"
    conda: 'envs/sysbio_singlecell'
    shell:
        """
        python3 src/02_file_download.py
        touch {output}
        """

rule build_library_mapping:
    '''
    This rule is specific to this project / dataset. It must be modified depending on the format
    of input metadata in order to get a list of libraries.
    '''
    input: ancient(f"{data_dir}/{config.synapse.metadata_summary}")
    output: 'cellranger-libraries.yaml'
    conda: 'envs/sysbio_singlecell'
    run:
        def get_library_file(fn1, fn2):
            with open(fn1, 'r') as infile:
                for line in infile:
                    line = line.rstrip().split()
                    if line[1].endswith(fn2):
                        return(line[1])
            # After reaching end of file without a match
            raise FileNotFoundError(f"Could not find file matching {fn2} within column 2 of ${fn1}")
        # generate library ID file if it does not exist
        library_metadata_file = get_library_file(fn1=input[0],
                                                fn2=config.metadata.sample_libraries)
        libraries = []
        delim=','
        header_start = 'specimenID'
        with open(data_dir / library_metadata_file, 'r') as infile:
            for line in infile:
                if line.startswith(header_start):
                    continue
                libraries.append(line.split(delim)[1])
        libraries = sorted(list(set(libraries)))
        library_yaml = { x : x+'-GEX' for x in libraries}
        
        with open(output[0], 'w') as outfile:
            yaml.dump(library_yaml, outfile, default_flow_style=False, sort_keys=False)

rule run_cellranger:
    '''
    This rule is specific to this project / dataset. It must be modified depending on the format
    of input metadata in order to get a list of libraries.
    '''
    input: ancient('cellranger-libraries.yaml')
    output: 
        h5=f"{data_dir}/CELLRANGER/{{libraryID}}/raw_feature_bc_matrix.h5",
        bam=f"{data_dir}/CELLRANGER/{{libraryID}}/possorted_genome_bam.bam",
    conda: 'envs/sysbio_singlecell'
    params:
        PAR_SLURM_ID=os.environ['SLURM_JOB_ID'],
        OUTDIR=os.getcwd()
    envmodules: 'cellranger/10.1.0'
    shell:
        '''
        if [ -z ${{SLURM_JOB_ID:-}} ]; then export SLURM_JOB_ID={params.PAR_SLURM_ID}; fi; export TMPDIR=/lscratch/$SLURM_JOB_ID
        python3 src/cellranger.py --library {wildcards.libraryID} && touch {output.h5}
        '''

rule run_cellbender:
    '''
    This rule is specific to this project / dataset. It must be modified depending on the format
    of input metadata in order to get a list of libraries.
    '''
    input: h5=ancient(f"{data_dir}/CELLRANGER/{{libraryID}}/raw_feature_bc_matrix.h5")
    output: filtered_h5=f"{data_dir}/CELLBENDER/{{libraryID}}/output.h5"
    conda: 'envs/sysbio_singlecell'
    params:
        PAR_SLURM_ID=os.environ['SLURM_JOB_ID'],
        OUTDIR= f"{data_dir}/CELLBENDER/{{libraryID}}/"
    envmodules: 'cellbender/0.3.2'
    shell:
        '''
        #if [ -z ${{SLURM_JOB_ID:-}} ]; then export SLURM_JOB_ID={params.PAR_SLURM_ID}; fi; export TMPDIR=/lscratch/$SLURM_JOB_ID
        TMPDIR="LOCALTMP/{wildcards.libraryID}/"
        mkdir -p $TMPDIR && cd $TMPDIR
        TMP=$TMPDIR
        
        cellbender remove-background --cuda --input {input.h5} --output output.h5
        
        cp ./* {params.OUTDIR} && rm -rf $TMPDIR
        
        '''
