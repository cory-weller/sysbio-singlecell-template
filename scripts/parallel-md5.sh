#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition quick,norm
#SBATCH --time 1:00:00



get_md5sum() {
    fn=${1}
    if [ ! -f "${fn}" ]; then
        echo "${fn} does not exist"
        return 0
    fi
    if [ ! -f "${fn}.md5" ]; then
        md5sum ${fn} | awk '{print $1}' > ${fn}.md5
    fi
}

export -f get_md5sum

parallel -j 32 --colsep '\t' get_md5sum {2} :::: ${1}