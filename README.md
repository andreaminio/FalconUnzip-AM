# FalconUnzip-DClab
FALCON Unzip based pipeline integrating DAmasker and boosting Unzip speed

To perform the pipeline, along with a complete FalconUnzip installation, that we assume is placed in `/path/to/falcon-verXX` , a proper smrtlink distribution (we assume installed ) is also necessary according to the sequencing chemistry used.

> **Note:** The tool best performs on a cluster, where single processes of the different steps may be run in parallel. Nevertheless, a single machine can run multiple process using [parallel](https://www.gnu.org/software/parallel/parallel_tutorial.html)<sup>[1](#parallel)</sup>. Such processe that di benefit from this will be marked and approximate RAM and core usage will be given for a single process. 

## 0 - Environment Setup
### 0.1 - Environmental variables  
Remeber to set the paths to Falcon libararies and environments
```
export PYTHONUSERBASE=/path/to/falcon-verXX/
export LD_LIBRARY_PATH=/path/to/falcon-verXX/lib:${LD_LIBRARY_PATH}
export PATH=/path/to/smrtlink/install/smrtlink-release_X.XX/bundles/smrttools/install/smrttools-release_X.XX/smrtcmds/bin/:/path/to/falcon-verXX/bin:${PATH}
```

### 0.2 - BAX.H5 file conversion to BAM
This step is necessary for Arrow polishing as Blasr mapping files with other input file formats will not be accepted
```
for name in $(find -name "*bax.h5" | sed 's:\..\.bax\.h5$::' | sort -u | less -S); do echo $name ; bax2bam $name.1.bax.h5 $name.2.bax.h5 $name.3.bax.h5; done
```

## 1 - Raw reads repeat marking
### 1.1 - Scripts and Folders Preparation
#### Damasker DB creation
```
mkdir 0-rawreads
cd 0-rawreads
fasta2DB -v raw_reads $reads
DBsplit -x500 -s100 raw_reads
LB=$(cat raw_reads.db | LD_LIBRARY_PATH= awk '$1 == "blocks" {print $3}')
echo $LB
```
> **Note:** As Daligner merging accepts 250 chunks or less, tune `-s` parameter in order to have a $LB <= 500 as this pipeline is thought for a maximum of 2x250 batches.


#### Create scripts
```bash
CUTOFF=3000
echo -n $CUTOFF > length_cutoff
 
mkdir scripts
cd scripts
ln -s ../raw_reads.db
HPC.TANmask -T8 -frun_jobs.01.TANmask raw_reads 1-$LB 
HPC.REPmask -T8 -g2 -c15 -mtan -frun_jobs.02.REPmask raw_reads 1-$LB 
HPC.daligner -d -v -B128 -D512 -t30 -e0.7 -M40 -l1000 -k16 -h64 -w7 -s1000 -mtan -mrep2 -H$CUTOFF -frun_jobs.03.Daligner -T8 raw_reads 1-$LB 
for block_id in $(seq 1 $LB); do echo "LA4Falcon -H$CUTOFF -fso raw_reads.db raw_reads.${block_id}.las | fc_consensus --n_core 8 --output_multi --min_idt .70 --min_cov 4 --max_n_read 400 > cns_${block_id}.fasta "; done > run_jobs.04.Call_consensus.01.call
cd ..
mkdir logs
```

### 1.2 Run TANmasker

1. `bash scripts/run_jobs.01.TANmask.01.OVL` | **Parallelizable:** task_mem="10G", task_cores="8"
2. 






## References

<a name="parallel">1</a>: Tange (2011): GNU Parallel - The Command-Line Power Tool, ;login: The USENIX Magazine, February 2011:42-47.
