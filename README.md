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
> **Note:** As DAligner merging accepts 250 chunks or less, tune `-s` parameter in order to have a $LB <= 500 as this pipeline is thought for a maximum of 2x250 batches.


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
```

### 1.2 Run TANmasker

1. `bash scripts/run_jobs.01.TANmask.01.OVL`
   * **Parallelizable:** task_mem="10G", task_cores="8"
2. `bash scripts/run_jobs.01.TANmask.02.SORT`
   * **Parallelizable:** task_mem="10G", task_cores="1"
1. [Optional Checkpoint] `bash scripts/run_jobs.01.TANmask.03.CHECK.OPT`
4. `bash scripts/run_jobs.01.TANmask.04.RM`
5. `bash scripts/run_jobs.01.TANmask.05.MASK`
   * **Parallelizable:** task_mem="10G", task_cores="8"
6. `Catrack -v raw_reads.db tan`
7. `bash scripts/run_jobs.01.TANmask.06.RM`



### 1.3 REPmask
1. `bash scripts/run_jobs.02.REPmask.01.OVL`
   * **Parallelizable:** task_mem="35G", task_cores="8"
1. `bash scripts/run_jobs.02.REPmask.02.SORT`
   * **Parallelizable:** task_mem="10G", task_cores="1"
1. [Optional Checkpoint] `bash scripts/bash scripts/run_jobs.02.REPmask.03.CHECK.OPT
1. `bash scripts/run_jobs.02.REPmask.04.RM`
1. `bash scripts/run_jobs.02.REPmask.05.MERGE`
   * **Parallelizable:** task_mem="15G", task_cores="1"
1. [Optional Checkpoint] `bash scripts/run_jobs.02.REPmask.06.CHECK.OPT`
1. `bash scripts/run_jobs.02.REPmask.07.RM`
1. `bash scripts/run_jobs.02.REPmask.08.MASK`
   * **Parallelizable:** task_mem="15G", task_cores="8"
1. `Catrack -v raw_reads.db rep2
1. `bash scripts/run_jobs.02.REPmask.09.RM`



### 1.4 - Daligner
1. `bash scripts/run_jobs.03.Daligner.00.MKDIR`
1. `bash scripts/run_jobs.03.Daligner.01.OVL
   * **Parallelizable:** task_mem="35G", task_cores="8"
1. `bash scripts/run_jobs.03.Daligner.02.SORT
   * **Parallelizable:** task_mem="10G", task_cores="1"
1. [Optional Checkpoint] `bash scripts/run_jobs.03.Daligner.03.CHECK.OPT
1. `bash scripts/run_jobs.03.Daligner.04.RM
1. (OPTIONAL) Split merging processes if the 250< # Chunks ($LB) <= 500
   * `sed -i 's:\(.*raw_reads\.\)\([0-9]*\) \(.*\)250 \(.*\):\1\2\.1 \3250; \1\2\.2 \4; \1\2 raw_reads\.\2\.1 raw_reads\.\2\.2:;s: L[0-9]*.[0-9]*.;:;:' scripts/run_jobs.03.Daligner.05.MERGE `
1. `bash scripts/run_jobs.03.Daligner.05.MERGE"
   * **Parallelizable:** task_mem="15G", task_cores="1"
1. [Optional Checkpoint] `bash scripts/run_jobs.03.Daligner.06.CHECK.OPT
1. `bash scripts/run_jobs.03.Daligner.07.RM.OPT



### 1.5 - Correct
1. `bash scripts/run_jobs.04.Call_consensus.01.call`
   * **Parallelizable:** task_mem="35G", task_cores="8"


## 2.
### 2.1 - Create new database for corrected reads
```bash
ls -1 $(pwd)/cns_*.fasta > input_preads.fofn
#mkdir safe_copy
#for file in $(cat input_preads.fofn); do echo $file; cp $file safe_copy/$file; done
touch rdb_build_done
mkdir ../1-preads_ovl
while read fn; do fasta2DB -v ../1-preads_ovl/preads $fn; done < input_preads.fofn
cd ../1-preads_ovl
DBsplit -x500 -s100 preads
LB=$(cat preads.db | LD_LIBRARY_PATH= awk '$1 == "blocks" {print $3}')
mkdir scripts
cd scripts/
ln -s ../preads.db
```


### 2.2 - Create Damasker scripts for corrected reads
```
HPC.TANmask -T8 -frun_jobs.06.TANmask_corr preads 1-$LB
HPC.REPmask -T8 -g2 -c15 -mtan -frun_jobs.07.REPmask_corr preads 1-$LB
cd ..
```



### 2.3 - TANmask
1. `bash scripts/run_jobs.06.TANmask_corr.01.OVL
   * **Parallelizable:** task_mem="16G", task_cores="8"
1. `bash scripts/run_jobs.06.TANmask_corr.02.SORT`
   * **Parallelizable:** task_mem="10G", task_cores="1"
1. [Optional Checkpoint] `bash scripts/run_jobs.06.TANmask_corr.03.CHECK.OPT`
1. `bash scripts/run_jobs.06.TANmask_corr.04.RM`
1. `bash scripts/run_jobs.06.TANmask_corr.05.MASK`
   * **Parallelizable:** task_mem="10G", task_cores="8"
1. `Catrack -v preads.db tan
1. `bash scripts/run_jobs.06.TANmask_corr.06.RM`




### 2.4 - REPmask
1. `bash scripts/run_jobs.07.REPmask_corr.01.OVL
   * **Parallelizable:** task_mem="35G", task_cores="8"
1. `bash scripts/run_jobs.07.REPmask_corr.02.SORT`
   * **Parallelizable:** task_mem="10G", task_cores="1"
1. [Optional Checkpoint] `bash scripts/run_jobs.07.REPmask_corr.03.CHECK.OPT`
1. `bash scripts/run_jobs.07.REPmask_corr.04.RM`
1. `bash scripts/run_jobs.07.REPmask_corr.05.MERGE`
   * **Parallelizable:** task_mem="15G", task_cores="1"
1. [Optional Checkpoint] `bash scripts/run_jobs.07.REPmask_corr.06.CHECK.OPT`
1. `bash scripts/run_jobs.07.REPmask_corr.07.RM`
1. `bash scripts/run_jobs.07.REPmask_corr.08.MASK`
   * **Parallelizable:** task_mem="15G", task_cores="8"
1. `Catrack -v preads.db rep2`
1. `bash scripts/run_jobs.07.REPmask_corr.09.RM`




### 2.5 - Generate check-point files
```bash
HPC.daligner -v -B128 -D512 -M40 -t60 -k20 -h256 -e.96 -l1000 -s100 -mtan -mrep2 -T8 -H10000 preads 1-$LB >| run_jobs.sh
mkdir preads-fofn-abs
echo "echo NOOP preads" > preads-fofn-abs/noop.sh
touch preads-fofn-abs/run.sh.done
cp ../0-rawreads/input_preads.fofn preads-fofn-abs/
touch prepare_pdb.sh
touch pdb_build_done
touch ../0-rawreads/rdb_build_done
```
>**Note:** Remeber to add `-mtan -mrep2` tags to DAligner options as to take advantage of repeats marking



## 3 - Assembly
### 3.1 - Falcon
  * **Note**: Falcon must be run starting from the corrected reads. In General tab set:
```
input_fofn = 0-rawreads/input_preads.fofn
input_type = preads
```
  * **Note**: Add `-mtan -mrep2` tags to DAligner options to use repeats marking information

```bash
$PYTHONUSERBASE/bin/fc_run.py fc_run.cfg
```
