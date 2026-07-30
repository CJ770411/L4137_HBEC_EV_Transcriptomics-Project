## Quick use

Test data has been included to allow functionality testing of the workflow. N.B. This will not produce the same results produced by the study.

**1. Create environments:**

```
conda env create -f environments/L4137.yaml
conda env create -f environments/L4137_cutadapt.yaml
conda env create -f environments/L4137_mirdeep2.yaml
```

**2. Run the high-performance computer (HPC) component of the workflow:**

```
cd scripts/workflow

sbatch run.sh
```

**3. Sequentially run the R scripts from the directory containing the script:**

```
cd scripts/workflow/*   
```