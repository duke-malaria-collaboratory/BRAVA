# Running the haplotype calling pipeline with Snakemake

[Snakemake](https://snakemake.readthedocs.io/en/stable) is a workflow manager
that enables massively parallel and reproducible
analyses.
Snakemake is a suitable tool to use when you can break a workflow down into
discrete steps, with each step having input and output files.

We provide this example workflow as a template to get started running the pipeline with Snakemake.
To adjust to your specific data, you can customize the config.yaml file.

Overview: Using human or mosquito samples from Webuye, Kenya, evaluate the extent to which there are unique haplotypes among two polymorphic gene targets: AMA, CSP.

Method: Targeted amplicon deep sequencing which produces forward and reverse fastq files for each sample.

For more details on Snakemake, see the
[Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

## The Workflow

The [`Snakefile`](Snakefile) contains rules which define the output files we want and how to make them.
Snakemake automatically builds a directed acyclic graph (DAG) of jobs to figure
out the dependencies of each of the rules and what order to run them in.
This workflow processes the data set by cleaning the sequencing reads, performing haplotype calling, and censoring the haplotypes to render numerous data analysis outputs, resulting in a final table that delineates the resulting haplotypes after censoring.

`clean_sequencing_reads` cleans, filters, and maps the raw reads. It uses BBmap to map all reads from the reference sequences to differentiate between the two targets, CutAdapt to trim the primers and adapter sequences from sequencing reads, and uses Trimmomatic to quality filter reads if average of every 4 nucleotides had a Phred Quality Score < 15 or was less than 80 nucleotides long. This is the first step in read processing on the cluster.

`organize_folders` creates an all_samples folder within the out/fastq/{target} folder and moves all forward and reverse sequences into that folder.

`call_haplotypes` call haplotypes for your target using the DADA2 program. This is the second step in the read processing on the cluster.

`censor_haplotypes` censors falsely detected haplotypes. Censoring criteria is applied in this order:
1. Haplotypes that occur in < 250 of the sample’s reads are removed.
2. Haplotypes that occur in < 3% of the sample’s reads are removed.
3. Haplotypes that are a different length than the majority of haplotypes (300 nucleotides for pfama1, 288 nucleotides for pfcsp) are removed.
4. For haplotypes that have 1 SNP difference, occur in the same sample, and have a >8 times read depth difference between them within that sample, removed the hapltoype with the lower read depth from that sample.
5. If a haplotype is defined by a single variant position that is only variable within that haplotype, then it is removed.

## Quick Start

1. Clone or download this repo.

    ``` sh
    git clone https://github.com/kathiehuang/haplotype_calling_pipeline.git
    ```
    Alternatively, if you're viewing this on GitHub,
    you can click the green `Use this template` button to create
    your own version of the repo on GitHub, then clone it.

1. Install the dependencies.

    1. If you don't have conda yet, we recommend installing 
       [miniconda](https://docs.conda.io/en/latest/miniconda.html).
       
    1. Next, install [mamba](https://mamba.readthedocs.io/en/latest/), 
       a fast drop-in replacement for conda:
       
       ``` sh
       conda install mamba -n base -c conda-forge
       ```
       
    1. Finally, create the environment and activate it:
    
       ``` sh
       mamba env create -f environment.yml
       conda activate haplotype_calling
       ```
       
    - Alternatively, you can install the dependencies listed in
    [`environment.yml`](environment.yml) however you like.

1. Edit the configuration file [`config.yml`](config.yml).
    - `ncores`: the number of cores to use. Do not exceed the number of cores you have available.
    - `target`: the polymorphic gene target.
    - `refs`: the path to the folder containing reference sequences for the polymorphic gene target that will be used to map the raw reads to the appropriate gene targets of interest
    - `pair1`: the path to the folder containing the forward reads.
    - `pair2`: the path to the folder containing the reverse reads.
    - `forward`: the path to the file with the list of forward primers.
    - `rev`: the path to the file with the list of reverse primers.
    - `out`: the name of desired output folder.
    - `length`: the length of the majority of haplotypes for the target.

    You can leave these options as-is if you'd like to first make sure the
    workflow runs without error on your machine before using your own dataset
    and custom parameters.

    The default config file is suitable for initial testing,
    but we recommend using more cores if available.

1. Do a dry run to make sure the snakemake workflow is valid.

    ``` sh
    snakemake -n
    ```

1. Run the workflow.

    Run it **locally** with:
    ``` sh
    snakemake
    ```

    To run the workflow on an **HPC with Slurm**:

    1. Edit your email (`YOUR_EMAIL_HERE`), Slurm account (`YOUR_ACCOUNT_HERE`), and other Slurm parameters as needed in:

        - [`code/submit_slurm.sh`](code/submit_slurm.sh)
        - [`config/cluster.json`](config/cluster.json)

    1. Submit the snakemake workflow with:

        ``` sh
        sbatch code/submit_slurm.sh
        ```

        The main job will then submit all other snakemake jobs, allowing
        independent steps of the workflow to run on different nodes in parallel.
        Slurm output files will be written to `log/hpc/`.


## More resources

- [Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)
- [conda user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html)