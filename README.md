# IPOD
----

This repository contains code for performing IPOD data analysis.
To promote efficient code distribution and analysis reproducibility,
we providea singularity container. For instructions on singularity
container use, see the [Singularity Use](#singularity-use) section below.

## Singularity Use

Download and install singularity to your system.

To use our singularity containers, email schroedj@umich.edu for access to
our Google Shared Drive. Once you are given access to the Shared Drive
you can follow [this link][singularity-link]. The version of the singularity
container with the most up-to-date code is in the "current" folder. Older
versions of the container can be found in the "archive" folder.

Download the sif file. These files are quite large for IPOD analysis (>1.5 GB).
This file contains all the necessary components to run our IPOD pipeline.
We recommend you also download the example data present in the Google Shared
Drive containing the singularity containers.

Right-click the "test_data" folder, then click "Download".
Once the download is complete, unzip the data to a new directory called "test"
by running `unzip <zipfile> -d test` from the command line, where you'll
substitute the name of your zip file for `<zipfile>`.

To test the downloaded singularity file on our provided dataset, enter the
newly-created "test" directory by running `cd test`. Now run the singularity
container by running the following at the command line:

```bash
singularity run ipod_<version>.sif
```

Substitute the appropriate version number for `<version>`.
If you need access to remote file systems within your singularity container
once it's run, you'll need to run something like the following:

```bash
singularity run -B <src>:<dest> ipod_<version>.sif
```

Here, `<src>` must be replaced by the remote location to be accessed by
the singularity container. Within the container's environment, `<dest>`
will the the location of the remote file system. If multiple locations
must be bound, simply add more `-B <src>:<dest>` arguments to your call
to the singularity command.

At this point, if your singularity container ran properly, your command prompt
should look something like `[ipod_<version>]$`. Once the singularity container
is running, you'll need to activate the conda environment by running the
following:

```bash
conda activate ipod
```

At this point, your command prompt should look something like
`(ipod) [ipod_<version>]$`. To test whether your environment is propery set
up, enter the "test" directory, and run the following command:

```bash
python ~/src/ipod/src_for_distrib/drivers/run_all_driver.py main.conf
```

This will run all steps in the IPOD analysis pipeline, through quantifying
IPOD signals. Once the pipeline has finished running, you should see the
following files in **insert file location here**: **insert file names here**.

[singularity-link]: https://drive.google.com/drive/folders/1k85Ew32F2Ek36yjEVvLss4OWv_EE7rUv?usp=sharing
