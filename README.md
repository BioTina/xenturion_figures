# xenturion_figures
Snakemake pipeline to reproduce Xenturion bioinformatics main figures panels.

# XENTURION is a population-level multidimensional resource of xenografts and tumoroids from metastatic colorectal cancer patients

Snakemake pipeline to reproduce Xenturion bioinformatics main figures panels.
Keep in mind that final figures were assembled manually (e.g. for
legend positioning) and that some rules produce more than one plot
(with/without legend) to ease
those manual fixes.

A Dockerfile with the same R version, R packages and overall linux
environment used to produce the Figures is available here:
https://github.com/vodkatad/snakemake_docker/blob/master/Dockerfiles/godot/Dockerfile

# Instructions to obtain figure panels on your pc

Install docker and git, set (or substitute) the $WDIR bash variable to a path on
your local system where you want to clone the repo and save the figure panels (/tmp here), then:

```
$ WDIR=/tmp
$ cd $WDIR 
$ git clone git@github.com:BioTina/xenturion_figures.git 
$ docker pull egrassi/godot:05_2024
$ docker run -u $(id -u $(whoami)) -e WDIR=$WDIR --rm -it --volume $WDIR:$WDIR egrassi/godot:05_2024
# inside the container:
$ cd ${WDIR}/xenturion_figures/dataset
$ snakemake oncoprint_0.05.pdf

```
