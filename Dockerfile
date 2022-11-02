FROM condaforge/mambaforge:22.9.0-1

RUN mamba install \
    --strict-channel-priority  \
    --channel conda-forge  \
    --channel bioconda \
    --yes \
    perl \
    perl-statistics-r \
    perl-path-class \
    perl-vcftools-vcf \
    perl-list-moreutils \
    perl-math-round \
    bioconductor-quantsmooth \
    r-ggplot2

# Should be replaced by the conda installation perl-iterator-simple once merged.
# https://github.com/conda-forge/staged-recipes/pull/20963
RUN cpan Iterator::Simple

ADD version_1.0 /usr/share/UPDio

RUN ln -s /usr/share/UPDio/UPDio.pl /usr/bin/updio