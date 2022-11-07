FROM condaforge/mambaforge:22.9.0-1

RUN mamba install \
    --copy \
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
    perl-iterator-simple \
    bioconductor-quantsmooth \
    r-ggplot2 && \
    mamba clean --all --yes

ADD version_1.0 /usr/share/UPDio

RUN ln -s /usr/share/UPDio/UPDio.pl /usr/bin/updio