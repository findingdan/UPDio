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

RUN mkdir /usr/share/UPDio
ADD ./UPDio.pl /usr/share/UPDio/UPDio.pl 
ADD ./scripts /usr/share/UPDio/scripts 
ADD ./exome_designs /usr/share/UPDio/exome_designs
ADD ./pre_processing /usr/share/UPDio/pre_processing
ADD ./sample_data /usr/share/UPDio/sample_data

RUN ln -s /usr/share/UPDio/UPDio.pl /usr/bin/updio