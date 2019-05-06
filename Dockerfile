FROM rocker/verse:3.5.3

RUN R -e 'devtools::install_dev_deps()'

