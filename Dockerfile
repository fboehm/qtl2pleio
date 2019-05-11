FROM rocker/binder:3.5.3

USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}

## Become normal user again
USER ${NB_USER}
RUN wget https://github.com/fboehm/qtl2pleio/raw/master/install.R && \
R -e 'source("install.R")'
