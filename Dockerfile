FROM rocker/binder:3.5.3

USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}

## Become normal user again
USER ${NB_USER}
RUN R -e 'source("install.R")'
