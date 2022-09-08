FROM rocker/verse:4.0.2

ENV RENV_VERSION 0.11.0
RUN R -e "install.packages('remotes',repos=c(CRAN='https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

WORKDIR /project
COPY renv.lock renv.lock
COPY environment.yml environment.yml

ENV RENV_PATHS_LIBRARY renv/library

RUN mkdir -p renv
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.dcf renv/settings.dcf

RUN R --vanilla -s -e 'renv::restore()'
