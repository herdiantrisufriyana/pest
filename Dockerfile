FROM rocker/verse:4.0.2

ENV RENV_VERSION 0.11.0
RUN R -e "install.packages('remotes',repos=c(CRAN='https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

WORKDIR /home/pest
COPY renv.lock renv.lock
COPY environment.yml environment.yml

ENV RENV_PATHS_LIBRARY renv/library

RUN mkdir -p renv
RUN mkdir -p R
COPY .Rprofile .Rprofile
COPY .gitignore .gitignore
COPY Dockerfile Dockerfile
COPY LICENSE LICENSE
COPY README.md README.md
COPY pest.Rmd pest.Rmd
COPY renv/activate.R renv/activate.R
COPY renv/settings.dcf renv/settings.dcf
COPY R/check_install_load-function.R R/check_install_load-function.R
COPY R/conduct_dea-function.R R/conduct_dea-function.R
COPY R/divnn_calibrator_evaluator-function.R R/divnn_calibrator_evaluator-function.R
COPY R/download_annotation-function.R R/download_annotation-function.R
COPY R/eval_divnn-function.R R/eval_divnn-function.R
COPY R/refresh_session-function.R R/refresh_session-function.R
COPY R/suspected_outliers-function.R R/suspected_outliers-function.R
COPY R/take_common_genes-function.R R/take_common_genes-function.R
COPY R/test_transformer-function.R R/test_transformer-function.R
COPY R/trainer_generator-function.R R/trainer_generator-function.R

# RUN R --vanilla -s -e 'renv::restore()'
