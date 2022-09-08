FROM rocker/verse:4.0.2

ENV RENV_VERSION 0.11.0
RUN R -e "install.packages('remotes',repos=c(CRAN='https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

COPY renv.lock .
COPY environment.yml .
COPY .Rprofile .
COPY .gitignore .
COPY Dockerfile .
COPY LICENSE .
COPY README.md .
COPY pest.Rmd .

ENV RENV_PATHS_LIBRARY renv/library

RUN mkdir -p renv
RUN mkdir -p R
COPY renv/activate.R ./renv/
COPY renv/settings.dcf ./renv/
COPY R/check_install_load-function.R ./R/
COPY R/conduct_dea-function.R ./R/
COPY R/divnn_calibrator_evaluator-function.R ./R/
COPY R/download_annotation-function.R ./R/
COPY R/eval_divnn-function.R ./R/
COPY R/refresh_session-function.R ./R/
COPY R/suspected_outliers-function.R .//R/
COPY R/take_common_genes-function.R ./R/
COPY R/test_transformer-function.R ./R/
COPY R/trainer_generator-function.R ./R/

# RUN R --vanilla -s -e 'renv::restore()'
