FROM rocker/verse:4.0.2

ENV RENV_VERSION 0.11.0
RUN R -e "install.packages('remotes',repos=c(CRAN='https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

COPY renv.lock /home/
COPY environment.yml /home/

ENV RENV_PATHS_LIBRARY renv/library

RUN mkdir -p renv
RUN mkdir -p R
COPY .Rprofile /home/
COPY .gitignore /home/
COPY Dockerfile /home/
COPY LICENSE /home/
COPY README.md /home/
COPY pest.Rmd /home/
COPY renv/activate.R /home/renv/
COPY renv/settings.dcf /home/renv/
COPY R/check_install_load-function.R /home/R/
COPY R/conduct_dea-function.R /home/R/
COPY R/divnn_calibrator_evaluator-function.R /home/R/
COPY R/download_annotation-function.R /home/R/
COPY R/eval_divnn-function.R /home/R/
COPY R/refresh_session-function.R /home/R/
COPY R/suspected_outliers-function.R /home/R/
COPY R/take_common_genes-function.R /home/R/
COPY R/test_transformer-function.R /home/R/
COPY R/trainer_generator-function.R /home/R/

# RUN R --vanilla -s -e 'renv::restore()'
