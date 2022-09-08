FROM rocker/verse:4.0.2

ENV RENV_VERSION 0.11.0
RUN R -e "install.packages('remotes',repos=c(CRAN='https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

COPY renv.lock /home/rstudio/
COPY environment.yml /home/rstudio/
COPY .Rprofile /home/rstudio/
COPY .gitignore /home/rstudio/
COPY Dockerfile /home/rstudio/
COPY LICENSE /home/rstudio/
COPY README.md /home/rstudio/
COPY pest.Rmd /home/rstudio/

ENV RENV_PATHS_LIBRARY renv/library

RUN mkdir -p renv
RUN mkdir -p R
COPY renv/activate.R /home/rstudio/renv/
COPY renv/settings.dcf /home/rstudio/renv/
COPY R/check_install_load-function.R /home/rstudio/R/
COPY R/conduct_dea-function.R /home/rstudio/R/
COPY R/divnn_calibrator_evaluator-function.R /home/rstudio/R/
COPY R/download_annotation-function.R /home/rstudio/R/
COPY R/eval_divnn-function.R /home/rstudio/R/
COPY R/refresh_session-function.R /home/rstudio/R/
COPY R/suspected_outliers-function.R /home/rstudio/R/
COPY R/take_common_genes-function.R /home/rstudio/R/
COPY R/test_transformer-function.R /home/rstudio/R/
COPY R/trainer_generator-function.R /home/rstudio/R/

# RUN R --vanilla -s -e 'renv::restore()'
