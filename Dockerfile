FROM rocker/ml-verse:4.0.2

ENV RENV_VERSION 0.14.0
RUN R -e "install.packages('remotes',repos=c(CRAN='https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

COPY .gitignore /home/rstudio/
COPY Dockerfile /home/rstudio/
COPY LICENSE /home/rstudio/
COPY pest.Rmd /home/rstudio/
COPY pest.R /home/rstudio/
COPY README.md /home/rstudio/

RUN mkdir -p R
COPY R/calib_comparison-function.R /home/rstudio/R/
COPY R/conduct_dea-function.R /home/rstudio/R/
COPY R/dea_qn-function.R /home/rstudio/R/
COPY R/deg_train-function.R /home/rstudio/R/
COPY R/deg-function.R /home/rstudio/R/
COPY R/disease_train-function.R /home/rstudio/R/
COPY R/divnn_calibrator_evaluator-function.R /home/rstudio/R/
COPY R/divnn_pre_object-function.R /home/rstudio/R/
COPY R/download_annotation-function.R /home/rstudio/R/
COPY R/eval_divnn-function.R /home/rstudio/R/
COPY R/get_auroc-function.R /home/rstudio/R/
COPY R/mb_predictor-function.R /home/rstudio/R/
COPY R/model_component-function.R /home/rstudio/R/
COPY R/opt_feat_disease_fn-function.R /home/rstudio/R/
COPY R/opt_pred_disease_fn-function.R /home/rstudio/R/
COPY R/re_evalm-function.R /home/rstudio/R/
COPY R/refresh_session-function.R /home/rstudio/R/
COPY R/result_comparison-function.R /home/rstudio/R/
COPY R/show_shortest_path-function.R /home/rstudio/R/
COPY R/std_no_outlier-function.R /home/rstudio/R/
COPY R/surrogate_deg-function.R /home/rstudio/R/
COPY R/suspected_outliers-function.R /home/rstudio/R/
COPY R/take_common_genes-function.R /home/rstudio/R/
COPY R/test_transformer-function.R /home/rstudio/R/
COPY R/trainer_generator-function.R /home/rstudio/R/
COPY R/tree_eval-function.R /home/rstudio/R/