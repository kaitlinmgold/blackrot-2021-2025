#!/usr/bin/env Rscript
# Fits additional reduced-data ML models for the exact archived interaction tests.
suppressPackageStartupMessages({library(glmmTMB); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
saved <- readRDS(args[[1]])
d <- saved$data
stopifnot(nrow(d)==256L, nlevels(d$Public_Program_ID)==64L)
out <- args[[2]]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
stopifnot(!file.exists(file.path(out, 'interaction_models.rds')))
ctrl <- glmmTMBControl(optCtrl=list(iter.max=20000L, eval.max=20000L))
ctrl_bfgs <- glmmTMBControl(optimizer=optim, optArgs=list(method='BFGS'), optCtrl=list(maxit=20000L))
usable <- function(m) !inherits(m,'try-error') && isTRUE(m$sdr$pdHess) && identical(unname(m$fit$convergence),0L) && is.finite(as.numeric(logLik(m)))
fit_sev <- function(formula) {
  m <- try(glmmTMB(formula, data=d, family=gaussian(), REML=FALSE, dispformula=~Year+Program_Class, control=ctrl),silent=TRUE)
  if(!usable(m)) m <- try(glmmTMB(formula, data=d, family=gaussian(), REML=FALSE, dispformula=~Year+Program_Class, control=ctrl_bfgs),silent=TRUE)
  stopifnot(usable(m)); m
}
inc_formula <- cbind(Cluster_Positive_Count, Cluster_Negative_Count) ~ Year * Program_Class + (1|BlockYear) + (1|Public_Program_ID)
# The unrestricted incidence model may be nonregular, as in the original analysis;
# only its full-model objective evaluated at the additive null is used by the score test.
inc_full <- glmmTMB(inc_formula, data=d, family=betabinomial(link='logit'), control=ctrl)
sev_add <- fit_sev(HA_Logit_Severity ~ Year + Program_Class + (1|BlockYear) + (1|Public_Program_ID))
sev_full <- fit_sev(HA_Logit_Severity ~ Year * Program_Class + (1|BlockYear) + (1|Public_Program_ID))
saved$incidence$full <- inc_full
saved$severity$additive_ml <- sev_add
saved$severity$full_ml <- sev_full
saveRDS(saved, file.path(out, 'interaction_models.rds'))
write.csv(d, file.path(out, 'clean_plot_data_public.csv'),row.names=FALSE)
info <- function(m) list(convergence=unname(m$fit$convergence), pd_hessian=isTRUE(m$sdr$pdHess), logLik=as.numeric(logLik(m)))
write_json(list(incidence_full=info(inc_full),severity_additive_ml=info(sev_add),severity_full_ml=info(sev_full)),file.path(out,'interaction_model_fit_info.json'),pretty=TRUE,auto_unbox=TRUE,digits=15)
writeLines(capture.output(sessionInfo()),file.path(out,'R_session_info.txt'))
