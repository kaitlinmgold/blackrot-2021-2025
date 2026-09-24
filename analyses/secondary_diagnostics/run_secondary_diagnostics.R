#!/usr/bin/env Rscript
# Executed diagnostic audit of the original within-year program+block OLS
# screens. Original specifications and Holm families are retained. HC3 and
# leave-one-plot-out results are sensitivity analyses, not replacement tests.
args <- commandArgs(trailingOnly=TRUE)
input <- if(length(args)>0) args[1] else '.'
out <- if(length(args)>1) args[2] else 'rerun/secondary_diagnostics'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
d <- read.csv(file.path(input,'data/clean_plot_data_public.csv'), stringsAsFactors=FALSE)
original <- read.csv(file.path(input,'results/descriptive/program_vs_utc_contrasts.csv'), stringsAsFactors=FALSE)
unresolved_ids <- unique(d$Public_Program_ID[grepl('2022_T2_T3_block_assignment_unresolved',d$QC_Flags, fixed=TRUE)])
d <- d[!d$Public_Program_ID %in% unresolved_ids,]
stopifnot(nrow(d)==256L, setequal(unresolved_ids,c('Y2022_P02','Y2022_P03')))
rows <- list(); diags <- list(); residuals_all <- list(); groups <- list(); loo <- list(); fits <- list()
model_contrasts <- function(fit, ids) {
  X <- model.matrix(fit); bread <- solve(crossprod(X)); e <- residuals(fit); h <- hatvalues(fit)
  V3 <- bread %*% crossprod(X, X * as.numeric((e/(1-h))^2)) %*% bread
  # Independently check the direct sandwich calculation against the package.
  stopifnot(max(abs(V3-sandwich::vcovHC(fit,type='HC3')))<1e-10)
  j <- match(paste0('Program_Model',ids), names(coef(fit)))
  est <- coef(fit)[j]; se <- sqrt(diag(vcov(fit)))[j]; se3 <- sqrt(diag(V3))[j]
  p <- 2*pt(-abs(est/se),df.residual(fit)); p3 <- 2*pt(-abs(est/se3),df.residual(fit))
  data.frame(Public_Program_ID=ids, Difference_vs_UTC_Logit=as.numeric(est), SE=se,
    t=est/se, Residual_df=df.residual(fit), p_raw=p,p_Holm=p.adjust(p,'holm'),
    HC3_SE=se3,HC3_t=est/se3,HC3_p_raw=p3,HC3_p_Holm=p.adjust(p3,'holm'),
    HC3_lower95=est-qt(.975,df.residual(fit))*se3,
    HC3_upper95=est+qt(.975,df.residual(fit))*se3, row.names=NULL)
}
for(endpoint in c('Incidence','Severity')) for(yr in 2021:2025) {
  dy <- d[d$Year==yr,]; utc <- unique(dy$Public_Program_ID[dy$Program_Class=='Untreated'])
  lev <- c(utc,setdiff(sort(unique(dy$Public_Program_ID)),utc))
  dy$Program_Model <- factor(dy$Public_Program_ID,levels=lev); dy$Block_Model <- factor(dy$Block)
  response <- paste0('HA_Logit_',endpoint)
  fit <- lm(reformulate(c('Program_Model','Block_Model'),response=response),dy)
  key <- paste(endpoint,yr,sep='_'); fits[[key]] <- list(fit=fit,data=dy,levels=lev)
  stopifnot(qr(model.matrix(fit))$rank==ncol(model.matrix(fit)))
  ct <- model_contrasts(fit,lev[-1]); ct$Endpoint <- endpoint; ct$Year <- yr
  md <- dy[match(ct$Public_Program_ID,dy$Public_Program_ID),]
  ct$Program_Label <- md$Program_Label; ct$Program_Class <- md$Program_Class
  rows[[key]] <- ct
  e <- residuals(fit); rs <- rstandard(fit); rt <- rstudent(fit); h <- hatvalues(fit); cd <- cooks.distance(fit)
  aux <- lm(I(e^2)~fitted(fit)); bp <- nrow(dy)*summary(aux)$r.squared
  # These p-values are descriptive flags. Residuals are correlated after fitting,
  # and program-level residual groups contain only four plots.
  bf_y <- abs(e-ave(e,dy$Program_Model,FUN=median)); bf <- anova(lm(bf_y~dy$Program_Model))
  diags[[key]] <- data.frame(Endpoint=endpoint,Year=yr,N=nrow(dy),Programs=length(lev),
    Residual_df=df.residual(fit),Residual_SD=summary(fit)$sigma,
    Shapiro_W=as.numeric(shapiro.test(e)$statistic),Shapiro_p=shapiro.test(e)$p.value,
    BP_fitted_LM=bp,BP_fitted_p=pchisq(bp,1,lower.tail=FALSE),
    Brown_Forsythe_program_F=bf[1,'F value'],Brown_Forsythe_program_p=bf[1,'Pr(>F)'],
    Max_abs_external_studentized=max(abs(rt)),Max_Cooks_distance=max(cd),
    N_Cook_gt_4_over_N=sum(cd>4/nrow(dy)),Max_hat=max(h),
    N_zero_raw=sum(dy[[paste0('Cluster_',endpoint,'_Pct')]]==0),
    N_hundred_raw=sum(dy[[paste0('Cluster_',endpoint,'_Pct')]]==100),
    Original_significant_lower=sum(ct$p_Holm<.05 & ct$Difference_vs_UTC_Logit<0),
    HC3_significant_lower=sum(ct$HC3_p_Holm<.05 & ct$Difference_vs_UTC_Logit<0))
  residuals_all[[key]] <- data.frame(Endpoint=endpoint,Year=yr,dy[,c('Public_Plot_ID','Public_Program_ID','Program_Label','Program_Class','Block')],
    Observed=dy[[response]],Fitted=fitted(fit),Residual=e,Standardized_residual=rs,
    Externally_studentized_residual=rt,Hat=h,Cooks_distance=cd)
  groups[[key]] <- do.call(rbind,lapply(lev,function(id) {
    ix <- dy$Public_Program_ID==id
    data.frame(Endpoint=endpoint,Year=yr,Public_Program_ID=id,Program_Label=dy$Program_Label[which(ix)[1]],
      N=sum(ix),Residual_SD=sd(e[ix]),Residual_min=min(e[ix]),Residual_max=max(e[ix]),
      Outcome_min=min(dy[[response]][ix]),Outcome_max=max(dy[[response]][ix]))
  }))
  for(i in seq_len(nrow(dy))) {
    fi <- update(fit,data=dy[-i,]); cti <- model_contrasts(fi,lev[-1])
    cti$Endpoint <- endpoint; cti$Year <- yr; cti$Deleted_Plot_ID <- dy$Public_Plot_ID[i]
    cti$Deleted_Program_ID <- dy$Public_Program_ID[i]; loo[[length(loo)+1]] <- cti
  }
  draw <- function() {
    par(mfrow=c(2,2),mar=c(4,4,2.6,1),oma=c(2,0,2,0))
    cols <- c('Untreated'='#333333','Biological only'='#188E8A','Mixed'='#C28132','Conventional only'='#7B65AE')[dy$Program_Class]
    plot(fitted(fit),rs,pch=19,col=cols,xlab='Fitted boundary-adjusted logit',ylab='Standardized residual',main='Residuals against fitted values'); abline(h=0,lty=2,col='grey50'); lines(lowess(fitted(fit),rs),col='#BD3E46',lwd=2)
    qqnorm(rs,pch=19,col=cols,main='Normal quantile plot'); qqline(rs,col='#BD3E46',lwd=2)
    plot(fitted(fit),sqrt(abs(rs)),pch=19,col=cols,xlab='Fitted boundary-adjusted logit',ylab='Square root |standardized residual|',main='Spread against fitted values'); lines(lowess(fitted(fit),sqrt(abs(rs))),col='#BD3E46',lwd=2)
    plot(seq_along(cd),cd,type='h',ylim=c(0,max(cd)*1.25),xlab='Plot (data order)',ylab="Cook's distance",main='Influence of individual plots'); abline(h=4/nrow(dy),lty=2,col='#BD3E46'); lab <- which(cd>4/nrow(dy)); text(lab,cd[lab],sub(paste0('Y',yr,'_'),'',dy$Public_Plot_ID[lab]),pos=3,cex=.6)
    mtext(paste(yr,endpoint,'- secondary program + block model'),outer=TRUE,font=2,cex=1.15)
    mtext('Black: untreated   Teal: biological only   Gold: mixed   Purple: conventional only',side=1,outer=TRUE,cex=.8)
  }
  pdf(file.path(out,paste0(key,'_diagnostics.pdf')),width=10,height=8); draw(); dev.off()
  png(file.path(out,paste0(key,'_diagnostics.png')),width=1800,height=1440,res=180,type='cairo'); draw(); dev.off()
}
contrasts <- do.call(rbind,rows); diag <- do.call(rbind,diags); residuals <- do.call(rbind,residuals_all); group_stats <- do.call(rbind,groups); leave_one_out <- do.call(rbind,loo)
check <- merge(contrasts,original,by=c('Endpoint','Year','Public_Program_ID'),suffixes=c('_rerun','_saved'))
fields <- c('Difference_vs_UTC_Logit','SE','t','Residual_df','p_raw','p_Holm')
agreement <- do.call(rbind,lapply(fields,function(nm) data.frame(Field=nm,Max_absolute_difference=max(abs(check[[paste0(nm,'_rerun')]]-check[[paste0(nm,'_saved')]])))))
stopifnot(nrow(check)==118,all(agreement$Max_absolute_difference<1e-10))
write.csv(agreement,file.path(out,'original_output_agreement.csv'),row.names=FALSE)
write.csv(contrasts,file.path(out,'original_vs_HC3_contrasts.csv'),row.names=FALSE)
write.csv(diag,file.path(out,'diagnostic_summary.csv'),row.names=FALSE)
write.csv(residuals,file.path(out,'plot_residuals_and_influence.csv'),row.names=FALSE)
write.csv(group_stats,file.path(out,'program_residual_spread.csv'),row.names=FALSE)
write.csv(leave_one_out,file.path(out,'leave_one_plot_out_contrasts.csv'),row.names=FALSE)
loo_groups <- split(leave_one_out,list(leave_one_out$Endpoint,leave_one_out$Year,leave_one_out$Deleted_Plot_ID),drop=TRUE)
loo_summary <- do.call(rbind,lapply(loo_groups,function(x) data.frame(Endpoint=x$Endpoint[1],Year=x$Year[1],Deleted_Plot_ID=x$Deleted_Plot_ID[1],
  N_significant_lower=sum(x$p_Holm<.05 & x$Difference_vs_UTC_Logit<0),
  Max_p_Holm=max(x$p_Holm),N_HC3_significant_lower=sum(x$HC3_p_Holm<.05 & x$Difference_vs_UTC_Logit<0))))
write.csv(loo_summary,file.path(out,'leave_one_plot_out_summary.csv'),row.names=FALSE)
screen_summary <- do.call(rbind,lapply(split(diag,list(diag$Endpoint,diag$Year),drop=TRUE),function(x) {
  q <- loo_summary[loo_summary$Endpoint==x$Endpoint & loo_summary$Year==x$Year,]
  data.frame(Endpoint=x$Endpoint,Year=x$Year,Contrasts=x$Programs-1,
    Original_significant_lower=x$Original_significant_lower,
    HC3_significant_lower=x$HC3_significant_lower,
    Leave_one_out_min_significant_lower=min(q$N_significant_lower),
    Leave_one_out_max_significant_lower=max(q$N_significant_lower),
    Leave_one_out_refits=nrow(q))
}))
write.csv(screen_summary,file.path(out,'screen_sensitivity_summary.csv'),row.names=FALSE)
jsonlite::write_json(list(audit_date='2026-09-11',scope='Separate within-year fixed program + block linear models',
  response='Previously specified boundary-adjusted logits of plot-level incidence or severity',
  original_contrasts_reproduced=118L,max_absolute_difference=max(agreement$Max_absolute_difference),
  Holm_family='All eligible program-versus-untreated comparisons within endpoint and year, retained in every sensitivity',
  HC3_inference='HC3 sandwich standard errors with original model residual degrees of freedom; approximate sensitivity, not an exact finite-sample test',
  exclusion='Eight unresolved 2022 plots excluded from these program-level models as originally specified',
  models=screen_summary),file.path(out,'screen_sensitivity_summary.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
changed <- contrasts[(contrasts$p_Holm<.05)!=(contrasts$HC3_p_Holm<.05),]
write.csv(changed,file.path(out,'comparisons_changing_Holm_threshold_HC3.csv'),row.names=FALSE)
saveRDS(fits,file.path(out,'within_year_program_models.rds'))
session_lines <- capture.output(sessionInfo())
session_lines <- gsub(normalizePath('.'),'<PROJECT>',session_lines,fixed=TRUE)
session_lines <- gsub(Sys.getenv('HOME'),'<HOME>',session_lines,fixed=TRUE)
writeLines(session_lines,file.path(out,'R_session_info.txt'))
print(agreement,row.names=FALSE); print(diag,row.names=FALSE)
