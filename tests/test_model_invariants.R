#!/usr/bin/env Rscript
# Scientific invariants, using synthetic output only; does not need study data.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1)stop('Usage: Rscript tests/test_model_invariants.R RESULTS_DIR')
file_arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
repo<-dirname(dirname(normalizePath(sub('^--file=','',file_arg[1]))))
source(file.path(repo,'R/model_core.R'))
run<-normalizePath(args[1]);dat<-readRDS(file.path(run,'analysis_data.rds'))
folds<-readRDS(file.path(run,'outer_folds.rds'));cfg<-readRDS(file.path(run,'run_signature.rds'))$config
paths<-list.files(file.path(run,'models'),pattern='[.]rds$',full.names=TRUE)
for(path in paths){
  z<-readRDS(path)
  stopifnot(!any(z$predictions$SUBJID %in% z$model$training_ids))
}
avg<-read.csv(file.path(run,'participant_predictions.csv'))
stopifnot(all(avg$n_repeats==cfg$repeats))

# Perturb every outer-test outcome/assay and all follow-up outcomes/assays,
# then repeat the entire first outer job. Fitted base and integration parameters
# must be unchanged. This includes inner tuning and integration cross-fitting.
perturbed<-dat
test_ids<-rownames(folds)[folds[,1]==1]
alter<-perturbed$data$SUBJID %in% test_ids | perturbed$data$Week==24
perturbed$y[alter]<-perturbed$y[alter]+1000
perturbed$x[alter,]<-perturbed$x[alter,]+10000
tmp<-tempfile('model_invariance_');dir.create(tmp)
cfg$output_dir<-tmp
fit_outer_job(1L,1L,perturbed,folds,cfg)
for(s in c('MetIR','ProtIR','PhenoIR')){
  nm<-paste0(s,'_direct_r01_f01.rds')
  old<-readRDS(file.path(run,'models',nm));new<-readRDS(file.path(tmp,'models',nm))
  if(s=='PhenoIR')stopifnot(isTRUE(all.equal(old$model$weights,new$model$weights,tolerance=1e-12)))
  else stopifnot(isTRUE(all.equal(old$model$beta,new$model$beta,tolerance=1e-12)),
    identical(old$model$prep,new$model$prep),identical(old$model$selected,new$model$selected))
}

# Only changing W24 outcome labels must not change frozen predictions at all.
changed<-dat;changed$y[changed$data$Week==24]<-changed$y[changed$data$Week==24]+1000
for(s in c('MetIR','ProtIR')){
  z<-readRDS(file.path(run,'models',paste0(s,'_direct_r01_f01.rds')))
  rows<-which(dat$data$SUBJID %in% test_ids)
  a<-ir_predict(z$model,dat,rows);b<-ir_predict(z$model,changed,rows)
  stopifnot(identical(a$predicted_log_homa,b$predicted_log_homa))
}
unlink(tmp,recursive=TRUE)
cat('PASS: held-out separation, repeated prediction coverage, no outer-test/follow-up training leakage, frozen follow-up label invariance.\n')
