# Core numerical functions extracted without mathematical changes from the study pipeline.
# Source: irscore_helpers.R in the archived study workflow; see repository documentation.
required_packages <- c('glmnet', 'dplyr', 'matrixStats', 'jsonlite')

read_model_config <- function(path, repo) {
  for (p in required_packages) if (!requireNamespace(p, quietly = TRUE)) stop('Install R package: ', p)
  if (!file.exists(path)) path <- file.path(repo, path)
  cfg <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  resolve <- function(x) {
    if (is.null(x) || !nzchar(x)) return(NULL)
    if (grepl('^(/|[A-Za-z]:)', x)) x else file.path(repo, x)
  }
  cfg$inputs <- lapply(cfg$inputs, resolve)
  cfg$output_dir <- resolve(cfg$output_dir)
  cfg$code_files <- c(file.path(repo,'R','model_core.R'),file.path(repo,'scripts','train_models.R'))
  for (p in c('phenotype', 'metabolites', 'proteins', 'features')) {
    if (is.null(cfg$inputs[[p]]) || !file.exists(cfg$inputs[[p]])) stop('Missing input ', p, ': ', cfg$inputs[[p]])
  }
  if (!is.null(cfg$inputs$folds) && !file.exists(cfg$inputs$folds)) stop('Missing supplied folds file')
  cfg$lambdas <- 10^seq(log10(cfg$lambda_max), log10(cfg$lambda_min), length.out = cfg$n_lambdas)
  cfg$scores <- c('MetIR', 'ProtIR', 'PhenoIR')
  cfg$variants <- 'direct'
  if (!identical(cfg$rules, c('lambda.1se', 'lambda.min'))) stop('Retain rules lambda.1se and lambda.min, in that order')
  ints <- c('seed', 'repeats', 'outer_folds', 'inner_folds', 'calibration_folds', 'minimum_observations', 'workers', 'n_lambdas')
  for (k in ints) {
    if (length(cfg[[k]]) != 1 || !is.finite(cfg[[k]]) || cfg[[k]] != as.integer(cfg[[k]]) || cfg[[k]] < 1) stop('Invalid config: ', k)
    cfg[[k]] <- as.integer(cfg[[k]])
  }
  stopifnot(cfg$outer_folds >= 2, cfg$inner_folds >= 2, cfg$calibration_folds >= 2,
    all(is.finite(cfg$alphas)), all(cfg$alphas > 0 & cfg$alphas <= 1), cfg$lambda_max > cfg$lambda_min,
    cfg$lambda_min > 0, cfg$n_lambdas >= 2, all(c(0, 24) %in% cfg$weeks))
  if (is.null(cfg$age_balance_swaps)) cfg$age_balance_swaps <- 20000L
  if(length(cfg$age_balance_swaps)!=1 || !is.finite(cfg$age_balance_swaps) || cfg$age_balance_swaps<0 ||
     cfg$age_balance_swaps!=as.integer(cfg$age_balance_swaps))stop('age_balance_swaps must be a nonnegative integer')
  cfg$age_balance_swaps<-as.integer(cfg$age_balance_swaps)
  cfg
}

read_table <- function(path, required) {
  key_columns <- intersect(c('sample_id','feature_id'),required)
  classes <- if(length(key_columns)) setNames(rep('character',length(key_columns)),key_columns) else NA
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = classes)
  if (!all(required %in% names(x))) stop('Missing columns in ', basename(path), ': ', paste(setdiff(required, names(x)), collapse = ', '))
  if (anyDuplicated(names(x))) stop('Duplicate column names in ', basename(path))
  x
}

visit_key <- function(x) paste(x$sample_id, x$week, sep = '\r')

numeric_column <- function(x, name) {
  y <- suppressWarnings(as.numeric(x))
  bad <- !is.na(x) & nzchar(trimws(as.character(x))) & is.na(y)
  if (any(bad)) stop('Non-numeric values in ', name)
  y
}

load_model_data <- function(cfg) {
  p <- read_table(cfg$inputs$phenotype, c('sample_id','week','trial','arm','sex','age','homa_ir','bmi'))
  for (n in c('week','sex','age','homa_ir','bmi')) p[[n]] <- numeric_column(p[[n]], n)
  if (anyNA(p$sample_id) || any(!nzchar(p$sample_id)) || anyDuplicated(visit_key(p))) stop('Participant/visit keys must be unique and nonempty')
  p <- p[p$week %in% cfg$weeks, , drop = FALSE]
  p <- p[order(p$sample_id, p$week), , drop = FALSE]
  for (field in c('trial','arm','sex')) {
    inconsistent <- vapply(split(p[[field]],p$sample_id),function(v)length(unique(v[!is.na(v)]))>1L,logical(1))
    if (any(inconsistent)) stop(field, ' must be consistent across visits within a participant')
  }
  f <- read_table(cfg$inputs$features, c('feature_id','omics','label'))
  if (anyDuplicated(f$feature_id) || anyNA(f$feature_id) || any(!nzchar(f$feature_id))) stop('Feature keys must be unique and nonempty')
  if(anyNA(f$label) || any(!nzchar(trimws(f$label))))stop('Every feature needs a nonempty display label')
  if (!all(f$omics %in% c('Metabolomics','Proteomics'))) stop('Feature omics must be Metabolomics or Proteomics')
  excluded <- f[f$label %in% cfg$excluded_feature_labels, , drop = FALSE]
  f <- f[!f$feature_id %in% excluded$feature_id, , drop = FALSE]
  blocks <- list()
  available <- list()
  for (nm in c('metabolites','proteins')) {
    z <- read_table(cfg$inputs[[nm]], c('sample_id','week'))
    z$week <- numeric_column(z$week, paste0(nm, '.week'))
    if (anyDuplicated(visit_key(z))) stop('Duplicate molecular participant/visit rows in ', nm)
    tp <- if (nm == 'metabolites') 'Metabolomics' else 'Proteomics'
    ids <- f$feature_id[f$omics == tp]
    if (length(ids) < 2 || !all(ids %in% names(z))) stop('Need at least two annotated features with matching input columns for ', nm)
    z[ids] <- lapply(ids, function(id) numeric_column(z[[id]], id))
    mat <- as.matrix(z[match(visit_key(p), visit_key(z)), ids, drop = FALSE])
    storage.mode(mat) <- 'double'
    mat[!is.finite(mat)] <- NA_real_
    blocks[[nm]] <- mat
    available[[nm]] <- rowSums(is.finite(mat)) > 0
  }
  baseline <- which(p$week == 0)
  if (!length(baseline) || anyDuplicated(p$sample_id[baseline])) stop('One baseline record per participant is required')
  valid_meta <- is.finite(p$homa_ir) & p$homa_ir > 0 & is.finite(p$bmi) & p$bmi > 0 &
    is.finite(p$age) & p$sex %in% c(1,2) & !is.na(p$trial) & nzchar(p$trial) & !is.na(p$arm) & nzchar(p$arm)
  eligible <- valid_meta & available$metabolites & available$proteins
  if (any(!eligible[baseline])) stop('Baseline inputs must be prefiltered to eligible participants: finite positive HOMA-IR/BMI, age/sex/trial/arm and both molecular platforms. No clinical imputation is performed.')
  eligible <- eligible & p$sample_id %in% p$sample_id[baseline]
  audit <- data.frame(sample_id=p$sample_id, week=p$week, clinical_complete=valid_meta,
    metabolomics_available=available$metabolites, proteomics_available=available$proteins, included=eligible)
  p <- p[eligible, , drop = FALSE]
  x <- cbind(blocks$metabolites[eligible,,drop=FALSE], blocks$proteins[eligible,,drop=FALSE])
  rownames(x) <- NULL
  baseline <- which(p$week == 0)
  if (length(baseline) < cfg$outer_folds * 4) stop('Too few participants for requested nested folds')
  d <- data.frame(SUBJID=p$sample_id, Week=p$week, Trial=p$trial, Arm=p$arm,
    Sex=as.character(p$sex), Age=p$age, BMI=p$bmi, HOMA.IR=p$homa_ir, stringsAsFactors=FALSE)
  cmat <- cbind('(Intercept)'=1, BMI=p$bmi, SexFemale=as.numeric(p$sex==2), Age=p$age)
  ir <- ifelse(p$homa_ir >= ifelse(p$sex==1,4.31,4.51), 'IR','IS')
  bmi_cat <- as.character(cut(p$bmi, c(-Inf,18.5,23,25,Inf), right=FALSE,
    labels=c('Underweight','Normal','Overweight','Obese')))
  if ('ir_status' %in% names(p)) ir <- as.character(p$ir_status)
  if ('bmi_category' %in% names(p)) bmi_cat <- as.character(p$bmi_category)
  if (anyNA(ir) || !all(ir %in% c('IR','IS'))) stop('ir_status, when supplied, must use original split labels IR and IS')
  if (anyNA(bmi_cat) || !all(bmi_cat %in% c('Underweight','Normal','Overweight','Obese'))) stop('Invalid bmi_category')
  # Original split labels were numeric trial codes. Preserve ordering for exact supplied-plan reuse.
  trial_split <- ifelse(p$trial=='CMAP','301',ifelse(p$trial=='CMAS','302',p$trial))
  strata <- data.frame(SUBJID=p$sample_id[baseline], Trial=trial_split[baseline], Sex=as.character(p$sex[baseline]),
    BMI=p$bmi[baseline], Age=p$age[baseline], log_homa=log(p$homa_ir[baseline]), IR=ir[baseline], BMI_category=bmi_cat[baseline], stringsAsFactors=FALSE)
  if (anyNA(strata)) stop('Baseline split variables must be complete')
  features <- data.frame(id=f$feature_id,name=f$label,type=ifelse(f$omics=='Metabolomics','MTB','PRT'),stringsAsFactors=FALSE)
  sets <- list(MetIR=features$id[features$type=='MTB'],ProtIR=features$id[features$type=='PRT'],PhenoIR=features$id)
  list(data=d, phenotype=p, x=x, c=cmat, y=log(p$homa_ir), baseline=baseline, strata=strata,
    features=features, annotations=f, excluded=excluded, sets=sets, audit=audit,
    available=matrix(TRUE,nrow(p),3,dimnames=list(NULL,c('MetIR','ProtIR','PhenoIR'))),cov_ok=rep(TRUE,nrow(p)))
}

balance_outer_age <- function(f, strata, k, seed, swaps) {
  age <- strata$Age
  if (swaps <= 0 || sd(age) < 1e-10) return(f)
  breaks <- unique(quantile(age,0:4/4))
  if (length(breaks) < 3) return(f)
  q <- as.integer(cut(age, breaks, include.lowest=TRUE)); nq <- max(q)
  a <- (age-mean(age))/sd(age)
  joint <- interaction(strata[c('Trial','IR','Sex','BMI_category')],drop=TRUE)
  groups <- split(seq_along(age),joint); groups <- groups[lengths(groups)>1]
  if (!length(groups)) return(f)
  set.seed(seed); original <- f; cap <- tabulate(f,k)
  sums <- as.numeric(rowsum(a,f)); counts <- table(factor(f,levels=seq_len(k)),factor(q,levels=seq_len(nq)))
  target <- outer(cap,tabulate(q,nq)/length(q))
  for (it in seq_len(swaps)) {
    group <- groups[[sample.int(length(groups),1)]]; pair <- sample(group,2)
    i <- pair[1]; j <- pair[2]; u <- f[i]; v <- f[j]
    if (u==v) next
    old <- sum(sums[c(u,v)]^2)+sum((counts[c(u,v),,drop=FALSE]-target[c(u,v),,drop=FALSE])^2)
    ss <- sums[c(u,v)]+c(a[j]-a[i],a[i]-a[j]); cc <- counts[c(u,v),,drop=FALSE]
    cc[1,q[i]] <- cc[1,q[i]]-1; cc[1,q[j]] <- cc[1,q[j]]+1
    cc[2,q[j]] <- cc[2,q[j]]-1; cc[2,q[i]] <- cc[2,q[i]]+1
    cost <- sum(ss^2)+sum((cc-target[c(u,v),,drop=FALSE])^2)
    if (cost < old-1e-12) {f[i]<-v;f[j]<-u;sums[c(u,v)]<-ss;counts[c(u,v),]<-cc}
  }
  stopifnot(identical(unname(table(joint,original)),unname(table(joint,f))),identical(tabulate(f,k),cap))
  as.integer(f)
}

make_outer_folds <- function(dat, cfg) {
  if (!is.null(cfg$inputs$folds)) {
    z <- read_table(cfg$inputs$folds,c('sample_id','repeat','fold'))
    z[['repeat']] <- numeric_column(z[['repeat']], 'repeat'); z$fold <- numeric_column(z$fold, 'fold')
    if (anyDuplicated(z[c('sample_id','repeat')]) || !setequal(z$sample_id,dat$strata$SUBJID) ||
        nrow(z)!=nrow(dat$strata)*cfg$repeats || !setequal(z[['repeat']],seq_len(cfg$repeats))) stop('Supplied outer assignments must cover each baseline participant exactly once per repetition')
    folds <- vapply(seq_len(cfg$repeats),function(r) {
      rr <- z[z[['repeat']]==r,]; f <- rr$fold[match(dat$strata$SUBJID,rr$sample_id)]
      if (anyNA(f) || !setequal(f,seq_len(cfg$outer_folds)) || any(f!=as.integer(f))) stop('Invalid supplied fold numbers')
      as.integer(f)
    }, integer(nrow(dat$strata)))
  } else {
    message('Generating new age-balanced outer assignments. Exact manuscript splits require the original supplied plan.')
    folds <- vapply(seq_len(cfg$repeats),function(r) {
      f <- ir_folds(dat$strata,cfg$outer_folds,cfg$seed+r)
      balance_outer_age(f,dat$strata,cfg$outer_folds,cfg$seed+800000L+r,cfg$age_balance_swaps)
    },integer(nrow(dat$strata)))
  }
  rownames(folds) <- dat$strata$SUBJID; colnames(folds) <- paste0('repeat_',seq_len(cfg$repeats))
  folds
}

fit_stack_weights <- function(train) {
  tolerance <- sqrt(.Machine$double.eps)*max(1,sd(train$target))
  dm <- sd(train$MetIR)<tolerance; dp <- sd(train$ProtIR)<tolerance
  if (dm) train$MetIR <- rep(0,nrow(train))
  if (dp) train$ProtIR <- rep(0,nrow(train))
  fit <- lm.fit(cbind(intercept=1,MetIR=train$MetIR,ProtIR=train$ProtIR),train$target)
  b <- fit$coefficients; b[is.na(b)] <- 0
  data.frame(intercept=unname(b[1]),beta_met=unname(b[2]),beta_prot=unname(b[3]),
    design_rank=fit$rank,n_meta_train=nrow(train),degenerate_met_score=dm,degenerate_prot_score=dp)
}

model_coefficients <- function(model, score, repetition, fold) {
  ix <- which(model$beta[-1,,drop=FALSE]!=0,arr.ind=TRUE)
  if (!nrow(ix)) return(data.frame(id=character(),rule=character(),beta_standardized_x=numeric(),
    beta_input_scale=numeric(),score=character(),variant=character(),repetition=integer(),fold=integer()))
  z <- data.frame(id=model$prep$ids[ix[,1]],rule=colnames(model$beta)[ix[,2]],
    beta_standardized_x=model$beta[-1,,drop=FALSE][ix]*model$prep$y_sd,
    score=score,variant='direct',repetition=repetition,fold=fold)
  z$beta_input_scale <- z$beta_standardized_x/model$prep$scale[match(z$id,model$prep$ids)]
  z
}

fit_outer_job <- function(r, f, dat, folds, cfg) {
  stem <- sprintf('direct_r%02d_f%02d.rds',r,f)
  out <- file.path(cfg$output_dir,'models'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
  final <- file.path(out,paste0('PhenoIR_',stem))
  if (file.exists(final)) return(final)
  tr <- dat$baseline[folds[,r]!=f]
  test_ids <- rownames(folds)[folds[,r]==f]
  te <- which(dat$data$SUBJID %in% test_ids)
  b <- dat$strata[match(dat$data$SUBJID[tr],dat$strata$SUBJID),,drop=FALSE]
  cf <- ir_folds(b,cfg$calibration_folds,cfg$seed+200000L+1000L*r+f)
  bases <- cal <- list()
  for (s in c('MetIR','ProtIR')) {
    path <- file.path(out,paste0(s,'_',stem))
    if (file.exists(path)) bases[[s]] <- readRDS(path) else {
      tuning <- ir_tune(dat,tr,dat$sets[[s]],'direct',cfg,cfg$seed+1000L*r+f)
      m <- ir_fit(dat,tr,dat$sets[[s]],'direct',cfg,tuning)
      pred <- ir_predict(m,dat,te); pred$score<-s;pred$variant<-'direct';pred$repetition<-r;pred$fold<-f
      bases[[s]] <- list(model=m,predictions=pred,tuning=tuning$selected,curve=tuning$curve,
        inner_folds=tuning$folds,coefficients=model_coefficients(m,s,r,f))
      saveRDS(bases[[s]],path)
    }
    parts <- lapply(seq_len(cfg$calibration_folds),function(k) {
      a<-tr[cf!=k];v<-tr[cf==k]
      cp<-file.path(cfg$output_dir,'crossfit',sprintf('%s_direct_r%02d_f%02d_c%02d.rds',s,r,f,k))
      if (file.exists(cp)) return(readRDS(cp)$predictions)
      tune<-ir_tune(dat,a,dat$sets[[s]],'direct',cfg,cfg$seed+300000L+10000L*r+100L*f+k)
      m<-ir_fit(dat,a,dat$sets[[s]],'direct',cfg,tune)
      p<-ir_predict(m,dat,v);p$calibration_fold<-k
      dir.create(dirname(cp),recursive=TRUE,showWarnings=FALSE)
      saveRDS(list(predictions=p,training_ids=dat$data$SUBJID[a],validation_ids=dat$data$SUBJID[v],selected=tune$selected),cp)
      p
    })
    cal[[s]] <- dplyr::bind_rows(parts)
  }
  weights <- preds <- eff <- list()
  for (rule in cfg$rules) {
    a<-cal$MetIR[cal$MetIR$rule==rule,];b<-cal$ProtIR[cal$ProtIR$rule==rule,]
    ix<-match(a$SUBJID,b$SUBJID)
    stopifnot(!anyNA(ix),!anyDuplicated(a$SUBJID),setequal(a$SUBJID,dat$data$SUBJID[tr]),
      !any(a$SUBJID%in%test_ids),max(abs(a$target-b$target[ix]))<1e-10)
    w<-fit_stack_weights(data.frame(MetIR=a$score_value,ProtIR=b$score_value[ix],target=a$target))
    w$rule<-rule;w$repetition<-r;w$fold<-f;w$score<-'PhenoIR';weights[[rule]]<-w
    a<-bases$MetIR$predictions; a<-a[a$rule==rule,]
    b<-bases$ProtIR$predictions; b<-b[b$rule==rule,]
    ix<-match(paste(a$SUBJID,a$Week),paste(b$SUBJID,b$Week));stopifnot(!anyNA(ix))
    p<-a;p$score<-'PhenoIR'
    p$score_value<-w$intercept+w$beta_met*a$score_value+w$beta_prot*b$score_value[ix]
    p$predicted_log_homa<-p$score_value;preds[[rule]]<-p
    eff[[rule]]<-dplyr::bind_rows(lapply(c('MetIR','ProtIR'),function(s) {
      z<-bases[[s]]$coefficients;z<-z[z$rule==rule,]
      z$base_score<-rep(s,nrow(z));z$base_beta_standardized_x<-z$beta_standardized_x
      weight<-if(s=='MetIR')w$beta_met else w$beta_prot
      z$meta_weight<-rep(weight,nrow(z));z$beta_standardized_x<-z$beta_standardized_x*weight
      z$beta_input_scale<-z$beta_input_scale*weight;z$score<-rep('PhenoIR',nrow(z))
      z[z$beta_standardized_x!=0,]
    }))
  }
  saveRDS(list(model=list(type='linear_score_stack',training_ids=dat$data$SUBJID[tr],weights=dplyr::bind_rows(weights)),
    predictions=dplyr::bind_rows(preds),coefficients=dplyr::bind_rows(eff),
    audit=data.frame(repetition=r,fold=f,n_train=length(tr),n_test_W0=sum(dat$data$Week[te]==0),
      n_test_W24=sum(dat$data$Week[te]==24),participant_overlap=0L)),final)
  final
}

run_model_pipeline <- function(cfg) {
  out<-cfg$output_dir;dir.create(out,recursive=TRUE,showWarnings=FALSE)
  paths<-unlist(cfg$inputs);hashes<-tools::md5sum(paths)
  code_hashes<-tools::md5sum(cfg$code_files)
  signature<-list(config=cfg,input_md5=hashes,code_md5=code_hashes)
  state<-file.path(out,'run_signature.rds')
  if(file.exists(state) && !identical(readRDS(state),signature)) stop('Configuration or input content changed. Use a new output_dir; old checkpoints cannot be reused.')
  dat<-load_model_data(cfg);folds<-make_outer_folds(dat,cfg)
  saveRDS(signature,state);saveRDS(dat,file.path(out,'analysis_data.rds'));saveRDS(folds,file.path(out,'outer_folds.rds'))
  write.csv(dat$annotations,file.path(out,'feature_annotations.csv'),row.names=FALSE)
  write.csv(dat$phenotype,file.path(out,'phenotype.csv'),row.names=FALSE)
  write.csv(dat$audit,file.path(out,'eligibility_audit.csv'),row.names=FALSE)
  write.csv(dat$excluded,file.path(out,'excluded_features.csv'),row.names=FALSE)
  fold_table<-dplyr::bind_rows(lapply(seq_len(cfg$repeats),function(r)data.frame(sample_id=rownames(folds),repetition=r,fold=folds[,r])))
  write.csv(fold_table,file.path(out,'outer_fold_assignments.csv'),row.names=FALSE)
  manifest<-list(config=cfg,input_md5=as.list(hashes),code_md5=as.list(code_hashes),n_baseline=length(dat$baseline),n_week24=sum(dat$data$Week==24),
    split_source=if(is.null(cfg$inputs$folds))'newly generated' else 'supplied',R=R.version.string,
    packages=as.list(setNames(vapply(required_packages,function(p)as.character(utils::packageVersion(p)),character(1)),required_packages)),
    frozen_followup=TRUE,genetic_PCs_used=FALSE,model_variant='direct; no molecular or outcome residualization',
    glmnet_call='Historical control=list(thresh,maxit) preserved for numerical compatibility; these values are not claimed to override glmnet top-level defaults.')
  jsonlite::write_json(manifest,file.path(out,'run_manifest.json'),pretty=TRUE,auto_unbox=TRUE,null='null')
  jobs<-expand.grid(fold=seq_len(cfg$outer_folds),repetition=seq_len(cfg$repeats))
  # Job-level seeds make ordering/worker count irrelevant. Windows uses the sequential fallback.
  workers<-if(.Platform$OS.type=='windows')1L else cfg$workers
  for(start in seq(1,nrow(jobs),by=workers)) {
    batch<-seq.int(start,min(start+workers-1L,nrow(jobs)))
    ans<-parallel::mclapply(batch,function(i)fit_outer_job(jobs$repetition[i],jobs$fold[i],dat,folds,cfg),mc.cores=workers,mc.set.seed=FALSE)
    if(any(vapply(ans,inherits,logical(1),'try-error'))) stop('One or more outer jobs failed')
    message('Completed outer jobs ',max(batch),' / ',nrow(jobs))
  }
  bundles<-lapply(list.files(file.path(out,'models'),pattern='^(MetIR|ProtIR|PhenoIR)_direct_r[0-9]+_f[0-9]+[.]rds$',full.names=TRUE),readRDS)
  p<-dplyr::bind_rows(lapply(bundles,`[[`,'predictions'))
  if(anyDuplicated(p[c('SUBJID','Week','score','rule','repetition')])) stop('Duplicate held-out predictions')
  avg<-p|>dplyr::group_by(SUBJID,Week,score,rule)|>dplyr::summarise(predicted_log_homa=mean(predicted_log_homa),observed_log_homa=dplyr::first(observed_log_homa),n_repeats=dplyr::n(),.groups='drop')
  stopifnot(all(avg$n_repeats==cfg$repeats))
  names(avg)[names(avg)=='SUBJID']<-'sample_id';names(avg)[names(avg)=='Week']<-'week'
  write.csv(avg,file.path(out,'participant_predictions.csv'),row.names=FALSE)
  write.csv(p,file.path(out,'predictions_by_repeat.csv'),row.names=FALSE)
  write.csv(dplyr::bind_rows(lapply(bundles,`[[`,'coefficients')),file.path(out,'effective_coefficients.csv'),row.names=FALSE)
  metrics<-p|>dplyr::group_by(score,rule,Week,repetition)|>dplyr::group_modify(function(d,k)ir_metrics(d$observed_log_homa,d$predicted_log_homa))|>dplyr::ungroup()
  write.csv(metrics,file.path(out,'performance_by_repeat.csv'),row.names=FALSE)
  writeLines(c('Completed all requested outer models.',paste('Baseline:',length(dat$baseline),'Week 24:',sum(dat$data$Week==24)),
    'Averaged predictions are on the natural-log HOMA-IR scale.','Participant-keyed files are local analysis outputs; do not commit private study data.'),file.path(out,'COMPLETE.txt'))
  invisible(out)
}

ir_folds <- function (strata, k, seed) 
{
    set.seed(seed)
    n <- nrow(strata)
    stopifnot(n >= k, !anyDuplicated(strata$SUBJID))
    cols <- c("Trial", "IR", "Sex", "BMI_category")
    labels <- lapply(strata[cols], as.character)
    joint <- do.call(paste, c(labels, sep = "|"))
    groups <- split(seq_len(n), joint)
    order_g <- order(lengths(groups), runif(length(groups)))
    capacity <- rep(n%/%k, k)
    if (n%%k) 
        capacity[sample.int(k, n%%k)] <- n%/%k + 1L
    total <- integer(k)
    counts <- lapply(labels, function(z) matrix(0L, length(unique(z)), k, dimnames = list(unique(z), NULL)))
    joint_counts <- matrix(0L, length(groups), k, dimnames = list(names(groups), NULL))
    result <- integer(n)
    for (g in names(groups)[order_g]) {
        indices <- groups[[g]]
        indices <- indices[sample.int(length(indices))]
        for (i in indices) {
            candidates <- which(total < capacity)
            cost <- 12 * joint_counts[g, candidates] + 2 * total[candidates]/pmax(capacity[candidates], 1)
            for (v in seq_along(labels)) {
                z <- labels[[v]][i]
                expected <- sum(labels[[v]] == z)/k
                cost <- cost + counts[[v]][z, candidates]/max(expected, 1)
            }
            best <- candidates[cost <= min(cost) + 1e-12]
            f <- best[sample.int(length(best), 1)]
            result[i] <- f
            total[f] <- total[f] + 1L
            joint_counts[g, f] <- joint_counts[g, f] + 1L
            for (v in seq_along(labels)) counts[[v]][labels[[v]][i], f] <- counts[[v]][labels[[v]][i], f] + 1L
        }
    }
    stopifnot(all(result %in% seq_len(k)), max(tabulate(result, k)) - min(tabulate(result, k)) <= 1L)
    result
}

ir_prepare <- function (x, y, cmat, variant, cfg) 
{
    stopifnot(all(is.finite(y)), all(is.finite(cmat)))
    qr_c <- qr(cmat)
    if (qr_c$rank != ncol(cmat)) 
        stop("Training covariate design is rank deficient")
    y_beta <- qr.coef(qr_c, y)
    y_cov <- as.numeric(cmat %*% y_beta)
    target <- if (variant == "adjusted") 
        y - y_cov
    else y
    keep <- colSums(is.finite(x)) >= cfg$minimum_observations
    x <- x[, keep, drop = FALSE]
    if (ncol(x) < 2L) 
        stop("Fewer than two observed training features")
    medians <- matrixStats::colMedians(x, na.rm = TRUE)
    names(medians) <- colnames(x)
    missing <- which(!is.finite(x), arr.ind = TRUE)
    if (nrow(missing)) 
        x[missing] <- medians[missing[, 2]]
    x_beta <- if (variant == "adjusted") 
        qr.coef(qr_c, x)
    else NULL
    if (variant == "adjusted") 
        x <- x - cmat %*% x_beta
    center <- colMeans(x)
    centered <- sweep(x, 2, center, "-")
    scales <- sqrt(colSums(centered^2)/(nrow(x) - 1))
    varying <- is.finite(scales) & scales > 1e-10
    if (sum(varying) < 2L) 
        stop("Fewer than two variable training features")
    ids <- colnames(x)[varying]
    y_mean <- mean(target)
    y_sd <- sd(target)
    stopifnot(is.finite(y_sd), y_sd > 1e-10)
    list(ids = ids, medians = medians[varying], center = center[varying], scale = scales[varying], x_beta = if (!is.null(x_beta)) x_beta[, 
        varying, drop = FALSE] else NULL, y_beta = y_beta, y_mean = y_mean, y_sd = y_sd, variant = variant, train_mean_y = mean(y), 
        train_n = length(y), x = sweep(centered[, varying, drop = FALSE], 2, scales[varying], "/"), target = target, 
        y = (target - y_mean)/y_sd)
}

ir_transform <- function (prep, x, cmat) 
{
    stopifnot(all(prep$ids %in% colnames(x)), all(is.finite(cmat)))
    x <- x[, prep$ids, drop = FALSE]
    missing <- which(!is.finite(x), arr.ind = TRUE)
    if (nrow(missing)) 
        x[missing] <- prep$medians[missing[, 2]]
    if (prep$variant == "adjusted") 
        x <- x - cmat %*% prep$x_beta
    x <- sweep(sweep(x, 2, prep$center, "-"), 2, prep$scale, "/")
    stopifnot(all(is.finite(x)))
    list(x = x, y_cov = as.numeric(cmat %*% prep$y_beta))
}

ir_path <- function (prep, alpha, cfg) 
{
    fit <- glmnet::glmnet(x = prep$x, y = prep$y, family = "gaussian", alpha = alpha, lambda = cfg$lambdas, intercept = TRUE, 
        standardize = FALSE, control = list(thresh = cfg$glmnet_thresh, maxit = cfg$glmnet_maxit))
    if (fit$jerr != 0L || length(fit$lambda) != length(cfg$lambdas)) {
        stop("Incomplete or nonconverged elastic-net lambda path: ", fit$jerr)
    }
    fit
}

ir_tune <- function (dat, rows, feature_ids, variant, cfg, seed) 
{
    b <- dat$strata[match(dat$data$SUBJID[rows], dat$strata$SUBJID), ]
    stopifnot(!anyNA(b$SUBJID))
    folds <- ir_folds(b, cfg$inner_folds, seed)
    errors <- array(NA_real_, c(cfg$inner_folds, length(cfg$alphas), length(cfg$lambdas)))
    for (f in seq_len(cfg$inner_folds)) {
        tr <- rows[folds != f]
        va <- rows[folds == f]
        prep <- ir_prepare(dat$x[tr, feature_ids, drop = FALSE], dat$y[tr], dat$c[tr, , drop = FALSE], variant, 
            cfg)
        test <- ir_transform(prep, dat$x[va, feature_ids, drop = FALSE], dat$c[va, , drop = FALSE])
        target <- if (variant == "adjusted") 
            dat$y[va] - test$y_cov
        else dat$y[va]
        for (a in seq_along(cfg$alphas)) {
            fit <- ir_path(prep, cfg$alphas[a], cfg)
            pred <- as.matrix(predict(fit, newx = test$x, s = cfg$lambdas)) * prep$y_sd + prep$y_mean
            errors[f, a, ] <- colMeans((pred - target)^2)
        }
    }
    mse <- apply(errors, c(2, 3), mean)
    se <- apply(errors, c(2, 3), sd)/sqrt(cfg$inner_folds)
    best <- which(mse == min(mse), arr.ind = TRUE)[1, ]
    ai <- best[1]
    li <- best[2]
    eligible <- which(mse[ai, ] <= mse[ai, li] + se[ai, li])
    one_se <- eligible[which.max(cfg$lambdas[eligible])]
    selected <- data.frame(rule = cfg$rules, alpha = cfg$alphas[ai], lambda = cfg$lambdas[c(one_se, li)], inner_mse = mse[ai, 
        c(one_se, li)], inner_se = se[ai, c(one_se, li)], lambda_index = c(one_se, li), lambda_boundary = c(one_se, 
        li) %in% c(1L, length(cfg$lambdas)), minimum_at_lower_boundary = li == length(cfg$lambdas))
    curve <- expand.grid(alpha = cfg$alphas, lambda = cfg$lambdas)
    curve$mse <- as.vector(mse)
    curve$se <- as.vector(se)
    list(selected = selected, curve = curve, folds = data.frame(SUBJID = b$SUBJID, fold = folds))
}

ir_fit <- function (dat, rows, feature_ids, variant, cfg, tuning) 
{
    prep <- ir_prepare(dat$x[rows, feature_ids, drop = FALSE], dat$y[rows], dat$c[rows, , drop = FALSE], variant, 
        cfg)
    fit <- ir_path(prep, tuning$selected$alpha[1], cfg)
    beta <- as.matrix(coef(fit, s = tuning$selected$lambda))
    colnames(beta) <- tuning$selected$rule
    prep$x <- prep$target <- prep$y <- NULL
    list(prep = prep, beta = beta, selected = tuning$selected, training_ids = dat$data$SUBJID[rows])
}

ir_predict <- function (model, dat, rows) 
{
    stopifnot(!any(dat$data$SUBJID[rows] %in% model$training_ids))
    transformed <- ir_transform(model$prep, dat$x[rows, , drop = FALSE], dat$c[rows, , drop = FALSE])
    pred <- cbind(1, transformed$x) %*% model$beta
    pred <- pred * model$prep$y_sd + model$prep$y_mean
    target <- if (model$prep$variant == "adjusted") 
        dat$y[rows] - transformed$y_cov
    else dat$y[rows]
    dplyr::bind_rows(lapply(seq_len(ncol(pred)), function(j) {
        data.frame(SUBJID = dat$data$SUBJID[rows], Week = dat$data$Week[rows], Trial = dat$data$Trial[rows], Arm = dat$data$Arm[rows], 
            rule = colnames(pred)[j], observed_log_homa = dat$y[rows], target = target, score_value = pred[, j], 
            predicted_log_homa = if (model$prep$variant == "adjusted") 
                pred[, j] + transformed$y_cov
            else pred[, j], covariate_prediction = transformed$y_cov, training_mean = model$prep$train_mean_y, 
            zero_residual_reference = 0)
    }))
}

ir_metrics <- function (y, pred) 
{
    ok <- is.finite(y) & is.finite(pred)
    y <- y[ok]
    pred <- pred[ok]
    n <- length(y)
    sst <- sum((y - mean(y))^2)
    varied <- n > 2 && sd(pred) > 1e-10 && sd(y) > 1e-10
    cal <- if (varied) 
        coef(lm(y ~ pred))
    else c(NA_real_, NA_real_)
    data.frame(n = n, r2 = if (sst > 0) 
        1 - sum((y - pred)^2)/sst
    else NA_real_, rmse = sqrt(mean((y - pred)^2)), mae = mean(abs(y - pred)), pearson = if (varied) 
        cor(y, pred)
    else NA_real_, spearman = if (varied) 
        cor(y, pred, method = "spearman")
    else NA_real_, calibration_intercept = unname(cal[1]), calibration_slope = unname(cal[2]))
}

