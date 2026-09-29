#!/usr/bin/env Rscript
# Fit five-year clinical progression from each selected CR MRI.
.this_file <- function() {
  frames <- vapply(sys.frames(),function(x) if(is.null(x$ofile)) "" else as.character(x$ofile),character(1))
  files <- frames[basename(frames)=='figure5_progression.R']
  if(length(files))return(normalizePath(tail(files,1L)))
  normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L]))
}
source(file.path(dirname(.this_file()),"common.R"))
f5p_endpoints <- c('mci_to_dementia','any_worsening')
f5p_identity <- function(x) {
  # Hash ordered values and factor levels, excluding runtime-only table pointers.
  if(is.data.table(x)) return(lapply(as.list(x),identity))
  if(is.list(x)) return(lapply(x,f5p_identity))
  x
}
f5p_hash <- function(x) digest::digest(f5p_identity(x),algo='sha256')
f5p_date <- function(x) as.Date(as.character(x))
f5p_finite <- function(d,cols) {
  answer <- complete.cases(d[,..cols])
  for(nm in cols) {
    x <- d[[nm]]
    if(is.numeric(x)) answer <- answer & is.finite(x) else
      answer <- answer & !is.na(text(x))
  }
  answer
}
f5p_event <- function(endpoint,baseline,future) {
  if(identical(endpoint,'mci_to_dementia')) return(future %in% c(3,4))
  if(identical(endpoint,'any_worsening')) return((baseline==1 & future %in% c(2,3,4)) | (baseline==2 & future %in% c(3,4)))
  stop('Unsupported progression endpoint')
}
f5p_followup <- function(history,origin,baseline_code,endpoint,guard_days=30) {
  # Follow-up starts at the selected CR MRI; missing follow-up is not a non-event.
  ans <- list(event_date=as.Date(NA),censor_date=as.Date(NA),event_dx_num=NA_real_,
    n_future_dx=0L,event=NA_integer_,time_to_event_days=NA_real_,time_to_event_years=NA_real_,
    followup_days=NA_real_,event_5y=NA_integer_,landmark_supported=FALSE,
    outcome_reason='no_future_diagnosis_after_guard')
  if(is.null(history)||!nrow(history)||is.na(origin)) return(ans)
  h <- copy(as.data.table(history)); need(h,c('diagnosis_date','diagnosis_code'))
  h[,diagnosis_date:=f5p_date(diagnosis_date)];h[,diagnosis_code:=number(diagnosis_code)]
  h <- h[!is.na(diagnosis_date)&diagnosis_date>origin+guard_days&is.finite(diagnosis_code)]
  if(!nrow(h)) return(ans)
  setorder(h,diagnosis_date)
  ev <- h[f5p_event(endpoint,baseline_code,diagnosis_code)]
  ans$censor_date <- max(h$diagnosis_date);ans$n_future_dx <- nrow(h)
  ans$event <- as.integer(nrow(ev)>0L)
  if(nrow(ev)) { ans$event_date <- ev$diagnosis_date[[1L]]; ans$event_dx_num <- ev$diagnosis_code[[1L]] }
  last <- if(ans$event==1L) ans$event_date else ans$censor_date
  ans$time_to_event_days <- as.numeric(last-origin)
  ans$time_to_event_years <- ans$time_to_event_days/365.25
  ans$followup_days <- as.numeric(ans$censor_date-origin)
  ans$landmark_supported <- is.finite(ans$time_to_event_days) && ans$time_to_event_days>0 &&
    (ans$event==1L || ans$time_to_event_days>=183)
  if(ans$event==1L&&ans$time_to_event_years<=5) ans$event_5y <- 1L else
    if(ans$time_to_event_years>=5) ans$event_5y <- 0L
  ans$outcome_reason <- if(!ans$landmark_supported) 'non_event_followup_below_183_days' else
    if(is.na(ans$event_5y)) 'five_year_outcome_unobserved' else 'observed'
  ans
}
f5p_normalize <- function(candidates) {
  d <- copy(as.data.table(candidates))
  need(d,c('rid','endpoint','scan_id','scan_set_id','index_id','index_date','mri_date',
    'diagnosis_mri','index_diagnosis_code','age_index','sex','education_years','apoe4_count',
    'field_strength','site_id','LC_z','SNVTA_z','LC_adjusted','SNVTA_adjusted'))
  for(nm in c('LC_z','SNVTA_z','LC_adjusted','SNVTA_adjusted','age_index',
              'education_years','apoe4_count','index_diagnosis_code')) set(d,j=nm,value=number(d[[nm]]))
  for(nm in c('index_date','mri_date')) set(d,j=nm,value=f5p_date(d[[nm]]))
  if(any(is.na(d$index_date))||any(d$index_date!=d$mri_date)) stop('Progression index must be the actual MRI date')
  if(any(!d$endpoint %in% f5p_endpoints)) stop('Unexpected progression endpoint')
  d[,`:=`(RID=as.character(rid),outcome=endpoint,baseline_dx_num=index_diagnosis_code,
    baseline_dx_label=as.character(diagnosis_mri),age_at_scan=age_index,
    education_match=education_years,APOE4=apoe4_count)]
  if(any(!((d$baseline_dx_num==1 & d$baseline_dx_label=='CN') |
           (d$baseline_dx_num==2 & d$baseline_dx_label=='MCI')))) stop('MRI diagnosis coding mismatch')
  if(any(d$outcome=='mci_to_dementia' & d$baseline_dx_num!=2)) stop('MCI endpoint contains a non-MCI MRI')
  if(any(!is.finite(d$age_at_scan)|d$age_at_scan<40|d$age_at_scan>110))
    stop('Missing or implausible true MRI age; original age support must be resolved')
  site_raw<-text(d$site_id)
  if(anyNA(site_raw)||any(!grepl('^[0-9]+$',site_raw)))
    stop('Original progression site codes must be nonmissing integer labels')
  site_original<-sub('^0+(?=[0-9])','',site_raw,perl=TRUE)
  site_map<-unique(data.table(prepared=site_raw,original=site_original))
  if(anyDuplicated(site_map$original)) stop('Original site label restoration would merge prepared site partitions')
  d[,`:=`(prepared_site_id=site_raw,site_id=site_original)]
  unique_keys(d,c('endpoint','scan_set_id'),'progression candidates')
  setorder(d,outcome,RID,index_date,scan_id,scan_set_id)
  d
}
f5p_center_landmarks <- function(d,population) {
  d <- copy(d)
  mu_age <- mean(d$age_at_scan,na.rm=TRUE);mu_edu <- mean(d$education_match,na.rm=TRUE)
  d[,`:=`(age_at_scan_c=age_at_scan-mu_age,education_c=education_match-mu_edu,
    age_center_reference=mu_age,education_center_reference=mu_edu,
    covariate_center_population=population)]
  list(data=d,scales=data.table(variable=c('age_at_scan','education_match'),center=c(mu_age,mu_edu),
    scale=1,n_reference=c(sum(is.finite(d$age_at_scan)),sum(is.finite(d$education_match))),
    reference_population=population,operation='mean_center_no_rescale'))
}
f5p_prepare_frames <- function(candidates,diagnoses) {
  d <- f5p_normalize(candidates)
  dx <- copy(as.data.table(diagnoses));need(dx,c('rid','diagnosis_date','diagnosis_code'))
  dx[,`:=`(rid=as.character(rid),diagnosis_date=f5p_date(diagnosis_date),diagnosis_code=number(diagnosis_code))]
  histories <- split(dx,by='rid',keep.by=TRUE)
  selections <- list();audit <- list();j <- 0L
  for(endpoint_value in f5p_endpoints) for(context in c('LC','SNVTA')) {
    eligible <- paste0(context,c('_z','_adjusted'))
    z <- d[outcome==endpoint_value]
    z <- z[f5p_finite(z,eligible)]
    setorder(z,RID,index_date,scan_id,scan_set_id)
    selected <- unique(z,by='RID') # BEFORE access to event/follow-up/CC values.
    j <- j+1L
    selected[,`:=`(selection_context=context,nucleus=context,region=context,
      landmark_strategy='first_scan_no_future')]
    selected[,`:=`(cr_resid_tech_z=get(paste0(context,'_z')),
      cr_pred_tech_adjusted=get(paste0(context,'_adjusted')))]
    if(nrow(selected)) {
      follow <- rbindlist(lapply(seq_len(nrow(selected)),function(i)
        as.data.table(f5p_followup(histories[[selected$RID[[i]]]],selected$index_date[[i]],
          selected$baseline_dx_num[[i]],endpoint_value))),fill=TRUE)
      # Compare the rebuilt clinical outcome with the prepared history fields.
      checks <- list()
      for(nm in intersect(c('event_date','event_5y','followup_days'),names(selected))) {
        a <- selected[[nm]];b <- follow[[nm]]
        if(nm=='event_date') {a<-f5p_date(a);b<-f5p_date(b)} else {a<-number(a);b<-number(b)}
        same <- (is.na(a)&is.na(b)) | (!is.na(a)&!is.na(b)&a==b)
        checks[[nm]] <- data.table(endpoint=endpoint_value,context,field=nm,n_rows=nrow(selected),
          n_equal=sum(same),n_different=sum(!same))
        selected[, (paste0('prepared_',nm)):=as.character(get(nm))]
      }
      selected[,intersect(names(follow),names(selected)):=NULL]
      selected <- cbind(selected,follow)
      audit[[j]] <- rbindlist(checks,fill=TRUE)
    } else {
      template <- as.data.table(f5p_followup(NULL,as.Date(NA),NA,'any_worsening'))[0]
      selected <- cbind(selected,template)
    }
    selections[[j]] <- selected
  }
  selected <- rbindlist(selections,fill=TRUE)
  primary <- selected[landmark_supported %in% TRUE]
  # The four progression pools contribute their original participant weights.
  primary <- f5p_center_landmarks(primary,'all_four_primary_landmark_pools_after_183d_nonevent_support_before_5y_or_CC')
  obj<-primary$data
  obj[,`:=`(sex=factor(text(sex)),field_strength_f=factor(text(field_strength)),
    baseline_dx_f=factor(baseline_dx_label,levels=c('CN','MCI')),site_id=factor(text(site_id)))]
  setorder(obj,outcome,nucleus,RID)
  list(primary=obj,selection=selected,scales=primary$scales,prepared_outcome_audit=rbindlist(audit,fill=TRUE))
}
f5p_tasks <- function() {
  ans <- rbindlist(lapply(f5p_endpoints,function(ep)data.table(endpoint=ep,
    context='primary_individual_index',nucleus=c('LC','SNVTA'),frame_pool='primary')))
  ans[,scientific_task:=paste(endpoint,context,nucleus,sep='__')]
  ans
}
f5p_covariates <- function(endpoint) c('age_at_scan_c','sex','education_c','APOE4','field_strength_f',
  if(endpoint=='any_worsening') 'baseline_dx_f')
f5p_model_frame <- function(spec,frames) {
  ep <- spec$endpoint[[1L]];target<-spec$nucleus[[1L]]
  d <- copy(frames[[spec$frame_pool[[1L]]]][outcome==ep])
  if(spec$frame_pool[[1L]]=='primary') d<-d[nucleus==target]
  covs <- f5p_covariates(ep)
  predictors_cc <- c('cr_resid_tech_z','cr_pred_tech_adjusted')
  cols <- c('RID','event_5y',predictors_cc,covs,'site_id')
  d <- d[f5p_finite(d,cols)]
  for(nm in c('sex','field_strength_f','baseline_dx_f','site_id')) d[,(nm):=droplevels(get(nm))]
  if(anyDuplicated(d$RID)) stop('Progression model frame has repeated participants')
  setorder(d,RID)
  d
}
f5p_estimable_term <- function(x) {
  if(is.factor(x)||is.character(x)) return(length(unique(na.omit(as.character(x))))>=2L)
  y<-x[is.finite(x)];length(y)>=2L&&is.finite(sd(y))&&sd(y)>0
}
f5p_fit_task <- function(spec,prepared) {
  d<-f5p_model_frame(spec,prepared);ep<-spec$endpoint[[1L]];rg<-spec$nucleus[[1L]]
  predictors<-'cr_resid_tech_z'
  covs<-f5p_covariates(ep);terms<-c(predictors,covs)
  formula_fixed<-paste('event_5y ~',paste(terms,collapse=' + '))
  formula_text<-paste(formula_fixed,'+ (1 | site_id)')
  model_id<-spec$scientific_task[[1L]];n<-nrow(d);nev<-sum(d$event_5y==1L);nno<-sum(d$event_5y==0L)
  base<-cbind(spec,data.table(model_id,n_subjects=n,n_events=nev,n_non_events=nno,
    formula=formula_text,random_effects='(1 | site_id)',estimator='binomial_logit_Laplace_nAGQ1',
    n_sites=nlevels(d$site_id),frame_sha256=f5p_hash(d),p_adjustment='none',
    multiplicity_correction_applied=FALSE,raw_p_is_primary=TRUE))
  unsupported<-function(reason) list(status='not_estimable',results=cbind(base,
    data.table(status='NOT_ESTIMABLE',reason)),frame=d,coefficients=data.table(),
    diagnostics=cbind(base,data.table(status='NOT_ESTIMABLE',reason)))
  if(n<50L||nev<20L||nno<20L) return(unsupported('requires_N50_events20_nonevents20'))
  if(nlevels(d$site_id)<2L) return(unsupported('mandatory_site_random_intercept_requires_two_levels'))
  bad<-terms[!vapply(terms,function(term)f5p_estimable_term(d[[term]]),logical(1L))]
  if(length(bad)) return(unsupported(paste('nonestimable_fixed_terms',paste(bad,collapse=','))))
  X<-model.matrix(as.formula(formula_fixed),d);rank<-qr(X)$rank
  if(rank<ncol(X)) return(unsupported('rank_deficient_fixed_design_no_term_deletion'))
  caught<-capture_fit(lme4::glmer(as.formula(formula_text),data=d,family=binomial(link='logit'),
    control=lme4::glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000L),
      check.rankX='stop.deficient',check.conv.singular='ignore'),nAGQ=1L))
  fit<-caught$fit;b<-lme4::fixef(fit);V<-as.matrix(vcov(fit));se<-sqrt(diag(V))
  conv<-fit@optinfo$conv$lme4$messages;conv_text<-paste(unique(as.character(conv)),collapse=' | ')
  opt<-unlist(fit@optinfo$conv$opt);singular<-lme4::isSingular(fit,tol=1e-4)
  good<-!nzchar(conv_text)&&(!length(opt)||all(opt==0))&&all(is.finite(b))&&all(is.finite(se)&se>0)
  status<-if(!good) 'UNRESOLVED_INFERENCE' else if(singular) 'OK_SINGULAR' else 'OK'
  co<-data.table(term=names(b),beta=unname(b),se=unname(se))
  co[,`:=`(z=beta/se,p_raw=2*pnorm(abs(beta/se),lower.tail=FALSE),
    beta_ci_low=beta-1.96*se,beta_ci_high=beta+1.96*se,
    odds_ratio=exp(beta),ci_low=exp(beta-1.96*se),ci_high=exp(beta+1.96*se),
    inference_usable=good,status=status,effect_scale='log_odds_and_exponentiated_odds_ratio')]
  focal<-co[term %in% predictors]
  focal[,`:=`(focal_nucleus=rg,
    predictor_scale='all_bilateral_technical_residual_rows_within_nucleus_before_diagnosis_or_index_selection')]
  result<-cbind(base[rep(1L,nrow(focal))],focal)
  result[,`:=`(converged=good,singular=singular,reason=conv_text)]
  coefficients<-cbind(base[rep(1L,nrow(co))],co)
  diag<-cbind(base,data.table(status=status,converged=good,singular=singular,
    optimizer='bobyqa',maxfun=100000L,nAGQ=1L,optimizer_code=paste(opt,collapse=','),
    convergence_messages=conv_text,warnings=paste(caught$warnings,collapse=' | '),
    design_rank=rank,design_columns=ncol(X),site_variance=as.numeric(lme4::VarCorr(fit)$site_id[1L,1L]),
    residual_deviance=deviance(fit),logLik=as.numeric(logLik(fit)),AIC=AIC(fit),
    iterations=as.integer(fit@optinfo$feval),fallback_used=FALSE,
    R_version=as.character(getRversion()),lme4_version=as.character(packageVersion('lme4'))))
  d[,model_id:=model_id]
  list(status='complete',results=result,coefficients=coefficients,
    diagnostics=diag,frame=d,fit=fit,covariance=V)
}

f5s_technical_bilateral <- function(tech) {
  need(tech,c('rid','scan_id','scan_set_id','nucleus','hemisphere','cr_residual','raw_stratum_mean'))
  unique_keys(tech,c('scan_set_id','nucleus','hemisphere'),'technical residuals')
  if(any(!is.finite(tech$cr_residual))||any(!is.finite(tech$raw_stratum_mean))) stop('Invalid technical residual/anchor')
  bi<-tech[,.(residual=mean(cr_residual),adjusted=mean(cr_residual+raw_stratum_mean),n_hemi=.N),by=.(rid,scan_id,scan_set_id,nucleus)]
  excluded<-sum(bi$n_hemi!=2L)
  # Require both hemispheres before within-region standardization.
  bi<-bi[n_hemi==2L]
  scales<-bi[,.(center=mean(residual),scale=sd(residual),n_reference=.N),by=nucleus]
  if(any(!is.finite(scales$scale)|scales$scale<=0)) stop('Nonestimable technical nucleus scale')
  bi<-merge(bi,scales,by='nucleus',all.x=TRUE,sort=FALSE);bi[,z:=(residual-center)/scale]
  wide<-dcast(bi,rid+scan_id+scan_set_id~nucleus,value.var=c('z','adjusted','residual'))
  for(n in c('LC','SNVTA')) for(v in c('z','adjusted','residual')) setnames(wide,paste(v,n,sep='_'),paste(n,v,sep='_'))
  scales[,`:=`(reference_population='all_finite_bilateral_technical_rows_within_nucleus_before_DX_index_outcome_filters',adjusted_definition='bilateral_mean(cr_residual + raw_stratum_mean); distinct_from_display_model')]
  list(wide=wide,scales=scales,bilateral=bi,excluded_unilateral=excluded)
}

f5p_prepare <- function() {
  transformed<-f5s_technical_bilateral(read_results("private/shared/technical_cr.csv"))
  candidates<-join_one(input("progression.csv"),transformed$wide,c("rid","scan_id","scan_set_id"),"technical progression linkage")
  frames<-f5p_prepare_frames(candidates,input("diagnoses.csv"))
  if(nrow(frames$prepared_outcome_audit)&&any(frames$prepared_outcome_audit$n_different>0))stop("Selected-index event reconstruction differs from prepared history")
  write_plan("progression",f5p_tasks(),frames)
  atomic_csv(frames$primary,"private/figure5/primary_landmarks.csv")
  atomic_csv(frames$selection,"private/figure5/selected_index_audit.csv")
  publish_results(frames$scales,"figure5/progression_covariate_scales.csv")
  publish_results(transformed$scales,"shared/technical_bilateral_scales.csv")
  record_diagnostics("progression_preparation",frames$prepared_outcome_audit)
}
f5p_finalize <- function() {
  tasks<-collect_tasks("progression")
  for(component in c("results","coefficients"))publish_results(bind_component(tasks,component),paste0("figure5/progression_",component,".csv"))
  record_diagnostics("progression",bind_component(tasks,"diagnostics"))
  atomic_csv(bind_component(tasks,"frame"),"private/figure5/progression_model_frames.csv")
}
if(sys.nframe()==0L){a<-cli();if(a$action=="prepare")f5p_prepare()else if(a$action=="task")run_task("progression",a$task,f5p_fit_task)else if(a$action=="finalize")f5p_finalize()else stop("Expected prepare/task/finalize")}
