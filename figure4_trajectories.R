#!/usr/bin/env Rscript
# Estimate CR, regional-volume and TPV trajectory intercepts and slopes.
.f4_file <- local({x<-grep("^--file=",commandArgs(FALSE),value=TRUE);if(length(x)) sub("^--file=","",x[1]) else "figure4_trajectories.R"})
if (!exists("need",mode="function")) source(file.path(dirname(.f4_file),"common.R"))
min_trajectory_followup_y <- 1
diagnosis_col <- "dx_bl"
diagnosis_stage_levels <- c("CN","MCI","AD")
diagnosis_reference_level <- "CN"
parse_numeric_plain <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  y0 <- trimws(as.character(x))
  y0[y0 %in% c('', '.', 'NA', 'N/A', 'na', 'n/a', 'NULL', 'null', 'NaN', 'nan')] <- NA_character_
  y0 <- gsub(',', '', y0, fixed=TRUE)
  suppressWarnings(as.numeric(y0))
}

as_num_na <- function(x) {
  y <- parse_numeric_plain(x)
  y[y %in% c(-4,-1,-9)] <- NA_real_
  y
}
fit_lme_trajectory_parameters <- function(dt, y_col, cov_terms, random_scan = FALSE, site_random = TRUE, time_col = 'time_c', min_n = 3L, label = 'y') {
  d <- copy(dt)
  if (!'RID' %in% names(d)) stop('trajectory table lacks RID')
  d[, y_traj := as_num_na(get(y_col))]
  d[, time_y := as_num_na(get(time_col))]
  d <- d[is.finite(y_traj) & is.finite(time_y)]
  if (!nrow(d)) return(list(values=data.table(), qc=data.table(label=label, method='lme', status='failed_no_finite_rows')))
  elig <- d[, .(n_scans=uniqueN(scan_id), followup_span_y=max(time_y, na.rm=TRUE)-min(time_y, na.rm=TRUE), n_time=uniqueN(time_y)), by=RID]
  if (min_n <= 1L) {
    elig <- elig[n_scans >= 1L]
  } else {
    elig <- elig[n_scans >= min_n & n_time >= 2L & is.finite(followup_span_y) & followup_span_y >= min_trajectory_followup_y]
  }
  d <- d[RID %in% elig$RID]
  if (uniqueN(d$RID) < 20L) return(list(values=data.table(), qc=data.table(label=label, method='lme', status='failed_too_few_eligible_subjects', n_subjects=uniqueN(d$RID))))

  # Use the first positive finite participant ICV from the complete MUSE source,
  # then standardize within the eligible outcome/model sample.
  icv_scale_n_rows <- NA_integer_; icv_scale_n_subjects <- NA_integer_
  icv_scale_mean_raw <- NA_real_; icv_scale_sd_raw <- NA_real_
  if ('ICV_z' %in% cov_terms) {
    if (!'ICV_raw' %in% names(d)) {
      return(list(values=data.table(), qc=data.table(label=label, method='lme', status='failed_missing_locked_ICV_raw')))
    }
    d[, ICV_raw := as_num_na(ICV_raw)]
    icv_ok <- is.finite(d$ICV_raw) & d$ICV_raw > 0
    icv_scale_n_rows <- sum(icv_ok)
    icv_conflict <- d[icv_ok, .(n_icv_values=uniqueN(ICV_raw)), by=RID][n_icv_values > 1L]
    if (nrow(icv_conflict)) {
      stop('Participant-baseline ICV map has multiple raw values for one or more RIDs inside eligible model sample for ', label)
    }
    icv_subject_scale <- unique(d[icv_ok, .(RID, ICV_raw)], by='RID')
    icv_scale_n_subjects <- nrow(icv_subject_scale)
    icv_scale_mean_raw <- if (icv_scale_n_subjects >= 1L) mean(icv_subject_scale$ICV_raw) else NA_real_
    icv_scale_sd_raw <- if (icv_scale_n_subjects >= 2L) sd(icv_subject_scale$ICV_raw) else NA_real_
    if (!is.finite(icv_scale_sd_raw) || icv_scale_sd_raw <= 0) {
      return(list(values=data.table(), qc=data.table(
        label=label, method='lme', status='failed_nonestimable_ICV_scale',
        icv_rule='first_positive_full_MUSE_source_then_unique_subject_model_sample_zscore',
        icv_scale_n_rows=icv_scale_n_rows, icv_scale_n_subjects=icv_scale_n_subjects,
        icv_scale_mean_raw=icv_scale_mean_raw, icv_scale_sd_raw=icv_scale_sd_raw
      )))
    }
    d[, ICV_z := (ICV_raw - icv_scale_mean_raw) / icv_scale_sd_raw]
  }

  # Site enters through the model-specific random effect, not the fixed covariates.
  covs <- setdiff(intersect(cov_terms, names(d)), 'site_id')
  if (length(covs)) {
    ok <- rep(TRUE, nrow(d))
    for (cc in covs) {
      if (is.numeric(d[[cc]]) || is.integer(d[[cc]])) ok <- ok & is.finite(as_num_na(d[[cc]]))
      else { vv <- as.character(d[[cc]]); ok <- ok & !is.na(vv) & trimws(vv) != '' }
    }
    d <- d[ok]
  }
  if (uniqueN(d$RID) < 20L) return(list(values=data.table(), qc=data.table(label=label, method='lme', status='failed_too_few_subjects_after_cov_complete', n_subjects=uniqueN(d$RID))))

  keep_cov <- character()
  for (cc in covs) {
    if (is.numeric(d[[cc]]) || is.integer(d[[cc]])) {
      vv <- as_num_na(d[[cc]]); d[[cc]] <- vv
      if (length(unique(vv[is.finite(vv)])) >= 2L) keep_cov <- c(keep_cov, cc)
    } else {
      if (is.factor(d[[cc]])) {
        d[[cc]] <- droplevels(d[[cc]])
      } else if (identical(cc, diagnosis_col)) {
        d[[cc]] <- factor(as.character(d[[cc]]), levels = diagnosis_stage_levels)
        d[[cc]] <- droplevels(d[[cc]])
      } else {
        d[[cc]] <- factor(as.character(d[[cc]]))
      }
      if (identical(cc, diagnosis_col) && !diagnosis_reference_level %in% levels(d[[cc]])) {
        stop('Baseline diagnosis reference level is absent from LME complete-case sample for ', label, ': ', diagnosis_reference_level)
      }
      if (nlevels(d[[cc]]) >= 2L) keep_cov <- c(keep_cov, cc)
    }
  }
  d[, RID_factor := factor(RID)]
  n_sites <- NA_integer_
  if (isTRUE(site_random)) {
    if (!'site_id' %in% names(d)) return(list(values=data.table(), qc=data.table(label=label, method='lme', status='failed_missing_mandatory_site_id', site_random=TRUE)))
    site_text <- trimws(as.character(d$site_id))
    if (any(is.na(site_text) | site_text == '' | toupper(site_text) == 'SITE_MISSING')) {
      return(list(values=data.table(), qc=data.table(label=label, method='lme', status='failed_missing_mandatory_site_id', site_random=TRUE)))
    }
    d[, site_id := droplevels(factor(site_text))]
    n_sites <- nlevels(d$site_id)
    if (n_sites < 2L) return(list(values=data.table(), qc=data.table(label=label, method='lme', status='failed_nonestimable_mandatory_site_random_effect', n_sites=n_sites, site_random=TRUE)))
  }
  if (random_scan && 'scan_id' %in% names(d)) d[, scan_id_factor := factor(scan_id)]

  rhs <- c('time_y', keep_cov)
  rand <- '(1 + time_y || RID_factor)'
  if (isTRUE(site_random)) rand <- paste(rand, '+ (1 | site_id)')
  if (isTRUE(random_scan) && 'scan_id_factor' %in% names(d) && nlevels(d$scan_id_factor) >= 2L) rand <- paste(rand, '+ (1 | scan_id_factor)')
  form <- as.formula(paste('y_traj ~', paste(rhs, collapse=' + '), '+', rand))
  fit_warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(lme4::lmer(form, data=d, REML=TRUE,
               control=lme4::lmerControl(check.conv.singular='ignore', optimizer='bobyqa')),
      warning=function(w) { fit_warnings <<- c(fit_warnings, conditionMessage(w)); invokeRestart('muffleWarning') }),
    error=function(e) e
  )
  if (inherits(fit, 'error')) {
    return(list(
      values=data.table(),
      qc=data.table(
        label=label, method='lme', status='failed_lmer_error', n_subjects=0L, n_rows=nrow(d),
        singular=NA, formula=paste(deparse(form), collapse=' '), covariates_fixed=paste(keep_cov, collapse='|'),
        site_random=isTRUE(site_random), n_sites=n_sites, lme4_available=TRUE,
        lme_message=conditionMessage(fit), ols_fallback_used=FALSE,
        icv_rule=if ('ICV_z' %in% keep_cov) 'first_positive_full_MUSE_source_then_unique_subject_model_sample_zscore' else 'not_used',
        icv_scale_n_rows=icv_scale_n_rows, icv_scale_n_subjects=icv_scale_n_subjects,
        icv_scale_mean_raw=icv_scale_mean_raw, icv_scale_sd_raw=icv_scale_sd_raw
      )
    ))
  }
  extracted <- tryCatch({
    fe <- lme4::fixef(fit)
    re <- lme4::ranef(fit)[['RID_factor']]
    if (!all(c('(Intercept)','time_y') %in% names(fe))) stop('required fixed intercept/time terms are absent')
    if (any(!is.finite(as.numeric(fe[c('(Intercept)','time_y')])))) stop('required fixed intercept/time terms are nonfinite')
    if (is.null(re) || !nrow(re)) stop('participant random effects are absent')
    ids <- rownames(re)
    ri <- if ('(Intercept)' %in% colnames(re)) re[,'(Intercept)'] else rep(0, length(ids))
    rs <- if ('time_y' %in% colnames(re)) re[,'time_y'] else rep(0, length(ids))
    data.table(
      RID=ids,
      intercept=as.numeric(fe['(Intercept)'] + ri),
      slope=as.numeric(fe['time_y'] + rs),
      n=as.integer(elig[match(ids, RID), n_scans]),
      n_time=as.integer(elig[match(ids, RID), n_time]),
      followup_span_y=as.numeric(elig[match(ids, RID), followup_span_y]),
      method='lme'
    )
  }, error=function(e) e)
  if (inherits(extracted, 'error')) {
    return(list(
      values=data.table(),
      qc=data.table(
        label=label, method='lme', status='failed_lme_parameter_extraction', n_subjects=0L, n_rows=nrow(d),
        singular=lme4::isSingular(fit), formula=paste(deparse(form), collapse=' '), covariates_fixed=paste(keep_cov, collapse='|'),
        site_random=isTRUE(site_random), n_sites=n_sites, lme4_available=TRUE,
        lme_message=conditionMessage(extracted), ols_fallback_used=FALSE,
        icv_rule=if ('ICV_z' %in% keep_cov) 'first_positive_full_MUSE_source_then_unique_subject_model_sample_zscore' else 'not_used',
        icv_scale_n_rows=icv_scale_n_rows, icv_scale_n_subjects=icv_scale_n_subjects,
        icv_scale_mean_raw=icv_scale_mean_raw, icv_scale_sd_raw=icv_scale_sd_raw
      )
    ))
  }
  qc <- data.table(
    label=label, method='lme', status='ok', n_subjects=nrow(extracted), n_rows=nrow(d),
    singular=lme4::isSingular(fit), formula=paste(deparse(form), collapse=' '),
    covariates_fixed=paste(keep_cov, collapse='|'), site_random=isTRUE(site_random), n_sites=n_sites,
    lme4_available=TRUE, lme_message='', ols_fallback_used=FALSE,
    icv_rule=if ('ICV_z' %in% keep_cov) 'first_positive_full_MUSE_source_then_unique_subject_model_sample_zscore' else 'not_used',
    icv_scale_n_rows=icv_scale_n_rows, icv_scale_n_subjects=icv_scale_n_subjects,
    icv_scale_mean_raw=icv_scale_mean_raw, icv_scale_sd_raw=icv_scale_sd_raw
  )
  qc[, `:=`(warnings=paste(unique(fit_warnings), collapse=' | '),
    convergence=paste(fit@optinfo$conv$lme4$messages, collapse=' | '),
    optimizer_code=paste(fit@optinfo$conv$opt, collapse=' | '),
    REML=lme4::isREML(fit), objective=as.numeric(-2*logLik(fit)),
    n_distinct_scans=uniqueN(d$scan_id), n_fixed=length(lme4::fixef(fit)),
    dropped_fixed_columns=paste(attr(lme4::getME(fit,'X'),'col.dropped'), collapse=' | '))]
  # Retain compact row identities alongside the model's stored frame.
  row_identity <- d[, .(RID, scan_id)]
  qc[, model_frame_sha256 := digest::digest(list(frame=fit@frame,
    X=lme4::getME(fit,'X'), rows=row_identity), algo='sha256')]
  list(values=extracted, qc=qc, fit=fit, row_identity=row_identity,
       eligibility=elig, warnings=unique(fit_warnings))
}


f4_identity <- function(x) {
  # Compare ordered values, types and levels without runtime-only table pointers.
  if (is.data.table(x)) return(lapply(as.list(x),identity))
  if (is.list(x)) return(lapply(x,f4_identity))
  x
}

f4_original_site_labels <- function(x) {
  value <- trimws(as.character(x))
  if (anyNA(value) || any(!grepl("^[0-9]+$",value))) stop("Original Figure4 site codes must be finite integer labels")
  restored <- sub("^0+(?=[0-9])","",value,perl=TRUE)
  # Leading-zero normalization must not merge two distinct prepared groups.
  mapping <- unique(data.table(prepared=value,restored=restored))
  if (anyDuplicated(mapping$restored)) stop("Original site restoration would merge prepared site groups")
  restored
}

f4_apply_baseline_reference <- function(d, reference) {
  need(d,c("nucleus","age_baseline","education_years"))
  centers <- reference$centers
  if(nrow(centers)!=2L || !setequal(centers$nucleus,c("LC","SNVTA")) || anyDuplicated(centers$nucleus))
    stop("One original baseline center per nucleus is required")
  index <- match(d$nucleus,centers$nucleus)
  if(anyNA(index) || any(!is.finite(centers$age_center)) || any(!is.finite(centers$education_center))) stop("Invalid baseline centers")
  d[,age_bl_c:=as_num_na(age_baseline)-centers$age_center[index]]
  d[,educ_c:=as_num_na(education_years)-centers$education_center[index]]
  d
}

f4_normalize <- function(d, type, baseline_reference) {
  d <- copy(d)
  aliases <- c(RID="rid", time_c="mri_time_years_c", age_bl_c="age_baseline_c",
    educ_c="education_years_c", apoe4="apoe4_count", dx_bl="diagnosis_baseline",
    echotime_f="echo_time_category", ICV_raw="baseline_icv")
  need(d, c(unname(aliases), "scan_id", "nucleus", "sex", "field_strength", "site_id", "efc_c"))
  for (nm in names(aliases)) d[, (nm) := get(aliases[[nm]])]
  for (nm in c("time_c", "age_bl_c", "educ_c", "apoe4", "ICV_raw", "efc_c"))
    set(d, j=nm, value=number(d[[nm]]))
  d[, `:=`(RID=as.character(RID), scan_id=as.character(scan_id), sex=as.character(sex),
           dx_bl=factor(dx_bl, levels=diagnosis_stage_levels),
           field_strength=as.character(field_strength), echotime_f=as.character(echotime_f),
           site_id=as.character(site_id))]
  if (type == "cr") {
    need(d, c("cr", "hemisphere")); d[, lc_value := number(cr)]
    unique_keys(d, c("nucleus", "scan_id", "hemisphere"), "Figure4 CR")
  } else unique_keys(d, c("nucleus", "scan_id"), "Figure4 MUSE")
  if (any(is.na(d$dx_bl)) || any(!d$nucleus %in% c("LC", "SNVTA")))
    stop("Unrecognized Figure4 nucleus or baseline diagnosis")
  d <- f4_apply_baseline_reference(d,baseline_reference)
  d[,site_id_prepared:=site_id]
  d[,site_id:=f4_original_site_labels(site_id)]
  d
}

f4_preparation_checks <- function(cr, muse, regions) {
  records <- list(); k <- 0L
  add <- function(nucleus, feature, check, passed, value=NA_real_) {
    k <<- k+1L; records[[k]] <<- data.table(nucleus, feature, check, passed, value)
  }
  scales <- list(); s <- 0L
  for (target in c("LC", "SNVTA")) {
    c <- cr[nucleus==target]; m <- muse[nucleus==target]
    keys <- unique(c[, .(scan_id)])
    add(target, "all", "matching_cr_and_structural_scan_sets",
        setequal(keys$scan_id, m$scan_id))
    counts <- c[, .(count=.N), by=scan_id]
    count <- counts$count[match(m$scan_id, counts$scan_id)]
    add(target, "all", "one_or_two_hemispheres_per_scan", all(count %in% c(1L,2L)))
    shared <- intersect(c("RID","time_c","age_bl_c","sex","educ_c","apoe4","dx_bl",
                          "field_strength","efc_c","echotime_f","site_id","ICV_raw"), names(m))
    for (nm in shared) {
      idx <- match(c$scan_id,m$scan_id)
      a <- c[[nm]]; b <- m[[nm]][idx]
      good <- if (is.numeric(a)) all(is.na(a)==is.na(b)) &&
        all(abs(a[is.finite(a)]-b[is.finite(a)]) <= 1e-10) else
        identical(as.character(a),as.character(b))
      add(target, nm, "scan_covariates_match_hemisphere_rows", good)
    }
    duplicate_keys <- character()
    for (feature in regions$variable_column) {
      id <- sub("^h_muse_volume_", "", feature)
      raw <- number(m[[paste0("muse_raw_", id)]])
      prepared <- number(m[[paste0("muse_", id)]])
      weighted <- rep(raw, count)
      mu <- mean(weighted,na.rm=TRUE); sigma <- sd(weighted,na.rm=TRUE)
      expected <- (raw-mu)/sigma
      finite <- is.finite(expected) & is.finite(prepared)
      delta <- if (any(finite)) max(abs(expected[finite]-prepared[finite])) else Inf
      add(target, feature, "matched_hemisphere_preparation_scale",
          identical(is.finite(expected),is.finite(prepared)) && delta <= 1e-10, delta)
      add(target, feature, "finite_scan_fraction", mean(is.finite(prepared)) >= .95,
          mean(is.finite(prepared)))
      duplicate_keys <- c(duplicate_keys,
        paste(ifelse(is.na(prepared),"NA",sprintf("%.17g",prepared)),collapse="\r"))
      s <- s+1L; scales[[s]] <- data.table(nucleus=target, feature,
        center=mu, scale=sigma, n_reference=sum(is.finite(weighted)),
        reference_population="finite matched hemisphere rows before scan collapse",
        prepared_column=paste0("muse_",id), raw_column=paste0("muse_raw_",id))
    }
    add(target,"all","all145_distinct_canonical_features",!anyDuplicated(duplicate_keys))
  }
  list(checks=rbindlist(records), scales=rbindlist(scales))
}

f4_fit_task <- function(spec, prepared) {
  target <- spec$nucleus[[1L]]; measure <- spec$measurement[[1L]]; feature <- spec$feature[[1L]]
  covariates <- c("age_bl_c","sex","educ_c","apoe4","dx_bl")
  if (measure == "CR") {
    d <- prepared$cr[nucleus==target & hemisphere==feature]
    covariates <- c(covariates,"field_strength","efc_c","echotime_f")
    y <- "lc_value"
  } else {
    d <- copy(prepared$muse[nucleus==target]); y <- "muse_y"
    d[, muse_y := get(paste0("muse_",sub("^h_muse_volume_","",feature)))]
    covariates <- c("ICV_z",covariates)
  }
  need(d,setdiff(covariates,"ICV_z"),"trajectory model")
  ans <- fit_lme_trajectory_parameters(d,y,covariates,random_scan=FALSE,
    site_random=measure=="CR",min_n=1L,label=spec$scientific_task[[1L]])
  ans$results <- cbind(spec,ans$qc)
  # Failed trajectories retain their reason; there is no model fallback.
  if (ans$qc$status[[1L]] != "ok") {
    ans$status <- "not_estimable"
    return(ans)
  }
  ans$parameters <- cbind(spec,ans$values)
  co <- as.data.table(coef(summary(ans$fit)),keep.rownames="term")
  ans$coefficients <- cbind(spec,co)
  ans$variance_components <- cbind(spec,as.data.table(as.data.frame(lme4::VarCorr(ans$fit))))
  ans$fixed_covariance <- as.matrix(vcov(ans$fit))
  ans$values <- NULL; ans$qc <- NULL
  ans
}

f4_bilateral <- function(left,right) {
  need(left,c("RID","intercept","slope","n","n_time","followup_span_y"))
  need(right,c("RID","intercept","slope","n","n_time","followup_span_y"))
  cols <- c("RID","intercept","slope","n","n_time","followup_span_y")
  d <- merge(left[,..cols],right[,..cols],by="RID",all=TRUE,sort=FALSE,suffixes=c("_left","_right"))
  d[, intercept_eligible := is.finite(intercept_left) & is.finite(intercept_right) &
      is.finite(n_left) & n_left>=1 & is.finite(n_right) & n_right>=1]
  d[, slope_eligible := is.finite(slope_left) & is.finite(slope_right) &
      is.finite(n_left) & n_left>=3 & is.finite(n_right) & n_right>=3 &
      is.finite(n_time_left) & n_time_left>=2 & is.finite(n_time_right) & n_time_right>=2 &
      is.finite(followup_span_y_left) & followup_span_y_left>=1 &
      is.finite(followup_span_y_right) & followup_span_y_right>=1]
  d[, `:=`(intercept=fifelse(intercept_eligible,(intercept_left+intercept_right)/2,NA_real_),
           slope=fifelse(slope_eligible,(slope_left+slope_right)/2,NA_real_))]
  d
}

f4_standardize <- function(X,Y) {
  if (!identical(as.character(X$RID),as.character(Y$RID))) stop("PLSC X/Y row identity mismatch")
  fx <- setdiff(names(X),"RID"); fy <- setdiff(names(Y),"RID")
  xm <- as.matrix(X[,..fx]); ym <- as.matrix(Y[,..fy])
  ok <- rowSums(!is.finite(xm))==0 & rowSums(!is.finite(ym))==0
  X <- X[ok]; Y <- Y[ok]; xm <- xm[ok,,drop=FALSE]; ym <- ym[ok,,drop=FALSE]
  if (nrow(X)<20L) stop("Fewer than20 complete participants; PLSC bundle unavailable")
  xx <- scale(xm,center=TRUE,scale=TRUE); yy <- scale(ym,center=TRUE,scale=TRUE)
  if (any(!is.finite(xx)) || any(!is.finite(yy)))
    stop("Nonestimable final regional/behavior scale; do not silently remove a region")
  scaling <- rbind(data.table(matrix="X",feature=fx,center=attr(xx,"scaled:center"),scale=attr(xx,"scaled:scale")),
                   data.table(matrix="Y",feature=fy,center=attr(yy,"scaled:center"),scale=attr(yy,"scaled:scale")))
  scaling[, `:=`(n_reference=nrow(X),reference_population="all regional and bilateral raw parameters finite, nucleus and parameter specific")]
  list(X=cbind(X[,.(RID)],as.data.table(xx)),Y=cbind(Y[,.(RID)],as.data.table(yy)),
       raw_X=X,raw_Y=Y,scaling=scaling,n_before=length(ok))
}

f5s_weighted_scale <- function(x,w) {
  x<-number(x);w<-number(w);ok<-is.finite(x)&is.finite(w)&w>0
  if(sum(ok)<2||sum(w[ok])<=1) stop('Invalid weighted reference')
  mu<-sum(w[ok]*x[ok])/sum(w[ok]);s<-sqrt(sum(w[ok]*(x[ok]-mu)^2)/(sum(w[ok])-1))
  if(!is.finite(s)||s<=0) stop('Nonestimable weighted reference')
  value<-rep(NA_real_,length(x));value[ok]<-(x[ok]-mu)/s
  list(value=value,center=mu,scale=s,n_rows=sum(ok),n_weight=sum(w[ok]))
}

f5s_frame <- function(spec,prepared) {
  target<-spec$nucleus[[1L]];d<-copy(prepared$muse[nucleus==target])
  need(d,c('RID','scan_id','tpv','n_hemispheres','time_c','ICV_raw','age_bl_c','sex','educ_c','apoe4','dx_bl','mri_date','muse_date'))
  # Scale all matched hemisphere-equivalent rows before CR support and covariate filtering.
  zz<-f5s_weighted_scale(d$tpv,pmax(1,number(d$n_hemispheres)))
  d[,`:=`(y_traj=zz$value,time_y=number(time_c))]
  cr_ids<-unique(prepared$cr[nucleus==target & is.finite(number(time_c)) & is.finite(number(cr)),RID])
  d<-d[is.finite(y_traj)&is.finite(time_y)&RID %in% cr_ids]
  eligibility<-d[,.(n_scans=uniqueN(scan_id),n_time=uniqueN(time_y),followup_span_y=max(time_y)-min(time_y)),by=RID]
  scientific_assert(uniqueN(d$RID)>=20,'Too few eligible TPV participants')
  d[,ICV_raw:=number(ICV_raw)]
  if(nrow(d[is.finite(ICV_raw)&ICV_raw>0,.(n=uniqueN(ICV_raw)),by=RID][n>1])) stop('Static ICV conflict')
  icv<-unique(d[is.finite(ICV_raw)&ICV_raw>0,.(RID,ICV_raw)],by='RID')
  im<-mean(icv$ICV_raw);iss<-sd(icv$ICV_raw)
  scientific_assert(nrow(icv)>=2&&is.finite(iss)&&iss>0,'Nonestimable ICV scale')
  d[,ICV_z:=(ICV_raw-im)/iss]
  # Apply intercept-map covariate availability after forming the full TPV/ICV references.
  map<-prepared$covariate_map_ids[[target]]
  if(is.null(map)||anyNA(map)||anyDuplicated(map)) stop('Missing/invalid original intercept covariate-map support')
  unavailable<-which(!as.character(d$RID) %in% as.character(map))
  if(length(unavailable)) for(k in c('age_bl_c','sex','educ_c','apoe4')) set(d,i=unavailable,j=k,value=NA)
  for(k in c('age_bl_c','educ_c','apoe4')) set(d,j=k,value=number(d[[k]]))
  d[,sex:=factor(toupper(substr(as.character(sex),1,1)))];d[,dx_bl:=factor(as.character(dx_bl),levels=c('CN','MCI','AD'))]
  requested<-c('ICV_z','age_bl_c','sex','educ_c','apoe4','dx_bl')
  d<-d[complete.cases(d[,..requested])]
  scientific_assert(uniqueN(d$RID)>=20,'Too few complete TPV participants')
  # Numeric/text/numeric round-trips, including ICV_z, preserve model-input precision.
  for(k in requested) if(is.numeric(d[[k]])||is.integer(d[[k]])) set(d,j=k,value=number(d[[k]]))
  keep<-requested[vapply(requested,function(k) uniqueN(d[[k]])>=2,logical(1))]
  scientific_assert('ICV_z' %in% keep,'ICV became nonestimable')
  for(k in requested) if(is.factor(d[[k]])) set(d,j=k,value=droplevels(d[[k]]))
  d[,RID_factor:=factor(RID)]
  form<-as.formula(paste('y_traj ~',paste(c('time_y',keep),collapse=' + '),'+ (1 + time_y || RID_factor)'))
  scales<-data.table(nucleus=target,variable=c('TPV','ICV'),center=c(zz$center,im),scale=c(zz$scale,iss),n_reference=c(zz$n_weight,nrow(icv)),reference_population=c('all_matched_hemisphere_equivalent_rows_before_eligibility_and_covariate_CC','unique_eligible_participants_before_covariate_CC'))
  list(frame=d,formula=form,scales=scales,eligibility=eligibility,requested=requested,keep=keep,
       covariate_map_excluded_rows=length(unavailable))
}

f5s_fit <- function(spec,prepared) {
  target<-spec$nucleus[[1L]];built<-f5s_frame(spec,prepared)
  d<-built$frame;form<-built$formula;eligibility<-built$eligibility
  requested<-built$requested;keep<-built$keep;scales<-built$scales
  cap<-capture_fit(lme4::lmer(form,data=d,REML=TRUE,control=lme4::lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000L),check.rankX='stop.deficient',check.conv.singular='ignore')))
  fit<-cap$fit
  if(!identical(names(lme4::getME(fit,'flist')),'RID_factor')) stop('Unexpected TPV random grouping')
  fx<-lme4::fixef(fit);re<-lme4::ranef(fit)[['RID_factor']]
  need(as.data.table(re),c('(Intercept)','time_y'),'TPV random parameters')
  values<-data.table(RID=rownames(re),intercept=as.numeric(fx['(Intercept)']+re[,'(Intercept)']),slope=as.numeric(fx['time_y']+re[,'time_y']))
  values<-merge(values,eligibility,by='RID',all.x=TRUE,sort=FALSE)
  support<-d[,.(first_mri_date=min(as.Date(mri_date)),last_mri_date=max(as.Date(mri_date)),first_structural_date=min(as.Date(muse_date)),last_structural_date=max(as.Date(muse_date))),by=RID]
  values<-merge(values,support,by='RID',all.x=TRUE,sort=FALSE);values[,nucleus:=target]
  info<-cbind(spec,fit_info(fit,cap$warnings))
  info[,`:=`(status='ok',requested_formula=paste(deparse(form),collapse=' '),n_subjects=uniqueN(d$RID),REML=TRUE,optimizer='bobyqa',site_policy='none_already_harmonized_MUSE',requested_covariates=paste(requested,collapse='|'),retained_covariates=paste(keep,collapse='|'))]
  opt<-unlist(fit@optinfo$conv$opt);msgs<-fit@optinfo$conv$lme4$messages
  numerical_ok<-(!length(opt)||all(opt==0))&&!length(msgs[!grepl('singular|boundary',msgs,ignore.case=TRUE)])
  info[,`:=`(converged=numerical_ok,current_calculation_status=if(numerical_ok)'estimated' else 'qualified_numerical_convergence',optimization_code=paste(opt,collapse='|'))]
  info[,covariate_map_excluded_rows:=built$covariate_map_excluded_rows]
  list(results=info,parameters=values,scales=scales,fit=fit,frame=d,coefficients=cbind(spec,as.data.table(coef(summary(fit)),keep.rownames='term')),variance_components=cbind(spec,as.data.table(as.data.frame(lme4::VarCorr(fit)))))
}

# Preparation constants define reference populations, not fitted parameters.
f4_reference <- function() {
  r <- input('reference_values.csv'); need(r,c('scope','nucleus','name','value','population'))
  r <- r[scope=='figure4_baseline']
  centers <- rbindlist(lapply(c('LC','SNVTA'),function(n) {
    a<-r[nucleus==n & name=='age_center'];e<-r[nucleus==n & name=='education_center']
    if(nrow(a)!=1L||nrow(e)!=1L||any(!is.finite(number(c(a$value,e$value))))) stop('Missing original full-population baseline references')
    data.table(nucleus=n,age_center=number(a$value),education_center=number(e$value),reference_population=a$population)
  }))
  list(centers=centers)
}
f4_prepare <- function() {
  cr<-f4_normalize(scan_members('figure4'),'cr',f4_reference())
  muse<-f4_normalize(input('muse.csv'),'muse',f4_reference())
  regions<-input('regions.csv');need(regions,c('variable_column','volume_index','datadic_text'))
  if(nrow(regions)!=145L||anyDuplicated(regions$variable_column)) stop('Exactly 145 ordered unique regions required')
  for(f in regions$variable_column) {
    id<-sub('^h_muse_volume_','',f);need(muse,c(paste0('muse_',id),paste0('muse_raw_',id)))
    set(muse,j=paste0('muse_',id),value=number(muse[[paste0('muse_',id)]]))
  }
  audit<-f4_preparation_checks(cr,muse,regions)
  atomic_csv(audit$checks,'.work/trajectories/input_checks.csv')
  if(any(!audit$checks$passed)) stop('Ordered input/scaling checks failed; no trajectory model prepared')
  tasks<-rbindlist(lapply(c('LC','SNVTA'),function(target) rbind(data.table(nucleus=target,measurement='CR',feature=c('L','R')),data.table(nucleus=target,measurement='MUSE',feature=regions$variable_column))))
  tasks[,scientific_task:=paste(nucleus,measurement,feature,sep='_')]
  write_plan('trajectories',tasks,list(cr=cr,muse=muse,regions=regions,preparation_scales=audit$scales))
}
f4_finalize <- function() {
  p<-load_plan('trajectories');pars<-coefs<-diag<-vars<-list()
  for(id in seq_len(nrow(p$tasks))) {
    f<-work_path('trajectories',sprintf('task_%04d.rds',id))
    if(!file.exists(f)) stop('Missing trajectory task ',id)
    x<-readRDS(f)
    if(!identical(x$signature,run_signature())) stop('Trajectory signature differs')
    diag[[id]]<-x$results;pars[[id]]<-x$parameters;coefs[[id]]<-x$coefficients;vars[[id]]<-x$variance_components
  }
  d<-rbindlist(diag,fill=TRUE);par<-rbindlist(pars,fill=TRUE)
  record_diagnostics('trajectories',d)
  atomic_csv(par,'private/figure4/trajectory_parameters.csv')
  atomic_csv(rbindlist(coefs,fill=TRUE),'figure4/trajectory_coefficients.csv')
  atomic_csv(rbindlist(vars,fill=TRUE),'figure4/trajectory_variance_components.csv')
  if(nrow(d)!=294L||any(d$status!='ok')) stop('Complete 294-task trajectory dependencies unavailable; diagnostics retained')
  matrices<-bilateral<-scales<-rows<-list();k<-0L
  for(target in c('LC','SNVTA')) {
    lr<-f4_bilateral(par[nucleus==target & measurement=='CR' & feature=='L'],par[nucleus==target & measurement=='CR' & feature=='R'])
    bilateral[[target]]<-cbind(nucleus=target,lr)
    for(parameter in c('intercept','slope')) {
      ids<-lr[get(paste0(parameter,'_eligible'))==TRUE,.(RID)];X<-copy(ids)
      for(f in p$data$regions$variable_column) {
        pp<-par[nucleus==target & measurement=='MUSE' & feature==f]
        valid<-is.finite(pp[[parameter]]) & pp$n>=1
        if(parameter=='slope') valid<-valid & pp$n>=3 & pp$n_time>=2 & is.finite(pp$followup_span_y) & pp$followup_span_y>=1
        pp<-pp[valid,.(RID,value=get(parameter))];setnames(pp,'value',f)
        X<-merge(X,pp,by='RID',all.x=TRUE,sort=FALSE)
      }
      Y<-lr[match(X$RID,RID),.(RID,value=get(parameter))];setnames(Y,'value',paste0(tolower(target),'_value'))
      m<-f4_standardize(X,Y);k<-k+1L
      dest<-file.path('.work','plsc_inputs',target,paste0('lme_',parameter))
      for(n in c('X','Y','scaling')) atomic_csv(m[[n]],file.path(dest,paste0(n,'.csv')))
      atomic_csv(m$X[,.(RID)],file.path(dest,'rows.csv'))
      matrices[[k]]<-cbind(nucleus=target,parameter,m$X)
      rows[[k]]<-cbind(nucleus=target,parameter,m$Y)
      scales[[k]]<-cbind(nucleus=target,parameter,m$scaling)
    }
  }
  atomic_csv(rbindlist(matrices,fill=TRUE),'private/figure4/matrix_X.csv')
  atomic_csv(rbindlist(rows,fill=TRUE),'private/figure4/matrix_Y.csv')
  atomic_csv(rbindlist(rows,fill=TRUE)[,.(nucleus,parameter,RID)],'private/figure4/matrix_rows.csv')
  atomic_csv(rbindlist(bilateral,fill=TRUE),'private/figure4/bilateral_parameters.csv')
  atomic_csv(rbindlist(scales),'figure4/matrix_scaling.csv')
  atomic_csv(p$data$preparation_scales,'figure4/preparation_scaling.csv')
  atomic_csv(p$data$regions,'.work/plsc_inputs/regions.csv')
}
f4_tpv_prepare <- function() {
  p<-load_plan('trajectories');d<-p$data
  fp<-file.path(output_root(),'private/figure4/matrix_rows.csv')
  if(!file.exists(fp)) stop('TPV requires this run completed trajectory intercept support')
  r<-fread(fp,colClasses='character',showProgress=FALSE)
  d$covariate_map_ids<-setNames(lapply(c('LC','SNVTA'),function(n) r[nucleus==n & parameter=='intercept',RID]),c('LC','SNVTA'))
  if(any(abs(number(d$muse$tpv)-number(d$muse$gm_volume)-number(d$muse$wm_volume))>1e-8,na.rm=TRUE)) stop('TPV differs from prepared GM+WM')
  for(k in c('tpv','n_hemispheres')) {need(d$muse,k);set(d$muse,j=k,value=number(d$muse[[k]]))}
  write_plan('tpv',data.table(nucleus=c('LC','SNVTA')),d)
}
f4_tpv_finalize <- function() {
  p<-load_plan('tpv');ans<-lapply(seq_len(nrow(p$tasks)),function(id) {
    f<-work_path('tpv',sprintf('task_%04d.rds',id));if(!file.exists(f)) stop('Missing TPV task ',id)
    x<-readRDS(f);if(!identical(x$signature,run_signature())||identical(x$status,'failed')) stop('TPV task failed or signature differs');x
  })
  for(k in c('parameters','scales','coefficients','variance_components')) {
    d<-rbindlist(lapply(ans,`[[`,k),fill=TRUE)
    atomic_csv(d,file.path(if(k=='parameters')'private/figure4' else 'figure4',paste0('tpv_',k,'.csv')))
  }
  record_diagnostics('tpv',rbindlist(lapply(ans,`[[`,'results'),fill=TRUE))
}
f4_primary_parameters <- function() {
  pars<-fread(file.path(output_root(),'private/figure4/tpv_parameters.csv'),colClasses=c(RID='character'),showProgress=FALSE)
  scores<-fread(file.path(output_root(),'private/figure4/plsc_participant_scores.csv'),colClasses=c(RID='character'),showProgress=FALSE)
  scores<-scores[task_id %in% c(1L,3L)];unique_keys(scores,c('branch','RID'),'direct-intercept scores')
  if(!setequal(unique(scores$model),'lme_intercept')) stop('Nonintercept primary score')
  if(any(abs(scores$brain_score_oriented-scores$brain_score_raw*scores$orientation_factor)>1e-12)) stop('Score orientation identity failure')
  parameters<-merge(scores,pars,by.x=c('branch','RID'),by.y=c('nucleus','RID'),all.x=TRUE,sort=FALSE)
  if(any(!is.finite(parameters$intercept))) stop('Missing TPV parameter for direct-intercept score')
  parameters[,`:=`(rid=RID,nucleus=branch,M_raw=brain_score_oriented,G_raw=intercept,score_scale='oriented_native_direct_intercept_all145_weights',support_policy='full_longitudinal_trajectory_may_include_post_clinical_index_observations')]
  prepared<-load_plan('tpv')$data
  pattern_support<-prepared$muse[,.(trajectory_first_date=min(as.Date(mri_date)),trajectory_last_date=max(as.Date(mri_date))),by=.(rid=as.character(RID),nucleus)]
  parameters<-merge(parameters,pattern_support,by=c('rid','nucleus'),all.x=TRUE,sort=FALSE)
  parameters[,score_definition:='oriented_native_direct_intercept; weights_and_participant_trajectories_estimated_using_full_longitudinal_cohort_not_index_only']
  parameters[,trajectory_support_definition:='bounds_of_prepared_matched_MRI_rows; potential_own_trajectory_support_not_exact_union_of_145_complete_model_frames']
  atomic_csv(parameters,'private/figure4/primary_parameters.csv')
}
f4_main <- function(a=cli()) {
  if(!requireNamespace('lme4',quietly=TRUE)||!requireNamespace('digest',quietly=TRUE)) stop('lme4 and digest required')
  unit<-if(is.null(a$stage))'trajectories' else a$stage
  action<-if(is.null(a$mode))a$action else a$mode
  if(action=='primary-parameters') return(f4_primary_parameters())
  if(!unit %in% c('trajectories','tpv')) stop('Stage must be trajectories or tpv')
  if(action=='prepare') {if(unit=='trajectories')f4_prepare() else f4_tpv_prepare()}
  else if(action=='task') run_task(unit,a[['task-id']],if(unit=='trajectories')f4_fit_task else f5s_fit)
  else if(action=='finalize') {if(unit=='trajectories')f4_finalize() else f4_tpv_finalize()}
  else stop('Use prepare/task/finalize')
}
if(sys.nframe()==0L) f4_main()
