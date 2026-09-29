# Primary phenotype associations and shared lag-model frame preparation.
.F3_FILE <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
if (is.null(.F3_FILE)) .F3_FILE <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1L])
F3_DIR <- dirname(normalizePath(.F3_FILE,mustWork=TRUE))
source(file.path(F3_DIR,"common.R"),local=environment())
F3_STAGES <- c("Stable_CN_or_Stable_MCI", "Future_any_worse_5y", "CurrentAD")
F3_FAMILIES <- c(PRIMARY_COGNITION_COMPOSITES=4L, MECHANISTIC_COGNITION=14L,
                 ATN_CORE12=12L, BEHAVIOR_FUNCTIONS=8L)


f3_registry <- function() {
  r <- input("outcomes.csv")
  need(r, c("outcome", "variable_name", "fdr_family", "requires_icv", "raw_direction_code"))
  unique_keys(r, "outcome", "Figure 3 registry")
  if (nrow(r) != 38L) stop("Figure 3 requires the fixed 38-outcome registry")
  r[, source_family := fdr_family]
  r[fdr_family == "GLOBAL_CLINICAL_FUNCTION", fdr_family := "BEHAVIOR_FUNCTIONS"]
  r[, `:=`(requires_icv = truth(requires_icv),
            cognitive = fdr_family %in% c("PRIMARY_COGNITION_COMPOSITES", "MECHANISTIC_COGNITION"))]
  counts <- r[, .N, by = fdr_family]
  if (any(!counts$fdr_family %in% names(F3_FAMILIES)) ||
      !all(counts$N == F3_FAMILIES[counts$fdr_family])) stop("Outcome family counts differ from 4/14/12/8")
  r[, plot_order := seq_len(.N)]
  r
}

f3_date <- function(x) as.Date(as.character(x),format="%Y-%m-%d")

f3_read <- function(unit) {
  if (!unit %in% c("concurrent", "prospective", "lag")) stop("Unknown Figure 3 unit")
  d <- input(paste0(unit, ".csv"))
  common <- c("rid", "participant_id", "scan_set_id", "mri_date", "outcome", "field_strength", "site_id",
              "lc_l", "lc_r", "snvta_l", "snvta_r", "lc_bilateral", "snvta_bilateral",
              "sex", "apoe4_count", "echo_time_category", "baseline_icv",
              "primary_diagnosis_stage", "primary_eligible", "primary_lag_eligible")
  centers <- unlist(lapply(c("lc", "snvta"), function(n) paste0(n, "_", c("age_baseline_c", "mri_time_years_c", "education_years_c", "efc_c"))))
  pair <- if (unit == "concurrent") c("assessment_id", "assessment_date", "value") else
    c("baseline_assessment_id", "baseline_assessment_date", "baseline_value",
      "future_assessment_id", "future_assessment_date", "future_value", "lag_days")
  extra <- if (unit == "prospective") c("window", "pairing") else character()
  need(d, c(common, centers, pair, extra), paste("Figure 3", unit))
  # Retain only model/identity inputs. Prepared outcome values are never re-scaled here.
  d <- d[, unique(c(common, centers, pair, extra)), with = FALSE]
  for (nm in c("lc_l", "lc_r", "snvta_l", "snvta_r", "lc_bilateral", "snvta_bilateral",
               "apoe4_count", "baseline_icv", centers,
               intersect(c("value", "baseline_value", "future_value", "lag_days"), names(d))))
    set(d, j = nm, value = number(d[[nm]]))
  for (nm in c("primary_eligible", "primary_lag_eligible"))
    set(d, j = nm, value = truth(d[[nm]]))
  for (nm in intersect(c("mri_date", "assessment_date", "baseline_assessment_date", "future_assessment_date"), names(d)))
    set(d, j = nm, value = f3_date(d[[nm]]))
  d[, field_strength := field(field_strength)]
  for(prefix in c("lc","snvta")) {
    a<-d[[paste0(prefix,"_l")]];b<-d[[paste0(prefix,"_r")]];bil<-d[[paste0(prefix,"_bilateral")]]
    ok<-is.finite(bil)
    if(any(ok & (!is.finite(a)|!is.finite(b))) ||
       any(abs(bil[ok]-(a[ok]+b[ok])/2)>1e-8*(1+abs(bil[ok]))))
      stop("Prepared bilateral CR is not the two-retained-hemisphere mean: ",prefix)
  }
  r <- f3_registry()
  if (any(!d$outcome %in% r$outcome)) stop("Prepared Figure 3 file contains an unregistered outcome")
  # Essential identities are not imputed, grouped as UNKNOWN, or matched NA to NA.
  idcols <- c("outcome", "scan_set_id", if (unit == "concurrent") "assessment_id" else c("baseline_assessment_id", "future_assessment_id"), extra)
  unique_keys(d, idcols, paste(unit, "pair identities"))
  if (unit == "concurrent") unique_keys(d, c("outcome", "scan_set_id"), "one concurrent phenotype per MRI")
  if (unit == "prospective") unique_keys(d, c("outcome", "scan_set_id", "window", "pairing"), "one selected future per window/pairing")
  if (anyNA(d$mri_date)) stop("Missing/invalid MRI dates in prepared Figure 3 data")
  if (unit == "concurrent") {
    if (anyNA(d$assessment_date) || any(abs(as.numeric(d$assessment_date - d$mri_date)) > 30)) stop("Invalid concurrent 30-day pairing")
  } else {
    if (anyNA(d$baseline_assessment_date) || anyNA(d$future_assessment_date) ||
        any(abs(as.numeric(d$baseline_assessment_date - d$mri_date)) > 30) ||
        any(d$future_assessment_date <= d$mri_date | d$future_assessment_date <= d$baseline_assessment_date) ||
        any(d$baseline_assessment_id == d$future_assessment_id) ||
        any(!is.finite(d$lag_days) | d$lag_days != as.numeric(d$future_assessment_date - d$mri_date))) stop("Invalid prepared future pairing or lag")
    if (unit == "prospective") {
      if (any(d$window != "180_1095") || any(d$pairing != "primary"))
        stop("Only primary 180-1095-day prospective pairs belong in this input")
      if (any(d$lag_days < 180 | d$lag_days > 1095)) stop("Prepared prospective window violated")
    }
  }
  scan_fields <- c("scan_set_id","rid","participant_id","mri_date","field_strength",
                   "lc_l","lc_r","snvta_l","snvta_r","lc_bilateral","snvta_bilateral",centers)
  scan_lookup <- unique(d[, ..scan_fields])
  unique_keys(scan_lookup,"scan_set_id",paste(unit,"repeated MRI measurements/covariates"))
  participants <- unique(d[, .(rid,participant_id)])
  unique_keys(participants,"rid",paste(unit,"RID mapping"))
  unique_keys(participants,"participant_id",paste(unit,"participant mapping"))
  d[, pair_key := do.call(paste, c(.SD, sep = "__")), .SDcols = idcols]
  setorder(d, outcome, participant_id, mri_date, pair_key)
  d
}

f3_tasks <- function(unit) {
  r <- f3_registry()
  if (!unit %in% c("concurrent","prospective","lag")) stop("Unknown Figure 3 unit")
  tasks <- CJ(outcome=r$outcome,nucleus=c("LC","SNVTA"),sorted=FALSE)
  tasks[, `:=`(model_type=unit,
    variant=switch(unit,concurrent="original",prospective="primary",lag="full_lag"),
    window=if(unit=="prospective")"180_1095" else "not_applicable",pairing="primary")]
  tasks <- merge(tasks,r[,.(outcome,variable_name,fdr_family,cognitive,requires_icv,raw_direction_code,plot_order)],by="outcome",sort=FALSE)
  setorder(tasks,variant,window,pairing,plot_order,nucleus)
  tasks
}
f3_prepare <- function(unit) {
  d <- f3_read(unit); tasks <- f3_tasks(unit)
  paths <- setNames(character(length(unique(tasks$outcome))),unique(tasks$outcome))
  baseline <- if(unit=="lag")f3_read("concurrent") else NULL
  for(j in seq_along(paths)) {
    id <- names(paths)[j]; path <- work_path(unit,sprintf("outcome_%02d.rds",j))
    item <- list(data=d[d$outcome==id],baseline=if(is.null(baseline))NULL else baseline[baseline$outcome==id])
    atomic_rds(item,path); paths[j] <- path
  }
  write_plan(unit,tasks,list(outcomes=paths))
}


f3_load_outcome <- function(spec, data) {
  path <- data$outcomes[[as.character(spec$outcome)]]
  if (is.null(path) || !file.exists(path)) stop("Prepared outcome shard is missing")
  readRDS(path)
}

f3_select <- function(d,spec) {
  d <- copy(d)
  if(spec$model_type=="prospective")d <- d[d$window==spec$window & d$pairing==spec$pairing]
  eligible <- if(spec$model_type=="lag")d$primary_lag_eligible else d$primary_eligible
  d[eligible %in% TRUE]
}


f3_scale_fields <- function(x, f, rule = "prospective") {
  z <- rep(NA_real_, length(x)); report <- list()
  for (g in c("1.5T", "3T")) {
    ix <- which(!is.na(f) & f == g & is.finite(x))
    mu <- if (length(ix)) mean(x[ix]) else NA_real_; sig <- if (length(ix)>1L) sd(x[ix]) else NA_real_
    if (is.finite(sig) && sig > 0) z[ix] <- (x[ix] - mu) / sig
    report[[g]] <- list(n = length(ix), mean = mu, sd = sig)
  }
  # Prospective models fill failed field strata from pooled scaling; lag models
  # use a pooled fallback only when fewer than two field-standardized values exist.
  pooled <- finite_z(x)
  fallback <- if (rule == "lag") sum(is.finite(z)) < 2L else TRUE
  use <- if (fallback) which(!is.finite(z) & is.finite(pooled)) else integer()
  z[use] <- pooled[use]
  list(value = z, reference = report, pooled_fallback_rows = length(use), rule = rule,
       pooled_mean = mean(x, na.rm = TRUE), pooled_sd = sd(x, na.rm = TRUE))
}

f3_scale_string <- function(scales) {
  if (!length(scales)) return("not_applicable")
  paste(vapply(names(scales), function(n) {
    s <- scales[[n]]
    paste0(n, ":", paste(vapply(names(s$reference), function(f) {
      a <- s$reference[[f]]; paste0(f, "(n=", a$n, ",mean=", format(a$mean,digits=17), ",sd=", format(a$sd,digits=17), ")")
    }, character(1)),collapse=";"), ";pooled_fallback=",s$pooled_fallback_rows)
  }, character(1)), collapse = " | ")
}

f3_frame <- function(item, spec) {
  d <- f3_select(item$data, spec)
  scientific_assert(nrow(d)>0L, "No eligible prepared pairs for this definition")
  original <- spec$model_type == "concurrent"
  nuclei <- as.character(spec$nucleus)
  center_nucleus <- tolower(nuclei)
  for (nm in c("age_baseline_c", "mri_time_years_c", "education_years_c", "efc_c"))
    d[, (nm) := get(paste0(center_nucleus, "_", nm))]
  scales <- list()
  if (original) {
    parts <- lapply(c("L", "R"), function(h) {
      x <- copy(d); x[, cr := get(paste0(tolower(nuclei), "_", tolower(h)))]
      x[, `:=`(hemisphere = h, hemisphere_c = if (h == "L") -0.5 else 0.5)]
      x[is.finite(cr) & is.finite(value)]
    })
    d <- rbindlist(parts); d[, response := cr]
  } else {
    for (n in nuclei) d <- d[is.finite(get(paste0(tolower(n), "_bilateral")))]
    d <- d[is.finite(baseline_value) & is.finite(future_value)]
    for (n in nuclei) {
      col <- paste0(tolower(n), "_z")
      s <- f3_scale_fields(d[[paste0(tolower(n), "_bilateral")]], d$field_strength,
                           if (spec$model_type == "lag") "lag" else "prospective")
      d[, (col) := s$value]; scales[[n]] <- s
    }
    for (n in nuclei) d <- d[is.finite(get(paste0(tolower(n), "_z")))]
    d[, response := future_value]
  }
  scientific_assert(nrow(d)>0L,"No finite response/exposure pairs")
  lag_center <- NA_real_; icv_reference <- list(mean=NA_real_,sd=NA_real_,n=0L)
  if (spec$model_type != "concurrent") {
    d[, lag_years := lag_days / 365.25]; lag_center <- mean(d$lag_years)
    d[, lag_c := lag_years - lag_center]
  }
  if (isTRUE(spec$requires_icv)) {
    icv <- number(d$baseline_icv)
    if (spec$model_type == "lag") {
      # Standardize ICV using baseline-paired MRIs, not repeated future assessments.
      bspec <- copy(spec); bspec[, `:=`(model_type = "lag", variant = "full_lag")]
      b <- f3_select(item$baseline, bspec)
      for (n in nuclei) b <- b[is.finite(get(paste0(tolower(n),"_bilateral")))]
      b <- b[is.finite(value)]
      refs <- number(b$baseline_icv)
    } else refs <- icv
    refs <- refs[is.finite(refs) & refs>0]
    icv_reference <- list(mean=mean(refs),sd=sd(refs),n=length(refs))
    scientific_assert(length(refs)>=(if (run_profile()=="selftest") 20L else 80L), "Insufficient baseline ICV support")
    scientific_assert(is.finite(icv_reference$sd) && icv_reference$sd>0, "Baseline ICV has no finite scale")
    d[, baseline_icv_z := fifelse(is.finite(baseline_icv) & baseline_icv>0,
                                  (baseline_icv-icv_reference$mean)/icv_reference$sd, NA_real_)]
  }
  d[, `:=`(participant=factor(text(participant_id)), site=factor(text(site_id)),
            scan_set=factor(text(scan_set_id)), sex=factor(text(sex),levels=c("F","M")),
            field_strength=factor(field_strength,levels=c("1.5T","3T")), echo_f=factor(text(echo_time_category)))]
  levels_dx <- if (spec$model_type=="lag") F3_STAGES[1:2] else F3_STAGES
  d[, diagnosis_stage := factor(primary_diagnosis_stage,levels=levels_dx)]
  focal <- if (original) "value" else paste0(tolower(nuclei),"_z")
  required <- unique(c("response",focal,"diagnosis_stage","participant","site","scan_set","rid",
                        if (original) "hemisphere_c",if (spec$model_type!="concurrent") c("baseline_value","lag_c"),
                        if (spec$requires_icv) "baseline_icv_z"))
  optional_numeric <- c("age_baseline_c","mri_time_years_c","education_years_c","apoe4_count","efc_c")
  optional_factor <- c("sex","field_strength","echo_f")
  supported <- function(x,kind) {
    if (kind=="numeric") sum(is.finite(x))>=(if (run_profile()=="selftest") 8L else 20L) else sum(!is.na(x))>0L
  }
  opts <- c(optional_numeric[vapply(optional_numeric,function(n) supported(d[[n]],"numeric") && uniqueN(d[[n]][is.finite(d[[n]])])>1L,logical(1))],
            optional_factor[vapply(optional_factor,function(n) {tt<-table(d[[n]]);tt<-tt[tt>0]; length(tt)>1L && all(tt>=3L)},logical(1))])
  if (spec$model_type == "lag") {
    # The original GAM forms complete cases before dropping constant covariates;
    # an entirely missing adjustment variable does not become an implicit waiver.
    opts <- c(optional_numeric, optional_factor)
    for (lev in levels_dx) {
      g <- d[as.character(diagnosis_stage) == lev]
      minimum <- if (run_profile()=="selftest") 3L else 10L
      scientific_assert(nrow(g)>=minimum && uniqueN(g$participant)>=minimum,
                        paste("Insufficient pre-frame lag stage support:",lev))
    }
  }
  # The source LME repeatedly re-evaluates support and rebuilds from its parent rows.
  # No mandatory stage/exposure/ICV/site term is removed to make a fit succeed.
  for (iter in 1:8) {
    dd <- d[complete.cases(d[,unique(c(required,opts)),with=FALSE])]
    for (n in intersect(c("participant","site","scan_set","sex","field_strength","echo_f"),names(dd))) dd[, (n):=droplevels(get(n))]
    if (!nrow(dd)) not_estimable("No complete model rows")
    bad <- opts[vapply(opts,function(n) {
      x<-dd[[n]]; if (is.factor(x)) {tt<-table(x);length(tt)<2L || (spec$model_type!="lag" && any(tt<3L))} else uniqueN(x)<2L
    },logical(1))]
    if (!length(bad)) break
    opts<-setdiff(opts,bad)
  }
  dd <- d[complete.cases(d[,unique(c(required,opts)),with=FALSE])]
  for (n in intersect(c("participant","site","scan_set","sex","field_strength","echo_f"),names(dd)))
    dd[, (n):=droplevels(get(n))]
  scientific_assert(nrow(dd)>0L,"No complete model rows")
  minrows<-if(run_profile()=="selftest") 40L else 160L
  minsubjects<-if(run_profile()=="selftest") 20L else 80L
  scientific_assert(nrow(dd)>=minrows && uniqueN(dd$participant)>=minsubjects,"Insufficient complete rows/participants")
  if(spec$model_type!="lag") scientific_assert(uniqueN(dd$scan_set)>=(if(run_profile()=="selftest")30L else 120L),"Insufficient independent MRI anchors")
  for (lev in levels_dx) {
    g<-dd[as.character(diagnosis_stage)==lev]
    minimum <- if(run_profile()=="selftest")3L else 10L
    if(spec$model_type=="lag") {
      scientific_assert(nrow(g)>=minimum,paste("Insufficient complete lag stage support:",lev))
    } else {
      scientific_assert(uniqueN(g$scan_set)>=minimum,paste("Insufficient stage MRI support:",lev))
      scientific_assert(uniqueN(g$participant)>=minimum,paste("Insufficient stage participant support:",lev))
    }
  }
  scientific_assert(sum(table(dd$site)>=2L)>=2L,"Mandatory site effect lacks repeated sites")
  for(n in unique(c("response",focal,if(spec$model_type!="concurrent") c("baseline_value","lag_c"),
                         if(spec$requires_icv)"baseline_icv_z")))
    scientific_assert(uniqueN(dd[[n]])>=2L,paste("Constant required variable:",n))
  if(spec$model_type=="lag") scientific_assert(uniqueN(dd$lag_days)>=8L,"Fewer than eight distinct lag days")
  if(spec$model_type=="lag") {
    repeat {
      fixed<-c("baseline_value",focal,"diagnosis_stage",opts,if(spec$requires_icv)"baseline_icv_z")
      mm<-model.matrix(as.formula(paste("response ~",paste(fixed,collapse=" + "))),data=dd)
      if(qr(mm)$rank==ncol(mm))break
      pref<-unique(c(intersect(c("echo_f","efc_c","apoe4_count","education_years_c"),opts),opts))
      scientific_assert(length(pref)>0L,"GAM fixed design rank deficient without removable optional covariate")
      opts<-setdiff(opts,tail(pref,1L))
    }
  }
  repeated_time <- if ("mri_time_years_c" %in% opts) dd[, .(nt=uniqueN(mri_time_years_c)),by=participant][nt>=2,.N] else 0L
  slope <- spec$model_type!="lag" && "mri_time_years_c" %in% opts &&
    uniqueN(dd$scan_set)>=(if(run_profile()=="selftest")60L else 500L) &&
    uniqueN(dd$participant)>=30L && repeated_time>=15L
  list(data=dd,nuclei=nuclei,original=original,focal=focal,optional=opts,
       dropped=setdiff(c(optional_numeric,optional_factor),opts),scales=scales,
       lag_center=lag_center,icv_reference=icv_reference,slope=isTRUE(slope),
       center_origin=paste(nuclei,"prepared hemisphere-reference centers"),
       n_parent_rows=nrow(d))
}

f3_counts <- function(d) {
  data.table(n_rows=nrow(d),n_participants=uniqueN(d$participant_id),n_scans=uniqueN(d$scan_set_id),
    n_assessments=if("future_assessment_id"%in%names(d))uniqueN(d$future_assessment_id) else uniqueN(d$assessment_id),
    n_baseline_assessments=if("baseline_assessment_id"%in%names(d))uniqueN(d$baseline_assessment_id) else NA_integer_,
    n_sites=uniqueN(d$site_id))
}

f3_frame_id <- function(d) {
  f<-tempfile();on.exit(unlink(f),add=TRUE)
  keys<-paste(d$pair_key,if("hemisphere"%in%names(d))d$hemisphere else "bilateral",sep="__")
  writeLines(sort(keys),f,useBytes=TRUE);unname(tools::md5sum(f))
}

f3_meta <- function(spec,frame,model,cr_region) {
  out<-copy(spec);out[,c("task_id"):=NULL]
  out[, `:=`(nucleus=cr_region,comparison_model=model,
    comparison_sample="nucleus_specific",
    response_scale=if(frame$original)"raw hemisphere CR" else "unchanged prepared phenotype",
    predictor_scale=if(frame$original)"unchanged prepared phenotype" else "within-field SD of raw bilateral CR",
    scale_reference=f3_scale_string(frame$scales),center_origin=frame$center_origin,
    lag_center_years=frame$lag_center,icv_reference_mean=frame$icv_reference$mean,
    icv_reference_sd=frame$icv_reference$sd,icv_reference_rows=frame$icv_reference$n,
    diagnosis_definition="prepared primary stage; no diagnosis rematching",
    stage_rows=paste(names(table(frame$data$diagnosis_stage)),as.integer(table(frame$data$diagnosis_stage)),sep="=",collapse=";"),
    n_parent_rows=frame$n_parent_rows,frame_id=f3_frame_id(frame$data),
    dropped_unsupported_covariates=paste(frame$dropped,collapse=";"),run_profile=run_profile())]
  cbind(out,f3_counts(frame$data))
}

f3_add_meta <- function(d,meta) {
  d<-copy(d)
  for(n in names(meta)) if(!n%in%names(d))d[, (n):=meta[[n]][1L]]
  d
}

f3_missing <- function(spec,reason) {
  ans <- copy(spec)
  if ("task_id" %in% names(ans)) ans[, task_id := NULL]
  ans[, `:=`(nucleus = as.character(nucleus), comparison_model = "single",
    status = "not_estimable", reason = reason, beta = NA_real_, se = NA_real_,
    statistic = NA_real_, df = NA_real_, p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    n_rows = NA_integer_, n_participants = NA_integer_, n_scans = NA_integer_,
    n_assessments = NA_integer_, run_profile = run_profile())]
  list(results=ans)
}

f3_call_task <- function(unit,id,fit_fun) {
  run_task(unit,id,function(spec,data) {
    item<-f3_load_outcome(spec,data);fr<-NULL
    tryCatch({fr<-f3_frame(item,spec);ans<-fit_fun(spec,fr)
      # Typed objects remain temporary; the finalizer emits participant frames and
      # factor references as private CSVs, never as public coefficient outputs.
      ans},
      model_not_estimable=function(e) {
        ans<-f3_missing(spec,conditionMessage(e))
        candidate<-f3_select(item$data,spec)
        ans$results[, `:=`(n_candidate_pairs=nrow(candidate),
          n_candidate_participants=uniqueN(candidate$participant_id),
          n_candidate_scans=uniqueN(candidate$scan_set_id))]
        if(!is.null(fr)) {
          counts<-f3_counts(fr$data)
          for(n in names(counts))ans$results[, (n):=counts[[n]][1L]]
          ans$frame<-fr
        }
        ans
      })
  })
}

f3_lme_family <- function(spec,fr,formula_builder) {
  combinations<-setNames(list(fr$focal),"single")
  fits<-list();warnings<-list();fallback_reason<-""
  for (attempt in 1:2) {
    failed<-NULL
    for(m in names(combinations)) {
      f<-formula_builder(fr,combinations[[m]],fr$slope && attempt==1L)
      mm<-model.matrix(lme4::nobars(f),data=as.data.frame(fr$data))
      scientific_assert(qr(mm)$rank==ncol(mm),paste("Rank-deficient fixed design:",m))
      cap<-tryCatch(capture_fit(lmerTest::lmer(f,data=as.data.frame(fr$data),REML=TRUE,na.action=na.fail,
        control=lme4::lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=200000L),check.rankX="stop.deficient"))),error=function(e)e)
      if(inherits(cap,"error")){failed<-cap;break}
      check_convergence(cap$fit)
      fits[[m]]<-cap$fit;warnings[[m]]<-cap$warnings
    }
    if(is.null(failed))break
    if(!fr$slope || attempt==2L)stop(failed)
    # Only a thrown random-time-slope fit permits the participant-intercept retry.
    fallback_reason<-conditionMessage(failed);fits<-list();warnings<-list()
  }
  results<-coefs<-list()
  for(m in names(fits)) {
    fit<-fits[[m]];c<-coef_lmm(fit)
    for(focal in combinations[[m]]) {
      nucleus<-if(focal=="lc_z")"LC" else if(focal=="snvta_z")"SNVTA" else as.character(spec$nucleus)
      meta<-f3_meta(spec,fr,m,nucleus)
      row<-c[term==focal]
      if(nrow(row)!=1L)stop("Focal coefficient absent from fitted model: ",focal)
      row[, `:=`(effect_scope=if(spec$cognitive)"reference-stage slope: Stable_CN_or_Stable_MCI" else "common adjusted slope",
        status="estimated",reason=NA_character_,formula=paste(deparse(formula(fit)),collapse=" "),
        singular=lme4::isSingular(fit,tol=1e-4),warnings=paste(warnings[[m]],collapse=" | "),
        random_slope_fallback=fallback_reason)]
      results[[length(results)+1L]]<-f3_add_meta(row,meta)
    }
    meta<-f3_meta(spec,fr,m,as.character(spec$nucleus))
    c[, `:=`(status="estimated",formula=paste(deparse(formula(fit)),collapse=" "),random_slope_fallback=fallback_reason)]
    coefs[[m]]<-f3_add_meta(c,meta)
  }
  list(results=rbindlist(results,fill=TRUE),coefficients=rbindlist(coefs,fill=TRUE),fits=fits,frame=fr)
}


fit_concurrent <- function(spec,frame) {
  builder <- function(fr,focal,slope) {
    interest <- if(isTRUE(spec$cognitive)) paste(paste0(focal," * diagnosis_stage"),collapse=" + ") else
      paste(c(focal,"diagnosis_stage"),collapse=" + ")
    covariates <- c(if(fr$original)"hemisphere_c",fr$optional,if(spec$requires_icv)"baseline_icv_z")
    random <- c(if(slope)"(1 + mri_time_years_c || participant)" else "(1 | participant)",
                "(1 | site)",if(fr$original)"(1 | scan_set)")
    as.formula(paste("response ~",paste(c(interest,covariates,random),collapse=" + ")))
  }
  f3_lme_family(spec,frame,builder)
}

fit_prospective <- function(spec,frame) {
  builder <- function(fr,focal,slope) {
    use_interaction <- isTRUE(spec$cognitive)
    interest <- if(use_interaction) paste(c("baseline_value * diagnosis_stage",
                       paste0(focal," * diagnosis_stage")),collapse=" + ") else
                       paste(c("baseline_value",focal,"diagnosis_stage"),collapse=" + ")
    covariates <- c("lag_c",fr$optional,if(spec$requires_icv)"baseline_icv_z")
    random <- c(if(slope)"(1 + mri_time_years_c || participant)" else "(1 | participant)","(1 | site)")
    as.formula(paste("response ~",paste(c(interest,covariates,random),collapse=" + ")))
  }
  f3_lme_family(spec,frame,builder)
}

f3_fdr <- function(d) {
  d <- copy(d)
  if(!nrow(d))return(d)
  d[, `:=`(p_fdr=NA_real_,fdr_n_planned=as.integer(F3_FAMILIES[fdr_family]),
    multiplicity_scope="BH within complete outcome family / analysis / nucleus")]
  groups <- c("variant","window","pairing","comparison_model","nucleus","fdr_family")
  if(any(d[,.N,by=c(groups,"outcome")]$N>1L))stop("More than one focal test per declared FDR cell")
  d[,p_fdr:={
    n <- unique(fdr_n_planned)
    if(length(n)!=1L || is.na(n) || .N!=n)stop("Incomplete multiplicity family")
    p.adjust(p_value,method="BH",n=n)
  },by=groups]
  d
}
f3_public <- function(d) {
  x <- copy(d)
  # Native invalid probabilities remain in run_diagnostics; never clamp to zero.
  for(n in intersect(c("p_value","p-value","p_fdr"),names(x))) {
    bad <- is.finite(x[[n]]) & (x[[n]]<0 | x[[n]]>1)
    if(any(bad))set(x,i=which(bad),j=n,value=NA_real_)
  }
  remove <- intersect(c("status","reason","warnings","convergence_message","singular",
    "random_slope_fallback","run_profile","frame_id","dropped_unsupported_covariates"),names(x))
  if(length(remove))x[,(remove):=NULL]
  if("p-value" %in% names(x)) {
    if(!"p_value" %in% names(x))x[,p_value:=NA_real_]
    rows <- which(is.na(x$p_value) & !is.na(x[["p-value"]]))
    if(length(rows))set(x,i=rows,j="p_value",value=x[["p-value"]][rows])
    x[,`p-value`:=NULL]
  }
  map <- c(p_value="P",p_fdr="P_fdr",se="SE")
  for(n in intersect(names(map),names(x)))setnames(x,n,map[[n]])
  # All-unavailable components still have an explicit, header-only CSV schema.
  if(!ncol(x))x<-data.table(model_id=character(),term=character(),estimate=numeric(),SE=numeric(),ci_low=numeric(),ci_high=numeric(),P=numeric(),P_fdr=numeric())
  x
}
f3_diagnostics <- function(unit,components) {
  rows <- list()
  fields <- c("outcome","nucleus","model_type","variant","term","record_type","status","reason",
    "warnings","convergence_message","singular","random_slope_fallback","formula","frame_id",
    "dropped_unsupported_covariates","run_profile","n_rows","n_participants","n_scans")
  for(name in names(components)) {
    x <- copy(components[[name]])
    if(!nrow(x))next
    for(p in intersect(c("p_value","p-value"),names(x))) {
      bad <- which(is.finite(x[[p]]) & (x[[p]]<0 | x[[p]]>1))
      if(length(bad)) {
        a <- x[bad,intersect(fields,names(x)),with=FALSE]
        a[, `:=`(record_type="invalid_auxiliary_probability",native_probability=x[[p]][bad],
          publication_probability=NA_real_,probability_column=p,
          reason="Native probability outside [0,1]; unavailable in scientific CSV; focal inference unchanged")]
        rows[[length(rows)+1L]] <- a
      }
    }
    if(name=="results")rows[[length(rows)+1L]] <- x[,intersect(fields,names(x)),with=FALSE]
  }
  record_diagnostics(paste0("figure3_",unit),rbindlist(rows,fill=TRUE))
}
f3_display_metadata <- function() {
  path <- file.path(output_root(),"figure3","prospective_associations.csv")
  r <- f3_registry()[,.(outcome,variable_name,plot_order)]
  r[, `:=`(primary_model_required=TRUE,main_lag_display_eligible=NA,
    main_lag_display_rule="primary 180-1095-day prospective P_fdr < 0.05 in either nucleus; both curves; registry order")]
  if(file.exists(path)) {
    a <- fread(path)
    if(nrow(a)!=76L || uniqueN(a[,.(outcome,nucleus)])!=76L)stop("Incomplete finalized primary prospective results")
    selected <- a[is.finite(P_fdr) & P_fdr<.05,outcome]
    r[,main_lag_display_eligible:=outcome %in% selected]
  }
  atomic_csv(r,"figure3/display_metadata.csv")
}
f3_collect_tasks <- function(unit) {
  plan <- load_plan(unit); frames <- factor_rows <- references <- list()
  tasks <- vector("list",nrow(plan$tasks))
  for(id in seq_len(nrow(plan$tasks))) {
    task <- readRDS(work_path(unit,sprintf("task_%04d.rds",id)))
    if(!identical(task$signature,run_signature()) || !identical(task$status,"complete"))
      stop("Incomplete task while exporting private model frames")
    tasks[[id]] <- drop_task_objects(task)
    if(is.null(task$frame))next
    spec <- task$task; d <- copy(task$frame$data)
    for(n in c("model_type","outcome","nucleus","variant","window","pairing"))d[,(n):=spec[[n]]]
    frames[[length(frames)+1L]] <- d
    factors <- names(d)[vapply(d,is.factor,logical(1))]
    for(n in factors)factor_rows[[length(factor_rows)+1L]] <- data.table(
      model_type=unit,outcome=spec$outcome,nucleus=spec$nucleus,variable=n,
      level_order=seq_along(levels(d[[n]])),level=levels(d[[n]]),
      ordered=is.ordered(d[[n]]),contrast=if(is.ordered(d[[n]]))"contr.poly" else "contr.treatment")
    fr <- task$frame
    for(n in names(fr$scales))for(g in names(fr$scales[[n]]$reference)) {
      x <- fr$scales[[n]]; a <- x$reference[[g]]
      references[[length(references)+1L]] <- data.table(model_type=unit,outcome=spec$outcome,
        nucleus=n,quantity="raw bilateral CR",reference_group=g,n=a$n,center=a$mean,scale=a$sd,
        pooled_center=x$pooled_mean,pooled_scale=x$pooled_sd,pooled_fallback_rows=x$pooled_fallback_rows,
        reference_rule=x$rule)
    }
    references[[length(references)+1L]] <- data.table(model_type=unit,outcome=spec$outcome,
      nucleus=spec$nucleus,quantity=c("lag_years","baseline_icv"),reference_group="parent_before_covariate_complete_cases",
      n=c(fr$n_parent_rows,fr$icv_reference$n),center=c(fr$lag_center,fr$icv_reference$mean),
      scale=c(NA_real_,fr$icv_reference$sd),
      reference_rule=c("mean of finite exposure/outcome parent pairs",
        if(unit=="lag")"finite positive ICV in eligible concurrent baseline-paired MRI anchors" else "finite positive ICV in finite exposure/outcome parent"))
  }
  if(length(frames))atomic_csv(rbindlist(frames,fill=TRUE),file.path("private","figure3",paste0(unit,"_model_frames.csv")))
  if(length(factor_rows))atomic_csv(rbindlist(factor_rows),file.path("private","figure3",paste0(unit,"_factor_references.csv")))
  if(length(references))atomic_csv(rbindlist(references,fill=TRUE),file.path("figure3",paste0(unit,"_scale_references.csv")))
  tasks
}
f3_finalize_lme <- function(unit) {
  tasks <- f3_collect_tasks(unit)
  a <- f3_fdr(bind_component(tasks,"results")); cc <- bind_component(tasks,"coefficients")
  if(!nrow(cc))cc <- data.table(term=character(),beta=numeric(),status=character())
  setorderv(a,c("variant","window","pairing","plot_order","nucleus","comparison_model"))
  f3_diagnostics(unit,list(results=a,coefficients=cc))
  atomic_csv(f3_public(cc),file.path("figure3",paste0(unit,"_coefficients.csv")))
  atomic_csv(f3_public(a),file.path("figure3",paste0(unit,"_associations.csv")))
  if(unit=="prospective")f3_display_metadata()
}
f3_association_main <- function() {
  args <- cli()
  if(is.null(args$unit) || !args$unit %in% c("concurrent","prospective"))stop("Use --unit concurrent or prospective")
  unit <- args$unit
  if(args$action=="prepare")f3_prepare(unit) else
  if(args$action=="task")f3_call_task(unit,args$task,if(unit=="concurrent")fit_concurrent else fit_prospective) else
  if(args$action=="finalize")f3_finalize_lme(unit) else stop("Use --action prepare, task or finalize")
}
if(sys.nframe()==0L)f3_association_main()
