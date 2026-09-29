#!/usr/bin/env Rscript
# Fit the site-clustered ordered mediation models.
.this_file <- function() {
  frames <- vapply(sys.frames(),function(x) if(is.null(x$ofile)) "" else as.character(x$ofile),character(1))
  files <- frames[basename(frames)=='figure5_mediation.R']
  if(length(files))return(normalizePath(tail(files,1L)))
  normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L]))
}
source(file.path(dirname(.this_file()),"common.R"))
f5m_sha <- function(p) {
  if (!file.exists(p) || dir.exists(p)) stop("Missing regular dependency: ", p)
  strsplit(system2("sha256sum", shQuote(p), stdout=TRUE), " +")[[1L]][1L]
}
f5m_runtime <- function(load=TRUE) {
  explicit <- Sys.getenv("FIGURE5_SEM_LIBRARY", "")
  before <- .libPaths()
  if (nzchar(explicit)) {
    if (!dir.exists(explicit)) stop("SEM_RUNTIME_UNAVAILABLE: explicit FIGURE5_SEM_LIBRARY does not exist")
    explicit <- normalizePath(explicit, mustWork=TRUE)
    path <- tryCatch(find.package("lavaan", lib.loc=explicit), error=function(e) "")
    if (!nzchar(path)) stop("SEM_RUNTIME_UNAVAILABLE: explicit FIGURE5_SEM_LIBRARY contains no usable lavaan package")
  } else path <- tryCatch(find.package("lavaan"), error=function(e) "")
  if (!nzchar(path)) stop("SEM_RUNTIME_UNAVAILABLE: approved lavaan 0.7.2 is required")
  path <- normalizePath(path, mustWork=TRUE)
  version <- as.character(utils::packageVersion("lavaan", lib.loc=dirname(path)))
  if (!identical(version, "0.7.2"))
    stop("SEM_VERSION_NOT_APPROVED: resolved lavaan ", version, "; exactly 0.7.2 is required")
  if ("lavaan" %in% loadedNamespaces() &&
      !identical(normalizePath(getNamespaceInfo(asNamespace("lavaan"), "path")), path))
    stop("SEM_NAMESPACE_PATH_MISMATCH: start a fresh R process with the approved library")
  if (nzchar(explicit)) .libPaths(c(explicit, before))
  if (load) {
    ns <- tryCatch(loadNamespace("lavaan", lib.loc=dirname(path)),
                   error=function(e) stop("SEM_RUNTIME_UNAVAILABLE: approved package cannot load: ", conditionMessage(e)))
    if (!identical(normalizePath(getNamespaceInfo(ns, "path")), path)) stop("SEM_NAMESPACE_PATH_MISMATCH")
  }
  files <- list.files(path, recursive=TRUE, full.names=TRUE)
  files <- files[!file.info(files)$isdir]
  list(available=TRUE, reason="", version=version, package_path=path,
       documented_historical_library=explicit,
       default_library_paths=before, actual_library_paths=.libPaths(),
       package_files=lapply(files, function(p) list(path=substring(p,nchar(path)+2L),sha256=f5m_sha(p))),
       estimator="WLSMV", parameterization="theta", ordered="event_5y", cluster="site_id",
       fixed_x=TRUE, meanstructure=TRUE, missing="listwise")
}

f5m_covariates <- function(outcome, field_included=TRUE) {
  c("age_at_scan_c", "sex", "education_c", "APOE4", if(field_included) "field_strength_f",
    if (outcome == "any_worsening") "baseline_dx_f")
}
f5m_frame <- function(d, outcome) {
  # SEM fields X_z, M_z and G_z encode standardized CR, PLSC intercept score and TPV.
  # Base and TPV models share the maximal complete-case reference population.
  d <- copy(as.data.table(d))
  req <- c("rid","event_5y","cr_resid_tech_z","M_raw","G_raw","site_id",f5m_covariates(outcome,FALSE))
  need(d, req, "mediation maximal complete-case frame")
  for (cc in intersect(c("rid","site_id","sex","field_strength_f","baseline_dx_f"),names(d)))
    set(d,j=cc,value=text(d[[cc]]))
  numeric_cols <- c("event_5y","cr_resid_tech_z","M_raw","G_raw","age_at_scan_c","education_c","APOE4")
  for (cc in numeric_cols) set(d,j=cc,value=suppressWarnings(as.numeric(as.character(d[[cc]]))))
  d <- d[complete.cases(d[,..req])]
  field_levels <- if("field_strength_f" %in% names(d)) sort(unique(na.omit(d$field_strength_f))) else character()
  include_field <- length(field_levels)>=2L
  if(include_field) d <- d[!is.na(field_strength_f)]
  d <- d[is.finite(event_5y)&is.finite(cr_resid_tech_z)&is.finite(M_raw)&is.finite(G_raw)]
  field_levels <- if("field_strength_f" %in% names(d)) sort(unique(na.omit(d$field_strength_f))) else character()
  include_field <- length(field_levels)>=2L
  if(any(!d$event_5y %in% c(0,1))) stop("Non-binary ordered mediation outcome")
  unique_keys(d,"rid","mediation pathway frame")
  scales <- list(); reason <- ""
  for (pair in list(c("cr_resid_tech_z","X_z"),c("M_raw","M_z"),c("G_raw","G_z"))) {
    x <- d[[pair[1L]]]; m <- mean(x); s <- sd(x)
    if (!is.finite(s)||s<=0) reason <- paste0("NONESTIMABLE_SCALE_",pair[1L])
    set(d,j=pair[2L],value=if(is.finite(s)&&s>0) (x-m)/s else rep(NA_real_,nrow(d)))
    scales[[length(scales)+1L]] <- data.table(variable=pair[1L],output_variable=pair[2L],mean=m,sd=s,n_reference=length(x),reference_population="maximal_complete_case_pathway_frame",reused_base_tpv=TRUE)
  }
  list(rows=d,scales=rbindlist(scales),covariates=f5m_covariates(outcome,include_field),
       field_strength_included=include_field,field_levels=field_levels,reason=reason)
}
f5m_support <- function(d) {
  data.table(n=nrow(d),n_events=sum(d$event_5y==1L),n_non_events=sum(d$event_5y==0L),
    n_sites=uniqueN(d$site_id),min_site_n=if(nrow(d))min(table(d$site_id))else NA_integer_,
    max_site_n=if(nrow(d))max(table(d$site_id))else NA_integer_,
    sites_no_event_variation=sum(vapply(split(d$event_5y,d$site_id),function(z)uniqueN(z)<2L,logical(1))))
}
f5m_support_reason <- function(d) {
  s <- f5m_support(d)
  if(s$n<50L) "n_lt_50" else if(s$n_events<20L) "events_lt_20" else
    if(s$n_non_events<20L) "non_events_lt_20" else if(s$n_sites<2L) "fewer_than_two_sites" else ""
}
f5m_build <- function(d,covariates,use_tpv) {
  req <- c("X_z","M_z","event_5y","site_id",covariates,if(use_tpv)"G_z")
  need(d,req,"SEM input")
  if(any(!complete.cases(d[,..req])))stop("SEM would change the locked complete-case frame")
  dd <- as.data.frame(d)
  for(cc in covariates) {
    if(cc %in% c("sex","field_strength_f","baseline_dx_f")) dd[[cc]]<-droplevels(factor(text(dd[[cc]])))
    else dd[[cc]]<-number(dd[[cc]])
  }
  mm<-model.matrix(reformulate(covariates),dd)
  mm<-mm[,colnames(mm)!="(Intercept)",drop=FALSE]
  colnames(mm)<-make.names(colnames(mm),unique=TRUE)
  if(ncol(mm)&&qr(mm)$rank<ncol(mm))stop("rank-deficient covariate design")
  out<-data.frame(X_z=dd$X_z,M_z=dd$M_z,event_5y=ordered(as.integer(dd$event_5y),levels=c(0L,1L)),site_id=factor(dd$site_id))
  if(use_tpv)out$G_z<-dd$G_z
  out<-cbind(out,mm)
  if(nlevels(out$site_id)<2L)stop("two-level mediation requires at least two site levels")
  list(data=out,predictor_names=colnames(mm),factor_levels=lapply(dd[intersect(covariates,c("sex","field_strength_f","baseline_dx_f"))],levels),contrasts=getOption("contrasts"))
}
f5m_syntax <- function(covariates,use_tpv=FALSE) paste(
  "level: 1",paste0("  M_z ~ ",paste(c("a*X_z",if(use_tpv)"gM*G_z",covariates),collapse=" + ")),
  paste0("  event_5y ~ ",paste(c("cprime*X_z","b*M_z",if(use_tpv)"gY*G_z",covariates),collapse=" + ")),
  "level: 2","  M_z ~~ M_z","  event_5y ~~ event_5y","  M_z ~~ event_5y",
  "indirect := a*b","direct := cprime","total := cprime + a*b",sep="\n")
f5m_fit <- function(d,covariates,use_tpv=FALSE) {
  built<-f5m_build(d,covariates,use_tpv);syntax<-f5m_syntax(built$predictor_names,use_tpv)
  warnings<-character()
  fit<-withCallingHandlers(tryCatch(lavaan::sem(model=syntax,data=built$data,cluster="site_id",ordered="event_5y",
     estimator="WLSMV",parameterization="theta",fixed.x=TRUE,meanstructure=TRUE,missing="listwise"),error=function(e)e),
     warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  if(inherits(fit,"error"))return(list(fit=NULL,error=conditionMessage(fit),warnings=warnings,syntax=syntax,built=built,converged=FALSE,post_check=FALSE))
  converged<-isTRUE(lavaan::lavInspect(fit,"converged"))
  post_check<-tryCatch(isTRUE(lavaan::lavInspect(fit,"post.check")),error=function(e)FALSE)
  list(fit=fit,error=if(!converged||!post_check)"lavaan convergence/post-check failure"else"",warnings=warnings,syntax=syntax,built=built,converged=converged,post_check=post_check)
}
f5m_effects <- function(pe) {
  pe<-as.data.table(pe)
  need(pe,c("op","label","est","se","z","pvalue","ci.lower","ci.upper","std.all"),"joint SEM parameter table")
  d<-pe[op==":=" & label %in% c("indirect","direct","total")]
  if(nrow(d)!=3L || anyDuplicated(d$label))stop("Joint SEM must return exactly indirect/direct/total")
  d[,.(effect=label,estimate=est,se,z,p_value=pvalue,ci_low=ci.lower,ci_high=ci.upper,
       standardized_all=std.all,effect_scale="latent_probit",uncertainty="joint_SEM_delta_method",odds_ratio=NA_real_)]
}
f5m_public_parameters <- function(pe, estimated) {
  d<-copy(as.data.table(pe))
  if(!estimated)for(nm in intersect(c("est","se","z","pvalue","ci.lower","ci.upper","std.lv","std.all","std.nox"),names(d)))
    set(d,j=nm,value=NA_real_)
  d
}
f5m_prepare <- function() {
  lm<-read_results("private/figure5/primary_landmarks.csv")
  st<-read_results("private/figure4/primary_parameters.csv")
  unique_keys(st,c("rid","nucleus"),"primary pattern and TPV")
  keep<-c("rid","nucleus","M_raw","G_raw",intersect(c("trajectory_first_date","trajectory_last_date","first_structural_date","last_structural_date","score_definition"),names(st)))
  joined<-merge(lm[!is.na(event_5y)],st[,..keep],by=c("rid","nucleus"),sort=FALSE,all=FALSE)
  tasks <- frames <- scales <- list()
  k <- 0L
  for (selected_region in c("LC", "SNVTA")) {
    for (selected_outcome in c("mci_to_dementia", "any_worsening")) {
      k <- k + 1L
      selected_rows <- which(
        joined[["nucleus"]] == selected_region &
        joined[["outcome"]] == selected_outcome
      )
      mediation_data <- joined[selected_rows]
      setorder(mediation_data, rid)
      frame <- f5m_frame(mediation_data, selected_outcome)
      frames[[k]] <- frame
      tasks[[k]] <- data.table(
        nucleus = selected_region, outcome = selected_outcome,
        definition = "primary_trajectory",
        pathway_id = paste("main", selected_region, "intercept", selected_outcome, sep = "__"),
        covariates = paste(frame$covariates, collapse = "|"),
        field_strength_included = frame$field_strength_included,
        n = nrow(frame$rows), unestimated_reason = frame$reason
      )
      scales[[k]] <- cbind(
        data.table(task_id = k, nucleus = selected_region, outcome = selected_outcome),
        frame$scales
      )
    }
  }
  paths<-file.path(output_root(),c("private/figure5/primary_landmarks.csv","private/figure4/primary_parameters.csv"))
  write_plan("mediation",rbindlist(tasks),list(frames=frames,runtime=f5m_runtime(),dependencies=lapply(paths,function(p)list(path=p,sha256=f5m_sha(p)))))
  publish_results(rbindlist(scales),"figure5/mediation_scales.csv")
  atomic_csv(rbindlist(lapply(seq_along(frames),function(i)cbind(data.table(task_id=rep(i,nrow(frames[[i]]$rows))),frames[[i]]$rows)),fill=TRUE),"private/figure5/mediation_model_frames.csv")
}
f5m_task <- function(spec, data) {
  for (dependency in data$dependencies) {
    if (!identical(f5m_sha(dependency$path), dependency$sha256)) {
      stop("Mediation preparation dependency changed")
    }
  }
  frame <- data$frames[[spec$task_id]]
  mediation_data <- frame$rows
  reason <- spec$unestimated_reason
  if (!nzchar(reason)) reason <- f5m_support_reason(mediation_data)
  runtime <- f5m_runtime()
  if (!identical(runtime$version, data$runtime$version) ||
      !identical(runtime$package_path, data$runtime$package_path) ||
      !identical(runtime$package_files, data$runtime$package_files)) {
    stop("SEM runtime identity changed")
  }
  if (!nzchar(reason) && !runtime$available) reason <- runtime$reason
  variants <- c("same_frame_base", "same_frame_tpv")
  status <- effects <- parameters <- diagnostics <- fit_measures <- models <- list()
  support <- f5m_support(mediation_data)
  for (variant in variants) {
    use_tpv <- variant == "same_frame_tpv"
    answer <- NULL
    current_reason <- reason
    if (!nzchar(current_reason)) {
      answer <- tryCatch(
        f5m_fit(mediation_data, frame$covariates, use_tpv),
        error = function(e) list(
          fit = NULL, error = conditionMessage(e), warnings = character(),
          converged = FALSE, post_check = FALSE
        )
      )
      current_reason <- answer$error
      if (nzchar(current_reason) && !grepl(
          "non-conformable|converg|post-check|positive definite|singular|variance|sigma|empty cluster|no within-cluster|not estimable|rank-deficient",
          current_reason, ignore.case = TRUE)) {
        stop("Unexpected SEM error: ", current_reason)
      }
    }
    estimated <- !nzchar(current_reason)
    metadata <- cbind(
      spec,
      data.table(model_variant = variant,
        status = if (estimated) "ESTIMATED" else "NOT_ESTIMABLE", reason = current_reason),
      support[, !"n"]
    )
    status[[variant]] <- metadata
    if (!is.null(answer$fit)) {
      models[[variant]] <- list(
        fit = answer$fit, model_data = answer$built$data,
        factor_levels = answer$built$factor_levels, contrasts = answer$built$contrasts
      )
      parameter_estimates <- as.data.table(
        lavaan::parameterEstimates(answer$fit, standardized = TRUE, ci = TRUE)
      )
      parameters[[variant]] <- cbind(
        metadata, f5m_public_parameters(parameter_estimates, estimated)
      )
      if (estimated) effects[[variant]] <- cbind(metadata, f5m_effects(parameter_estimates))
      metrics <- c("chisq.scaled", "df.scaled", "pvalue.scaled", "cfi.scaled",
        "tli.scaled", "rmsea.scaled", "srmr")
      values <- tryCatch(
        lavaan::fitMeasures(answer$fit, metrics),
        error = function(e) setNames(rep(NA_real_, 7), metrics)
      )
      if (!estimated) values[] <- NA_real_
      fit_measures[[variant]] <- cbind(
        metadata, data.table(metric = names(values), value = as.numeric(values))
      )
    }
    covariance <- if (!is.null(answer$fit)) {
      tryCatch(lavaan::lavInspect(answer$fit, "vcov"), error = function(e) NULL)
    } else NULL
    eigenvalues <- if (!is.null(covariance)) {
      tryCatch(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values,
        error = function(e) NA_real_)
    } else NA_real_
    diagnostics[[variant]] <- cbind(metadata, data.table(
      converged = isTRUE(answer$converged), post_check = isTRUE(answer$post_check),
      warnings = paste(answer$warnings, collapse = " | "),
      syntax = if (!is.null(answer$syntax)) answer$syntax else
        f5m_syntax(frame$covariates, use_tpv),
      vcov_min_eigenvalue = min(eigenvalues), vcov_max_eigenvalue = max(eigenvalues),
      estimator = "WLSMV", parameterization = "theta", ordered_outcome = "event_5y",
      effect_scale = "latent_probit",
      random_structure = "two-level site cluster; between-site mediator/outcome variances and covariance",
      fallback = "none", lavaan_version = runtime$version
    ))
  }
  list(
    results = rbindlist(status, fill = TRUE), effects = rbindlist(effects, fill = TRUE),
    parameters = rbindlist(parameters, fill = TRUE),
    diagnostics = rbindlist(diagnostics, fill = TRUE),
    fit_measures = rbindlist(fit_measures, fill = TRUE), models = models
  )
}

f5m_finalize <- function() {
  tasks <- collect_tasks("mediation")
  status <- bind_component(tasks, "results")
  need(status, c("task_id", "nucleus", "outcome", "definition", "pathway_id",
    "model_variant", "status", "reason", "n"), "completed mediation status")
  if (nrow(status) != 8L || uniqueN(status$task_id) != 4L ||
      anyDuplicated(status[, .(task_id, model_variant)]) ||
      any(status$definition != "primary_trajectory")) {
    stop("Incomplete eight-model primary SEM coverage")
  }
  if (anyDuplicated(status[, .(nucleus, outcome, model_variant)]) ||
      !setequal(status$nucleus, c("LC", "SNVTA")) ||
      !setequal(status$outcome, c("mci_to_dementia", "any_worsening")) ||
      !setequal(status$model_variant, c("same_frame_base", "same_frame_tpv"))) {
    stop("Wrong primary mediation families")
  }
  if (any(status$pathway_id !=
      paste("main", status$nucleus, "intercept", status$outcome, sep = "__"))) {
    stop("Primary mediation requires direct-intercept pathways")
  }
  effects <- bind_component(tasks, "effects")
  families <- list()
  if (!nrow(effects)) effects <- data.table(
    task_id = integer(), nucleus = character(), outcome = character(),
    definition = character(), model_variant = character(), effect = character(),
    estimate = numeric(), se = numeric(), p_value = numeric(), ci_low = numeric(),
    ci_high = numeric(), effect_scale = character()
  )
  effects[, `:=`(q_value = NA_real_, correction_family = NA_character_,
    family_status = "raw_effect_no_adjustment")]
  for (selected_outcome in c("mci_to_dementia", "any_worsening")) {
    for (variant in c("same_frame_base", "same_frame_tpv")) {
      selected_rows <- which(
        status$outcome == selected_outcome & status$model_variant == variant
      )
      family_status <- status[selected_rows]
      effect_rows <- which(effects$outcome == selected_outcome &
        effects$model_variant == variant & effects$effect == "indirect")
      complete <- nrow(family_status) == 2L &&
        all(family_status$status == "ESTIMATED") && length(effect_rows) == 2L &&
        all(is.finite(effects$p_value[effect_rows]))
      family_name <- paste("primary_indirect", selected_outcome, variant, sep = "__")
      if (length(effect_rows)) effects[effect_rows, `:=`(
        q_value = if (complete) p.adjust(p_value, method = "BH") else NA_real_,
        correction_family = family_name,
        family_status = if (complete) "complete" else "qualified_unestimated_member"
      )]
      families[[family_name]] <- data.table(
        correction_family = family_name, expected_nuclei = 2L,
        accounted_statuses = nrow(family_status),
        estimated = sum(family_status$status == "ESTIMATED"),
        status = if (complete) "complete" else "qualified_unestimated_member",
        method = "BH", effect = "indirect"
      )
    }
  }
  # Missing planned effects are NA rows, never substitute another estimator.
  missing <- copy(status)
  padding <- lapply(seq_len(nrow(missing)), function(i) {
    model_status <- missing[i]
    present <- effects[task_id == model_status$task_id &
      model_variant == model_status$model_variant]$effect
    absent <- setdiff(c("indirect", "direct", "total"), present)
    if (!length(absent)) return(NULL)
    cbind(model_status[rep(1L, length(absent))], data.table(
      effect = absent, estimate = NA_real_, se = NA_real_, z = NA_real_,
      p_value = NA_real_, q_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      effect_scale = "latent_probit", uncertainty = "joint_SEM_delta_method",
      odds_ratio = NA_real_
    ))
  })
  effects <- rbindlist(c(list(effects), padding), fill = TRUE)
  publish_results(effects, "figure5/mediation_effects.csv")
  for (component in c("parameters", "fit_measures")) {
    publish_results(bind_component(tasks, component),
      paste0("figure5/mediation_", component, ".csv"))
  }
  diagnostics <- rbindlist(list(
    bind_component(tasks, "diagnostics"), rbindlist(families)
  ), fill = TRUE)
  record_diagnostics("mediation", diagnostics)
}

if(sys.nframe()==0L){a<-cli();if(a$action=="runtime")cat(jsonlite::toJSON(f5m_runtime(),auto_unbox=TRUE),"\n")else if(a$action=="prepare")f5m_prepare()else if(a$action=="task")run_task("mediation",a$task,f5m_task)else if(a$action=="finalize")f5m_finalize()else stop("Expected prepare/task/finalize")}
