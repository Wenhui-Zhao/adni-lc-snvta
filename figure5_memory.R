#!/usr/bin/env Rscript
# Fit CR, PLSC and TPV interactions with amyloid and memory time.
.this_file <- function() {
  frames <- vapply(sys.frames(),function(x) if(is.null(x$ofile)) "" else as.character(x$ofile),character(1))
  files <- frames[basename(frames)=='figure5_memory.R']
  if(length(files))return(normalizePath(tail(files,1L)))
  normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L]))
}
source(file.path(dirname(.this_file()),"common.R"))
a5_spec <- list(
  schema = "figure5_amyloid_r2_v1", minimum_subjects = 20L,
  minimum_observations = 40L, minimum_binary_group_subjects = 10L,
  minimum_memory_dates = 2L, minimum_memory_span_days = 180L,
  source_fit_sha256 = "2a55ecc421a85aa04412a969aaa11d4c3cc902824188d5adb56f24a74972fac5",
  source_frame_sha256 = "2cf7a2eb63d21619fee3d1eaf7171836a58c73a4d94fd43127edce8f78309c3e",
  source_parent_sha256 = "8761723f97e28b17e7403673f4d8bde3faff8f80e6b72b23e1df7fd67eb768cf",
  source_collector_sha256 = "d84bf915c9d8f68aa5f757017f77d498032047ca408cf31d10aeab4e302eabe3",
  comparison = list(discrete = "exact", frame_atol = 1e-12, frame_rtol = 0,
    coefficient_atol = 1e-8, coefficient_rtol = 1e-6,
    policy = "Recorded before new fits. Changed upstream scores qualify coefficient comparisons; failed fields remain failed."),
  estimator = "lmerTest::lmer REML=TRUE; bobyqa maxfun=100000; check.rankX=stop.deficient",
  random_effects = "(1 | site_id) + (1 | SubjectID)",
  resampling = "none", model_variant = "R2_local_plsc_tpv",
  memory_clock = "(memory_date-max(original MRI date,original PET date))/365.25; center over all unique primary parent RID/memory_date rows before span/structural/model selection",
  standardization = "LC/SNVTA parent z and oriented native direct-intercept M_raw/TPV G_raw: participant-static sample z on nucleus/encoding maximal primary frame",
  multiplicity = "BH planned m=2 nuclei separately for scope/encoding/local_pathology_time or plsc_pathology_time; TPV no BH; finalize only after all planned task statuses")

a5_write <- function(d,name) {
  private<-startsWith(name,"private/")||name=="site_level_mapping.tsv"
  name<-sub("^private/","",name);name<-sub("[.]tsv$",".csv",name)
  rel<-file.path(if(private)"private/figure5"else"figure5",paste0("memory_",name))
  if(private)atomic_csv(d,rel)else publish_results(d,rel)
}
a5_hash <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
a5_object_hash <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)
a5_nobars <- function(x) if (requireNamespace("reformulas", quietly = TRUE)) reformulas::nobars(x) else lme4::nobars(x)
a5_required_numeric <- c("memory", "pathology_model", "lc_parent_z", "snvta_parent_z",
  "memory_time_years", "memory_time_c", "age_baseline_parent_c", "education_years_parent_c", "apoe4_count", "efc_parent_c")
a5_flag <- function(x, label) {
  value <- tolower(trimws(as.character(x)))
  if (any(is.na(value) | !value %in% c("true", "false", "1", "0", "t", "f"))) stop("Invalid/missing Boolean: ", label)
  value %in% c("true", "1", "t")
}
a5_same <- function(actual, expected, label, atol = 1e-12) {
  if (length(actual) != length(expected) || any(!is.finite(actual)) ||
      any(!is.finite(expected)) || any(abs(actual - expected) > atol)) stop("Invalid prepared ", label)
  invisible(TRUE)
}
a5_parent <- function(d) {
  d <- copy(as.data.table(d))
  need(d, c("rid", "participant_id", "scan_id", "amyloid_encoding", "mri_date", "pet_date", "index_date",
    "memory_date", "memory_assessment_id", "site_id", "sex", "diagnosis_baseline", "echo_time_category",
    "field_strength", "amyloid_positive", "centiloid", "memory_time_parent_center", "centiloid_parent_mean",
    "centiloid_parent_sd", "lc_bilateral", "snvta_bilateral", "lc_parent_mean", "lc_parent_sd",
    "snvta_parent_mean", "snvta_parent_sd",
    "memory_sem", "memory_precise", "image_qc", "timing_qc", "pet_warning", a5_required_numeric), "amyloid candidate")
  unique_keys(d, c("amyloid_encoding", "rid", "memory_date"), "amyloid candidate")
  if (!all(d$amyloid_encoding %in% c("binary", "centiloid"))) stop("Unsupported amyloid encoding")
  for (nm in unique(c(a5_required_numeric, "amyloid_positive", "centiloid", "memory_time_parent_center",
                     "centiloid_parent_mean", "centiloid_parent_sd", "lc_bilateral", "snvta_bilateral",
                     "lc_parent_mean", "lc_parent_sd", "snvta_parent_mean", "snvta_parent_sd",
                     "memory_sem", "memory_precise"))) set(d, j = nm, value = number(d[[nm]]))
  for (nm in c("mri_date", "pet_date", "index_date", "memory_date")) {
    set(d, j = nm, value = as.Date(d[[nm]]))
    if (anyNA(d[[nm]])) stop("Invalid prepared amyloid date: ", nm)
  }
  # These optional provenance flags do not select the active primary cohort.
  legacy_flags <- intersect(c("ad_any_record", "lifetime_no_ad", "current_nonad_verified"), names(d))
  for (nm in c(legacy_flags, "image_qc", "timing_qc", "pet_warning"))
    set(d, j = nm, value = a5_flag(d[[nm]], nm))
  if (all(c("lifetime_no_ad", "ad_any_record") %in% names(d)) &&
      any(d$lifetime_no_ad != !d$ad_any_record)) stop("Lifetime AD flags disagree")
  if (any(!d$image_qc | !d$timing_qc | d$pet_warning) || any(!d$amyloid_positive %in% c(0, 1)) ||
      any(!is.finite(d$memory_sem) | d$memory_sem > 0.6 | d$memory_precise != 1)) stop("Prepared amyloid/memory QC contract violated")
  if (any(!d$diagnosis_baseline %in% c("CN", "MCI")) || any(abs(as.numeric(d$pet_date - d$mri_date)) > 90) ||
      any(d$memory_date < d$index_date)) stop("Prepared primary pair/diagnosis/clock eligibility violated")
  expected_index <- pmax(as.numeric(d$mri_date), as.numeric(d$pet_date))
  a5_same(as.numeric(d$index_date), expected_index, "index=max(MRI,PET)", 0)
  a5_same(d$memory_time_years, as.numeric(d$memory_date - d$index_date) / 365.25, "memory time")
  audit <- list()
  for (enc in unique(d$amyloid_encoding)) {
    z <- d[amyloid_encoding == enc]
    static_fields <- c("scan_id", "mri_date", "pet_date", "index_date", "pathology_model",
      "lc_parent_z", "snvta_parent_z", intersect(c("ad_any_record", "lifetime_no_ad"), names(d)))
    for (nm in static_fields) if (any(z[, uniqueN(get(nm)), by = rid]$V1 != 1L)) stop("Parent field not participant-static: ", nm)
    ref_time <- unique(z[, .(rid, memory_date, memory_time_years)])
    mu <- mean(ref_time$memory_time_years)
    a5_same(z$memory_time_parent_center, rep(mu, nrow(z)), "parent time center")
    a5_same(z$memory_time_c, z$memory_time_years - mu, "centered memory time")
    if (enc == "binary") a5_same(z$pathology_model, z$amyloid_positive, "binary encoding", 0) else {
      if (any(!is.finite(z$centiloid_parent_sd) | z$centiloid_parent_sd <= 0)) stop("Invalid primary Centiloid scale")
      a5_same(z$pathology_model, (z$centiloid - z$centiloid_parent_mean) / z$centiloid_parent_sd, "Centiloid source scale")
    }
    for (nuc in c("lc", "snvta"))
      a5_same(z[[paste0(nuc, "_parent_z")]],
        (z[[paste0(nuc, "_bilateral")]] - z[[paste0(nuc, "_parent_mean")]]) / z[[paste0(nuc, "_parent_sd")]],
        paste(nuc, "source scale"))
    audit[[enc]] <- data.table(encoding = enc, n_parent_rows = nrow(z), n_parent_subjects = uniqueN(z$rid),
      memory_time_center = mu, parent_memory_order_hash = a5_object_hash(ref_time),
      pair_clock_qc_verified = TRUE, no_rematching = TRUE)
  }
  # Use integer-text site levels without merging site groups.
  original_site <- text(d$site_id)
  if (anyNA(original_site) || any(!grepl("^[0-9]+$", original_site))) stop("Expected numeric original reserve site identifiers")
  d[, site_id := as.character(as.integer(original_site))]
  mapping <- unique(data.table(before = original_site, after = d$site_id))
  if (anyDuplicated(mapping$before) || anyDuplicated(mapping$after)) stop("Site relabeling changed site partitions")
  list(data = d, audit = rbindlist(audit), site_mapping = mapping)
}

a5_covariates <- function(d) {
  continuous <- c("Age_bl_c", "Educ_c", "APOE4", "EFC_c")
  categorical <- c("Sex", "DX_bl", "EchoTime_f")
  selected <- character(); omitted <- character()
  for (nm in continuous) {
    v <- d[[nm]][is.finite(d[[nm]])]
    if (length(v) >= 2L && uniqueN(v) >= 2L) selected <- c(selected, nm) else
      omitted <- c(omitted, paste0(nm, ":nonestimable_zero_or_missing_variance_in_locked_frame"))
  }
  for (nm in categorical) {
    if (uniqueN(na.omit(text(d[[nm]]))) >= 2L) selected <- c(selected, nm) else
      omitted <- c(omitted, paste0(nm, ":nonestimable_fewer_than_2_levels_in_locked_frame"))
  }
  field <- uniqueN(na.omit(text(d$field_strength))) >= 2L
  if (field) selected <- append(selected, "field_strength", after = min(3L, length(selected))) else
    omitted <- c(omitted, "field_strength:fewer_than_2_observed_levels_in_locked_frame")
  list(included = selected, omitted = omitted, field = field)
}
a5_make_frames <- function(parent, structure) {
  validated <- a5_parent(parent); d <- validated$data
  s <- copy(as.data.table(structure)); need(s, c("rid", "nucleus", "M_raw", "G_raw", "trajectory_first_date", "trajectory_last_date",
    "first_structural_date", "last_structural_date", "score_definition"), "Figure 5 primary structure")
  unique_keys(s, c("rid", "nucleus"), "Figure 5 primary structure")
  s[, `:=`(rid = text(rid), M_raw = number(M_raw), G_raw = number(G_raw))]
  for (nm in c("trajectory_first_date", "trajectory_last_date", "first_structural_date", "last_structural_date")) {
    set(s, j = nm, value = as.Date(s[[nm]]))
    if (anyNA(s[[nm]])) stop("Invalid/missing structural support date: ", nm)
  }
  if (any(s$trajectory_last_date < s$trajectory_first_date | s$last_structural_date < s$first_structural_date)) stop("Invalid structural support interval")
  d[, rid := text(rid)]
  tasks <- list(); frames <- list(); scales <- list(); flow <- list(); k <- 0L
  aliases <- c(RID = "rid", SubjectID = "participant_id", cognition_date = "memory_date", outcome_value = "memory",
    LC_z = "lc_parent_z", SNVTA_z = "snvta_parent_z", Time_yr = "memory_time_years", Time_c = "memory_time_c",
    Age_bl_c = "age_baseline_parent_c", Educ_c = "education_years_parent_c", APOE4 = "apoe4_count", EFC_c = "efc_parent_c",
    Sex = "sex", DX_bl = "diagnosis_baseline", EchoTime_f = "echo_time_category")
  scope <- "primary"
  for (encoding in c("binary", "centiloid")) for (selected_region in c("LC", "SNVTA")) {
    k <- k + 1L; key <- paste(encoding, selected_region, sep = "__")
    z <- d[amyloid_encoding == encoding]
    z <- merge(z, s[nucleus == selected_region], by = "rid", all = FALSE, sort = FALSE)
    for (nm in names(aliases)) set(z, j = nm, value = z[[aliases[[nm]]]])
    for (nm in c("RID", "SubjectID", "site_id", "Sex", "DX_bl", "EchoTime_f", "field_strength")) set(z, j = nm, value = text(z[[nm]]))
    required <- c("RID", "SubjectID", "site_id", "cognition_date", "outcome_value", "pathology_model", "LC_z", "SNVTA_z",
      "M_raw", "G_raw", "Time_yr", "Time_c", "Age_bl_c", "Sex", "DX_bl", "Educ_c", "APOE4", "EFC_c", "EchoTime_f")
    before_n <- nrow(z); before_subjects <- uniqueN(z$RID)
    z <- z[complete.cases(z[, ..required])]
    if (uniqueN(na.omit(z$field_strength)) >= 2L) z <- z[!is.na(field_strength)]
    support <- z[, .(n_dates = uniqueN(cognition_date), span_days = as.numeric(max(cognition_date) - min(cognition_date))), by = RID]
    z <- z[RID %in% support[n_dates >= a5_spec$minimum_memory_dates & span_days >= a5_spec$minimum_memory_span_days, RID]]
    setorder(z, RID, cognition_date)
    unique_keys(z, c("RID", "cognition_date"), "locked amyloid frame")
    constants <- list()
    for (nm in c("LC_z", "SNVTA_z", "M_raw", "G_raw")) {
      if (any(z[, uniqueN(get(nm)), by = RID]$V1 != 1L)) stop("Locked structural/local predictor is not participant-static")
      v <- unique(z, by = "RID")[[nm]]
      constants[[nm]] <- c(mean = mean(v), sd = sd(v), n = length(v))
      if (length(v) && (!is.finite(sd(v)) || sd(v) <= 0)) stop("Nonestimable primary amyloid standardization: ", key, "/", nm)
      target <- switch(nm, M_raw = "M_z", G_raw = "G_z", nm)
      if (nm %in% c("LC_z", "SNVTA_z")) set(z, j = paste0(nm, "_source_locked"), value = z[[nm]])
      set(z, j = target, value = (z[[nm]] - mean(v)) / sd(v))
      scales[[length(scales) + 1L]] <- data.table(encoding = encoding, nucleus = selected_region, variable = nm,
        output_variable = target, center = mean(v), scale = sd(v), n_reference = length(v),
        reference_population = "participant-static maximal complete-case primary R2 frame after original memory support",
        reference_row_hash = a5_object_hash(unique(z[, .(RID)], by = "RID")))
    }
    policy <- a5_covariates(z)
    if (length(policy$included)) z <- z[complete.cases(z[, policy$included, with = FALSE])]
    flow[[length(flow) + 1L]] <- data.table(encoding = encoding, nucleus = selected_region,
      scope = scope, n_before_complete_rows = before_n, n_before_complete_subjects = before_subjects,
      n_rows = nrow(z), n_subjects = uniqueN(z$RID), n_sites = uniqueN(z$site_id))
    frame_key <- paste(scope, key, sep = "__")
    frames[[frame_key]] <- z
    tasks[[k]] <- data.table(scope = scope, encoding = encoding, nucleus = selected_region, frame_key = frame_key,
      analysis_id = paste("f5_amyloid", scope, encoding, tolower(selected_region), sep = "_"), model_variant = a5_spec$model_variant,
      included_covariates = paste(policy$included, collapse = "|"), omitted_covariates = paste(policy$omitted, collapse = "|"),
      field_strength_included = policy$field, n_rows = nrow(z), n_subjects = uniqueN(z$RID),
      frame_hash = a5_object_hash(z), primary_constants_hash = a5_object_hash(constants))
  }
  flow <- rbindlist(flow, fill = TRUE)
  for (i in seq_len(nrow(flow))) {
    k <- paste(flow$scope[i], flow$encoding[i], flow$nucleus[i], sep = "__")
    dd <- frames[[k]]
    set(flow, i = i, j = "n_subjects_own_trajectory_after_index", value = uniqueN(dd[trajectory_last_date > index_date, RID]))
    set(flow, i = i, j = "n_subjects_tpv_after_index", value = uniqueN(dd[last_structural_date > index_date, RID]))
  }
  flow[, `:=`(pattern_weights_reference = "Learned from full completed Figure4 cohort; not index-only prospective measurements",
    information_timing_qualification = "Counts reflect own trajectory/TPV support after original max(MRI,PET) index; learned weights additionally use full cohort")]
  list(tasks = rbindlist(tasks), frames = frames,
    scaling = rbindlist(scales), flow = flow, parent_audit = validated$audit,
    site_mapping = validated$site_mapping)
}

a5_prepare <- function() {
  structure_path <- file.path(output_root(), "private", "figure4", "primary_parameters.csv")
  if (!file.exists(structure_path)) stop("Completed Figure 5 structure dependency unavailable")
  z <- a5_make_frames(input("amyloid_memory.csv"), fread(structure_path, colClasses = "character"))
  z$structure_sha256 <- a5_hash(structure_path); z$input_sha256 <- a5_hash(file.path(Sys.getenv("ADNI_MODEL_DATA_DIR"), "amyloid_memory.csv"))
  z$spec <- a5_spec
  if (file.exists(work_path("amyloid", "prepared.rds"))) {
    old <- load_plan("amyloid")
    expected_tasks <- copy(z$tasks); expected_tasks[, task_id := seq_len(.N)]
    setcolorder(expected_tasks, c("task_id", setdiff(names(expected_tasks), "task_id")))
    if (!identical(old$data$structure_sha256, z$structure_sha256) || !identical(old$data$input_sha256, z$input_sha256) ||
        !identical(old$data$spec, z$spec) ||
        !identical(a5_object_hash(as.data.frame(old$tasks)), a5_object_hash(as.data.frame(expected_tasks))))
      stop("Amyloid completed/prepared dependency, specification or frame identity changed; use a fresh signed run")
  } else write_plan("amyloid", z$tasks, z)
  a5_write(z$scaling, "scaling.tsv"); a5_write(z$flow, "sample_flow.tsv")
  a5_write(z$parent_audit, "parent_validation.tsv")
  a5_write(z$site_mapping, "site_level_mapping.tsv")
  a5_write(rbindlist(lapply(names(z$frames), function(k) cbind(data.table(frame_key = k), z$frames[[k]])), fill = TRUE), "private/model_frames.tsv")
  atomic_rds(a5_spec,work_path("amyloid","analysis_spec.rds"))
}

a5_formula <- function(cr_region, covariates) {
  cr_predictor <- if (cr_region == "LC") "LC_z" else if (cr_region == "SNVTA") "SNVTA_z" else stop("Invalid CR region")
  rhs <- c(paste0("pathology_model * ", cr_predictor, " * Time_c"),
    "pathology_model * M_z * Time_c", "pathology_model * G_z * Time_c", covariates)
  as.formula(paste("outcome_value ~", paste(rhs, collapse = " + "), "+ (1 | site_id) + (1 | SubjectID)"))
}
a5_factors <- function(d) {
  d <- copy(d)
  for (nm in c("Sex", "DX_bl", "field_strength", "EchoTime_f", "site_id", "SubjectID")) {
    preferred <- switch(nm, Sex = c("F", "M"), DX_bl = c("CN", "MCI", "AD"), field_strength = c("1.5T", "3T"), character())
    observed <- sort(unique(na.omit(text(d[[nm]]))))
    set(d, j = nm, value = factor(d[[nm]], levels = c(intersect(preferred, observed), setdiff(observed, preferred))))
  }
  d
}
a5_extract <- function(fit) {
  c <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
  setnames(c, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"), c("estimate", "se", "statistic", "p_value"))
  crit <- qt(0.975, df = pmax(c$df, 1))
  c[, `:=`(ci_low = estimate - crit * se, ci_high = estimate + crit * se,
    uncertainty_method = "Satterthwaite t; historical df lower bound 1")]
  c
}
a5_focal <- function(coefs, cr_region) {
  cr_predictor <- if (cr_region == "LC") "LC_z" else "SNVTA_z"
  roles <- c(local_pathology_time = cr_predictor, plsc_pathology_time = "M_z", tpv_pathology_time = "G_z")
  rbindlist(lapply(names(roles), function(role) {
    want <- sort(c("pathology_model", roles[[role]], "Time_c"))
    hits <- which(vapply(strsplit(coefs$term, ":", fixed = TRUE), function(x) identical(sort(x), want), logical(1)))
    if (length(hits) != 1L) stop("Missing/ambiguous R2 three-way coefficient: ", role)
    cbind(copy(coefs[hits]), data.table(effect_role = role))
  }))
}
a5_fit <- function(spec, prepared) {
  if (!spec$scope %in% "primary") stop("Unsupported amyloid scope")
  d <- a5_factors(prepared$frames[[spec$frame_key]])
  if (!identical(a5_object_hash(prepared$frames[[spec$frame_key]]), spec$frame_hash)) stop("Prepared amyloid frame identity changed")
  covariates <- strsplit(spec$included_covariates, "|", fixed = TRUE)[[1L]]; covariates <- covariates[nzchar(covariates)]
  formula <- a5_formula(spec$nucleus, covariates)
  ns <- uniqueN(d$RID); n0 <- uniqueN(d[pathology_model == 0, RID]); n1 <- uniqueN(d[pathology_model == 1, RID])
  reason <- if (ns < a5_spec$minimum_subjects) "fewer_than_20_subjects" else if (nrow(d) < a5_spec$minimum_observations) "fewer_than_40_observations" else ""
  if (!nzchar(reason) && (nlevels(d$site_id) < 2L || nlevels(d$SubjectID) < 2L)) reason <- "fewer_than_two_mandatory_group_levels"
  if (!nzchar(reason) && (!is.finite(sd(d$pathology_model)) || sd(d$pathology_model) <= 0)) reason <- "zero_pathology_variance"
  if (!nzchar(reason) && spec$encoding == "binary" && min(n0, n1) < a5_spec$minimum_binary_group_subjects) reason <- "fewer_than_10_subjects_in_binary_group"
  X <- tryCatch(model.matrix(a5_nobars(formula), d), error = identity)
  if (!nzchar(reason) && inherits(X, "error")) reason <- paste0("model_matrix_error:", conditionMessage(X))
  if (!nzchar(reason) && qr(X)$rank < ncol(X)) reason <- "rank_deficient_fixed_effect_design"
  meta <- cbind(spec, data.table(formula = paste(deparse(formula), collapse = " "), random_effects = a5_spec$random_effects,
    REML = TRUE, optimizer = "bobyqa", fallback = "none", n_sites = nlevels(d$site_id),
    n_binary_negative = if (spec$encoding == "binary") n0 else NA_integer_,
    n_binary_positive = if (spec$encoding == "binary") n1 else NA_integer_,
    model_status = if (nzchar(reason)) "NOT_ESTIMABLE" else "ESTIMATED", reason = reason,
    design_rank = if (inherits(X, "error")) NA_integer_ else qr(X)$rank,
    design_columns = if (inherits(X, "error")) NA_integer_ else ncol(X),
    condition_number = if (inherits(X, "error")) NA_real_ else kappa(X)))
  if (nzchar(reason)) return(list(results = meta, diagnostics = copy(meta), model = NULL))
  caught <- NULL
  fit_result <- tryCatch(capture_fit(lmerTest::lmer(formula, data = d, REML = TRUE,
    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000), check.rankX = "stop.deficient"))),
    error = function(e) { caught <<- e; NULL })
  if (!is.null(caught)) {
    meta[, `:=`(model_status = "NOT_ESTIMABLE", reason = paste0("lmer_error:", conditionMessage(caught)))]
    return(list(results = meta, diagnostics = copy(meta), model = NULL))
  }
  fit <- fit_result$fit
  if (!setequal(names(lme4::getME(fit, "flist")), c("site_id", "SubjectID")) || nrow(model.frame(fit)) != nrow(d)) stop("R2 fit changed random terms or locked rows")
  coefs <- a5_extract(fit); focal <- a5_focal(coefs, spec$nucleus)
  conv <- fit@optinfo$conv$lme4$messages
  vc <- as.data.table(as.data.frame(lme4::VarCorr(fit)))
  info <- data.table(singular = lme4::isSingular(fit, tol = 1e-4),
    convergence_messages = paste(conv, collapse = " | "), warnings = paste(fit_result$warnings, collapse = " | "),
    optimizer_code = paste(unlist(fit@optinfo$conv$opt), collapse = "|"),
    objective = unname(lme4::getME(fit, "devcomp")$cmp[["REML"]]), logLik = as.numeric(logLik(fit)), sigma = sigma(fit),
    same_frame_verified = TRUE, participant_time_random_slope = FALSE)
  info[, diagnostic_status := if (singular || nzchar(convergence_messages) || nzchar(warnings)) "ESTIMATED_REVIEW" else "ESTIMATED"]
  list(results = meta, coefficients = cbind(spec, coefs), focal = cbind(spec, focal),
    diagnostics = cbind(meta, info), variance_components = cbind(spec, vc),
    model = fit,
    model_frame_hash = a5_object_hash(model.frame(fit)),
    contrast_levels = lapply(d[, names(d)[vapply(d, is.factor, logical(1))], with = FALSE], levels))
}

a5_adjust <- function(focal, statuses, partial = FALSE) {
  focal <- copy(focal)
  if (any(!focal$scope %in% "primary") || any(!statuses$scope %in% "primary"))
    stop("Unsupported amyloid correction-family scope")
  if (!nrow(focal)) return(focal)
  focal[, `:=`(q_bh = NA_real_, family_n_planned = NA_integer_, correction_status = "not_applicable_TPV")]
  scope <- "primary"
  for (encoding in c("binary", "centiloid")) for (role in c("local_pathology_time", "plsc_pathology_time")) {
    ii <- which(focal$scope == scope & focal$encoding == encoding & focal$effect_role == role)
    selected_status <- which(statuses$scope == scope & statuses$encoding == encoding)
    ss <- statuses[selected_status]
    accounted <- nrow(ss) == 2L && setequal(ss$nucleus, c("LC", "SNVTA")) && all(ss$model_status %in% c("ESTIMATED", "NOT_ESTIMABLE"))
    if (length(ii)) {
      set(focal, i = ii, j = "family_n_planned", value = 2L)
      set(focal, i = ii, j = "correction_status", value = if (accounted && !partial) "final_all_planned_statuses_accounted" else "incomplete_not_final")
      if (accounted && !partial) {
        if (any(!is.finite(focal$p_value[ii]) | focal$p_value[ii] < 0 | focal$p_value[ii] > 1)) stop("Invalid focal probability in complete amyloid family")
        set(focal, i = ii, j = "q_bh", value = p.adjust(focal$p_value[ii], "BH", n = 2L))
      }
    }
  }
  focal
}


a5_finalize <- function() {
  tasks<-collect_tasks("amyloid");statuses<-bind_component(tasks,"results")
  if(nrow(statuses)!=4L)stop("Incomplete four-model memory scope")
  focal<-a5_adjust(bind_component(tasks,"focal"),statuses,FALSE)
  for(component in c("coefficients","variance_components"))a5_write(bind_component(tasks,component),paste0(component,".tsv"))
  # Preserve planned unavailable focal rows after the unchanged correction step.
  padding<-lapply(seq_len(nrow(statuses)),function(i){
    st<-statuses[i];roles<-c("local_pathology_time","plsc_pathology_time","tpv_pathology_time")
    have<-if(nrow(focal))focal[task_id==st$task_id]$effect_role else character()
    absent<-setdiff(roles,have)
    if(!length(absent))return(NULL)
    cbind(st[rep(1L,length(absent))],data.table(effect_role=absent,term=NA_character_,estimate=NA_real_,se=NA_real_,df=NA_real_,statistic=NA_real_,p_value=NA_real_,q_bh=NA_real_,ci_low=NA_real_,ci_high=NA_real_))
  })
  focal<-rbindlist(c(list(focal),padding),fill=TRUE)
  a5_write(focal,"focal_effects.tsv")
  record_diagnostics("memory",rbindlist(list(bind_component(tasks,"diagnostics"),statuses,focal),fill=TRUE))
}
if(sys.nframe()==0L){a<-cli();if(a$action=="prepare")a5_prepare()else if(a$action=="task")run_task("amyloid",a$task,a5_fit)else if(a$action=="finalize")a5_finalize()else stop("Expected prepare/task/finalize")}
