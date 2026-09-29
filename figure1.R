#!/usr/bin/env Rscript
# Technical CR adjustment, raw-CR repeatability and gross hypopigmentation models.
.local_file <- function() {
  ss <- vapply(sys.frames(), function(x) if(is.null(x$ofile)) "" else as.character(x$ofile), character(1))
  z <- ss[basename(ss)=="figure1.R"]
  if(length(z)) return(normalizePath(tail(z,1L)))
  z <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])
  normalizePath(z)
}
source(file.path(dirname(.local_file()), "common.R"))
suppressPackageStartupMessages(library(lme4))

# Recorded order columns define the prepared acquisition-series frame.
icc_from_scans <- function(s) {
  need(s, c("rid","participant_id","scan_id","scan_set_id","mri_date","diagnosis_mri","field_strength"))
  unique_keys(s,"scan_set_id","portable scan registry")
  out <- rbindlist(lapply(c("LC","SNVTA"), function(rg) rbindlist(lapply(c("L","R"), function(hm) {
    member <- paste("icc",rg,hm,"order",sep="_")
    value <- paste(tolower(rg),tolower(hm),sep="_")
    need(s,c(member,value)); ord <- number(s[[member]])
    if(any(is.finite(ord) & (ord <= 0 | ord != as.integer(ord)))) stop("Invalid prepared ICC row order")
    keep <- which(is.finite(ord))
    data.table(rid=s$rid[keep],participant_id=s$participant_id[keep],scan_id=s$scan_id[keep],
      scan_set_id=s$scan_set_id[keep],mri_date=s$mri_date[keep],diagnosis_mri=s$diagnosis_mri[keep],
      field_strength=s$field_strength[keep],nucleus=rg,hemisphere=hm,cr=s[[value]][keep],.source_order=ord[keep])
  }))))
  unique_keys(out,".source_order","prepared ICC order")
  if(!identical(sort(out$.source_order),as.numeric(seq_len(nrow(out))))) stop("Prepared ICC order is not contiguous")
  setorder(out,.source_order); out[, .source_order := NULL]; out
}

technical_prepare <- function() {
  s <- input("scans.csv")
  need(s, c("rid", "participant_id", "scan_set_id", "scan_id", "field_strength", "qc_pass",
            "lc_l", "lc_r", "snvta_l", "snvta_r", "efc", "fber", "echo_time_category", "site_id"))
  unique_keys(s, "scan_set_id", "scan registry")
  d <- rbindlist(lapply(c("LC", "SNVTA"), function(rg) {
    rbindlist(lapply(c("L", "R"), function(hm) {
      y <- number(s[[paste0(tolower(rg), "_", tolower(hm))]])
      keep <- truth(s$qc_pass) & is.finite(y)
      data.table(rid = s$rid[keep], participant_id = s$participant_id[keep],
                 scan_set_id = s$scan_set_id[keep], scan_id = s$scan_id[keep],
                 nucleus = rg, hemisphere = hm, field_strength = field(s$field_strength[keep]),
                 cr = y[keep], efc = number(s$efc[keep]), fber = number(s$fber[keep]),
                 echo = text(s$echo_time_category[keep]), site_id = text(s$site_id[keep]))
    }))
  }))
  unique_keys(d, c("scan_set_id", "nucleus", "hemisphere"), "technical observations")
  if (anyNA(d$field_strength) || anyNA(d$participant_id)) stop("Unresolved technical-model keys")
  setorder(d, nucleus, field_strength, hemisphere, participant_id, scan_id)
  tasks <- CJ(nucleus = c("LC", "SNVTA"), field_strength = c("1.5T", "3T"), hemisphere = c("L", "R"))
  write_plan("technical", tasks, d)
}

technical_fit <- function(spec, all) {
  d <- copy(all[nucleus == spec$nucleus & field_strength == spec$field_strength & hemisphere == spec$hemisphere])
  if (nrow(d) < 100L || uniqueN(d$participant_id) < 40L) stop("Technical stratum does not meet original 100-row/40-participant support")
  # Reference statistics precede model-specific complete-case selection.
  efc_mean <- mean(d$efc, na.rm = TRUE); efc_sd <- sd(d$efc, na.rm = TRUE)
  fber_mean <- mean(d$fber, na.rm = TRUE); fber_sd <- sd(d$fber, na.rm = TRUE)
  raw_anchor <- mean(d$cr)
  d[, `:=`(efc_z = finite_z(efc), fber_z = finite_z(fber), echo_f = factor(echo),
            site = factor(site_id), participant = factor(participant_id))]
  technical_terms <- function(x, fber_min) {
    terms <- character()
    if (sum(is.finite(x$efc_z)) >= 20L && is.finite(sd(x$efc_z, na.rm = TRUE)) && sd(x$efc_z, na.rm = TRUE) > 0) terms <- c(terms, "efc_z")
    if (sum(is.finite(x$fber_z)) >= fber_min && is.finite(sd(x$fber_z, na.rm = TRUE)) && sd(x$fber_z, na.rm = TRUE) > 0) terms <- c(terms, "fber_z")
    if (all(c("efc_z", "fber_z") %in% terms)) terms <- c(terms, "efc_z:fber_z")
    if (nlevels(droplevels(x$echo_f)) > 1L) terms <- c(terms, "echo_f")
    terms
  }
  display_terms <- technical_terms(d, 20L)
  if (!length(display_terms) || nlevels(d$site) < 2L) stop("Technical display model cannot retain its mandatory site effect")
  subject_display <- nlevels(d$participant) > 1L && anyDuplicated(d$participant) > 0L
  display_formula <- as.formula(paste("cr ~", paste(c(display_terms, "(1 | site)", if (subject_display) "(1 | participant)"), collapse = " + ")))
  if (!all(complete.cases(d[, all.vars(display_formula), with = FALSE]))) stop("Technical display inputs are incomplete; no raw fallback is permitted")
  fit_display <- capture_fit(lme4::lmer(display_formula, data = d, REML = TRUE, na.action = na.fail,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000L), check.rankX = "stop.deficient", check.conv.singular = "ignore")))
  check_convergence(fit_display$fit)
  pred_display <- as.numeric(predict(fit_display$fit, newdata = d, re.form = ~(1 | site), allow.new.levels = FALSE))
  if (any(!is.finite(pred_display))) stop("Nonfinite technical display predictions")

  inference_terms <- technical_terms(d, 50L)
  if (!length(inference_terms)) inference_terms <- "1"
  repeats <- d[, .N, by = participant_id]
  subject_inference <- uniqueN(d$participant_id) >= 40L && sum(repeats$N >= 2L) >= 10L && sum(pmax(0L, repeats$N - 1L)) >= 20L
  inference_formula <- as.formula(paste("cr ~", paste(c(inference_terms, "(1 | site)", if (subject_inference) "(1 | participant)"), collapse = " + ")))
  reuse <- identical(deparse(inference_formula), deparse(display_formula))
  fit_inference <- if (reuse) fit_display else capture_fit(lme4::lmer(inference_formula, data = d, REML = TRUE,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000L)), na.action = na.fail))
  check_convergence(fit_inference$fit)
  policy <- "site_and_participant"
  if (subject_inference && lme4::isSingular(fit_inference$fit)) {
    vc <- as.data.frame(VarCorr(fit_inference$fit)); vsub <- vc$vcov[vc$grp == "participant"][1L]
    if (!is.finite(vsub) || vsub <= 1e-8) {
      # A singular fit with negligible participant variance permits site-only inference.
      inference_formula <- as.formula(paste("cr ~", paste(c(inference_terms, "(1 | site)"), collapse = " + ")))
      fit_inference <- capture_fit(lme4::lmer(inference_formula, data = d, REML = TRUE,
        control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000L)), na.action = na.fail))
      check_convergence(fit_inference$fit); reuse <- FALSE; policy <- "site_only_original_zero_subject_variance_fallback"
    }
  }
  if (!subject_inference) policy <- "site_only_original_support_rule"
  pred <- as.numeric(predict(fit_inference$fit, newdata = d, re.form = ~(1 | site), allow.new.levels = FALSE))
  d[, `:=`(cr_residual = cr - pred, display_residual = cr - pred_display,
            cr_adjusted_anchored = cr - pred_display + raw_anchor,
            raw_stratum_mean = raw_anchor, fixed_site_removed = pred,
            display_fixed_site_removed = pred_display)]
  if (any(!is.finite(d$cr_residual)) || any(!is.finite(d$cr_adjusted_anchored))) stop("Nonfinite adjusted CR")
  values <- d[, .(rid, participant_id, scan_id, scan_set_id, nucleus, field_strength, hemisphere,
                   cr_raw = cr, cr_residual, display_residual, cr_adjusted_anchored, raw_stratum_mean,
                   fixed_site_removed, display_fixed_site_removed)]
  model_rows <- rbindlist(lapply(c("inference", "display"), function(role) {
    obj <- if (role == "inference") fit_inference else fit_display
    coeff <- coef_lmm(obj$fit)
    coeff[, `:=`(nucleus = spec$nucleus, field_strength = spec$field_strength, hemisphere = spec$hemisphere,
                  role = role, n_participants = uniqueN(d$rid), n_sites = uniqueN(d$site_id),
                  n_rows = nrow(d), raw_stratum_mean = raw_anchor, efc_mean = efc_mean, efc_sd = efc_sd,
                  fber_mean = fber_mean, fber_sd = fber_sd, fit_reused = reuse, participant_effect_subtracted = FALSE,
                  random_effect_policy = if (role == "inference") policy else "original_display_subject_support",
                  formula = paste(deparse(formula(obj$fit)), collapse = " "),
                  singular = isSingular(obj$fit), status = "estimated", warnings = paste(obj$warnings, collapse = " | "))]
    coeff
  }), fill = TRUE)
  list(values = values, results = model_rows,
       fits = list(inference = fit_inference$fit, display = if (reuse) NULL else fit_display$fit),
       frame = d, interpretation = "Remove fixed technical terms (including intercept) and site; retain participant variation; raw stratum mean anchors display only")
}
icc_prepare <- function() {
  d <- icc_from_scans(input("scans.csv"))
  need(d, c("rid", "participant_id", "scan_id", "scan_set_id", "nucleus", "hemisphere", "field_strength",
            "cr", "diagnosis_mri", "mri_date"))
  unique_keys(d, c("scan_set_id", "nucleus", "hemisphere"), "ICC observations")
  d[, cr := number(cr)]
  d[, field_strength := field(field_strength)]
  if (anyNA(d$cr) || anyNA(d$participant_id) || anyNA(d$field_strength) ||
      any(is.na(d$diagnosis_mri) | d$diagnosis_mri == "AD")) stop("Prepared ICC eligibility does not satisfy the declared non-AD/finite rule")
  counts <- d[, .N, by = .(nucleus, field_strength, hemisphere, participant_id)]
  if (any(counts$N < 2L)) stop("Prepared ICC input contains a single-observation participant stratum")
  setorder(d, nucleus, field_strength, hemisphere, participant_id, mri_date, scan_id)
  write_plan("icc", CJ(nucleus = c("LC", "SNVTA"), field_strength = c("1.5T", "3T"), hemisphere = c("L", "R")), d)
}
variance_icc <- function(fit) {
  v <- as.data.frame(VarCorr(fit)); a <- v$vcov[v$grp == "participant"][1L]; b <- v$vcov[v$grp == "Residual"][1L]
  c(icc = a / (a + b), variance_participant = a, variance_residual = b)
}
icc_variants <- function() "primary"
icc_fit <- function(spec, all) {
  d <- copy(all[nucleus == spec$nucleus & field_strength == spec$field_strength & hemisphere == spec$hemisphere])
  scientific_assert(nrow(d) >= 8L && uniqueN(d$rid) >= 4L, "Insufficient repeated observations for ICC")
  d[, participant := factor(participant_id)]
  models <- list(); results <- list(); draws <- list()
  for (variant in icc_variants()) {
    f <- cr ~ 1 + (1 | participant)
    fit <- capture_fit(lmer(f, data = d, REML = TRUE, na.action = na.fail,
                           control = lmerControl(check.conv.singular = "ignore", check.conv.grad = "ignore")))
    check_convergence(fit$fit)
    v <- variance_icc(fit$fit)
    if (!is.finite(v[["icc"]])) stop("Nonfinite ICC variance ratio")
    B <- resample_n(); seed <- 20260611L
    set.seed(seed)
    boot <- bootMer(fit$fit, FUN = function(m) unname(variance_icc(m)[["icc"]]), nsim = B,
                    use.u = FALSE, type = "parametric", parallel = "no", .progress = "none")
    vals <- as.numeric(boot$t[, 1L]); good <- vals[is.finite(vals)]
    if (length(good) < 10L) stop("Fewer than ten finite bootstrap ICC estimates")
    ci <- quantile(good, c(.025, .975), names = FALSE, type = 6)
    r <- cbind(spec, fit_info(fit$fit, fit$warnings))
    r[, `:=`(variant = variant, icc = unname(v[["icc"]]), ci_low = ci[1], ci_high = ci[2],
              variance_participant = unname(v[["variance_participant"]]), variance_residual = unname(v[["variance_residual"]]),
              within_subject_sd = sqrt(v[["variance_residual"]]), n_participants = uniqueN(d$rid), n_scans = uniqueN(d$scan_set_id),
              n_observations = nrow(d), mean_repeats = nrow(d) / uniqueN(d$rid), seed = seed,
              bootstrap_requested = B, bootstrap_valid = length(good), bootstrap_failed = B - length(good),
              bootstrap_failure_messages = paste(capture.output(str(attr(boot, "boot.fail.msgs"))), collapse = " "),
              quantile_type = 6L, use_u = FALSE, time_method = "none",
              status = if (length(good) == B) "estimated" else "estimated_incomplete_bootstrap",
              inference_scale = "raw retained CR; no technical residualization")]
    results[[variant]] <- r; models[[variant]] <- fit$fit; draws[[variant]] <- vals
  }
  list(results = rbindlist(results, fill = TRUE), fits = models, bootstrap_draws = draws,
       frame = d[, .(rid, participant_id, scan_set_id, cr)])
}

# Pathology-specific echo-time coding and four-decimal CR remain separate from
# the shared three-decimal technical-CR adjustment used by other analyses.
pathology_diagnosis <- function(x) {
  value <- text(x); code <- suppressWarnings(as.numeric(value))
  label <- toupper(gsub("[^A-Za-z0-9]+", "_", value))
  code[!is.finite(code) & label %in% c("CN", "NORMAL", "COGNITIVELY_NORMAL")] <- 1
  code[!is.finite(code) & grepl("MCI", label)] <- 2
  code[!is.finite(code) & grepl("ALZ|^AD$|AD_DEMENTIA", label)] <- 3
  code[!is.finite(code) & grepl("NON.?AD|OTHER_DEMENTIA|DEMENTIA_NON", label)] <- 4
  # The bare Other diagnosis label is unavailable under the pathology definition.
  factor(fifelse(code == 1, "CN", fifelse(code == 2, "MCI",
    fifelse(code == 3, "AD", fifelse(code == 4, "nonAD_dementia", NA_character_)))))
}
pathology_echo <- function(x) {
  value <- number(x); value[is.finite(value) & value > 1] <- value[is.finite(value) & value > 1] / 1000
  factor(ifelse(is.finite(value), sprintf("%.4fs", round(value, 4)), NA_character_))
}
pathology_technical_rows <- function(s) {
  need(s, c("rid", "participant_id", "scan_set_id", "scan_id", "field_strength", "qc_pass",
    "lc_l", "lc_r", "snvta_l", "snvta_r", "efc", "fber", "echo_time", "site_id"))
  unique_keys(s, "scan_set_id", "pathology scan registry")
  rows <- rbindlist(lapply(c("LC", "SNVTA"), function(rg) rbindlist(lapply(c("L", "R"), function(hm) {
    y <- number(s[[paste(tolower(rg), tolower(hm), sep = "_")]])
    keep <- which(truth(s$qc_pass) & is.finite(y))
    data.table(rid = text(s$rid[keep]), participant_id = text(s$participant_id[keep]),
      scan_set_id = text(s$scan_set_id[keep]), scan_id = text(s$scan_id[keep]),
      nucleus = rg, hemisphere = hm, field_strength = field(s$field_strength[keep]),
      CR = y[keep], EFC = number(s$efc[keep]), FBER = number(s$fber[keep]),
      EchoTime_raw = text(s$echo_time[keep]), site_id = text(s$site_id[keep]))
  }))))
  # Unpadded integer-text site labels preserve site groups and factor-level order.
  site <- rows$site_id; numeric_site <- !is.na(site) & grepl("^[0-9]+$", site)
  site[numeric_site] <- sub("^0+(?=[0-9])", "", site[numeric_site], perl = TRUE)
  pairs <- unique(data.table(original = rows$site_id, normalized = site))
  if (anyDuplicated(pairs$normalized)) stop("Pathology site normalization merges source groups")
  rows[, `:=`(site_id = site, EchoTime_f = pathology_echo(EchoTime_raw))]
  unique_keys(rows, c("scan_set_id", "nucleus", "hemisphere"), "pathology technical observations")
  if (anyNA(rows[, .(rid, participant_id, field_strength)])) stop("Unresolved pathology technical keys")
  # Use lexical participant/scan order within each stratum.
  setorder(rows, nucleus, field_strength, hemisphere, rid, scan_id)
  rows
}
pathology_technical_fit <- function(spec, all) {
  d <- copy(all[nucleus == spec$nucleus & field_strength == spec$field_strength & hemisphere == spec$hemisphere])
  raw_anchor <- mean(d$CR, na.rm = TRUE)
  refs <- list(efc_mean = mean(d$EFC, na.rm = TRUE), efc_sd = sd(d$EFC, na.rm = TRUE),
    fber_mean = mean(d$FBER, na.rm = TRUE), fber_sd = sd(d$FBER, na.rm = TRUE), raw_stratum_mean = raw_anchor)
  d[, `:=`(EFC_Z_fs = finite_z(EFC), FBER_Z_fs = finite_z(FBER),
    EchoTime_f = droplevels(factor(EchoTime_f)), site_id = factor(site_id), SubjectID = factor(participant_id))]
  fixed_terms <- function(x, echo_column, fber_min) {
    terms <- character()
    if (sum(is.finite(x$EFC_Z_fs)) >= 20L && is.finite(sd(x$EFC_Z_fs, na.rm = TRUE)) && sd(x$EFC_Z_fs, na.rm = TRUE) > 0) terms <- c(terms, "EFC_Z_fs")
    if (sum(is.finite(x$FBER_Z_fs)) >= fber_min && is.finite(sd(x$FBER_Z_fs, na.rm = TRUE)) && sd(x$FBER_Z_fs, na.rm = TRUE) > 0) terms <- c(terms, "FBER_Z_fs")
    if (all(c("EFC_Z_fs", "FBER_Z_fs") %in% terms)) terms <- c(terms, "EFC_Z_fs:FBER_Z_fs")
    if (nlevels(droplevels(factor(x[[echo_column]]))) > 1L) terms <- c(terms, echo_column)
    terms
  }
  subject_supported <- function(x) {
    repeats <- x[!is.na(participant_id), .N, by = participant_id]
    nrow(repeats) >= 40L && sum(repeats$N >= 2L) >= 10L && sum(pmax(0L, repeats$N - 1L)) >= 20L
  }
  # Raw-echo anchored values use the entire QC stratum before pathology selection.
  display <- copy(d)
  display[, `:=`(EchoTime_adjustment_f = factor(EchoTime_raw),
    site_adjustment_f = factor(site_id), subject_adjustment_f = factor(SubjectID))]
  display_terms <- fixed_terms(display, "EchoTime_adjustment_f", 20L)
  if (!length(display_terms) || nlevels(display$site_adjustment_f) < 2L) stop("Pathology display requires technical terms and its mandatory site effect")
  subject_display <- nlevels(display$subject_adjustment_f) > 1L && anyDuplicated(display$subject_adjustment_f) > 0L
  display_formula <- as.formula(paste("CR ~", paste(c(display_terms, "(1 | site_adjustment_f)", if (subject_display) "(1 | subject_adjustment_f)"), collapse = " + ")))
  fit_display <- capture_fit(lme4::lmer(display_formula, data = display, REML = TRUE, na.action = na.exclude,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000L),
      check.rankX = "stop.deficient", check.conv.singular = "ignore", check.conv.grad = "ignore")))
  pred_display <- as.numeric(predict(fit_display$fit, newdata = display, re.form = ~(1 | site_adjustment_f), allow.new.levels = FALSE))
  if (length(pred_display) != nrow(display) || any(!is.finite(pred_display)) || !is.finite(raw_anchor)) stop("Nonfinite pathology display adjustment")

  requested <- fixed_terms(d, "EchoTime_f", 50L); if (!length(requested)) requested <- "1"
  vars <- unique(c("CR", all.vars(as.formula(paste("~", paste(requested, collapse = "+")))),
    if (nlevels(d$site_id) >= 2L) "site_id", if (subject_supported(d)) "SubjectID"))
  formal <- droplevels(d[complete.cases(d[, ..vars])])
  terms <- fixed_terms(formal, "EchoTime_f", 50L); terms <- terms[terms %in% requested]
  if (!length(terms)) terms <- "1"
  if (nrow(formal) < 100L || uniqueN(formal$participant_id) < 40L) stop("Pathology technical stratum does not meet original 100-row/40-participant support")
  if (nlevels(formal$site_id) < 2L) stop("Pathology formal residual requires its mandatory site effect")
  supported <- subject_supported(formal); notes <- character(); fit_formal <- NULL
  for (with_subject in if (supported) c(TRUE, FALSE) else FALSE) {
    f <- as.formula(paste("CR ~", paste(c(terms, "(1 | site_id)", if (with_subject) "(1 | SubjectID)"), collapse = " + ")))
    attempt <- tryCatch(capture_fit(lme4::lmer(f, data = formal, REML = TRUE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000L)))), error = identity)
    if (inherits(attempt, "error")) { notes <- c(notes, conditionMessage(attempt)); next }
    if (with_subject && isSingular(attempt$fit)) {
      vc <- as.data.frame(VarCorr(attempt$fit)); vs <- vc$vcov[vc$grp == "SubjectID"][1L]
      if (!is.finite(vs) || vs <= 1e-8) { notes <- c(notes, "original_zero_subject_variance_site_only_fallback"); next }
    }
    fit_formal <- attempt; break
  }
  if (is.null(fit_formal)) stop("Pathology technical fit failed: ", paste(notes, collapse = " | "))
  pred <- as.numeric(predict(fit_formal$fit, newdata = formal, re.form = ~(1 | site_id), allow.new.levels = TRUE))
  if (any(!is.finite(pred))) stop("Nonfinite pathology formal residual")
  values <- display[, .(rid, participant_id, scan_id, scan_set_id, nucleus, field_strength, hemisphere,
    cr_raw = CR, display_residual = CR - pred_display, cr_adjusted_anchored = CR - pred_display + raw_anchor,
    raw_stratum_mean = raw_anchor, display_fixed_site_removed = pred_display)]
  residuals <- formal[, .(scan_set_id, cr_residual = CR - pred, fixed_site_removed = pred)]
  values <- join_one(values, residuals, "scan_set_id", "pathology formal technical residual")
  audit <- rbindlist(lapply(c("inference", "display"), function(role) {
    obj <- if (role == "inference") fit_formal else fit_display
    frame <- if (role == "inference") formal else display
    cbind(spec, data.table(role = role, formula = paste(deparse(formula(obj$fit)), collapse = " "),
      n_rows = nrow(frame), n_participants = uniqueN(frame$participant_id), n_sites = uniqueN(frame$site_id),
      echo_levels = nlevels(frame[[if (role == "inference") "EchoTime_f" else "EchoTime_adjustment_f"]]),
      echo_policy = if (role == "inference") "seconds_rounded_four_decimals" else "raw_echo_factor",
      participant_effect_subtracted = FALSE, singular = isSingular(obj$fit),
      warnings = paste(obj$warnings, collapse = " | "), optimizer_messages = paste(obj$fit@optinfo$conv$lme4$messages, collapse = " | "),
      original_fallback_notes = if (role == "inference") paste(notes, collapse = " | ") else "",
      status = "estimated"), as.data.table(refs))
  }), fill = TRUE)
  list(values = values, audit = audit, fits = list(inference = fit_formal$fit, display = fit_display$fit),
    frames = list(inference = formal, display = display), references = refs)
}

pathology_prepare <- function() {
  g <- input("gross_pathology.csv")
  # A signed completed plan is reused without recalculating technical fits.
  if (file.exists(work_path("pathology", "prepared.rds"))) return(invisible(load_plan("pathology")))
  rows <- pathology_technical_rows(input("scans.csv"))
  strata <- CJ(nucleus = c("LC", "SNVTA"), field_strength = c("1.5T", "3T"), hemisphere = c("L", "R"))
  technical <- lapply(seq_len(nrow(strata)), function(i) pathology_technical_fit(strata[i], rows))
  v <- rbindlist(lapply(technical, `[[`, "values"))
  need(g, c("rid", "nucleus", "scan_set_id", "mri_date", "predeath_order", "years_to_death", "death_age",
            "sex", "diagnosis_mri", "field_strength", "lc_hypopigmentation_any", "sn_hypopigmentation_any"))
  unique_keys(g, c("scan_set_id", "nucleus"), "predeath candidates")
  for (n in c("cr_residual", "cr_adjusted_anchored")) set(v, j = n, value = number(v[[n]]))
  unique_keys(v, c("scan_set_id", "nucleus", "hemisphere"), "technical CR")
  # Require both adjusted hemispheres, from this same scan. Never use a unilateral fallback.
  b <- v[, .(bilateral_residual = if (setequal(hemisphere, c("L", "R")) && .N == 2L && all(is.finite(cr_residual))) mean(cr_residual) else NA_real_,
              bilateral_anchored = if (setequal(hemisphere, c("L", "R")) && .N == 2L && all(is.finite(cr_adjusted_anchored))) mean(cr_adjusted_anchored) else NA_real_),
         by = .(scan_set_id, nucleus)]
  g <- join_one(g, b, c("scan_set_id", "nucleus"), "adjusted bilateral MRI")
  for (n in c("predeath_order", "years_to_death", "death_age", "lc_hypopigmentation_any", "sn_hypopigmentation_any")) set(g, j = n, value = number(g[[n]]))
  g <- g[is.finite(bilateral_residual) & is.finite(bilateral_anchored) & is.finite(predeath_order)]
  setorder(g, nucleus, rid, predeath_order, scan_set_id)
  selected <- unique(g, by = c("rid", "nucleus"))
  # Selection precedes pathology outcome/covariate availability and is shared by all four gross models.
  unique_keys(selected, c("rid", "nucleus"), "selected predeath MRI")
  write_plan("pathology", CJ(kind="gross", nucleus=c("LC","SNVTA"), rating=c("LC","SN")), list(selected=selected, pathology_technical=technical))
}

gross <- function(spec, data) {
  d <- copy(data$selected[nucleus == spec$nucleus])
  scientific_assert(nrow(d) >= 12L, "Fewer than 12 adjusted predeath MRI participants")
  d[, `:=`(outcome = bilateral_residual, display_cr = bilateral_anchored,
            pathology_var = if (spec$rating == "LC") lc_hypopigmentation_any else sn_hypopigmentation_any,
            age_at_death_c = center(death_age), years_to_death_c = center(years_to_death),
            Sex_f = factor(text(sex)), DX_scan_f = pathology_diagnosis(diagnosis_mri), field_strength_f = factor(field(field_strength)))]
  terms <- c("pathology_var", "age_at_death_c", "years_to_death_c")
  base <- c("outcome", "display_cr", terms)
  d <- d[complete.cases(d[, ..base])]
  dropped <- character()
  for (term in c("Sex_f", "DX_scan_f", "field_strength_f")) {
    counts <- table(droplevels(d[[term]]))
    if (length(counts) >= 2L && all(counts >= 3L)) terms <- c(terms, term) else dropped <- c(dropped, term)
  }
  vars <- c("outcome", terms); d <- droplevels(d[complete.cases(d[, ..vars])])
  scientific_assert(nrow(d) >= 12L && all(sort(unique(d$pathology_var)) == c(0, 1)) && min(table(d$pathology_var)) >= 3L,
                    "Insufficient gross-pathology sample or binary-group support")
  f <- as.formula(paste("outcome ~", paste(terms, collapse = " + ")))
  scientific_assert(qr(model.matrix(f, d))$rank == ncol(model.matrix(f, d)), "Rank-deficient gross model")
  fit <- lm(f, d); co <- hc3(fit)
  co[, `:=`(record_type = "gross_coefficient", nucleus = spec$nucleus, rating = spec$rating,
             n = nrow(d), n_absent = sum(d$pathology_var == 0), n_present = sum(d$pathology_var == 1),
             formula = paste(deparse(f), collapse = " "), uncertainty = "HC3_t",
             dropped_optional_terms = paste(dropped, collapse = ";"), status = "estimated", focal = term == "pathology_var")]
  list(results = co, fit = fit, frame = d)
}

technical_finalize <- function() {
  x <- collect_tasks("technical"); v <- bind_component(x,"values"); m <- bind_component(x,"results")
  if(any(vapply(x,function(a)is.null(a$values),logical(1)))) stop("All eight technical strata are required")
  unique_keys(v,c("scan_set_id","nucleus","hemisphere"),"technical output")
  setorder(v,nucleus,field_strength,hemisphere,participant_id,scan_id)
  atomic_csv(v,"private/shared/technical_cr.csv")
  record_diagnostics("technical", m)
  publish_results(m,"shared/technical_coefficients.csv")
  refs <- unique(m[,.(nucleus,field_strength,hemisphere,role,raw_stratum_mean,efc_mean,efc_sd,fber_mean,fber_sd,
    n_participants,n_sites,n_rows,participant_effect_subtracted,random_effect_policy,formula)])
  atomic_csv(refs,"shared/technical_references.csv")
}

icc_finalize <- function() {
  x <- collect_tasks("icc"); r <- bind_component(x,"results")
  if(any(!r$variant %in% icc_variants())) stop("Only primary raw ICC is supported")
  record_diagnostics("icc", r)
  publish_results(r,"figure1/icc.csv")

}

pathology_finalize <- function() {
  x <- collect_tasks("pathology")
  record_diagnostics("pathology", bind_component(x,"results"))
  publish_results(bind_component(x,"results"),"figure1/gross_coefficients.csv")
}

if(sys.nframe()==0L) {
  a <- cli(); unit <- a$unit
  if(is.null(unit) || !unit %in% c("technical","icc","pathology")) stop("--unit must be technical, icc, or pathology")
  if(a$action=="prepare") get(paste0(unit,"_prepare"))()
  else if(a$action=="task") run_task(unit,a$task,switch(unit,technical=technical_fit,icc=icc_fit,pathology=gross))
  else if(a$action=="finalize") get(paste0(unit,"_finalize"))()
  else stop("Unknown action")
}
