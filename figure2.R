#!/usr/bin/env Rscript
# Main Figure 2 raw-hemisphere CR models; one fit supplies all main outputs.
.local_file <- function() {
  ss <- vapply(sys.frames(), function(x) if(is.null(x$ofile)) "" else as.character(x$ofile), character(1))
  z <- ss[basename(ss)=="figure2.R"]
  if(length(z)) return(normalizePath(tail(z,1L)))
  z <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])
  normalizePath(z)
}
source(file.path(dirname(.local_file()), "common.R"))
suppressPackageStartupMessages({library(lme4); library(lmerTest)})

figure2_from_scans <- function(s) {
  base <- c("rid","participant_id","scan_set_id","scan_id","field_strength","site_id","age_baseline","mri_time_years",
    "sex","diagnosis_baseline","education_years","apoe4_count","efc","echo_time_category","lc_bilateral","snvta_bilateral")
  need(s,base); unique_keys(s,"scan_set_id","portable scan registry")
  centered <- c("age_baseline_c","mri_time_years_c","education_years_c","efc_c")
  out <- rbindlist(lapply(c("LC","SNVTA"),function(rg) rbindlist(lapply(c("L","R"),function(hm) {
    member <- paste("figure2",rg,hm,"order",sep="_")
    response <- paste(tolower(rg),tolower(hm),sep="_")
    saved_centers <- paste0("f2_",tolower(rg),"_",centered)
    need(s,c(member,response,saved_centers))
    ord <- number(s[[member]])
    if(any(is.finite(ord)&(ord<=0|ord!=as.integer(ord)))) stop("Invalid prepared Figure 2 order")
    keep <- which(is.finite(ord)); d <- copy(s[keep,..base])
    d[,`:=`(nucleus=rg,hemisphere=hm,hemisphere_c=if(hm=="L")-.5 else .5,
      cr=s[[response]][keep],.source_order=ord[keep])]
    for(j in seq_along(centered)) set(d,j=centered[j],value=s[[saved_centers[j]]][keep])
    d
  }))))
  unique_keys(out,".source_order","prepared Figure 2 order")
  if(!identical(sort(out$.source_order),as.numeric(seq_len(nrow(out))))) stop("Prepared Figure 2 order is not contiguous")
  setorder(out,.source_order); out[,.source_order:=NULL]; out
}

prepare <- function() {
  d <- figure2_from_scans(input("scans.csv"))
  unique_keys(d,c("scan_set_id","nucleus","hemisphere"),"Figure 2 observations")
  for(n in c("hemisphere_c","cr","age_baseline","mri_time_years","education_years","apoe4_count","efc",
    "age_baseline_c","mri_time_years_c","education_years_c","efc_c","lc_bilateral","snvta_bilateral")) set(d,j=n,value=number(d[[n]]))
  d[,field_strength:=field(field_strength)]
  v <- fread(file.path(output_root(),"private/shared/technical_cr.csv"),colClasses="character",na.strings=c("","NA"))
  v <- v[,.(scan_set_id,nucleus,hemisphere,anchored=number(cr_adjusted_anchored))]
  d <- join_one(d,v,c("scan_set_id","nucleus","hemisphere"),"anchored display CR")
  write_plan("figure2",data.table(kind="primary",nucleus=c("LC","SNVTA")),list(hemi=d))
}

factor_frame <- function(d) {
  d[, `:=`(sex = factor(text(sex), levels = c("F", "M")),
            diagnosis_baseline = factor(text(diagnosis_baseline), levels = c("CN", "MCI", "AD")),
            field_strength = factor(field(field_strength), levels = c("1.5T", "3T")),
            echo_f = factor(text(echo_time_category)), participant = factor(text(participant_id)),
            site = factor(text(site_id)), scan_set = factor(text(scan_set_id)))]
  d
}
full_covariates <- c("age_baseline_c", "mri_time_years_c", "sex", "diagnosis_baseline", "field_strength", "education_years_c", "apoe4_count", "efc_c", "echo_f")
fit_model <- function(f, d) {
  for (n in c("sex", "diagnosis_baseline", "field_strength", "echo_f")) scientific_assert(nlevels(d[[n]]) >= 2L, paste("Required Figure 2 factor has fewer than two levels:", n))
  scientific_assert(nlevels(d$site) >= 2L, "Figure 2 needs at least two sites")
  repeated <- d[, .(n_times = uniqueN(mri_time_years_c)), by = participant_id]
  scientific_assert(nrow(d) >= 121L && uniqueN(d$scan_set_id) >= 121L && uniqueN(d$rid) >= 40L && sum(repeated$n_times >= 2L) >= 30L,
                    "Figure 2 sample lacks the original longitudinal support")
  scientific_assert(qr(model.matrix(lme4::nobars(f), d))$rank == ncol(model.matrix(lme4::nobars(f), d)), "Rank-deficient Figure 2 fixed effects")
  z <- capture_fit(lmerTest::lmer(f, data = d, REML = TRUE, na.action = na.fail,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000L), check.rankX = "stop.deficient")))
  check_convergence(z$fit); z
}
primary <- function(spec, data) {
  d <- factor_frame(copy(data$hemi[nucleus == spec$nucleus]))
  vars <- c("rid", "participant", "site", "scan_set", "cr", "hemisphere_c", full_covariates)
  d <- droplevels(d[complete.cases(d[, ..vars])]); setorder(d, participant_id, scan_id, hemisphere)
  scientific_assert(nrow(d) > 0, "No complete Figure 2 frame")
  if (any(!is.finite(d$anchored))) stop("Primary Figure 2 display lacks technical-adjusted values on its fitted frame")
  f <- cr ~ age_baseline_c + mri_time_years_c + sex + diagnosis_baseline + hemisphere_c +
    field_strength + education_years_c + apoe4_count + efc_c + echo_f +
    (1 + mri_time_years_c || participant) + (1 | site) + (1 | scan_set)
  z <- fit_model(f, d); fit <- z$fit; co <- coef_lmm(fit)
  co[, `:=`(record_type = "fixed_coefficient", variant = "primary", nucleus = spec$nucleus,
             n_rows = nrow(d), n_scans = uniqueN(d$scan_set_id), n = uniqueN(d$rid),
             focal = term %in% c("age_baseline_c", "mri_time_years_c", "sexM", "hemisphere_c"),
             uncertainty = "Satterthwaite_t", formula = paste(deparse(f), collapse = " "), status = "estimated",
             singular = isSingular(fit), warnings = paste(z$warnings, collapse = " | "))]
  shift <- mean(d$anchored) - mean(d$cr)
  # Prepared centers belong to the original full reference population.
  for (variable in c("age_baseline", "mri_time_years")) {
    offset <- d[[variable]] - d[[paste0(variable, "_c")]]
    if (diff(range(offset)) > 1e-7) {
      stop("Prepared centering reference is not constant: ", variable)
    }
  }
  list(results = co, fit = fit, frame = d, display_shift = shift)
}

finalize <- function() {
  x <- collect_tasks("figure2")
  r <- bind_component(x,"results")
  record_diagnostics("figure2", r)
  publish_results(r[record_type=="fixed_coefficient"],"figure2/coefficients.csv")
  refs <- rbindlist(lapply(x,function(z) {
    if(is.null(z$frame))return(NULL)
    d <- z$frame
    data.table(nucleus=z$task$nucleus,reference=c("age_baseline","mri_time_years","education_years","efc","display_shift"),
      value=c(d$age_baseline[1]-d$age_baseline_c[1],d$mri_time_years[1]-d$mri_time_years_c[1],
              d$education_years[1]-d$education_years_c[1],d$efc[1]-d$efc_c[1],z$display_shift),
      population=c(rep("prepared full retained nucleus-specific hemisphere population",4),"fitted frame anchored mean minus raw CR mean"))
  }),fill=TRUE)
  atomic_csv(refs,"figure2/reference_values.csv")
}

if(sys.nframe()==0L) {
  a <- cli()
  if(a$action=="prepare") prepare()
  else if(a$action=="task") run_task("figure2",a$task,primary)
  else if(a$action=="finalize") finalize()
  else stop("Unknown action")
}
