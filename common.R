# Shared I/O, model summaries and deterministic task dispatch.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.poly"))

cli <- function() {
  x <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[i], "--") || i == length(x)) stop("Expected --name value arguments")
    key <- substring(x[i], 3L)
    if (key %in% names(out)) stop("Repeated argument: ", key)
    out[[key]] <- x[i + 1L]
    i <- i + 2L
  }
  if (is.null(out$action)) out$action <- "run"
  out
}
need <- function(d, cols, label = "input") {
  absent <- setdiff(cols, names(d))
  if (length(absent)) stop(label, " missing columns: ", paste(absent, collapse = ", "))
  invisible(d)
}
text <- function(x) {
  y <- trimws(as.character(x))
  y[is.na(y) | y %in% c("", "NA", "NaN", "NULL", "<NA>")] <- NA_character_
  y
}
number <- function(x) {
  y <- suppressWarnings(as.numeric(as.character(x)))
  y[!is.finite(y)] <- NA_real_
  y
}
truth <- function(x) !is.na(x) & tolower(trimws(as.character(x))) %in% c("true", "1", "t", "yes")
field <- function(x) {
  y <- toupper(text(x))
  out <- rep(NA_character_, length(y))
  out[y %in% c("1P5T", "1.5T", "1.5")] <- "1.5T"
  out[y %in% c("3T", "3.0T", "3", "3.0")] <- "3T"
  out
}
finite_z <- function(x) {
  y <- number(x)
  s <- sd(y, na.rm = TRUE)
  m <- mean(y, na.rm = TRUE)
  if (!is.finite(s) || s <= 0 || !is.finite(m)) return(rep(NA_real_, length(y)))
  (y - m) / s
}
center <- function(x) { y <- number(x); if (any(is.finite(y))) y - mean(y, na.rm = TRUE) else rep(NA_real_, length(y)) }
unique_keys <- function(d, keys, label = "input", allow_missing = FALSE) {
  need(d, keys, label)
  bad <- !complete.cases(d[, ..keys])
  for (k in keys) if (is.character(d[[k]])) bad <- bad | is.na(text(d[[k]]))
  if (!allow_missing && any(bad)) stop(label, " has missing identity keys")
  if (anyDuplicated(d[!bad, ..keys])) stop(label, " has duplicate keys: ", paste(keys, collapse = ", "))
  invisible(d)
}
join_one <- function(left, right, keys, label = "lookup") {
  unique_keys(right, keys, label)
  if (length(intersect(setdiff(names(right), keys), names(left)))) stop(label, ": overlapping non-key columns")
  d <- copy(left); d[, .join_order := .I]
  d <- merge(d, right, by = keys, all.x = TRUE, sort = FALSE, allow.cartesian = FALSE)
  setorder(d, .join_order); d[, .join_order := NULL]; d
}
input <- function(relative, allow_empty = FALSE) {
  root <- Sys.getenv("ADNI_MODEL_DATA_DIR")
  if (!nzchar(root)) stop("ADNI_MODEL_DATA_DIR must explicitly name the prepared candidate")
  p <- file.path(root, relative)
  if (!file.exists(p)) stop("Missing prepared input: ", p)
  d <- fread(p, colClasses = "character", na.strings = c("", "NA", "NaN"), showProgress = FALSE)
  if (!allow_empty && !nrow(d)) stop("Empty required input: ", relative)
  d
}
output_root <- function() {
  p <- Sys.getenv("ADNI_MODEL_RESULTS_DIR")
  if (!nzchar(p) || !dir.exists(p)) stop("Use run.sh to create an explicit results directory")
  p
}
atomic_csv <- function(d, relative) {
  p <- file.path(output_root(), relative); dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  t <- tempfile(pattern = ".write_", tmpdir = dirname(p)); on.exit(unlink(t), add = TRUE)
  encoded <- copy(as.data.table(d))
  for (nm in names(encoded)) {
    x <- encoded[[nm]]
    if (is.double(x) && !inherits(x,c("Date","POSIXt","integer64"))) {
      y <- sprintf("%.17g",x); y[is.na(x)] <- NA_character_;set(encoded,j=nm,value=y)
    }
  }
  fwrite(encoded, t, na = "NA", quote = "auto", nThread = 1L, dateTimeAs = "ISO")
  if (!file.rename(t, p)) stop("Unable to commit output: ", p)
  invisible(p)
}
atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  t <- tempfile(pattern = ".rds_", tmpdir = dirname(path)); on.exit(unlink(t), add = TRUE)
  saveRDS(x, t, compress = "gzip")
  if (!file.rename(t, path)) stop("Unable to commit RDS: ", path)
}
work_path <- function(unit, name) file.path(output_root(), ".work", unit, name)
run_signature <- function() {
  x <- Sys.getenv("MODEL_SIGNATURE")
  if (!nzchar(x)) stop("Missing run signature; use run.sh")
  x
}
canonical_value <- function(x) {
  if(is.data.frame(x))return(list(columns=names(x),data=lapply(x,canonical_value),rows=nrow(x)))
  if(is.factor(x))return(list(values=as.integer(x),levels=levels(x),ordered=is.ordered(x)))
  if(is.list(x))return(lapply(x,canonical_value))
  x
}
write_plan <- function(unit, tasks, data) {
  tasks <- as.data.table(tasks); tasks[, task_id := seq_len(.N)]
  setcolorder(tasks, c("task_id", setdiff(names(tasks), "task_id")))
  plan<-list(signature=run_signature(),execution_scope="main",tasks=tasks,data=data)
  path<-work_path(unit,"prepared.rds")
  if(file.exists(path)) {
    old<-load_plan(unit)
    if(!identical(canonical_value(old),canonical_value(plan)))stop("Prepared inputs/tasks changed on resume: ",unit)
    return(invisible(old))
  }
  atomic_rds(plan,path)
  atomic_csv(tasks, file.path(".work", unit, "tasks.csv"))
  writeLines(as.character(nrow(tasks)), work_path(unit, "n_tasks.txt"))
  cat(unit, ": prepared", nrow(tasks), "tasks\n")
}
load_plan <- function(unit) {
  x <- readRDS(work_path(unit, "prepared.rds"))
  if (!identical(x$signature, run_signature())) stop("Prepared work belongs to another run")
  if (!identical(x$execution_scope, "main")) stop("Prepared work belongs to another execution scope")
  x
}
not_estimable <- function(message) {
  stop(structure(list(message = message, call = NULL), class = c("model_not_estimable", "error", "condition")))
}
scientific_assert <- function(test, message) { if (!isTRUE(test)) not_estimable(message) }
run_task <- function(unit, task_id, fun) {
  p <- load_plan(unit); id <- as.integer(task_id)
  if (length(id) != 1L || is.na(id) || id < 1L || id > nrow(p$tasks)) stop("Invalid task ID")
  destination <- work_path(unit, sprintf("task_%04d.rds", id))
  if (file.exists(destination)) {
    old <- readRDS(destination)
    if (!identical(old$signature, run_signature())) stop("Task belongs to another run")
    if (!identical(old$status, "failed")) { cat(unit, id, "already complete\n"); return(invisible(old)) }
  }
  spec <- p$tasks[id]; caught <- NULL
  ans <- tryCatch(fun(spec, p$data), model_not_estimable = function(e) {
    list(results = cbind(spec, data.table(status = "not_estimable", reason = conditionMessage(e))))
  }, error = function(e) { caught <<- e; NULL })
  if (!is.null(caught)) {
    atomic_rds(list(signature = run_signature(), task = spec, status = "failed",
                    error = conditionMessage(caught)), destination)
    stop(caught)
  }
  if (!is.list(ans)) stop("Task must return a named list")
  ans$signature <- run_signature(); ans$task <- spec
  if (is.null(ans$status)) ans$status <- "complete"
  atomic_rds(ans, destination)
  cat(unit, "task", id, "completed\n")
  invisible(ans)
}
drop_task_objects <- function(x) {
  x[intersect(names(x),c("fit","fits","model","models","companion_models","frame","data","model_data","basis"))]<-NULL
  x
}
collect_tasks <- function(unit, drop_objects=FALSE) {
  p <- load_plan(unit)
  lapply(seq_len(nrow(p$tasks)), function(id) {
    f <- work_path(unit, sprintf("task_%04d.rds", id))
    if (!file.exists(f)) stop("Missing task result: ", f)
    x <- readRDS(f)
    if (!identical(x$signature, run_signature()) || identical(x$status, "failed")) stop("Unsuccessful task: ", f)
    if(drop_objects)x<-drop_task_objects(x)
    x
  })
}
bind_component <- function(tasks, name) rbindlist(lapply(tasks, function(x) x[[name]]), use.names = TRUE, fill = TRUE)
run_profile <- function() Sys.getenv("MODEL_PROFILE", "analysis")
resample_n <- function() { if (run_profile() == "selftest") 20L else 5000L }
capture_fit <- function(expr) {
  warnings <- character()
  fit <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  list(fit = fit, warnings = unique(warnings))
}
fit_info <- function(fit, warnings = character()) {
  lmm <- inherits(fit, "merMod")
  conv <- if (lmm) fit@optinfo$conv$lme4$messages else NULL
  data.table(formula = paste(deparse(formula(fit)), collapse = " "), n_rows = nobs(fit),
             singular = if (lmm) lme4::isSingular(fit, tol = 1e-4) else FALSE,
             convergence_message = paste(conv, collapse = " | "),
             warnings = paste(warnings, collapse = " | "), run_profile = run_profile())
}
check_convergence <- function(fit) {
  if (!inherits(fit, "merMod")) return(invisible(fit))
  code <- unlist(fit@optinfo$conv$opt); msgs <- fit@optinfo$conv$lme4$messages
  bad <- msgs[!grepl("singular|boundary", msgs, ignore.case = TRUE)]
  if ((length(code) && any(code != 0)) || length(bad)) stop("Mixed model failed convergence: ", paste(bad, collapse = "; "))
  invisible(fit)
}
coef_lmm <- function(fit) {
  d <- as.data.table(as.data.frame(summary(fit)$coefficients), keep.rownames = "term")
  mp <- c("Estimate" = "beta", "Std. Error" = "se", "t value" = "statistic", "Pr(>|t|)" = "p_value")
  for (n in names(mp)) if (n %in% names(d)) setnames(d, n, mp[[n]])
  if (!"df" %in% names(d)) d[, df := Inf]
  if (!"p_value" %in% names(d)) d[, p_value := NA_real_]
  crit <- ifelse(is.finite(d$df), qt(0.975, d$df), qnorm(0.975))
  d[, `:=`(ci_low = beta - crit * se, ci_high = beta + crit * se)]
  d
}
hc3 <- function(fit) {
  X <- model.matrix(fit); r <- residuals(fit); h <- hatvalues(fit)
  keep <- is.finite(r) & is.finite(h) & h < 0.999999
  X <- X[keep, , drop = FALSE]; r <- r[keep]; h <- h[keep]
  scientific_assert(qr(X)$rank == ncol(X) && nrow(X) > ncol(X), "HC3 design not estimable")
  bread <- solve(crossprod(X)); V <- bread %*% crossprod(X, X * (r / (1 - h))^2) %*% bread
  se <- sqrt(diag(V)); b <- coef(fit); df <- max(1, nrow(X) - ncol(X)); st <- b / se
  data.table(term = names(b), beta = unname(b), se = unname(se), statistic = unname(st), df = df,
             p_value = 2 * pt(abs(st), df, lower.tail = FALSE),
             ci_low = b - qt(.975, df) * se, ci_high = b + qt(.975, df) * se)
}
# Stage fragments are combined only by the coordinator after worker completion.
record_diagnostics <- function(unit,d) {
  d <- copy(as.data.table(d)); d[, unit := unit]
  atomic_csv(d,file.path(".work","diagnostics",paste0(unit,".csv")))
}
collect_diagnostics <- function() {
  files <- list.files(file.path(output_root(),".work","diagnostics"),pattern="[.]csv$",full.names=TRUE)
  external<-file.path(output_root(),".work",c("plsc/diagnostics.csv","molecular/diagnostics.csv","molecular/null_diagnostics.csv"))
  files<-c(files,external[file.exists(external)])
  d <- if(length(files)) rbindlist(lapply(files,fread),fill=TRUE) else data.table(unit=character(),status=character())
  atomic_csv(d,"run_diagnostics.csv")
}
publish_results <- function(d,relative) {
  d <- copy(as.data.table(d))
  if(!ncol(d))d<-data.table(model_id=character(),term=character(),estimate=numeric(),se=numeric(),ci_low=numeric(),ci_high=numeric(),P=numeric(),P_fdr=numeric())
  aliases <- c(p_raw="P",pvalue="P",p_value="P",p_fdr="P_fdr",q_value="P_fdr",q_bh="P_fdr")
  for (nm in names(aliases)) if(nm %in% names(d) && !aliases[[nm]] %in% names(d)) setnames(d,nm,aliases[[nm]])
  diagnostic <- grep("(^|_)(status|reason|warnings?|converged|convergence_message|singular|run_profile|bootstrap_failure_messages)$",names(d),value=TRUE)
  diagnostic <- setdiff(diagnostic, "pathology_status") # Scientific binary amyloid exposure, not an execution status.
  if(length(diagnostic))d[,(diagnostic):=NULL]
  atomic_csv(d,relative)
}
read_results <- function(relative,allow_empty=FALSE) {
  p <- file.path(output_root(),relative)
  if(!file.exists(p))stop("Missing completed within-run dependency: ",relative)
  d <- fread(p,colClasses=c(rid="character",RID="character",participant_id="character",site_id="character"),showProgress=FALSE)
  if(!allow_empty&&!nrow(d))stop("Empty within-run dependency: ",relative)
  d
}
# Reconstruct the recorded prepared frame in its source order.
scan_members <- function(stage) {
  if(stage!="figure4")stop("Unsupported scan membership")
  s<-input("scans.csv")
  base<-c("rid","participant_id","scan_id","scan_set_id","scan","mri_date","field_strength","site_id","age_baseline","mri_time_years","sex","diagnosis_baseline","diagnosis_mri","education_years","apoe4_count","efc","echo_time","echo_time_category")
  centers<-c("age_baseline_c","mri_time_years_c","education_years_c","efc_c")
  extras<-c("baseline_icv","baseline_icv_date","muse_visit_id","muse_image_id","muse_date","muse_gap_days","visit_icv")
  need(s,base);unique_keys(s,"scan_set_id","scan registry")
  rows<-rbindlist(lapply(c("LC","SNVTA"),function(rg)rbindlist(lapply(c("L","R"),function(hm){
    key<-paste(stage,rg,hm,"order",sep="_");need(s,key);ord<-number(s[[key]])
    if(any(is.finite(ord)&(ord<=0|ord!=as.integer(ord))))stop("Invalid prepared source order")
    ids<-which(is.finite(ord));d<-copy(s[ids,..base])
    value<-paste(tolower(rg),tolower(hm),sep="_");need(s,value)
    d[,`:=`(nucleus=rg,hemisphere=hm,hemisphere_c=if(hm=="L")-.5 else .5,cr=s[[value]][ids],.source_order=ord[ids])]
    for(nm in centers){col<-paste0("f2_",tolower(rg),"_",nm);need(s,col);set(d,j=nm,value=s[[col]][ids])}
    for(nm in extras){col<-paste0("f4_",tolower(rg),"_",nm);need(s,col);set(d,j=nm,value=s[[col]][ids])}
    d
  }))))
  unique_keys(rows,".source_order","prepared order")
  if(!identical(sort(rows$.source_order),as.numeric(seq_len(nrow(rows)))))stop("Prepared source order not contiguous")
  setorder(rows,.source_order);rows[, .source_order:=NULL];rows
}
if(sys.nframe()==0L) {
  a<-cli()
  if(a$action=="diagnostics")collect_diagnostics() else stop("common.R is a helper; use run.sh")
}
