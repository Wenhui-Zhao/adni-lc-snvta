# Full-lag GAMs for 38 outcomes and both CR regions.
.F3_LAG_FILE <- tryCatch(sys.frame(1)$ofile,error=function(e) NULL)
if(is.null(.F3_LAG_FILE)).F3_LAG_FILE <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1L])
source(file.path(dirname(normalizePath(.F3_LAG_FILE,mustWork=TRUE)),"figure3_associations.R"),local=environment())
if(!requireNamespace("mgcv",quietly=TRUE))stop("Figure 3 lag models require mgcv")


lag_seed <- function(...) {
  h<-20260729
  for(v in utf8ToInt(paste(...,collapse="|")))h<-(h*131+v)%%2147483646
  as.integer(max(1,floor(h)))
}

lag_covariance <- function(fit) {
  # mgcv may not provide a smoothing-parameter-corrected covariance for every fit.
  # Report what is actually available instead of labeling a conditional band unconditional.
  unconditional<-!is.null(fit$Vc)
  V<-as.matrix(vcov(fit,freq=FALSE,unconditional=unconditional))
  if(any(!is.finite(V)))stop("GAM coefficient covariance contains nonfinite values")
  list(V=(V+t(V))/2,method=if(unconditional)"unconditional_Bayesian" else "conditional_Bayesian_smoothing_correction_unavailable")
}

lag_design_contrast <- function(fit,d,focal,other_focals,grid,center) {
  new<-as.data.frame(d[rep(1L,length(grid))])
  new$lag_c<-grid-center
  for(n in unique(c(focal,other_focals)))new[[n]]<-0
  exclude<-vapply(Filter(function(s)inherits(s,"random.effect"),fit$smooth),function(s)s$label,character(1))
  X0<-predict(fit,newdata=new,type="lpmatrix",exclude=exclude)
  new[[focal]]<-1
  X1<-predict(fit,newdata=new,type="lpmatrix",exclude=exclude)
  X1-X0
}

lag_band <- function(L,b,V,seed,nsim) {
  estimate<-as.numeric(L%*%b)
  active<-which(colSums(abs(L))>1e-14)
  if(!length(active))stop("GAM nuclear contrast has no active coefficients")
  A<-L[,active,drop=FALSE];C<-V[active,active,drop=FALSE]
  variance<-rowSums((A%*%C)*A)
  tol<-max(1,max(abs(variance)))*1e-10
  if(any(variance < -tol))stop("Negative GAM contrast variance exceeds numerical tolerance")
  se<-sqrt(pmax(0,variance));positive<-is.finite(se)&se>sqrt(.Machine$double.eps)
  if(!any(positive))stop("GAM nuclear contrast has no positive standard error")
  e<-eigen((C+t(C))/2,symmetric=TRUE)
  etol<-max(abs(e$values))*1e-10
  if(any(e$values < -max(1,max(abs(e$values)))*1e-8))stop("GAM covariance is not positive semidefinite")
  keep<-e$values>etol
  if(!any(keep))stop("GAM covariance has no positive eigenvalues")
  # Deterministic Gaussian coefficient simulation, not a participant bootstrap.
  set.seed(seed)
  Z<-matrix(rnorm(sum(keep)*nsim),nrow=sum(keep),ncol=nsim)
  draws<-e$vectors[,keep,drop=FALSE]%*%sweep(Z,1L,sqrt(e$values[keep]),"*")
  errors<-A[positive,,drop=FALSE]%*%draws
  maxima<-apply(abs(sweep(errors,1L,se[positive],"/")),2L,max)
  critical<-as.numeric(quantile(maxima,.95,type=8,names=FALSE))
  list(estimate=estimate,se=se,critical=critical,
       low=estimate-critical*se,high=estimate+critical*se,
       point_low=estimate-qnorm(.975)*se,point_high=estimate+qnorm(.975)*se)
}

lag_intervals <- function(grid,low,high) {
  state<-ifelse(low>0,1L,ifelse(high<0,-1L,0L));runs<-rle(state)
  finish<-cumsum(runs$lengths);start<-finish-runs$lengths+1L;out<-list()
  root<-function(x0,y0,x1,y1) {
    if(!is.finite(y0-y1)||abs(y1-y0)<.Machine$double.eps)return(x1)
    x0-y0*(x1-x0)/(y1-y0)
  }
  for(i in which(runs$values!=0L)) {
    a<-start[i];z<-finish[i];sign<-runs$values[i];bound<-if(sign>0)low else high
    left<-grid[a];right<-grid[z]
    if(a>1L)left<-root(grid[a-1L],bound[a-1L],grid[a],bound[a])
    if(z<length(grid))right<-root(grid[z],bound[z],grid[z+1L],bound[z+1L])
    out[[length(out)+1L]]<-data.table(interval_id=length(out)+1L,direction=if(sign>0)"positive" else "negative",
      start_years=left,end_years=right,touches_lower_support=a==1L,touches_upper_support=z==length(grid),
      first_grid_index=a,last_grid_index=z,status="estimated")
  }
  if(!length(out))return(data.table(interval_id=0L,direction="None",start_years=NA_real_,end_years=NA_real_,
    touches_lower_support=FALSE,touches_upper_support=FALSE,first_grid_index=NA_integer_,last_grid_index=NA_integer_,status="estimated"))
  rbindlist(out)
}

lag_support <- function(d) {
  # Disjoint descriptive bins; no extra exclusion or curve extrapolation.
  x<-copy(d);breaks<-c(0,1,3,5,7,9,Inf)
  x[, support_bin:=cut(lag_years,breaks=breaks,right=FALSE,include.lowest=TRUE)]
  x[, .(n_pairs=.N,n_participants=uniqueN(participant_id),n_scans=uniqueN(scan_set_id),
        n_future_assessments=uniqueN(future_assessment_id),min_lag_years=min(lag_years),
        max_lag_years=max(lag_years)),by=support_bin]
}

fit_lag <- function(spec,fr) {
  d<-fr$data
  combinations<-setNames(list(fr$focal),"single")
  results<-curves<-intervals<-support<-coefs<-smooths<-fits<-list()
  for(model in names(combinations)) {
    focal<-combinations[[model]]
    terms<-c("baseline_value",focal,"diagnosis_stage",fr$optional,if(spec$requires_icv)"baseline_icv_z",
      's(lag_c, bs = "tp", k = 5)',paste0('s(lag_c, by = ',focal,', bs = "tp", k = 5, pc = 0)'))
    if(sum(table(d$participant)>=2L)>=2L)terms<-c(terms,'s(participant, bs = "re")')
    if(sum(table(d$scan_set)>=2L)>=2L)terms<-c(terms,'s(scan_set, bs = "re")')
    terms<-c(terms,'s(site, bs = "re")')
    f<-as.formula(paste("response ~",paste(terms,collapse=" + ")))
    cap<-capture_fit(mgcv::bam(f,data=as.data.frame(d),family=gaussian(),method="fREML",
                              discrete=TRUE,nthreads=1L,na.action=na.fail))
    fit<-cap$fit
    if(!is.null(fit$converged)&&!isTRUE(fit$converged))stop("GAM convergence failed")
    if(any(!is.finite(coef(fit))))stop("GAM coefficients are nonfinite")
    fits[[model]]<-fit;cov<-lag_covariance(fit)
    limits<-as.numeric(quantile(d$lag_years,c(.05,.95),type=7,names=FALSE))
    scientific_assert(all(is.finite(limits))&&limits[2L]>limits[1L],"No supported lag range")
    grid<-seq(limits[1L],limits[2L],length.out=201L)
    nsim<-if(run_profile()=="selftest")200L else 10000L
    su<-summary(fit)
    pt<-as.data.table(as.data.frame(su$p.table),keep.rownames="term")
    setnames(pt,names(pt)[2:5],c("beta","se","statistic","p_value"))
    st<-as.data.table(as.data.frame(su$s.table),keep.rownames="term")
    for(nm in focal) {
      nucleus<-if(nm=="lc_z")"LC" else "SNVTA"
      meta<-f3_meta(spec,fr,model,nucleus)
      L<-lag_design_contrast(fit,d,nm,focal,grid,fr$lag_center)
      Lcenter<-lag_design_contrast(fit,d,nm,focal,fr$lag_center,fr$lag_center)
      center_est<-as.numeric(Lcenter%*%coef(fit));main<-unname(coef(fit)[nm])
      if(!is.finite(main)||abs(center_est-main)>1e-7*(1+abs(main)))stop("GAM point constraint does not anchor the nuclear curve at centered lag")
      seed<-lag_seed("CLAG_GAM",nucleus,"effect",spec$outcome)
      band<-lag_band(L,coef(fit),cov$V,seed,nsim)
      row<-pt[term==nm]
      if(nrow(row)!=1L)stop("Missing parametric nuclear coefficient")
      row[, `:=`(df=fit$df.residual,ci_low=NA_real_,ci_high=NA_real_,
        status="estimated",reason=NA_character_,formula=paste(deparse(f),collapse=" "),
        covariance_method=cov$method,parametric_covariance_method="conditional_model_summary",
        parametric_test_scope="mean-lag coefficient, not whole-curve significance",band_simulations=nsim,band_seed=seed,band_critical=band$critical,
        supported_min_years=limits[1L],supported_max_years=limits[2L],
        interval_correction="simultaneous within this curve only; not across outcomes",
        warnings=paste(cap$warnings,collapse=" | "),
        convergence_message=paste(if(!is.null(fit$outer.info)) fit$outer.info$conv else "not_supplied",collapse=" | "),
        rank=fit$rank,n_coefficients=length(coef(fit)),
        deviance_explained=su$dev.expl,adjusted_r_squared=su$r.sq)]
      # Mean-lag coefficient intervals use 1.96 SE, separately from curve bands.
      row[, `:=`(ci_low=beta-1.96*se,ci_high=beta+1.96*se)]
      results[[length(results)+1L]]<-f3_add_meta(row,meta)
      curve<-data.table(lag_years=grid,lag_days=grid*365.25,estimate=band$estimate,se=band$se,
        ci_low=band$low,ci_high=band$high,pointwise_low=band$point_low,pointwise_high=band$point_high,
        excludes_zero=band$low>0|band$high<0,band_critical=band$critical,covariance_method=cov$method,status="estimated")
      curves[[length(curves)+1L]]<-f3_add_meta(curve,meta)
      ints<-lag_intervals(grid,band$low,band$high)
      intervals[[length(intervals)+1L]]<-f3_add_meta(ints,meta)
      support[[length(support)+1L]]<-f3_add_meta(lag_support(d),meta)
    }
    meta<-f3_meta(spec,fr,model,spec$nucleus)
    pt[, `:=`(status="estimated",record_type="parametric_coefficient")]
    st[, `:=`(status="estimated",record_type="smooth_test")]
    coefs[[model]]<-f3_add_meta(pt,meta);smooths[[model]]<-f3_add_meta(st,meta)
  }
  list(results=rbindlist(results,fill=TRUE),curves=rbindlist(curves,fill=TRUE),intervals=rbindlist(intervals,fill=TRUE),
       support=rbindlist(support,fill=TRUE),coefficients=rbindlist(c(coefs,smooths),fill=TRUE),fits=fits,frame=fr)
}

finalize_lag <- function() {
  all<-f3_collect_tasks("lag");m<-bind_component(all,"results");c<-bind_component(all,"curves")
  i<-bind_component(all,"intervals");s<-bind_component(all,"support");p<-bind_component(all,"coefficients")
  # A curve-level p for the mean-lag coefficient is not the simultaneous-band result.
  m[, p_fdr:=NA_real_]
  m[, multiplicity_scope:="No across-outcome correction of simultaneous bands; parametric p is not a whole-curve p"]
  missing<-m[status!="estimated"]
  if(nrow(missing)) {
    missing[, `:=`(interval_id=NA_integer_,direction=NA_character_,start_years=NA_real_,end_years=NA_real_,
      touches_lower_support=NA,touches_upper_support=NA)]
    i<-rbindlist(list(i,missing),fill=TRUE)
  }
  key<-c("model_type","outcome","variable_name","variant","window","pairing","nucleus","comparison_model","fdr_family","status")
  summary<-i[, .(n_intervals=sum(interval_id>0,na.rm=TRUE),
    significant_intervals=if(all(status!="estimated"))NA_character_ else if(all(is.na(interval_id)|interval_id==0L))"None" else
      paste(sprintf("%s: %.6f to %.6f years",direction[interval_id>0],start_years[interval_id>0],end_years[interval_id>0]),collapse="; "),
    n_participants=first(n_participants),n_scans=first(n_scans),n_rows=first(n_rows)),by=key]
  f3_diagnostics("lag",list(results=m,coefficients=p))
  for(pair in list(list(m,"lag_models.csv"),list(c,"lag_curves.csv"),list(i,"lag_intervals.csv"),
                   list(summary,"lag_interval_summary.csv"),list(s,"lag_support.csv"),list(p,"lag_coefficients.csv"))) {
    d<-pair[[1L]];if(!ncol(d))d<-data.table(status=character(),reason=character())
    atomic_csv(f3_public(d),file.path("figure3",pair[[2L]]))
  }
}

if(sys.nframe()==0L) {
  args <- cli()
  if(args$action=="prepare")f3_prepare("lag") else
  if(args$action=="task")f3_call_task("lag",args$task,fit_lag) else
  if(args$action=="finalize")finalize_lag() else stop("Use --action prepare, task or finalize")
}
