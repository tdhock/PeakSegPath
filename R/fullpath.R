now <- function(){
  as.integer(Sys.time())
}

startPath <- structure(function
(prob.dir,
  n.penalties=1000,
  penalty.vec=exp(c(-Inf, Inf, sample(seq(-15, 15, l=n.penalties)))),
  ...
){
  submitPenalties(prob.dir, penalty.vec, ...)
}, ex=function(){
  prob.dir <- "/scratch/thocking/uci-chipseq-data/H3K27ac-H3K4me3_TDHAM_BP/samples/GR3_H3K4me3/S00K4G_WTSI/problems/chr3:60000-66170270"
  prob.dir <- "/scratch/thocking/uci-chipseq-data/ATAC_JV_adipose/samples/AC1/MSC77/problems/chr10:18024675-38818835"
  PeakSegPath::startPath(prob.dir)
  prob.dir.vec <- Sys.glob("/scratch/thocking/uci-chipseq-data/H3K27ac-H3K4me3_TDHAM_BP/samples/GR3_H3K4me3/*/problems/chr3:60000-66170270")
  prob.dir <- prob.dir.vec[3]
  for(prob.dir in prob.dir.vec){
    PeakSegPath::submitNewPenalties(prob.dir)
  }
  some.dir.vec <- grep("K4G|NBR", prob.dir.vec, value=TRUE, invert=TRUE)
  for(prob.dir in prob.dir.vec){
    PeakSegPath::startPath(prob.dir, seconds.per.penalty=60)
  }
})

getProbSuffix <- function(prob.dir){
  sub(".*uci-chipseq-data/", "", normalizePath(prob.dir, mustWork=TRUE))
}

initProblems <- function(){
  labels.bed.vec <- Sys.glob(
    "~/project/uci-chipseq-data/*/samples/*/*/problems/*/labels.bed")
  i <- 1
  prob.dir <- prob.dir.vec[i]
  prob.dir.vec <- dirname(labels.bed.vec)
  prob.dt <- data.table(
    prob.suffix=getProbSuffix(prob.dir.vec),
    prob.id=seq_along(labels.bed.vec))
  pen.str <- "Inf"
  computing.dt <- data.table(
    prob.id=prob.dt$prob.id[1],
    pen.str,
    job.id="foo",
    time.started=now(),
    max.job.seconds=Inf)
  oneTransaction(function(con){
    DBI::dbRemoveTable(con, "problems")
    DBI::dbRemoveTable(con, "problems_loss")
    DBI::dbRemoveTable(con, "problems_computing")
    DBI::dbWriteTable(con, "problems", prob.dt, row.names=FALSE)
    DBI::dbWriteTable(con, "problems_computing", computing.dt, row.names=FALSE)
  })
  onePenDB(prob.dir, pen.str)#creates problems_loss table in DB.
  oneTransaction(function(con){
    dbSendClear(con, 'create index on problems_loss ("prob.id")')
    dbSendClear(con, 'alter table problems_loss add constraint loss_pen_prob_unique unique ("pen.str", "prob.id")')
    dbSendClear(con, 'alter table problems_loss alter "total.loss" type double precision')
    dbSendClear(con, 'alter table problems_loss alter "mean.pen.cost" type double precision')
    dbSendClear(con, 'alter table problems_loss alter "penalty" type double precision')
    dbSendClear(con, 'create index on problems_computing ("prob.id")')
    dbSendClear(con, 'alter table problems_computing add constraint comp_pen_prob_unique unique ("pen.str", "prob.id")')
  })
}

dbSendClear <- function(...){
  res <- DBI::dbSendQuery(...)
  DBI::dbClearResult(res)
}

oneTransaction <- function
### Run one function on the database. The idea is to make these runs
### short as possible in order to limit the number of open
### connections. Must use function so that we can insert objects into
### the database that were present at function definition (this does
### not work with expressions).
(fun,
### expression containing "con" database connection object, which will
### be evaluated after creating a connection.
  drv=RPostgres::Postgres(),
### database driver to pass to dbConnect.
  sleep.seconds=10
### seconds to wait if dbConnect fails (typically because there are
### already too many connections.
){
  con <- NULL
  while(is.null(con)){
    tryCatch({
      con <- DBI::dbConnect(drv)
      on.exit(DBI::dbDisconnect(con))
    }, error=function(e){
      ## if all connections are busy, try waiting for one to open up.
      Sys.sleep(sleep.seconds)
    })
  }
  DBI::dbWithTransaction(con, fun(con))
}

## This object is created at build time.
problems.dt <- data.table(oneTransaction(function(con){
  DBI::dbReadTable(con, "problems")
}))
setkey(problems.dt, prob.suffix)

submitPenalties <- function
(prob.dir,
  pen.vec,
  max.tasks.per.array=9999,
  res.list=list(
    walltime = 8*60*60,#seconds
    memory = 1000,#megabytes per cpu
    ncpus=1,
    ntasks=1,
    chunks.as.arrayjobs=TRUE),
  seconds.per.penalty=res.list$walltime
){
  coverage.bedGraph.gz <- file.path(prob.dir, "coverage.bedGraph.gz")
  coverage.bedGraph <- file.path(prob.dir, "coverage.bedGraph")
  if(file.exists(coverage.bedGraph.gz)){
    system(paste("gunzip", coverage.bedGraph.gz))
  }
  if(!file.exists(coverage.bedGraph)){
    stop(coverage.bedGraph, " does not exist")
  }
  reg.dir <- file.path(
    prob.dir,
    "PeakSegPath",
    sub(" ", "-", Sys.time()))
  dir.create(dirname(reg.dir), recursive=TRUE, showWarnings=FALSE)
  reg <- batchtools::makeRegistry(reg.dir)
  pen.str.vec <- paste(pen.vec)
  penalties.per.job <- ceiling(res.list$walltime / seconds.per.penalty / 2)  
  n.jobs <- length(pen.str.vec) / penalties.per.job
  job.vec <- 1:n.jobs
  pen.dt <- data.table(
    pen.str=pen.str.vec,
    job.id=rep(job.vec, l=length(pen.str.vec)))
  batchtools::batchMap(function(job.id, prob.dir, pen.dt){
    ## Need PeakSegPath:: because this fun will be executed later by a
    ## worker node/batchtools.
    select.dt <- data.table(job.id)
    job.penalties <- pen.dt[select.dt, on=list(job.id)]
    for(penalty.i in 1:nrow(job.penalties)){
      pen.str <- job.penalties[penalty.i, pen.str]
      cat(sprintf("%4d / %4d penalty=%s\n", penalty.i, nrow(job.penalties), pen.str))
      PeakSegPath::onePenDB(prob.dir, pen.str)
    }
    PeakSegPath::submitNewPenalties(prob.dir)
  }, job.vec, reg=reg, more.args=list(
    pen.dt=pen.dt,
    prob.dir=prob.dir))
  jt <- batchtools::getJobTable(reg=reg)
  chunks <- data.table(jt, chunk=(1:nrow(jt)) %/% max.tasks.per.array)
  batchtools::submitJobs(chunks, resources=res.list, reg=reg)
  after <- batchtools::getJobTable(reg=reg)[, .(job.id, batch.id)]
  job.vec <- namedCapture::str_match_variable(after$batch.id, 
    job="[0-9]+", 
    "_")[, "job"]
  sacct.dt <- data.table()
  do.print <- TRUE
  while(nrow(sacct.dt) < nrow(after)){
    sacct.dt <- sacct.jobs(job.vec)
    if(do.print){
      cat("some jobs still not in sacct\n")
      print(sacct.dt)
    }
    do.print <- TRUE
    Sys.sleep(1)
  }
  join.dt <- pen.dt[after, on=list(job.id)]
  prob.id <- getProbID(prob.dir)
  data.table(
    prob.id,
    join.dt[, .(pen.str, job.id=batch.id)],
    time.started=-1,
    max.job.seconds=res.list$walltime)
}

getProbID <- function(prob.dir){
  suffix <- getProbSuffix(prob.dir)
  problems.dt[suffix, prob.id]
}

onePenDB <- function(prob.dir, pen.str){
  prob.id <- getProbID(prob.dir)
  on.exit({#in case of error exit:
    oneTransaction(function(con){
      dbSendClear(con, '
delete from problems_computing
where "prob.id"=$1 and "pen.str"=$2
', params=list(prob.id, pen.str))
    })
  })
  oneTransaction(function(con){
    dbSendClear(con, '
update problems_computing
set "time.started"=$1
where "prob.id"=$2 and "pen.str"=$3
', params=list(now(), prob.id, pen.str))
  })
  fit <- PeakSegDisk::problem.PeakSegFPOP(prob.dir, paste(pen.str))
  ord.segs <- fit$segments[.N:1]
  diff.vec <- diff(ord.segs$mean)
  ord.segs[, diff.after  := c(diff.vec, Inf)]
  ord.segs[, diff.before := c(-Inf, diff.vec)]
  change.before.and.after <- ord.segs[diff.before != 0 & diff.after != 0]
  bkg <- change.before.and.after[status=="background"]
  start.end <- function(dt){
    dt[status=="peak", data.table(chromStart, chromEnd)]
  }
  peak.list <- list(
    use=start.end(ord.segs),
    remove=start.end(change.before.and.after),
    join=bkg[, data.table(
      chromStart=chromEnd[-.N],
      chromEnd=chromStart[-1])])
  labels.bed <- file.path(prob.dir, "labels.bed")
  labels.dt <- fread(labels.bed, col.names=c(
    "chrom", "chromStart", "chromEnd", "annotation"))
  loss.dt <- data.table(
    pen.str,
    prob.id,
    time.computed=now(),
    fit$loss)    
  for(rule in names(peak.list)){
    peak.dt <- peak.list[[rule]]
    error.df <- tryCatch({
      PeakError::PeakErrorChrom(peak.dt, labels.dt)
    }, error=function(e){
      print(e)#typically peak bases not increasing, weird bug in solver.
      cat("error in PeakErrorChrom, returning NA\n")
      data.frame(fp=NA_integer_, fn=NA_integer_)
    })
    error.vec <- with(error.df, c(
      errors=sum(fp+fn),
      fp=sum(fp),
      fn=sum(fn)))
    for(error.type in names(error.vec)){
      error.col <- paste0(rule, ".", error.type)
      loss.dt[[error.col]] <- error.vec[[error.type]]
    }
  }
  loss.dt[, mean.intervals := as.numeric(mean.intervals)]
  tryCatch({
    oneTransaction(function(con){
      DBI::dbWriteTable(con, "problems_loss", loss.dt, append=TRUE, row.names=FALSE)
    })
  }, error=function(e){
    print(e)
    cat("error writing loss\n")
  })
  to.delete <- file.path(
    prob.dir,
    paste0(
      "coverage.bedGraph_penalty=",
      pen.str, "_",
      c("timing.tsv", "segments.bed", "loss.tsv")))
  unlink(to.delete)
}

sacct.jobs <- function(job.vec){
  u.vec <- unique(job.vec)
  jobs <- paste(u.vec, collapse=",")
  job.arg <- paste0("-j", jobs)
  tryCatch({
    slurm::sacct(job.arg)
  }, error=function(e){
    print(e)
    print(prob.dir)
    cat(sprintf(
      "slurm::sacct('%s') errored so assuming jobs are still working.\n",
      job.arg))
    data.table()
  })
}

### return either NULL (if there are no new penalties, or if writing
### to problems_computing failed due to unique constraint -- another
### process already wrote the same value) or a list of arguments to
### pass to submitPenalties (new penalties that are reserved but not
### yet submitted).
reservePenalties <- function(prob.dir){
  prob.id <- getProbID(prob.dir)
  tryCatch({
    oneTransaction(function(con){
      ##browser()
      computing.dt <- data.table(DBI::dbGetQuery(con, '
select "pen.str", "job.id", "time.started", "max.job.seconds"
from problems_computing
where "prob.id"=$1
', params=list(prob.id)))
      computing.dt[, started := 0 < time.started]
      computing.dt[, running := started & (now() < time.started+max.job.seconds)]
      check.dt <- computing.dt[running==FALSE]
      pen.not.completed <- c(computing.dt[running==TRUE, pen.str], if(nrow(check.dt)){
        job.dt <- namedCapture::df_match_variable(check.dt, job.id=list(
          job="[0-9]+", as.integer,
          "_",
          task="[0-9]+", as.integer))
        sacct.dt <- sacct.jobs(job.dt$job.id.job)
        if(nrow(sacct.dt)==0){
          sacct.dt <- data.table(
            JobID.job=integer(),
            task=integer(),
            State_blank=character())
        }
        join.dt <- sacct.dt[job.dt, on=list(JobID.job=job.id.job, task=job.id.task)]
        ## State_blank could be FAILED NODE_FAIL   TIMEOUT COMPLETED
        is.complete <- join.dt[, !State_blank %in% c("RUNNING", "PENDING")]
        completed <- join.dt[is.complete]
        if(nrow(completed)){
          print(completed)
          cat("deleting stale jobs from problems_computing table.\n")
          dbSendClear(con, '
delete from problems_computing
where "prob.id"=$1 and "pen.str"=$2
', params=list(rep(prob.id, nrow(completed)), completed$pen.str))
        }
        join.dt[!is.complete, pen.str]
      })
      loss.dt <- data.table(DBI::dbGetQuery(con, '
select "time.computed", "total.loss", "peaks", "seconds", "pen.str"
from problems_loss
where "prob.id"=$1
', params=list(prob.id)))
      (loss.inc <- loss.dt[, list(
        total.loss=total.loss[which.min(time.computed)]
      ), by=list(peaks)][order(peaks)])
      loss.inc[, cummin := cummin(total.loss)]
      loss.min <- loss.inc[total.loss==cummin]
      path.dt <- data.table(penaltyLearning::modelSelection(
        loss.min, "total.loss", "peaks"))
      cat(sprintf(
        "penalties=%d peaks=%d selected=%d computing=%d incomplete=%d\n",
        nrow(loss.dt), nrow(loss.inc), nrow(path.dt), nrow(computing.dt), length(pen.not.completed)))
      path.dt[, max.computed := max.lambda %in% loss.dt$pen.str]
      path.dt[, no.next := c(diff(peaks) == -1, NA)]
      path.dt[, done := max.computed | no.next]
      ##print(path.dt[, table(max.computed, no.next)])
      new.pen.vec <- unique(paste(path.dt[done==FALSE, max.lambda]))
      candidate.pen.vec <- new.pen.vec[!new.pen.vec %in% pen.not.completed]
      if(0 < length(candidate.pen.vec)){
        candidate.dt <- data.table(
          prob.id,
          pen.str=candidate.pen.vec,
          job.id="NONE",
          time.started=-1,
          max.job.seconds=-1)
        print(candidate.dt)
        DBI::dbWriteTable(
          con, "problems_computing", candidate.dt,
          append=TRUE, row.names=FALSE)
        list(prob.dir=prob.dir, pen.vec=candidate.pen.vec, seconds.per.penalty=mean(loss.dt$seconds))
      }else{
        cat("no candidate penalties\n")
      }
    })
  }, error=function(e){
    print(e)
    cat("error in reservePenalties\n")
    NULL
  })
}

stopComputing <- function(){
  system("psql -c 'delete from problems_computing'")
  system("squeue -u thocking|sed 's/ .*//'|xargs scancel")
}

submitNewPenalties <- function(prob.dir){
  submit.list <- reservePenalties(prob.dir)
  if(is.list(submit.list)){
    cat("submitting", length(submit.list$pen.vec), "new penalties\n")
    submitted.dt <- do.call(submitPenalties, submit.list)
    params.list <- with(submitted.dt, list(
      job.id, max.job.seconds, prob.id, pen.str))
    oneTransaction(function(con){
      dbSendClear(con, '
update problems_computing
set "job.id"=$1, "max.job.seconds"=$2
where "prob.id"=$3 and "pen.str"=$4
', params=params.list)
    })
  }
}

getPath <- function(prob.dir){
  prob.id <- getProbID(prob.dir)
  oneTransaction(function(con){
    loss.dt <- data.table(DBI::dbGetQuery(con, '
select "time.computed", "total.loss", "peaks", "seconds", "pen.str", penalty, "remove.errors"
from problems_loss
where "prob.id"=$1
', params=list(prob.id)))[order(-penalty)]
    (loss.inc <- loss.dt[, list(
      total.loss=total.loss[which.min(time.computed)]
    ), by=list(peaks)][order(peaks)])
    loss.inc[, cummin := cummin(total.loss)]
    loss.min <- loss.inc[total.loss==cummin]
    path.dt <- data.table(penaltyLearning::modelSelection(
      loss.min, "total.loss", "peaks"))
    cat(sprintf(
      "penalties=%d peaks=%d selected=%d\n",
      nrow(loss.dt), nrow(loss.inc), nrow(path.dt)))
    path.dt[, max.computed := max.lambda %in% loss.dt$pen.str]
    path.dt[, no.next := c(diff(peaks) == -1, NA)]
    path.dt[, done := max.computed | no.next]
    list(loss=loss.dt, path=path.dt)
  })
}
