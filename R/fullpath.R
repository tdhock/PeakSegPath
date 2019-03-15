now <- function(){
  as.integer(Sys.time())
}

startPath <- structure(function
(prob.dir,
  n.penalties=1000,
  penalty.vec=exp(c(-Inf, seq(-15, 15, l=n.penalties), Inf))
){
  submitPenalties(prob.dir, penalty.vec)
}, ex=function(){
  prob.dir <- "/scratch/thocking/uci-chipseq-data/ATAC_JV_adipose/samples/AC1/MSC77/problems/chr10:18024675-38818835"
  prob.dir <- "~/scratch/uci-chipseq-data/H3K27ac-H3K4me3_TDHAM_BP/samples/GR3_H3K4me3/S00K4G_WTSI/problems/chr3:60000-66170270"
  startPath(prob.dir, 2)
})

getProbSuffix <- function(prob.dir){
  sub(".*uci-chipseq-data/", "", prob.dir)
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
    chunks.as.arrayjobs=TRUE)
){
  reg.dir <- file.path(
    prob.dir,
    "PeakSegPath",
    sub(" ", "-", Sys.time()))
  dir.create(dirname(reg.dir), recursive=TRUE, showWarnings=FALSE)
  reg <- batchtools::makeRegistry(reg.dir)
  pen.str.vec <- paste(pen.vec)
  batchtools::batchMap(function(pen.str, prob.dir){
    ## Need PeakSegPath:: because this fun will be executed later by a
    ## worker node/batchtools.
    PeakSegPath::onePenDB(prob.dir, pen.str)
    PeakSegPath::submitNewPenalties(prob.dir)
  }, pen.str.vec, reg=reg, more.args=list(prob.dir=prob.dir))
  jt <- batchtools::getJobTable(reg=reg)
  chunks <- data.table(jt, chunk=(1:nrow(jt)) %/% max.tasks.per.array)
  batchtools::submitJobs(chunks, resources=res.list, reg=reg)
  after <- batchtools::getJobTable(reg=reg)
  prob.id <- getProbID(prob.dir)
  data.table(
    prob.id,
    pen.str=pen.str.vec,
    job.id=after$batch.id,
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
  coverage.bedGraph.gz <- file.path(prob.dir, "coverage.bedGraph.gz")
  if(file.exists(coverage.bedGraph.gz)){
    system(paste("gunzip", coverage.bedGraph.gz))
  }
  fit <- PeakSegDisk::problem.PeakSegFPOP(prob.dir, paste(pen.str))
  ord.segs <- fit$segments[order(chromStart)]
  diff.vec <- diff(ord.segs$mean)
  ord.segs[, diff.after  := c(diff.vec, Inf)]
  ord.segs[, diff.before := c(-Inf, diff.vec)]
  change.before.and.after <- ord.segs[diff.before != 0 & diff.after != 0]
  bkg <- change.before.and.after[status=="background"]
  start.end <- function(dt){
    dt[status=="peak", .(chromStart, chromEnd)]
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
    error.df <- PeakError::PeakErrorChrom(peak.dt, labels.dt)
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
  oneTransaction(function(con){
    DBI::dbWriteTable(con, "problems_loss", loss.dt, append=TRUE, row.names=FALSE)
  })
  to.delete <- file.path(
    prob.dir,
    paste0(
      "coverage.bedGraph_penalty=",
      pen.str, "_",
      c("timing.tsv", "segments.bed", "loss.tsv")))
  unlink(to.delete)
}

submitNewPenalties <- function(prob.dir){
  prob.id <- getProbID(prob.dir)
  oneTransaction(function(con){
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
      job.vec <- unique(job.dt$job.id.job)
      jobs <- paste(job.vec, collapse=",")
      job.arg <- paste0("-j", jobs)
      sacct.dt <- tryCatch({
        slurm::sacct(job.arg)
      }, error=function(e){
        cat("slurm un-available so assuming jobs are still working.\n")
        data.table()
      })
      if(nrow(sacct.dt)){
        join.dt <- sacct.dt[job.dt, on=list(JobID.job=job.id.job, task=job.id.task)]
        completed <- join.dt[State_blank=="COMPLETED"]
        if(nrow(completed)){
          print(completed)
          cat("deleting stale jobs from problems_computing table.\n")
          dbSendClear(con, '
delete from problems_computing
where "prob.id"=$1 and "pen.str"=$2
', params=list(prob.id, completed$pen.str))
        }
        join.dt[State_blank!="COMPLETED", pen.str]
      }
    })
    loss.dt <- data.table(DBI::dbGetQuery(con, '
select "time.computed", "total.loss", "peaks"
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
    path.dt[, max.computed := max.lambda %in% loss.dt$penalty]
    path.dt[, no.next := c(diff(peaks) == -1, NA)]
    path.dt[, done := max.computed | no.next]
    new.pen.vec <- path.dt[done==FALSE, max.lambda]
    candidate.pen.vec <- new.pen.vec[!new.pen.vec %in% pen.not.completed]
    if(0 < length(candidate.pen.vec)){
      submitted.dt <- submitPenalties(prob.dir, candidate.pen.vec)
      DBI::dbWriteTable(
        con, "problems_computing", submitted.dt,
        append=TRUE, row.names=FALSE)
    }
    cat("submitted", length(candidate.pen.vec), "new penalties\n")
  })
}

