# TODO: Add proper roxygen documentation

# R functions to automatically perform the workflow described in section 6 of the fgwas manual on fine-mapping input.
# (fgwas manual: https://github.com/joepickrell/fgwas/blob/master/man/fgwas_manual.pdf)

# These functions assume sumstat_file is formatted properly for "fgwas -fine", i.e.:
  # has "SNPID", "CHR", "POS", "Z", "F", "N", "SE", "SEGNUMBER", and anno cols.
    # CHR must named "chr#" not just "#".
    # F and N columns are overridden by SE, but must be present nonetheless for some reason.
    # Rows sorted by SEGNUMBER,CHR,POS


fgwasPruneInsignifAnnos <- function(sumstat_file, annos_to_try, fgwas_bin, silent=T) {
    runs <- data.table(
      anno_to_try = annos_to_try,
      model_to_run = annos_to_try,
      output_name = paste0("out/fgwas/prune_",annos_to_try)) # TODO: eventually, don't hardcode the output dir 
    fwrite(runs, "in/fgwas/runs.tsv", sep='\t', col.names=F)

    system(ignore.stdout=silent, ignore.stderr=silent, paste(
      "parallel -a in/fgwas/runs.tsv --colsep '\t'",
      fgwas_bin, "-fine",
        "-i", sumstat_file,
        "-w {2}",
        "-o {3}"))

    # If confidence interval for anno's effect estimate contains 0, not significant.
    cis <- do.call(rbind, lapply(paste0(runs$output_name,".params"), fread))
    cis[, `:=`(parameter = sub("_ln","",parameter),
                   CI_lo = as.numeric(sub( "<" ,"",CI_lo)),
                   CI_hi = as.numeric(sub( ">", "",CI_hi)))]

    cis[!(CI_lo<=0 & CI_hi>=0), parameter]
}

# Iteratively adds to a model the anno which most increases likelihood, until no more improvement. 
fgwasGetBestLikelihoodAnnos <- function(sumstat_file, annos_to_try, fgwas_bin, silent=T) {

  aux_fxn <- function(sumstat_file, annos_to_try, annos_acc, best_llk, fgwas_bin, silent) {
    # It may seem like it would be easier to parallelize a for loop in R than use GNU parallel.
    # However, every call to system() launches a new shell (source: https://ro-che.info/articles/2020-12-11-r-system2).
    # fgwas is run many times, so the overhead is significant. Hence why a single system() call w/ GNU parallel is better.
    runs <- data.table(
      anno_to_try = annos_to_try,
      model_to_run = sapply(annos_to_try, function(anno_to_try) paste(collapse='+', c(annos_acc, anno_to_try))),
      output_name  = paste0("out/fgwas/try_",length(annos_acc),"+",annos_to_try))
    fwrite(runs, "in/fgwas/runs.tsv", sep='\t', col.names=F)

    system(ignore.stdout=silent, ignore.stderr=silent, paste(
      "parallel -a in/fgwas/runs.tsv --colsep '\t'",
      fgwas_bin, "-fine",
        "-i", sumstat_file,
        "-w {2}",  # model, e.g. anno1+anno3+anno7
        "-o {3}")) # output name prefix

    llks <- setNames(sapply(paste0(runs$output_name,".llk"), function(f) read.table(f)[1,2]), runs$anno_to_try)

    if(max(llks) > best_llk) {
      new_best_llk <- max(llks)
      anno_to_add  <- names(llks)[which.max(llks)]
      new_anno_acc <- c(annos_acc, anno_to_add)
      new_annos_to_try <- annos_to_try[annos_to_try!=anno_to_add]

      if(length(new_annos_to_try)==0) return(new_anno_acc)

      return(aux_fxn(sumstat_file, new_annos_to_try, new_anno_acc, new_best_llk, fgwas_bin, silent) )
    } else return(annos_acc)
  }

  # The above function is tail-recursive.
  # tailr::loop_transform makes it slightly faster and makes sure it doesn't overflow the stack, because R doesn't implement tail recursion optimizations.
    # (Tail recursion is when you don't need to hold all recursive calls in memory because your recursive function returns only either a lone call to itself or a final value.)
  aux_fxn <- tailr::loop_transform(aux_fxn)

  aux_fxn(sumstat_file, annos_to_try, c(), 0, fgwas_bin, silent)
}



# Tests which penalty value has the best cross-validation likelihood, given the model fgwasGetBestLikelihoodModel()
fgwasGetBestXValidationPenalty <- function(sumstat_file, annos, fgwas_bin, silent=T) {
  # TODO: Does the function likelihood(model,penalty) have only one minimum? If so, could optimize by making more intelligent steps than just 0.05 intervals.
    # NOTE: can withhold from parallelizing this function until it's clear whether making more intelligent steps is possible
  # TODO: Go finer than steps of 0.05?
  penalties <- seq(0, 1, 0.05)
  model <- paste(collapse='+', annos)

  xv_llks <- sapply(penalties, function(penalty) {
    fgwas_out_prefix <- paste0("out/fgwas/penalty_",as.character(penalty))

    system(ignore.stdout=silent, ignore.stderr=silent, paste(
      fgwas_bin, "-fine -xv -print",
        "-i", sumstat_file,
        "-w", model,
        "-p", penalty,
        "-o", fgwas_out_prefix))

    # fgwas puts the x-validation likelihood result in the last line of the .ridgeparams output files.
    tmp <- readLines(paste0(fgwas_out_prefix,".ridgeparams"))
    xv_llk <- as.numeric(tail(strsplit(tmp[length(tmp)], ' ')[[1]], 1))
  })

  return(list(
    best_xv_llk = max(xv_llks),
    best_penalty = penalties[which.max(xv_llks)]
  ))
}



# Iteratively drops annos from the model until the cross-validation likelihood is maximized.
fgwasGetBestXValidatedAnnos <- function(sumstat_file, annos, xv_penalty, best_xv_llk, fgwas_bin, silent=T) {
  # TODO: Here I assume that dropping an anno later in the model *might* cause an anno earlier in the model to now also become worth dropping even if it wasn't before.
    # This might be overly cautious, in which case, this could be optimized by only needing to traverse the model once, dropping whatever increases the llk without having to look back.
    # NOTE: can withhold from parallelizing this function until the answer to this is clear.
  xv_llks <- sapply(annos, function(anno_to_drop) { # TODO: embarassingly parallel lapply
    fgwas_out_prefix <- paste0("out/fgwas/drop_",anno_to_drop)

    trimmed_annos <- annos[annos!=anno_to_drop]
    model <- paste(collapse='+', annos)
    #system(ignore.stdout=T, ignore.stderr=T, paste(
    system(ignore.stdout=silent, ignore.stderr=silent, paste(
      fgwas_bin, "-fine -xv -print",
        "-i", sumstat_file,
        "-w", model,
        "-p", xv_penalty,
        "-o", fgwas_out_prefix))

    tmp <- readLines(paste0(fgwas_out_prefix,".ridgeparams"))
    xv_llk <- as.numeric(tail(strsplit(tmp[length(tmp)], ' ')[[1]], 1))
  })

  best_new_xv_llk <- max(xv_llks)
  if(max(xv_llks) > best_xv_llk) {
    new_best_xv_llk <- max(xv_llks)
    anno_to_drop <- names(xv_llks)[which.max(xv_llks)]

    return(fgwasGetBestXValidatedAnnos(sumstat_file, xv_penalty, trimmed_annos, new_best_xv_llk))
  } else return(annos)

}
fgwasGetBestXValidatedAnnos <- tailr::loop_transform(fgwasGetBestXValidatedAnnos)
