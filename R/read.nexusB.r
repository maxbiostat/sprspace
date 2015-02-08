## Slight modification of the ape::read.nexus() function to store posterior densities from BEAST trees
read.nexusB <- function (file, type = "BEAST", tree.names = NULL) 
{
  X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
  if(type == "BEAST"){ # TODO: implement MrBayes option
    Y <- strsplit(X, "Begin trees;")
    first.half <- lapply(strsplit(X, "posterior="), function(x) x[2])
    posteriors <- na.omit(as.numeric(unlist(lapply(first.half, function (x) strsplit(x, "\\]")[[1]][1])))) # horrible, I know, but life is short. YOLO!
  } 
  print(mean(posteriors))
  LEFT <- grep("\\[", X)
  RIGHT <- grep("\\]", X)
  if (length(LEFT)) {
    w <- LEFT == RIGHT
    if (any(w)) {
      s <- LEFT[w]
      X[s] <- gsub("\\[[^]]*\\]", "", X[s])
    }
    w <- !w
    if (any(w)) {
      s <- LEFT[w]
      X[s] <- gsub("\\[.*", "", X[s])
      sb <- RIGHT[w]
      X[sb] <- gsub(".*\\]", "", X[sb])
      if (any(s < sb - 1)) 
        X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
    }
  }
  endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
  semico <- grep(";", X)
  i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
  i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
  translation <- if (length(i2) == 1 && i2 > i1) 
    TRUE
  else FALSE
  if (translation) {
    end <- semico[semico > i2][1]
    x <- X[(i2 + 1):end]
    x <- unlist(strsplit(x, "[,; \t]"))
    x <- x[nzchar(x)]
    TRANS <- matrix(x, ncol = 2, byrow = TRUE)
    TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
    n <- dim(TRANS)[1]
  }
  start <- if (translation) 
    semico[semico > i2][1] + 1
  else i1 + 1
  end <- endblock[endblock > i1][1] - 1
  tree <- X[start:end]
  rm(X)
  tree <- tree[tree != ""]
  semico <- grep(";", tree)
  Ntree <- length(semico)
  if (Ntree == 1 && length(tree) > 1) 
    STRING <- paste(tree, collapse = "")
  else {
    if (any(diff(semico) != 1)) {
      STRING <- character(Ntree)
      s <- c(1, semico[-Ntree] + 1)
      j <- mapply(":", s, semico)
      if (is.list(j)) {
        for (i in 1:Ntree) STRING[i] <- paste(tree[j[[i]]], 
                                              collapse = "")
      }
      else {
        for (i in 1:Ntree) STRING[i] <- paste(tree[j[, 
                                                     i]], collapse = "")
      }
    }
    else STRING <- tree
  }
  rm(tree)
  STRING <- STRING[grep("^[[:blank:]]*tree.*= *", STRING, ignore.case = TRUE)]
  Ntree <- length(STRING)
  nms.trees <- sub(" *= *.*", "", STRING)
  nms.trees <- sub("^[[:blank:]]*tree[[:blank:]\\*]*", "", 
                   nms.trees, ignore.case = TRUE)
  STRING <- sub("^.*= *", "", STRING)
  STRING <- gsub(" ", "", STRING)
  colon <- grep(":", STRING)
  if (!length(colon)) {
    trees <- lapply(STRING, clado.build)
  }
  else if (length(colon) == Ntree) {
    trees <- if (translation) 
      lapply(STRING, .treeBuildWithTokens)
    else lapply(STRING, tree.build)
  }
  else {
    trees <- vector("list", Ntree)
    trees[colon] <- lapply(STRING[colon], tree.build)
    nocolon <- (1:Ntree)[!1:Ntree %in% colon]
    trees[nocolon] <- lapply(STRING[nocolon], clado.build)
    if (translation) {
      for (i in 1:Ntree) {
        tr <- trees[[i]]
        for (j in 1:n) {
          ind <- which(tr$tip.label[j] == TRANS[, 1])
          tr$tip.label[j] <- TRANS[ind, 2]
        }
        if (!is.null(tr$node.label)) {
          for (j in 1:length(tr$node.label)) {
            ind <- which(tr$node.label[j] == TRANS[, 
                                                   1])
            tr$node.label[j] <- TRANS[ind, 2]
          }
        }
        trees[[i]] <- tr
      }
      translation <- FALSE
    }
  }
  for (i in 1:Ntree) {
    tr <- trees[[i]]
    if (!translation) 
      n <- length(tr$tip.label)
    ROOT <- n + 1
    if (sum(tr$edge[, 1] == ROOT) == 1 && dim(tr$edge)[1] > 
          1) {
      stop(paste("The tree apparently has singleton node(s): cannot read tree file.\n  Reading NEXUS file aborted at tree no.", 
                 i, sep = ""))
    }
  }
  if (Ntree == 1) {
    trees <- trees[[1]]
    if (translation) {
      trees$tip.label <- if (length(colon)) 
        TRANS[, 2]
      else TRANS[, 2][as.numeric(trees$tip.label)]
    }
  }
  else {
    if (!is.null(tree.names)) 
      names(trees) <- tree.names
    if (translation) {
      if (length(colon) == Ntree) 
        attr(trees, "TipLabel") <- TRANS[, 2]
      else {
        for (i in 1:Ntree) trees[[i]]$tip.label <- TRANS[, 
                                                         2][as.numeric(trees[[i]]$tip.label)]
        trees <- .compressTipLabel(trees)
      }
    }
    class(trees) <- "multiPhylo"
    if (!all(nms.trees == "")) 
      names(trees) <- nms.trees
  }
  for (i in 1:Ntree) {
    trees[[i]]$posterior <- posteriors[i]
  }
  trees
}
require(ape)
environment(read.nexusB) <- environment(read.nexus)
