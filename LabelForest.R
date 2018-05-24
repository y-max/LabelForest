
testCase <- function(filename, numPerLabel){
  dat <- read.csv(filename)
  dat <- apply(dat, 2, as.numeric)
  
  qid <- getGroundTruth(dat, numPerLabel)
  groundTruth <- dat[qid,]
  unlabeled <- dat[-qid,]
  
  lv <- LabelForest(groundTruth, unlabeled)
  if(length(lv) == 0){
    print("No data sample been labeled")
  }
  else if(is.matrix(lv)){
    y.estimate <- lv[,1]
    y.truth <- unlabeled[lv[,2],ncol(unlabeled)]
    ans <- labelAccuracy(y.truth, y.estimate)
    print(data.frame(f1_per_class = ans[3], labeling_accuracy = ans[4]))
  }
  else{
    ans <- sum(lv[1] == unlabeled[lv[2],ncol(unlabeled)])
    print(sprintf("Only 1 sample been labeled with an accuracy of %d", ans))
  }
}


# ------------------------ LabelForest Function -----------------------------
# Input:
# Q - initially labeled dataset WITH the last column as the labels
# udat - unlabeled dataset WITH the last column as the true labels
# Output:
# lv - newly labeled data samples in rows, column 1 is estimated labels,
#      column 2 is the corresponding index in udat
# ---------------------------------------------------------------------------
LabelForest <- function(Q, udat){
  dat <- rbind(Q[,-ncol(Q)], udat[,-ncol(udat)])
  label <- Q[,ncol(Q)]
  ntree <- length(unique(label))
  
  # S1. compute the center of each labeled tree and initialize nset
  vid <- c()
  rej <- c()
  Margin <- rep(0, ntree)
  if(ntree == nrow(Q)){
    TC <- Q
    nset <- rep(1, ntree)
  }
  else{
    nset <- rep(0, ntree)
    for(i in 1:ntree){
      nset[label[i]] <- sum(Q[,ncol(Q)] == label[i])
    }
    TC <- treeCenter(Q, ntree)
  }
  
  # S2. heapify distance of edges adjacent to at least one vertex in vid
  for(i in 1:nrow(Q)){
    vid <- c(vid, i)
    ds <- dist(rbind(Q[i, -ncol(Q)], dat))
    ds <- ds[1:nrow(dat)]
    if(i == 1){
      hp <- dist2heap(ds[2:nrow(dat)], 1, 2:nrow(dat))
    }
    else{
      other <- setdiff(1:nrow(dat), vid)
      es <- cbind(rep(i, length(other)), other, ds[other])
      hp <- add2heap(hp, es)
    }
  }
  
  # S3. build forest using the valid edge returned from the heap
  evad <- examineEdge(Margin, hp[1,], Q, c(), TC, dat)
  while(evad != -2){
    if(evad > 0){
      s1 <- Q[hp[1,1], ncol(Q)]
      s2 <- Q[hp[1,2], ncol(Q)]
      if(Margin[s1] == 0 || Margin[s1] > hp[1,3]){
        Margin[s1] <- hp[1,3]
      }
      if(Margin[s2] == 0 || Margin[s2] > hp[1,3]){
        Margin[s2] <- hp[1,3]
      }
    }
    else if(evad <= -3){
      rej <- c(rej, hp[1,(-evad - 2)])
    }
    hp <- removeMin(hp)
    evad <- examineEdge(Margin, hp[1,], Q, c(), TC, dat)
  }
  
  # S4. update tree center and the heap after adding new vertex to the tree
  e <- hp[1,]
  hp <- removeMin(hp)
  
  vi <- max(e[1:2])
  vset <- Q[min(e[1:2]), ncol(Q)]
  Forest <- c(vset, e)
  vid <- c(vid, vi)
  ds <- dist(rbind(dat[vi,], dat))
  ds <- ds[1:nrow(dat)]
  other <- setdiff(1:nrow(dat), vid)
  es <- cbind(rep(vi, length(other)), other, ds[other])
  hp <- add2heap(hp, es)
  
  vx <- dat[vi,]
  TC[vset,-ncol(TC)] <- (nset[vset] * TC[vset, -ncol(TC)] + vx)/(nset[vset] + 1)
  nset[vset] <- nset[vset] + 1
  
  # S5. loop process until all the vertices have been visited
  while(length(vid) + length(rej) < nrow(dat)){
    e <- hp[1,]
    if(e[1] %in% rej || e[2] %in% rej){
      hp <- removeMin(hp)
    }
    else{
      evad <- examineEdge(Margin, e, Q, Forest, TC, dat)
      if(evad == 0){
        hp <- removeMin(hp)
      }
      else if(evad > 0){
        s1 <- findSet(e[1], Q, Forest)
        s2 <- findSet(e[2], Q, Forest)
        if(Margin[s1] == 0 || e[3] < Margin[s1]){
          Margin[s1] <- e[3]
        }
        if(Margin[s2] == 0 || e[3] < Margin[s2]){
          Margin[s2] <- e[3]
        }
        
        hp <- removeMin(hp)
      }
      else if(evad <= -3){
        rej <- c(rej, e[(-evad - 2)])
        hp <- removeMin(hp)
      }
      else{
        hp <- removeMin(hp)
        vi <- e[(-evad)]
        vset <- findSet(setdiff(e[1:2], vi), Q, Forest)
        Forest <- rbind(Forest, c(vset, e))
        vid <- c(vid, vi)
        other <- setdiff(1:nrow(dat), c(vid, rej))
        if(length(other) == 0){
        }
        else if(length(other) == 1){
          ds <- dist(rbind(dat[vi,], dat[other,]))
          es <- c(vi, other, ds)
          hp <- add2heap(hp, es)
        }
        else{
          ds <- dist(rbind(dat[vi,], dat[other,]))
          ds <- ds[1:length(other)]
          es <- cbind(rep(vi, length(other)), other, ds)
          hp <- add2heap(hp, es)
        }
        
        vx <- dat[vi,]
        TC[vset,-ncol(TC)] <- (nset[vset] * TC[vset, -ncol(TC)] + vx)/(nset[vset] + 1)
        nset[vset] <- nset[vset] + 1
      }
    }
  }
  
  yhat <- batchCluster(rej, Forest, nrow(Q))
  rid <- which(yhat[,1] == 0)
  if(length(rid) > 0){
    ys <- yhat[-rid,]
  }
  else{
    ys <- yhat
  }

  if(length(rid) == nrow(udat)){
    print("ZERO new vertex after GSF")
    return(c())
  }
  
  # S6. sample filtering phase
  lv <- filterVx(udat, Q, TC, Margin, Forest, ys, ntree)
  if(length(lv) == 0){
    print("ZERO new vertex after filter")
    return(c())
  }
  else{
    return(lv)
  }
}


# ---------------------------- FILTER FUNCTIONs ---------------------------
filterVx <- function(udat, Q, TC, Margin, Forest, ys, ntree){
  vsel <- filterAuto(Q, TC, ntree, udat[ys[,2], -ncol(udat)], ys[,1], rbind(Q, TC))
  if(length(vsel) == 0){
    return(c())
  }
  return(ys[vsel,])
}

# --------------------------------------------------------
# vdat - WITHOUT label selected from udat
# vy - estimated labels for vdat
# dcomp - WITH labels Q alone or [Q, TC]
#
# OUTPUT: vsel - selected index of vdat
# --------------------------------------------------------
filterAuto <- function(Q, TC, ntree, vdat, vy, dcomp){
  sbase <- rep(0, ntree)
  vsel <- c()
  for(i in 1:ntree){
    if(length(TC) > 0){
      qix <- which(Q[,ncol(Q)] == i)
      
      if(length(qix) == 1){
        sbase[i] <- max(0, silOne(Q[qix,-ncol(Q)], i, dcomp[-qix,]))
      }
      else{
        ibase <- c()
        for(q in qix){
          ibase <- c(ibase, silOne(Q[q,-ncol(Q)], i, Q[-q,]))
        }
        sbase[i] <- max(c(0, ibase))
      }
      
    }
    else{
      if(sum(Q[,ncol(Q)] == i) > 1){
        qix <- which(Q[,ncol(Q)] == i)
        ibase <- c()
        for(q in qix){
          ibase <- c(ibase, silOne(Q[q,-ncol(Q)], i, dcomp[-q,]))
        }
        sbase[i] <- max(c(0, ibase))
      }
    }
    
    if(sum(vy == i) == 1){
      sx <- silOne(vdat[vy == i,], i, Q)
      if(sx >= sbase[i]){
        vsel <- c(vsel, which(vy == i))
      }
    }
    else if(sum(vy == i) > 1){
      tix <- which(vy == i)
      six <- apply(vdat[tix,], 1, function(x) silOne(x, i, Q))
      if(sum(six >= sbase[i]) > 0){
        vsel <- c(vsel, tix[which(six >= sbase[i])])
      }
    }
  }
  return(vsel)
}


# -------------------------------------------------------
# vx - WITHOUT label
# dcomp - WITHOUT data sample of vx
# -------------------------------------------------------
silOne <- function(vx, vset, dcomp){
  L <- unique(dcomp[,ncol(dcomp)])
  tself <- dist(rbind(vx, dcomp[dcomp[,ncol(dcomp)] == vset, -ncol(dcomp)]))
  tself <- tself[1:sum(dcomp[,ncol(dcomp)] == vset)]
  a <- mean(tself)
  bs <- c()
  for(ol in setdiff(L, vset)){
    ot <- dist(rbind(vx, dcomp[dcomp[,ncol(dcomp)] == ol, -ncol(dcomp)]))
    ot <- ot[1:sum(dcomp[,ncol(dcomp)] == ol)]
    bs <- c(bs, mean(ot))
  }
  b <- min(bs)
  sil <- (b - a)/max(b,a)
  return(sil)
}


# --------------------------- GSF Functions ------------------------
batchCluster <- function(rej, Forest, nq){
  if(length(rej) > 0){
    vlist <- cbind(rep(0, length(rej)), rej)
  }
  else{
    vlist <- c()
  }
  temp <- unique(rbind(Forest[,1:2], Forest[,c(1,3)]))
  if(length(unique(temp[,2])) < nrow(temp)){
    print("ERROR in forest - multiple sets for one vertex")
  }
  temp <- temp[-which(temp[,2] <= nq),]
  vlist <- rbind(vlist, temp)
  vlist[,2] <- vlist[,2] - nq
  vlist
}

examineEdge <- function(Margin, edge, Q, Forest, TC, dat){
  vset <- edge2set(edge, Forest, Q)
  if(vset[1] == vset[2]){ # case 1. irrelevant edge or circle - return 0
    return(0)
  }
  
  if(min(vset) > 0){ # case 2. conflicting edge
    return(edge[3])
  }
  
  vnew <- which(vset == -1)
  vx <- dat[edge[vnew],]
  k <- closestTC(TC, vx)
  vi <- vset[vset > 0]
  if(k == vi){
    if(Margin[vi] == 0 || edge[3] < Margin[vi]){  # case 3. add new vertex to vi
      return(-vnew)
    }
    else{
      return(-(2 + vnew))                   # case 4. reject new vertex
    }
  }
  else{
    return(0)                       # case 5. skip making decision for now
  }
  
}

# ---------------------------------------------------
# return label set for two vertices in one edge
# default - [-1, -1]
# ---------------------------------------------------
edge2set <- function(edge, Forest, Q){
  ans <- rep(-1, 2)
  if(max(edge[1:2]) <= nrow(Q)){
    return(Q[edge[1:2],ncol(Q)])
  }
  if(is.null(Forest)){
    if(min(edge[1:2]) <= nrow(Q)){
      t <- which.min(edge[1:2])
      ans[t] <- Q[edge[t],ncol(Q)]
    }
  }
  else{
    ans[1] <- findSet(edge[1], Q, Forest)
    ans[2] <- findSet(edge[2], Q, Forest)
  }
  return(ans)
}

findSet <- function(vx, Q, Forest){
  if(vx <= nrow(Q)){
    return(Q[vx, ncol(Q)])
  }
  if(is.matrix(Forest)){
    a <- c(which(Forest[,2] == vx), which(Forest[,3] == vx))
    if(length(a) == 0){
      return(-1)
    }
    else{
      return(Forest[a[1],1])
    }
  }
  else{
    if(vx %in% Forest[2:3]){
      return(Forest[1])
    }
    else{
      return(-1)
    }
  }
}

closestTC <- function(TC, vx){
  ds <- dist(rbind(vx, TC[,-ncol(TC)]))
  ds <- ds[1:nrow(TC)]
  k <- which.min(ds)
  return(k)
}

treeCenter <- function(Q, ntree){
  if(nrow(Q) == ntree){
    return(Q)
  }
  tc <- c()
  for(i in 1:ntree){
    if(length(which(Q[,ncol(Q)] == i)) > 1){
      t <- c(colMeans(Q[Q[,ncol(Q)] == i, -ncol(Q)]), i)
    }
    else{
      t <- Q[Q[,ncol(Q)] == i,]
    }
    tc <- rbind(tc, t)
  }
  return(tc)
}

# -------------------------- HEAP FUNCTIONS -------------------------
dist2heap <- function(ds, vi, others){
  tail <- 0
  Hp <- c()
  for(i in 1:length(ds)){
    Hp = rbind(Hp, c(vi, others[i], ds[i]))
    tail = tail + 1
    
    child = tail
    parent = round(child/2)
    while(parent > 0){
      if(Hp[child,3] < Hp[parent,3]){
        temp = Hp[parent,]
        Hp[parent,] = Hp[child,]
        Hp[child,] = temp
      }
      child = parent
      parent = round(child/2)
    }
  }
  Hp
}

add2heap <- function(Hp, es){
  if(is.matrix(es)){
    for(i in 1:nrow(es)){
      Hp <- insertOne(Hp, es[i,])
    }
  }
  else{
    Hp <- insertOne(Hp, es)
  }
  return(Hp)
}

removeMin <- function(Hp){
  if(!is.matrix(Hp)){
    return(NULL)
  }
  tail = nrow(Hp)
  Hp[1,] = Hp[tail,]
  Hp = Hp[-tail,]
  
  parent = 1
  left = 2
  right = 3
  while((left < tail) && ((Hp[parent,3] > Hp[left,3]) || ((right < tail) && (Hp[parent,3] > Hp[right,3])))){
    if(right < tail){
      if(Hp[left,3] <= Hp[right,3]){
        temp = Hp[left,]
        Hp[left,] = Hp[parent,]
        Hp[parent,] = temp
        parent = left
      }
      else{
        temp = Hp[right,]
        Hp[right,] = Hp[parent,]
        Hp[parent,] = temp
        parent = right
      }
    }
    else{
      temp = Hp[left,]
      Hp[left,] = Hp[parent,]
      Hp[parent,] = temp
      parent = left
    }
    left = parent * 2
    right = parent * 2 + 1
  }
  
  Hp
}

insertOne <- function(Hp, e){
  if(is.matrix(Hp)){
    child <- nrow(Hp) + 1
  }
  else{
    child <- 2
  }
  parent <- floor(child/2)
  Hp <- rbind(Hp, e)
  while(Hp[parent,3] > Hp[child,3]){
    temp <- Hp[parent,]
    Hp[parent,] <- Hp[child,]
    Hp[child,] <- temp
    child <- parent
    parent <- floor(child/2)
    if(parent < 1){
      break
    }
  }
  Hp
}


# ---------------------------- TEST Functions -------------------------------
getGroundTruth <- function(dat, inilabel){
  L <- sort(unique(dat[,ncol(dat)]))
  qid <- c()
  for(i in 1:length(L)){
    yi <- which(dat[,ncol(dat)] == L[i])
    if(length(yi) <= inilabel){
      qid <- c(qid, yi)
    }
    else{
      qid <- c(qid, sample(yi, inilabel))
    }
  }
  qid
}

labelAccuracy <- function(ytruth, yhat){
  L <- unique(ytruth)
  t <- rbind(rep(0, length(L)), rep(0, length(L)))
  for(i in 1:length(L)){
    tp <- sum(yhat[which(ytruth == L[i])] == L[i])
    t[1,i] <- tp/sum(ytruth == L[i])
    if(sum(yhat == L[i]) > 0){
      t[2,i] <- tp/sum(yhat == L[i])
    }
  }
  prec <- mean(t[1,])
  rec <- mean(t[2,])
  if((rec + prec) == 0){
    f1 <- 0
  }
  else{
    f1 <- (2 * prec * rec)/(prec + rec)
  }
  acc <- round(sum(ytruth == yhat)/length(ytruth), 3)
  return(c(prec, rec, f1, acc))
}

