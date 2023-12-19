
if(getRversion() >= "2.15.1")  utils::globalVariables(c("C_I", "G", "S1", "S2", "b0", "b1", "b2", "b3", "b4", "c0", "c1", "c2", "ddom", "dia", "hdom", "hojd", "plot.id", "t_l"))


#' Species-specific response calibration for Scots pine, Norway spruce, and Birch
#'
#' This function uses the plot- and revision-level random effects to predict missing tree height of Scots pine, Norway spruce and birch. It contains the species-specific height functions.
#' @param DATA data frame containing at least yta, rev, t_l, dia, QMD, G, hdom, hojd.
#' @returns data frame with estimated height
#' @note yta: plot number, rev: revision, t_l: species code (1: Scot pine; 2: Norway spruce; 3: Silver birch; 4: Downy birch), dia: diameter at breast height, QMD: quadratic mean diameter, G: basal area per ha, hdom: height of the tree with the largest diameter (ddom) regardless of the species (m), hojd: sample tree height with missing values. In cases where hdom is not present in the inventory data, site index (SI) can serve as an alternative, although the estimated height may exhibit slight variations.
#' @keywords species_specific
#' @author Ogana F.N. and Arias-Rodil M.
#' @seealso [main_species()], which estimate the tree height of the main species based on equation 16 i.e., generalized function with species as a covariate (dummy variable).
#' @references Ogana et al. (2023) https://doi.org/10.1016/j.foreco.2023.120843
#' @references Arias-Rodil et al. (2015) https://doi.org/10.1371/JOURNAL.PONE.0143521
#' @importFrom stats numericDeriv
#' @importFrom stats na.omit
#' @importFrom Matrix bdiag
#' @import magic
#' @import Deriv
#' @import dplyr
#' @import tidyr
#' @export
#' @examples
#' library(THREC)
#'
#' # sample data
#' data(Treeht)
#'
#' species_specific(Treeht)

species_specific <- function(DATA){

  if(any(names(DATA) != "hdom")){
    names(DATA)[names(DATA)=="SI"] <- "hdom"
  }else{
    names(DATA)[names(DATA)=="hdom"] <- "hdom"
  }

  DATA$plot.id = with(DATA, match(paste(yta), unique(paste(yta))))

  if(any(DATA$t_l %in% c(1:4))){

    if(any(DATA$t_l == 1)){

      DATA$C_I <- with(DATA, dia/QMD)

      Pine <- function(dia, hdom, G, C_I, parms, b, randparms, parmnames){

        b0 = 1.0108471
        b1 = 0.7921256
        b2 = -0.46148570
        b3 = -0.0001625
        b4=0.0015050

        rp <- parmnames %in% randparms
        bp <- bt <- rep(0, length(parmnames))
        bp[which(rp)] <- b[1:length(randparms)]
        bt[which(rp)] <- b[(length(randparms) + 1):(2 * length(randparms))]
        prms <- parms + rp * bp + rp * bt

        hojd <- 1.3+(dia / (prms[["b0"]] + (prms[["b1"]]*hdom^prms[["b2"]]+prms[["b3"]]*G+prms[["b4"]]*C_I)*dia))^2
        return(hojd)
      }

      PredictRandomEffect <- function(dfplot, randparms, fparms, Dplot, Drevision, sigma2, delta, tolerance = 1e-2){
        nrand <- length(randparms)
        revname <- unique(dfplot$rev)
        nrev <- length(revname)
        nobs <- nrow(dfplot)
        posrev <- tapply(1:nobs, dfplot$rev, min)
        lrev <- tapply(dfplot$rev, dfplot$rev, length)
        lrand <- nrand + nrand * nrev
        rp <- names(fparms) %in% randparms

        Dlist <- list()
        Dlist[[1]] <- Dplot
        for(i in 2:(nrev + 1)){
          Dlist[[i]] <- Drevision
        }
        D <- as.matrix(bdiag(Dlist))

        Mi <- diag(1, nrow = nobs, ncol = nobs)
        if(!missing(delta)){
          Mi.list <- list()
          for(i in 1:nrev){
            M.it <- Mi[posrev[i]:(posrev[i] + lrev[i] - 1), posrev[i]:(posrev[i] + lrev[i] - 1)]
            dia.it <- dfplot$dia[posrev[i]]
            Mi.list[[i]] <- M.it * dia.it ^ (2 * delta)
          }
          Mi <- as.matrix(bdiag(Mi.list))
        }
        Ri <- sigma2 * Mi

        yi <- dfplot$hojd

        b.0 <- rep(0, lrand)
        tol <- rep(1, lrand)

        while (sum(tol > tolerance) > 1e-2){

          Zlist <- list()
          fxiBblist <- list()
          posrand <- seq(nrand + 1, lrand, nrand)
          for(i in 1:nrev){
            posrand.i <- posrand[i]
            b.z <- c(b.0[1:nrand], b.0[posrand.i:(posrand.i + nrand - 1)])
            rev.obs <- dfplot[posrev[i]:(posrev[i] + lrev[i] - 1), ]
            Zit <- attr(numericDeriv(quote(Pine(dia = rev.obs$dia, hdom = rev.obs$hdom, G = rev.obs$G, C_I = rev.obs$C_I, parms = fparms, b = b.z, randparms = randparms, parmnames = names(fparms))), theta = "b.z"), "gradient")[, (nrand + 1):(2 * nrand)]
            fxiBb <- Pine(dia = rev.obs$dia, hdom = rev.obs$hdom, G = rev.obs$G, C_I = rev.obs$C_I, parms = fparms, b = b.z, randparms = randparms, parmnames = names(fparms))
            if(is.null(nrow(Zit))) Zit <- matrix(Zit, nrow = 1)
            Zlist[[i]] <- Zit
            fxiBblist[[i]] <- fxiBb
          }
          Zi.plot <- do.call(rbind, Zlist)
          Zit.rev <- as.matrix(bdiag(Zlist))
          Zi <- cbind(Zi.plot, Zit.rev)
          fxiBb <- do.call(c, fxiBblist)

          b <- D %*% t(Zi) %*% solve (Ri + Zi %*% D %*% t(Zi)) %*% ((yi - fxiBb) + Zi %*% b.0)
          if (all(b.0 == 0)) b.prev <- rep(1, lrand) else b.prev <- b.0
          tol <- abs((b - b.prev) / b.prev)
          b.0 <- b
        }
        bi <- split(b, ceiling(seq_along(b) / nrand), )
        names(bi) <- c("plot.id", paste(" ", revname, sep = ""))
        bi <- do.call(rbind, bi)
        colnames(bi) <- randparms
        return(bi)
      }

      fparms.sp <- c(b0 = 1.0108471, b1 = 0.7921256, b2 = -0.46148570, b3 = -0.0001625, b4=0.0015050)

      Dplot <- matrix(c(0.19008400^2, -0.744*0.19008400*0.01158396, -0.744*0.19008400*0.01158396, 0.01158396^2), nrow = 2, byrow = T)
      Drevision <- matrix(c(0.128720699^2, -0.676*0.128720699*0.009275845, -0.676*0.128720699*0.009275845, 0.009275845^2), nrow = 2, byrow = T)


      sigma2 <- 0.755024152^2

      delta <- 0.09244063

      val_pine_1 <- DATA[order(DATA$plot.id, DATA$rev), ] %>%
        select(plot.id, rev, hojd, dia, hdom, G, C_I) %>%
        na.omit()

      sp <- split(val_pine_1, val_pine_1$plot.id)
      results<-c()
      for (i in 1:length(sp)){
        caldat<-sp[[i]]
        bq <- PredictRandomEffect(dfplot = caldat, randparms = c("b0", "b2"), fparms = fparms.sp, Dplot = Dplot, Drevision = Drevision, sigma2 = sigma2, delta = delta)
        pl <- rep(i,nrow(bq)-1)
        bq2 <- bq[-1,]
        rev <- if(is.null(rownames(bq2))) {
          rownames(bq)[2]
        }
        else {rownames(bq2)
        }
        b0.u <- rep(bq[1,1],nrow(bq)-1)
        b2.u <- rep(bq[1,2],nrow(bq)-1)
        bqt <- bq[-1,]
        rownames(bqt)<-NULL
        b0.v=if(length(bqt)[1]<=2){
          bqt[1]
        }
        else {bqt[,1]
        }
        b2.v=if(length(bqt)[1]<=2){
          bqt[2]
        }
        else {bqt[,2]
        }
        results<-rbind(results, data.frame(plot.id=pl, rev=rev, b0.u=b0.u, b2.u=b2.u, b0.v=b0.v, b2.v=b2.v))
      }

      results$rev <- with(results, as.integer(rev))
      coef_p <- merge(results, t(fparms.sp))


      output <- merge(DATA, coef_p, by=c("plot.id", "rev")) %>%
        mutate(hojd=ifelse(is.na(hojd), round(1.3+(dia / ((b0+b0.u+b0.v) + (b1*hdom^(b2+b2.u+b2.v)+b3*G+b4*C_I)*dia))^2,1), hojd))%>%
        select(-c(plot.id,C_I,b0.u,b2.u,b0.v,b2.v,b0,b1,b2,b3,b4))

    }

    else if(any(DATA$t_l == 2)){

      DATA$C_I <- with(DATA, dia/QMD)

      Spruce <- function(dia, hdom, G, C_I, parms, b, randparms, parmnames){

        b0 = 1.5683142
        b1 = 0.5235892
        b2 = -0.3189100
        b3 = -0.0005743
        b4 = 0.0074741

        rp <- parmnames %in% randparms
        bp <- bt <- rep(0, length(parmnames))
        bp[which(rp)] <- b[1:length(randparms)]
        bt[which(rp)] <- b[(length(randparms) + 1):(2 * length(randparms))]
        prms <- parms + rp * bp + rp * bt

        hojd <- 1.3+(dia / (prms[["b0"]] + (prms[["b1"]]*sqrt(hdom)^prms[["b2"]]+prms[["b3"]]*G+prms[["b4"]]*C_I)*dia))^3
        return(hojd)
      }

      PredictRandomEffect <- function(dfplot, randparms, fparms, Dplot, Drevision, sigma2, delta, tolerance = 1e-2){
        nrand <- length(randparms)
        revname <- unique(dfplot$rev)
        nrev <- length(revname)
        nobs <- nrow(dfplot)
        posrev <- tapply(1:nobs, dfplot$rev, min)
        lrev <- tapply(dfplot$rev, dfplot$rev, length)
        lrand <- nrand + nrand * nrev
        rp <- names(fparms) %in% randparms

        Dlist <- list()
        Dlist[[1]] <- Dplot
        for(i in 2:(nrev + 1)){
          Dlist[[i]] <- Drevision
        }
        D <- as.matrix(bdiag(Dlist))

        Mi <- diag(1, nrow = nobs, ncol = nobs)
        if(!missing(delta)){
          Mi.list <- list()
          for(i in 1:nrev){
            M.it <- Mi[posrev[i]:(posrev[i] + lrev[i] - 1), posrev[i]:(posrev[i] + lrev[i] - 1)]
            dia.it <- dfplot$dia[posrev[i]]
            Mi.list[[i]] <- M.it * dia.it ^ (2 * delta)
          }
          Mi <- as.matrix(bdiag(Mi.list))
        }
        Ri <- sigma2 * Mi

        yi <- dfplot$hojd

        b.0 <- rep(0, lrand)
        tol <- rep(1, lrand)

        while (sum(tol > tolerance) > 1e-2){

          Zlist <- list()
          fxiBblist <- list()
          posrand <- seq(nrand + 1, lrand, nrand)
          for(i in 1:nrev){
            posrand.i <- posrand[i]
            b.z <- c(b.0[1:nrand], b.0[posrand.i:(posrand.i + nrand - 1)])
            rev.obs <- dfplot[posrev[i]:(posrev[i] + lrev[i] - 1), ]
            Zit <- attr(numericDeriv(quote(Spruce(dia = rev.obs$dia, hdom = rev.obs$hdom, G = rev.obs$G, C_I = rev.obs$C_I, parms = fparms, b = b.z, randparms = randparms, parmnames = names(fparms))), theta = "b.z"), "gradient")[, (nrand + 1):(2 * nrand)]
            fxiBb <- Spruce(dia = rev.obs$dia, hdom = rev.obs$hdom, G = rev.obs$G, C_I = rev.obs$C_I, parms = fparms, b = b.z, randparms = randparms, parmnames = names(fparms))
            if(is.null(nrow(Zit))) Zit <- matrix(Zit, nrow = 1)
            Zlist[[i]] <- Zit
            fxiBblist[[i]] <- fxiBb
          }
          Zi.plot <- do.call(rbind, Zlist)
          Zit.rev <- as.matrix(bdiag(Zlist))
          Zi <- cbind(Zi.plot, Zit.rev)
          fxiBb <- do.call(c, fxiBblist)

          b <- D %*% t(Zi) %*% solve (Ri + Zi %*% D %*% t(Zi)) %*% ((yi - fxiBb) + Zi %*% b.0)
          if (all(b.0 == 0)) b.prev <- rep(1, lrand) else b.prev <- b.0
          tol <- abs((b - b.prev) / b.prev)
          b.0 <- b
        }
        bi <- split(b, ceiling(seq_along(b) / nrand), )
        names(bi) <- c("plot.id", paste(" ", revname, sep = ""))
        bi <- do.call(rbind, bi)
        colnames(bi) <- randparms
        return(bi)
      }

      fparms.ns <- c(b0 = 1.5683142, b1 = 0.5235892, b2 = -0.3189100, b3 = -0.0005743, b4 = 0.0074741)

      Dplot <- matrix(c(0.200880919^2, -0.555*0.200880919*0.008432095, -0.555*0.200880919*0.008432095,0.008432095^2), nrow = 2, byrow = T)
      Drevision <- matrix(c(0.138460407^2, -0.383*0.138460407*0.004704052, -0.383*0.138460407*0.004704052, 0.004704052^2), nrow = 2, byrow = T)

      sigma2 <- 0.331556629^2

      delta <- 0.4955218


      val_spruce_1 <- DATA[order(DATA$plot.id, DATA$rev), ] %>%
        select(plot.id, rev, hojd, dia, hdom, G, C_I) %>%
        na.omit()

      ns <- split(val_spruce_1, val_spruce_1$plot.id)
      results<-c()
      for (i in 1:length(ns)){
        caldat<-ns[[i]]
        bq <- PredictRandomEffect(dfplot = caldat, randparms = c("b0", "b4"), fparms = fparms.ns, Dplot = Dplot, Drevision = Drevision, sigma2 = sigma2, delta = delta)
        pl <- rep(i,nrow(bq)-1)
        bq2 <- bq[-1,]
        rev <- if(is.null(rownames(bq2))) {
          rownames(bq)[2]
        }
        else {rownames(bq2)
        }
        b0.u <- rep(bq[1,1],nrow(bq)-1)
        b4.u <- rep(bq[1,2],nrow(bq)-1)
        bqt <- bq[-1,]
        rownames(bqt)<-NULL
        b0.v=if(length(bqt)[1]<=2){
          bqt[1]
        }
        else {bqt[,1]
        }
        b4.v=if(length(bqt)[1]<=2){
          bqt[2]
        }
        else {bqt[,2]
        }
        results<-rbind(results, data.frame(plot.id=pl, rev=rev, b0.u=b0.u, b4.u=b4.u, b0.v=b0.v, b4.v=b4.v))
      }

      results$rev <- with(results, as.integer(rev))
      coef_ns <- merge(results, t(fparms.ns))


      output <- merge(DATA, coef_ns, by=c("plot.id", "rev")) %>%
        mutate(hojd=ifelse(is.na(hojd), round(1.3+(dia / ((b0+b0.u+b0.v) + (b1*sqrt(hdom)^b2+b3*G+(b4+b4.u+b4.v)*C_I)*dia))^3,1), hojd))%>%
        select(-c(plot.id,C_I,b0.u,b4.u,b0.v,b4.v,b0,b1,b2,b3,b4))

    }

    else if(any(DATA$t_l == 3|DATA$t_l == 4)){

      DATA$C_I <- with(DATA, dia/QMD)

      Birch <- function(dia, hdom, G, C_I, parms, b, randparms, parmnames){

        b0 = -0.7999016
        b1 = -0.6412611
        b2 = -0.7338666
        b3 = 0.0003234
        b4 = -0.0012111

        rp <- parmnames %in% randparms
        bp <- bt <- rep(0, length(parmnames))
        bp[which(rp)] <- b[1:length(randparms)]
        bt[which(rp)] <- b[(length(randparms) + 1):(2 * length(randparms))]
        prms <- parms + rp * bp + rp * bt

        hojd <- 1.3+(dia / (prms[["b0"]] + (prms[["b1"]]*sqrt(hdom)^prms[["b2"]]+prms[["b3"]]*G+prms[["b4"]]*C_I)*dia))^2
        return(hojd)
      }



      PredictRandomEffect <- function(dfplot, randparms, fparms, Dplot, Drevision, sigma2, delta, tolerance = 1e-2){
        nrand <- length(randparms)
        revname <- unique(dfplot$rev)
        nrev <- length(revname)
        nobs <- nrow(dfplot)
        posrev <- tapply(1:nobs, dfplot$rev, min)
        lrev <- tapply(dfplot$rev, dfplot$rev, length)
        lrand <- nrand + nrand * nrev
        rp <- names(fparms) %in% randparms

        Dlist <- list()
        Dlist[[1]] <- Dplot
        for(i in 2:(nrev + 1)){
          Dlist[[i]] <- Drevision
        }
        D <- as.matrix(bdiag(Dlist))

        Mi <- diag(1, nrow = nobs, ncol = nobs)
        if(!missing(delta)){
          Mi.list <- list()
          for(i in 1:nrev){
            M.it <- Mi[posrev[i]:(posrev[i] + lrev[i] - 1), posrev[i]:(posrev[i] + lrev[i] - 1)]
            dia.it <- dfplot$dia[posrev[i]]
            Mi.list[[i]] <- M.it * dia.it ^ (2 * delta)
          }
          Mi <- as.matrix(bdiag(Mi.list))
        }
        Ri <- sigma2 * Mi

        yi <- dfplot$hojd

        b.0 <- rep(0, lrand)
        tol <- rep(1, lrand)

        while (sum(tol > tolerance) > 1e-2){
          Zlist <- list()
          fxiBblist <- list()
          posrand <- seq(nrand + 1, lrand, nrand)
          for(i in 1:nrev){
            posrand.i <- posrand[i]
            b.z <- c(b.0[1:nrand], b.0[posrand.i:(posrand.i + nrand - 1)])
            rev.obs <- dfplot[posrev[i]:(posrev[i] + lrev[i] - 1), ]
            Zit <- attr(numericDeriv(quote(Birch(dia = rev.obs$dia, hdom = rev.obs$hdom, G = rev.obs$G, C_I = rev.obs$C_I, parms = fparms, b = b.z, randparms = randparms, parmnames = names(fparms))), theta = "b.z"), "gradient")[, (nrand + 1):(2 * nrand)]
            fxiBb <- Birch(dia = rev.obs$dia, hdom = rev.obs$hdom, G = rev.obs$G, C_I = rev.obs$C_I, parms = fparms, b = b.z, randparms = randparms, parmnames = names(fparms))
            if(is.null(nrow(Zit))) Zit <- matrix(Zit, nrow = 1)
            Zlist[[i]] <- Zit
            fxiBblist[[i]] <- fxiBb
          }
          Zi.plot <- do.call(rbind, Zlist)
          Zit.rev <- as.matrix(bdiag(Zlist))
          Zi <- cbind(Zi.plot, Zit.rev)
          fxiBb <- do.call(c, fxiBblist)

          b <- D %*% t(Zi) %*% solve (Ri + Zi %*% D %*% t(Zi)) %*% ((yi - fxiBb) + Zi %*% b.0)
          if (all(b.0 == 0)) b.prev <- rep(1, lrand) else b.prev <- b.0
          tol <- abs((b - b.prev) / b.prev)
          b.0 <- b
        }
        bi <- split(b, ceiling(seq_along(b) / nrand), )
        names(bi) <- c("plot.id", paste(" ", revname, sep = ""))
        bi <- do.call(rbind, bi)
        colnames(bi) <- randparms
        return(bi)
      }


      fparms.birch <- c(b0 = -0.7999016, b1 = -0.6412611, b2 = -0.7338666, b3 = 0.0003234, b4 = -0.0012111)

      Dplot <- matrix(c(0.19177979^2, 0.633*0.19177979*0.03900942, 0.633*0.19177979*0.03900942, 0.03900942^2), nrow = 2, byrow = T)
      Drevision <- matrix(c(0.12223499^2, 0.595*0.12223499*0.01782558, 0.595*0.12223499*0.01782558, 0.01782558^2), nrow = 2, byrow = T)

      sigma2 <- 0.54450698^2

      delta <- 0.2398267


      val_birch_1 <- DATA[order(DATA$plot.id, DATA$rev), ] %>%
        select(plot.id, rev, hojd, dia, hdom, G, C_I) %>%
        na.omit()

      birch <- split(val_birch_1, val_birch_1$plot.id)
      results<-c()
      for (i in 1:length(birch)){
        caldat<-birch[[i]]
        bq <- PredictRandomEffect(dfplot = caldat, randparms = c("b0", "b2"), fparms = fparms.birch, Dplot = Dplot, Drevision = Drevision, sigma2 = sigma2, delta = delta)
        pl <- rep(i,nrow(bq)-1)
        bq2 <- bq[-1,]
        rev <- if(is.null(rownames(bq2))) {
          rownames(bq)[2]
        }
        else {rownames(bq2)
        }
        b0.u <- rep(bq[1,1],nrow(bq)-1)
        b2.u <- rep(bq[1,2],nrow(bq)-1)
        bqt <- bq[-1,]
        rownames(bqt)<-NULL
        b0.v=if(length(bqt)[1]<=2){
          bqt[1]
        }
        else {bqt[,1]
        }
        b2.v=if(length(bqt)[1]<=2){
          bqt[2]
        }
        else {bqt[,2]
        }
        results<-rbind(results, data.frame(plot.id=pl, rev=rev, b0.u=b0.u, b2.u=b2.u, b0.v=b0.v, b2.v=b2.v))
      }

      results$rev <- with(results, as.integer(rev))
      coef_birch <- merge(results, t(fparms.birch))

      output <- merge(DATA, coef_birch, by=c("plot.id", "rev")) %>%
        mutate(hojd=ifelse(is.na(hojd), round(1.3+(dia / ((b0+b0.u+b0.v) + (b1*sqrt(hdom)^(b2+b2.u+b2.v)+b3*G+b4*C_I)*dia))^2,1), hojd))%>%
        select(-c(plot.id,C_I,b0.u,b2.u,b0.v,b2.v,b0,b1,b2,b3,b4))


    }


    return(output)

  }

  else{
    stop("Function can only be applied to Scots pine, Norway spruce, Silver birch, and Downy birch", call. = FALSE)
  }


}



#' Main species generalized function response calibration
#'
#' This function performs a response calibration for main species at once. The main species include Scots pine, Norway spruce, and birch.
#' @param DATA data frame containing at least yta, rev, t_l, dia, QMD, G, hdom, and hojd.
#' @returns data frame with estimated height
#' @note yta: plot number, rev: revision, t_l: species code, dia: diameter at breast height, QMD: quadratic mean diameter, G: basal area per ha, hdom: height of the tree with the largest diameter (ddom) regardless of the species (m), hojd: sample tree height with missing values. In cases where hdom is not present in the inventory data, site index (SI) can serve as an alternative, although the estimated height may exhibit slight variations.
#' @keywords main_species
#' @author Ogana F.N. and Arias-Rodil M.
#' @seealso [species_specific()], which estimate the tree height of the main species based on species-specific height functions.
#' @references Ogana et al. (2023) https://doi.org/10.1016/j.foreco.2023.120843
#' @references Arias-Rodil et al. (2015) https://doi.org/10.1371/JOURNAL.PONE.0143521
#' @importFrom stats numericDeriv
#' @importFrom stats na.omit
#' @importFrom Matrix bdiag
#' @import magic
#' @import Deriv
#' @import dplyr
#' @import tidyr
#' @export
#' @examples
#' library(THREC)
#'
#' # sample data
#' data(Treeht)
#'
#' main_species(Treeht)

main_species <- function(DATA){

  if(any(names(DATA) != "hdom")){
    names(DATA)[names(DATA)=="SI"] <- "hdom"
  }else{
    names(DATA)[names(DATA)=="hdom"] <- "hdom"
  }

  DATA$plot.id = with(DATA, match(paste(yta), unique(paste(yta))))

  if(any(DATA$t_l %in% c(1:4))){
    DATA$C_I <- with(DATA, dia/QMD)

    FMS <- function(dia, S1, S2, hdom, G, C_I, parms, b, randparms, parmnames){

      c0 = 1.1169918
      c1 = 0.1651139
      c2 = -0.1455055
      b1 = 0.7307220
      b2 = -0.4483489
      b3 = -0.0002491
      b4 = 0.0041421

      rp <- parmnames %in% randparms
      bp <- bt <- rep(0, length(parmnames))
      bp[which(rp)] <- b[1:length(randparms)]
      bt[which(rp)] <- b[(length(randparms) + 1):(2 * length(randparms))]
      prms <- parms + rp * bp + rp * bt

      hojd <- 1.3+(dia / ((prms[["c0"]]+prms[["c1"]]*S1+prms[["c2"]]*S2) + (prms[["b1"]]*hdom^prms[["b2"]]+prms[["b3"]]*G+prms[["b4"]]*C_I)*dia))^2
      return(hojd)
    }

    PredictRandomEffect <- function(dfplot, randparms, fparms, Dplot, Drevision, sigma2, delta, tolerance = 1e-2){
      nrand <- length(randparms)
      revname <- unique(dfplot$rev)
      nrev <- length(revname)
      nobs <- nrow(dfplot)
      posrev <- tapply(1:nobs, dfplot$rev, min)
      lrev <- tapply(dfplot$rev, dfplot$rev, length)
      lrand <- nrand + nrand * nrev
      rp <- names(fparms) %in% randparms

      Dlist <- list()
      Dlist[[1]] <- Dplot
      for(i in 2:(nrev + 1)){
        Dlist[[i]] <- Drevision
      }
      D <- as.matrix(bdiag(Dlist))

      Mi <- diag(1, nrow = nobs, ncol = nobs)
      if(!missing(delta)){
        Mi.list <- list()
        for(i in 1:nrev){
          M.ij <- Mi[posrev[i]:(posrev[i] + lrev[i] - 1), posrev[i]:(posrev[i] + lrev[i] - 1)]
          dia.ij <- dfplot$dia[posrev[i]]
          Mi.list[[i]] <- M.ij * dia.ij ^ (2 * delta)
        }
        Mi <- as.matrix(bdiag(Mi.list))
      }
      Ri <- sigma2 * Mi

      yi <- dfplot$hojd

      b.0 <- rep(0, lrand)
      tol <- rep(1, lrand)

      while (sum(tol > tolerance) > 1e-2){

        Zlist <- list()
        fxiBblist <- list()
        posrand <- seq(nrand + 1, lrand, nrand)
        for(i in 1:nrev){
          posrand.i <- posrand[i]
          b.z <- c(b.0[1:nrand], b.0[posrand.i:(posrand.i + nrand - 1)])
          rev.obs <- dfplot[posrev[i]:(posrev[i] + lrev[i] - 1), ]
          Zij <- attr(numericDeriv(quote(FMS(dia = rev.obs$dia, S1 = rev.obs$S1, S2 = rev.obs$S2, hdom = rev.obs$hdom, G = rev.obs$G, C_I = rev.obs$C_I, parms = fparms, b = b.z, randparms = randparms, parmnames = names(fparms))), theta = "b.z"), "gradient")[, (nrand + 1):(2 * nrand)]
          fxiBb <- FMS(dia = rev.obs$dia, S1 = rev.obs$S1, S2 = rev.obs$S2, hdom = rev.obs$hdom, G = rev.obs$G, C_I = rev.obs$C_I, parms = fparms, b = b.z, randparms = randparms, parmnames = names(fparms))
          if(is.null(nrow(Zij))) Zij <- matrix(Zij, nrow = 1)
          Zlist[[i]] <- Zij
          fxiBblist[[i]] <- fxiBb
        }
        Zi.plot <- do.call(rbind, Zlist)
        Zij.rev <- as.matrix(bdiag(Zlist))
        Zi <- cbind(Zi.plot, Zij.rev)
        fxiBb <- do.call(c, fxiBblist)

        b <- D %*% t(Zi) %*% solve (Ri + Zi %*% D %*% t(Zi)) %*% ((yi - fxiBb) + Zi %*% b.0)
        if (all(b.0 == 0)) b.prev <- rep(1, lrand) else b.prev <- b.0
        tol <- abs((b - b.prev) / b.prev)
        b.0 <- b
      }
      bi <- split(b, ceiling(seq_along(b) / nrand), )
      names(bi) <- c("plot.id", paste(" ", revname, sep = ""))
      bi <- do.call(rbind, bi)
      colnames(bi) <- randparms
      return(bi)
    }


    fparms.fms <- c(c0=1.1169918, c1=0.1651139, c2=-0.1455055, b1 = 0.7307220, b2 = -0.4483489, b3 = -0.0002491, b4=0.0041421)

    Dplot <- matrix(c(0.19377930^2, -0.685*0.19377930*0.01449417, -0.685*0.19377930*0.01449417, 0.01449417^2), nrow = 2, byrow = T)
    Drevision <- matrix(c(0.19600558^2, -0.861*0.19600558*0.01553089, -0.861*0.19600558*0.01553089, 0.01553089^2), nrow = 2, byrow = T)

    sigma2 <- 0.91518389^2

    delta <- 0.05072372

    DATA <- DATA[order(DATA$plot.id, DATA$rev), ] %>%
      mutate(S1=ifelse(DATA$t_l=='2', 1, 0), S2= ifelse(DATA$t_l=='3' | DATA$t_l=='4', 1, 0))

    dbval_1 <- DATA %>%
      select(plot.id, rev, t_l, hojd, dia, hdom, G, C_I, S1, S2) %>%
      na.omit()

    ms <- split(dbval_1, dbval_1$plot.id)
    results<-c()
    for (i in 1:length(ms)){
      caldat<-ms[[i]]
      bq <- PredictRandomEffect(dfplot = caldat, randparms = c("c0", "b2"), fparms = fparms.fms, Dplot = Dplot, Drevision = Drevision, sigma2 = sigma2, delta = delta)
      pl <- rep(i,nrow(bq)-1)
      bq2 <- bq[-1,]
      rev <- if(is.null(rownames(bq2))) {
        rownames(bq)[2]
      }
      else {rownames(bq2)
      }
      c0.u <- rep(bq[1,1],nrow(bq)-1)
      b2.u <- rep(bq[1,2],nrow(bq)-1)
      bqt <- bq[-1,]
      rownames(bqt)<-NULL
      c0.v=if(length(bqt)[1]<=2){
        bqt[1]
      }
      else {bqt[,1]
      }
      b2.v=if(length(bqt)[1]<=2){
        bqt[2]
      }
      else {bqt[,2]
      }
      results<-rbind(results, data.frame(plot.id=pl, rev=rev, c0.u=c0.u, b2.u=b2.u, c0.v=c0.v, b2.v=b2.v))
    }

    results$rev <- with(results, as.integer(rev))
    coef_fms <- merge(results, t(fparms.fms))


    output <- merge(DATA, coef_fms, by=c("plot.id", "rev")) %>%
      mutate(hojd=ifelse(is.na(hojd), round(1.3+(dia / (((c0+c0.u+c0.v)+c1*S1+c2*S2) + (b1*hdom^(b2+b2.u+b2.v)+b3*G+b4*C_I)*dia))^2,1), hojd)) %>%
      select(-c(plot.id,C_I,S1,S2,c0.u,b2.u,c0.v,b2.v,c0,c1,c2,b1,b2,b3,b4))

    return(output)

  }
  else{
    stop("Function can only be applied to Scots pine, Norway spruce, Silver birch and, Downy birch", call. = FALSE)
  }


}




#' A response calibration function for the other conifer species
#'
#' This function performs a response calibration of the height function for other conifer species. The other conifer species include: Silver fir, Douglas fir, Sitka spruce, Siberian larch, European larch, other larch, Lodgepole pine, other fir, other spruces.
#' @param DATA data frame containing at least yta, rev, t_l, dia, hdom, ddom, and hojd.
#' @returns data frame with estimated height
#' @note yta: plot number, rev: revision, t_l: species code, dia: diameter at breast height, hdom: height of the tree with the largest diameter (ddom) regardless of the species (m), hojd: sample tree height with missing values. In cases where hdom is not present in the inventory data, site index (SI) can serve as an alternative, although the estimated height may exhibit slight variations.
#' @keywords conifers
#' @author Ogana F.N. and Arias-Rodil M.
#' @seealso [broadleaves()], which estimate the tree height of 'other broadleaves'.
#' @references Ogana et al. (2023) https://doi.org/10.1016/j.foreco.2023.120843
#' @references Arias-Rodil et al. (2015) https://doi.org/10.1371/JOURNAL.PONE.0143521
#' @importFrom stats numericDeriv
#' @importFrom stats na.omit
#' @importFrom Matrix bdiag
#' @import magic
#' @import Deriv
#' @import dplyr
#' @import tidyr
#' @export
#' @examples
#' library(THREC)
#'
#' # sample data
#' data(broad)
#'
#' conifers(broad)

conifers <- function(DATA){

  if(any(names(DATA) != "hdom")){
    names(DATA)[names(DATA)=="SI"] <- "hdom"
  }else{
    names(DATA)[names(DATA)=="hdom"] <- "hdom"
  }

  DATA$plot.id = with(DATA, match(paste(yta), unique(paste(yta))))

  if(any(DATA$t_l > 4)){

    Conifers <- function(dia, hdom, ddom, parms, b, randparms, parmnames){

      b0 = 1.0865300
      b1 = 0.7812784
      b2 = -0.5126774
      b3 = -0.0037750

      rp <- parmnames %in% randparms
      bp <- bt <- rep(0, length(parmnames))
      bp[which(rp)] <- b[1:length(randparms)]
      bt[which(rp)] <- b[(length(randparms) + 1):(2 * length(randparms))]
      prms <- parms + rp * bp + rp * bt

      hojd <- 1.3+(dia / (prms[["b0"]] + (prms[["b1"]]*sqrt(hdom)^prms[["b2"]]+prms[["b3"]]*sqrt(ddom))*dia))^3
      return(hojd)
    }

    PredictRandomEffect <- function(dfplot, randparms, fparms, Dplot, Drevision, sigma2, delta, tolerance = 1e-2){
      nrand <- length(randparms)
      revname <- unique(dfplot$rev)
      nrev <- length(revname)
      nobs <- nrow(dfplot)
      posrev <- tapply(1:nobs, dfplot$rev, min)
      lrev <- tapply(dfplot$rev, dfplot$rev, length)
      lrand <- nrand + nrand * nrev
      rp <- names(fparms) %in% randparms

      Dlist <- list()
      Dlist[[1]] <- Dplot
      for(i in 2:(nrev + 1)){
        Dlist[[i]] <- Drevision
      }
      D <- as.matrix(bdiag(Dlist))

      Mi <- diag(1, nrow = nobs, ncol = nobs)
      if(!missing(delta)){
        Mi.list <- list()
        for(i in 1:nrev){
          M.it <- Mi[posrev[i]:(posrev[i] + lrev[i] - 1), posrev[i]:(posrev[i] + lrev[i] - 1)]
          dia.it <- dfplot$dia[posrev[i]]
          Mi.list[[i]] <- M.it * dia.it ^ (2 * delta)
        }
        Mi <- as.matrix(bdiag(Mi.list))
      }
      Ri <- sigma2 * Mi

      yi <- dfplot$hojd

      b.0 <- rep(0, lrand)
      tol <- rep(1, lrand)

      while (sum(tol > tolerance) > 1e-2){

        Zlist <- list()
        fxiBblist <- list()
        posrand <- seq(nrand + 1, lrand, nrand)
        for(i in 1:nrev){
          posrand.i <- posrand[i]
          b.z <- c(b.0[1:nrand], b.0[posrand.i:(posrand.i + nrand - 1)])
          rev.obs <- dfplot[posrev[i]:(posrev[i] + lrev[i] - 1), ]
          Zit <- attr(numericDeriv(quote(Conifers(dia = rev.obs$dia, hdom = rev.obs$hdom, ddom = rev.obs$ddom, parms = fparms, b = b.z, randparms = randparms, parmnames = names(fparms))), theta = "b.z"), "gradient")[, (nrand + 1):(2 * nrand)]
          fxiBb <- Conifers(dia = rev.obs$dia, hdom = rev.obs$hdom, ddom = rev.obs$ddom, parms = fparms, b = b.z, randparms = randparms, parmnames = names(fparms))
          if(is.null(nrow(Zit))) Zit <- matrix(Zit, nrow = 1)
          Zlist[[i]] <- Zit
          fxiBblist[[i]] <- fxiBb
        }
        Zi.plot <- do.call(rbind, Zlist)
        Zit.rev <- as.matrix(bdiag(Zlist))
        Zi <- cbind(Zi.plot, Zit.rev)
        fxiBb <- do.call(c, fxiBblist)

        b <- D %*% t(Zi) %*% solve (Ri + Zi %*% D %*% t(Zi)) %*% ((yi - fxiBb) + Zi %*% b.0)
        if (all(b.0 == 0)) b.prev <- rep(1, lrand) else b.prev <- b.0
        tol <- abs((b - b.prev) / b.prev)
        b.0 <- b
      }
      bi <- split(b, ceiling(seq_along(b) / nrand), )
      names(bi) <- c("plot.id", paste(" ", revname, sep = ""))
      bi <- do.call(rbind, bi)
      colnames(bi) <- randparms
      return(bi)
    }


    fparms.con <- c(b0 = 1.0865300, b1 = 0.7812784, b2 = -0.5126774, b3 = -0.0037750)


    Dplot <- matrix(c(0.27445753^2, -0.729*0.27445753*0.02467243,  -0.729*0.27445753*0.02467243, 0.02467243^2), nrow = 2, byrow = T)
    Drevision <- matrix(c(0.116277177^2,-0.49*0.116277177*0.009885883, -0.49*0.116277177*0.009885883, 0.009885883^2), nrow = 2, byrow = T)

    sigma2 <- 0.607996858^2

    delta <- 0.179807

    val_conifer_1 <- DATA[order(DATA$plot.id, DATA$rev), ] %>%
      select(plot.id, rev, hojd, dia, hdom, ddom) %>%
      na.omit()

    conifer <- split(val_conifer_1, val_conifer_1$plot.id)
    results<-c()
    for (i in 1:length(conifer)){
      caldat<-conifer[[i]]
      bq <- PredictRandomEffect(dfplot = caldat, randparms = c("b0", "b1"), fparms = fparms.con, Dplot = Dplot, Drevision = Drevision, sigma2 = sigma2, delta = delta)
      pl <- rep(i,nrow(bq)-1)
      bq2 <- bq[-1,]
      rev <- if(is.null(rownames(bq2))) {
        rownames(bq)[2]
      }
      else {rownames(bq2)
      }
      b0.u <- rep(bq[1,1],nrow(bq)-1)
      b1.u <- rep(bq[1,2],nrow(bq)-1)
      bqt <- bq[-1,]
      rownames(bqt)<-NULL
      b0.v=if(length(bqt)[1]<=2){
        bqt[1]
      }
      else {bqt[,1]
      }
      b1.v=if(length(bqt)[1]<=2){
        bqt[2]
      }
      else {bqt[,2]
      }
      results<-rbind(results, data.frame(plot.id=pl, rev=rev, b0.u=b0.u, b1.u=b1.u, b0.v=b0.v, b1.v=b1.v))
    }

    results$rev <- with(results, as.integer(rev))
    coef_conifer <- merge(results, t(fparms.con))

    output <- merge(DATA, coef_conifer, by=c("plot.id", "rev")) %>%
      mutate(hojd=ifelse(is.na(hojd), round(1.3+(dia / ((b0+b0.u+b0.v) + ((b1+b1.u+b1.v)*sqrt(hdom)^b2 + b3*sqrt(ddom))*dia))^3,1), hojd))%>%
      select(-c(plot.id,b0.u,b1.u,b0.v,b1.v,b0,b1,b2,b3))

    return(output)
  }
  else if(any(DATA$t_l < 4)){
    stop("Function cannnot be applied to this species. Try species_specific or main_species functions", call. = FALSE)
  }

}




#' A response calibration function for the other broadleaves
#'
#' This function performs a response calibration of the height function for other broadleaves. The other broadleaves include Black alder, Grey alder, Trembling aspen, Beech, Small-leaved lime, Rowan, Goat willow, Poplars, and other deciduous species. This function is not recommended for Oak species.
#' @param DATA data frame containing at least yta, rev, t_l, dia, hdom, ddom, and hojd.
#' @returns data frame with estimated height
#' @note yta: plot number, rev: revision, t_l: species code, dia: diameter at breast height, hdom: height of the tree with the largest diameter (ddom) regardless of the species (m), hojd: sample tree height with missing values. In cases where hdom is not present in the inventory data, site index (SI) can serve as an alternative, although the estimated height may exhibit slight variations.
#' @keywords broadleaves
#' @author Ogana F.N. and Arias-Rodil M.
#' @seealso [conifers()], which estimate the tree height of other conifer species.
#' @references Ogana et al. (2023) https://doi.org/10.1016/j.foreco.2023.120843
#' @references Arias-Rodil et al. (2015) https://doi.org/10.1371/JOURNAL.PONE.0143521
#' @importFrom stats numericDeriv
#' @importFrom stats na.omit
#' @importFrom Matrix bdiag
#' @import magic
#' @import Deriv
#' @import dplyr
#' @import tidyr
#' @export
#' @examples
#' library(THREC)
#'
#' # sample data
#' data(broad)
#'
#' broadleaves(broad)


broadleaves <- function(DATA){

  if(any(names(DATA) != "hdom")){
    names(DATA)[names(DATA)=="SI"] <- "hdom"
  }else{
    names(DATA)[names(DATA)=="hdom"] <- "hdom"
  }

  DATA$plot.id = with(DATA, match(paste(yta), unique(paste(yta))))

  if(any(DATA$t_l > 4)){

    Broadleaves <- function(dia, hdom, ddom, parms, b, randparms, parmnames){

      b0 = 0.9305084
      b1 = 0.8895082
      b2 = -0.7246591
      b3 = 0.0088517

      rp <- parmnames %in% randparms
      bp <- bt <- rep(0, length(parmnames))
      bp[which(rp)] <- b[1:length(randparms)]
      bt[which(rp)] <- b[(length(randparms) + 1):(2 * length(randparms))]
      prms <- parms + rp * bp + rp * bt

      hojd <- 1.3+(dia / (prms[["b0"]] + (prms[["b1"]]*sqrt(hdom)^prms[["b2"]]+prms[["b3"]]*sqrt(ddom))*dia))^3
      return(hojd)
    }


    PredictRandomEffect <- function(dfplot, randparms, fparms, Dplot, Drevision, sigma2, delta, tolerance = 1e-2){
      nrand <- length(randparms)
      revname <- unique(dfplot$rev)
      nrev <- length(revname)
      nobs <- nrow(dfplot)
      posrev <- tapply(1:nobs, dfplot$rev, min)
      lrev <- tapply(dfplot$rev, dfplot$rev, length)
      lrand <- nrand + nrand * nrev
      rp <- names(fparms) %in% randparms

      Dlist <- list()
      Dlist[[1]] <- Dplot
      for(i in 2:(nrev + 1)){
        Dlist[[i]] <- Drevision
      }
      D <- as.matrix(bdiag(Dlist))

      Mi <- diag(1, nrow = nobs, ncol = nobs)
      if(!missing(delta)){
        Mi.list <- list()
        for(i in 1:nrev){
          M.it <- Mi[posrev[i]:(posrev[i] + lrev[i] - 1), posrev[i]:(posrev[i] + lrev[i] - 1)]
          dia.it <- dfplot$dia[posrev[i]]
          Mi.list[[i]] <- M.it * dia.it ^ (2 * delta)
        }
        Mi <- as.matrix(bdiag(Mi.list))
      }
      Ri <- sigma2 * Mi

      yi <- dfplot$hojd

      b.0 <- rep(0, lrand)
      tol <- rep(1, lrand)

      while (sum(tol > tolerance) > 1e-2){
        Zlist <- list()
        fxiBblist <- list()
        posrand <- seq(nrand + 1, lrand, nrand)
        for(i in 1:nrev){
          posrand.i <- posrand[i]
          b.z <- c(b.0[1:nrand], b.0[posrand.i:(posrand.i + nrand - 1)])
          rev.obs <- dfplot[posrev[i]:(posrev[i] + lrev[i] - 1), ]
          Zit <- attr(numericDeriv(quote(Broadleaves(dia = rev.obs$dia, hdom = rev.obs$hdom, ddom = rev.obs$ddom, parms = fparms, b = b.z, randparms = randparms, parmnames = names(fparms))), theta = "b.z"), "gradient")[, (nrand + 1):(2 * nrand)]
          fxiBb <- Broadleaves(dia = rev.obs$dia, hdom = rev.obs$hdom, ddom = rev.obs$ddom, parms = fparms, b = b.z, randparms = randparms, parmnames = names(fparms))
          if(is.null(nrow(Zit))) Zit <- matrix(Zit, nrow = 1)
          Zlist[[i]] <- Zit
          fxiBblist[[i]] <- fxiBb
        }
        Zi.plot <- do.call(rbind, Zlist)
        Zit.rev <- as.matrix(bdiag(Zlist))
        Zi <- cbind(Zi.plot, Zit.rev)
        fxiBb <- do.call(c, fxiBblist)

        b <- D %*% t(Zi) %*% solve (Ri + Zi %*% D %*% t(Zi)) %*% ((yi - fxiBb) + Zi %*% b.0)
        if (all(b.0 == 0)) b.prev <- rep(1, lrand) else b.prev <- b.0
        tol <- abs((b - b.prev) / b.prev)
        b.0 <- b
      }
      bi <- split(b, ceiling(seq_along(b) / nrand), )
      names(bi) <- c("plot.id", paste(" ", revname, sep = ""))
      bi <- do.call(rbind, bi)
      colnames(bi) <- randparms
      return(bi)
    }

    fparms.broad <- c(b0 = 0.9305084, b1 = 0.8895082, b2 = -0.7246591, b3 = 0.0088517)

    Dplot <- matrix(c(0.158579038^2, -0.445*0.158579038*0.002558484,  -0.445*0.158579038*0.002558484, 0.002558484^2), nrow = 2, byrow = T)
    Drevision <- matrix(c(0.212150220^2,-0.823*0.212150220*0.002141993, -0.823*0.212150220*0.002141993, 0.002141993^2), nrow = 2, byrow = T)

    sigma2 <- 0.767414630^2

    delta <- 0.2267002

    val_broad_1 <- DATA[order(DATA$plot.id, DATA$rev), ] %>%
      select(plot.id, rev, hojd, dia, hdom, ddom) %>%
      na.omit()

    broadleaf <- split(val_broad_1, val_broad_1$plot.id)
    results<-c()
    for (i in 1:length(broadleaf)){
      caldat<-broadleaf[[i]]
      bq <- PredictRandomEffect(dfplot = caldat, randparms = c("b0", "b3"), fparms = fparms.broad, Dplot = Dplot, Drevision = Drevision, sigma2 = sigma2, delta = delta)
      pl <- rep(i,nrow(bq)-1)
      bq2 <- bq[-1,]
      rev <- if(is.null(rownames(bq2))) {
        rownames(bq)[2]
      }
      else {rownames(bq2)
      }
      b0.u <- rep(bq[1,1],nrow(bq)-1)
      b3.u <- rep(bq[1,2],nrow(bq)-1)
      bqt <- bq[-1,]
      rownames(bqt)<-NULL
      b0.v=if(length(bqt)[1]<=2){
        bqt[1]
      }
      else {bqt[,1]
      }
      b3.v=if(length(bqt)[1]<=2){
        bqt[2]
      }
      else {bqt[,2]
      }
      results<-rbind(results, data.frame(plot.id=pl, rev=rev, b0.u=b0.u, b3.u=b3.u, b0.v=b0.v, b3.v=b3.v))
    }

    results$rev <- with(results, as.integer(rev))
    coef_broad <- merge(results, t(fparms.broad))

    output <- merge(DATA, coef_broad, by=c("plot.id", "rev")) %>%
      mutate(hojd=ifelse(is.na(hojd), round(1.3+(dia / ((b0+b0.u+b0.v) + (b1*sqrt(hdom)^b2 + (b3+b3.u+b3.v)*sqrt(ddom))*dia))^3,1), hojd))%>%
      select(-c(plot.id,b0.u,b3.u,b0.v,b3.v,b0,b1,b2,b3))

    return(output)

  }
  else if(any(DATA$t_l < 4)){
    stop("Function cannnot be applied to this species. Try species_specific or main_species functions", call. = FALSE)
  }


}



#' Tree data from the long-term forest experiments in Sweden
#'
#' A subset of the data used for the development of multilevel mixed-effect height functions for the Swedish long-term forest experiments
#' Report ...
#'
#' @docType data
#' @usage data(Treeht)
#'
#' @format ## `Treeht`
#' A data frame with 334 rows and 9 columns:
#' \describe{
#'   \item{yta}{Plot number}
#'   \item{rev}{Revision number}
#'   \item{t_l}{Species code}
#'   \item{hojd}{Sampled tree height in meters}
#'   \item{dia}{Diameter at breast height in centimeter}
#'   \item{hdom}{Height of the tree with the largest diameter per plot regardless of the species}
#'   \item{ddom}{Largest diameter per plot regardless of the species}
#'   \item{QMD}{Quadratic mean diameter in centimeter per plot}
#'   \item{G}{Basal area in meter square per hectare}
#' }
#' @source <doi:10.1016/j.foreco.2023.120843>
#' @examples
#' data(Treeht)
#'
"Treeht"



#' Tree data from the long-term forest experiments in Sweden
#'
#' A subset of the data used for the development of multilevel mixed-effect height functions for the Swedish long-term forest experiments. For illustration purpose, the species code was changed.
#' Report ...
#'
#' @docType data
#' @usage data(broad)
#'
#' @format ## `broad`
#' A data frame with 117 rows and 9 columns:
#' \describe{
#'   \item{yta}{Plot number}
#'   \item{rev}{Revision number}
#'   \item{t_l}{Species code}
#'   \item{hojd}{Sampled tree height in meters}
#'   \item{dia}{Diameter at breast height in centimeter}
#'   \item{hdom}{Height of the tree with the largest diameter per plot regardless of the species}
#'   \item{ddom}{Largest diameter per plot regardless of the species}
#'   \item{QMD}{Quadratic mean diameter in centimeter per plot}
#'   \item{G}{Basal area in meter square per hectare}
#' }
#' @source <doi:10.1016/j.foreco.2023.120843>
#' @examples
#' data(broad)
#'
"broad"
