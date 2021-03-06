# Error curves
get_bi <- function(data, model, outDir, prefix="", suffix="", saveplot=TRUE, test=FALSE, cv=TRUE, ...){
    data <- model[[1]]
    if(saveplot){
    dots <- list(...)
    if(is.null(dots$width)) w <- 1600 else w <- dots$width
    if(is.null(dots$height)) h <- 1600 else h <- dots$height
    if(is.null(dots$res)) r <- 200 else r <- dots$res
    if(prefix!="") prefix <- paste0(prefix, "_")
    if(suffix!="") suffix <- paste0("_", suffix)
    png(file.path(outDir, paste0(prefix, "gbm_ERR", suffix, ".png")), width=w, height=h, res=r)
    if(test) gbm.perf(data, method="test")
    if(cv) gbm.perf(data, method="cv")
    dev.off()
    }
    list(
        Test=gbm.perf(data, method="test", plot.it=FALSE),
        CV=gbm.perf(data, method="cv", plot.it=FALSE)
    )
}

# Relative influence
get_ri <- function(data, model, n.trees, outDir, prefix="", suffix="", saveplot=TRUE, ...){
    data <- model[[1]]
    n <- n.trees[[1]]
    cbpal <- c("#8B4500", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    ri.test <- if(length(n$Test)) cbind(Method="Test", summary(data, n.trees=n$Test, order=F, plotit=F)) else NULL
    ri.cv <- if(length(n$CV)) cbind(Method="CV", summary(data, n.trees=n$CV, order=F, plotit=F)) else NULL
    ri <- data.table(rbind(ri.test, ri.cv))
    setnames(ri, c("Method","Predictor","RI"))
    if(saveplot){
    dots <- list(...)
    if(is.null(dots$width)) w <- 1600 else w <- dots$width
    if(is.null(dots$height)) h <- 1600 else h <- dots$height
    if(is.null(dots$res)) r <- 200 else r <- dots$res
    if(prefix!="") prefix <- paste0(prefix, "_")
    if(suffix!="") suffix <- paste0("_", suffix)
    png(file.path(outDir, paste0(prefix, "gbm_RIcv", suffix, ".png")), width=w, height=h, res=r)
    g <- ggplot(filter(ri, Method=="CV"), aes(Predictor, RI)) + geom_bar(stat="identity", position="dodge") +
    theme_bw() + theme(legend.position="bottom") +
    ggtitle(paste("Predictor relative influence:\nStart of", prefix, "growing season"))
    print(g)
    dev.off()
    png(file.path(outDir, paste0(prefix, "gbm_RI", suffix, ".png")), width=w, height=h, res=r)
    g <- ggplot(ri, aes(Method, RI, fill=Predictor)) + geom_bar(stat="identity", position="dodge") +
    scale_fill_manual(values=cbpal[-2]) + theme_bw() + theme(legend.position="bottom") +
    ggtitle(paste("Predictor relative influence:\nStart of", prefix, "growing season"))
    print(g)
    dev.off()
    }
    ri
}

# Time series
get_preds <- function(data, model, newdata, n.trees, type.err="test", grp=NULL){
    m <- model[[1]]
    if(is.null(grp)) grp <- as.character(groups(data))
    for(i in 1:length(grp)) newdata <- filter_(newdata, .dots=list(paste0(grp[i], "==\'", data[[grp[i]]][1], "\'")))
    if(!is.numeric(n.trees)){
        n.trees <- if(type.err=="test") n.trees[[1]]$Test else if(type.err=="cv") n.trees[[1]]$CV else stop("type.err must be 'test' or 'cv'.")
    }
    predict(m, newdata=newdata, n.trees=n.trees)
}

# make data table of TDD map values
make_TDD_dt <- function(d, extractBy, y, keep.y=FALSE, years=NULL, dem, ...){
  inner_fun <- function(d, extractBy, y, keep.y=FALSE, years=NULL){
    if(is.null(years)) years <- 1:nlayers(d)
    d[d <= 1] <- NA
    r <- calc(y, mean)
    cells <- which(!is.na(r[]))
    d <- projectRaster(d, r) %>% resample(r, method="bilinear") %>% mask(r) %>% raster::extract(extractBy, cellnumbers=TRUE)
    s <- raster::extract(y, extractBy)
    d <- rbindlist(lapply(1:length(d),
      function(i, tdd, sos, years, eco, shp, mask){
        cells <- as.numeric(tdd[[i]][,1])
        xy <- xyFromCell(mask, cells)
        x <- data.table(Region=eco[i], Year=rep(years, each=nrow(tdd[[i]])), SOS=NA, TDD=as.numeric(tdd[[i]][,-1]), x=xy[,1], y=xy[,2], Cell=cells)
        idx <- which(x$Year %in% 1982:2010)
        x$SOS[idx] <- as.numeric(sos[[i]])
        x
      }, tdd=d, sos=s, years=years, eco=names(extractBy), shp=extractBy, mask=r))
    cells <- intersect(cells, (group_by(d, Cell) %>% summarise(n=n()) %>% filter(n==length(years)))$Cell)
    filter(d, Cell %in% cells)
  }
  d <- mclapply(d, inner_fun, extractBy=extractBy, y=y, keep.y=keep.y, years=years, ...)
  cells <- sort(unique(unlist(purrr::map(d, ~.x$Cell))))
  d <- d %>% purrr::map2(c("05",10,15,20), ~mutate(.x, Threshold=paste0(.y, "pct"))) %>% bind_rows() %>% filter(Cell %in% cells) %>% group_by(Region, Year, Threshold) %>%
    mutate(Obs=1:n()) %>% select(Region, Year, Threshold, Obs, x, y, Cell, SOS, TDD) %>% group_by %>% mutate(Elev=raster::extract(dem, cbind(x, y))) %>%
    dcast(Region + Year + Obs + x + y + Cell + Elev + SOS ~ Threshold, value.var="TDD") %>% data.table
  nam <- c("Region", "Year", "Obs", "x", "y", "Cell", "Elev", "SOS", paste0("DOY_TDD", c("05", 10, 15, 20)))
  if(!keep.y){
    d <- select(d, -SOS)
    nam <- nam[-which(nam=="SOS")]
  }
  setnames(d, nam)
  d
}

# Partial dependence
get_pd <- function(data=NULL, source_data, x, y=NULL, model, outDir, vars=NULL, n.vars=length(vars), order.by.ri=TRUE,
                   density.adjust=2, spline.df=12, prefix="", suffix="", saveplot=TRUE, grp=NULL, ...){
    if(is.null(data) & class(model)!="gbm") stop("Must pass a gbm object to 'model' if not using a gbm data table for 'data'.")
    if(!is.null(data)){
        model <- model[[1]]
        max.vars <- length(unique(data$RI[[1]]$Predictor))
        if(is.null(grp)) grp <- as.character(groups(data))
        for(i in 1:length(grp)) source_data <- filter_(source_data, .dots=list(paste0(grp[i], "==\'", data[[grp]][1], "\'")))
    }
    data <- source_data
    stopifnot(length(n.vars) > 0)
    #x <- substitute(x)
    if(is.character(x)) x <- as.name(x)
    #y <- substitute(y)
    if(is.character(y)) y <- as.name(y)
    if(is.null(vars)) vars <- 1:n.vars
    if(is.null(n.vars)) n.vars <- length(vars)
    stopifnot(n.vars <= max.vars)
    dots <- list(...)
    if(is.null(dots$width)) w <- 1600 else w <- dots$width
    if(is.null(dots$height)) h <- 1600 else h <- dots$height
    if(is.null(dots$res)) r <- 200 else r <- dots$res
    if(prefix!="") prefix <- paste0(prefix, "_")
    if(suffix!="") suffix <- paste0("_", suffix)
    facet.formula <- if(is.null(dots$facet.formula)) "~ Var" else dots$facet.formula
    scales <- if(is.null(dots$scales)) "free" else dots$scales
    cols <- if(is.null(dots$ncol)) ceiling(sqrt(n.vars)) else dots$ncol
    color.var <- if(is.null(dots$colour)) NULL else dots$colour
    xlb <- if(is.null(dots$xlab)) as.character(x) else dots$xlab
    ylb <- if(!is.null(dots$ylab)) dots$ylab else if(is.null(y)) "y" else as.character(y)

    # internal support functions
    dtDen <- function(x, n=1000, adj=0.1, out="vector", min.zero=TRUE, diversify=FALSE){
        b <- max(1, 0.05*diff(range(x)))
        z <- density(x, adjust=adj, n=n, from=min(x)-b, to=max(x)+b)
        if(min.zero && any(z$x < 0)) z <- density(x, adjust=adj, n=n, from=0, to=max(x)+b)
        if(out=="vector") return(as.numeric(c(z$x, z$y))) else if(out=="list") return(z)
    }
    get_pd_dt <- function(i, model, order=TRUE){
        n.trees <- gbm.perf(model, method="cv", plot.it=F)
        ri <- summary(model, n.trees=n.trees, order=order, plotit=F)
        if(i > nrow(ri)) stop("vars index exceeds the number of gbm predictor variables.")
        lev <- ri$var
        d <- data.table(plot(model, i.var=i, return.grid=T))
        d <- mutate(d, RI=ri$rel.inf[which(ri$var==names(d)[1])])
        d <- mutate(d, Var=factor(names(d)[1], levels=lev)) %>% setnames(c("x", "y", "RI", "Var"))
        d
    }
    map_range <- function(x, y=NULL){
        x <- (x - min(x))/diff(range(x))
        if(is.null(y)) x else x * diff(range(y)) + min(y)
    }

    # process
    d.pd <- rbindlist(lapply(1:max.vars, get_pd_dt, model=model, order=order.by.ri)) %>% group_by(Var, RI)
    spx <- function(a, b, d) smooth.spline(x=a, y=b, df=d)$x
    spy <- function(a, b, d) smooth.spline(x=a, y=b, df=d)$y
    d.pd <- d.pd %>%
      mutate_(
        x2=lazyeval::interp(~spx(x, y, spline.df), x=as.name("x"), y=as.name("y")),
        y2=lazyeval::interp(~spy(x, y, spline.df), x=as.name("x"), y=as.name("y"))) %>% group_by(RI, add=T) %>%
      do(x=approx(.$x2, .$y2, n=1000)$x, y=approx(.$x2, .$y2, n=1000)$y) %>% unnest %>% arrange(Var)
    lev <- levels(d.pd$Var)
    if(order.by.ri){
        lev <- paste0(unique(d.pd$Var), ": RI = ", round(unique(d.pd$RI), 2))
        d.pd <- mutate(d.pd, Var=factor(paste0(Var, ": RI = ", round(RI, 2)), levels=lev))
    }

    #data <- unnest(select(ungroup(data), -Region)) %>% group_by(Region, Var) %>%
    data <- group_by(data, Region, Var) %>%
      do(Val=dtDen(.$Val, adj=density.adjust, out="list")$x, Prob=dtDen(.$Val, adj=density.adjust, out="list")$y) %>% unnest %>%
      mutate(Var=factor(lev[pmatch(Var, lev, duplicates.ok=TRUE)], levels=lev)) %>% arrange(Var) %>%
      cbind(select(d.pd, x, y)) %>% group_by(Region, Var) %>% mutate(Prob=map_range(Prob, y)) %>% group_by(Var, add=T)
    if(saveplot){
    g <- ggplot(data %>% filter(Var %in% levels(Var)[1:n.vars]) %>% mutate(Ymin=min(Prob)), aes(x=Val)) +
      facet_wrap(as.formula(facet.formula), ncol=cols, scales=scales) + labs(x=xlb, y=ylb)
    if(is.null(color.var)) g <- g + geom_ribbon(aes(ymin=Ymin, ymax=Prob), fill="orange") + geom_line(aes(x=x, y=y), size=1)
    if(!is.null(color.var)) g <- g + geom_ribbon(aes_string(ymin="Ymin", ymax="Prob", colour=color.var, fill=color.var), alpha=0.5) +
        geom_line(aes_string(x="x", y="y", colour=color.var), size=1)
    png(file.path(outDir, paste0(prefix, "gbm_PD", suffix, ".png")), width=w, height=h, res=r)
    print(g)
    dev.off()
    }
    data
}

get_pd_OLD <- function(data, model, outDir, prefix="", suffix="", saveplot=TRUE, ...){
    data <- model[[1]]
    data$sa <- plot.pd.gbm(object=data, saData=NULL, nVar=4)
    if(saveplot){
    dots <- list(...)
    if(is.null(dots$width)) w <- 1600 else w <- dots$width
    if(is.null(dots$height)) h <- 1600 else h <- dots$height
    if(is.null(dots$res)) r <- 200 else r <- dots$res
    if(prefix!="") prefix <- paste0(prefix, "_")
    if(suffix!="") suffix <- paste0("_", suffix)
    png(file.path(outDir, paste0(prefix, "gbm_PD", suffix, ".png")), width=w, height=h, res=r)
    print(sa.pd.plot(data, saData=NULL, nVars=4, den.adj=1, spl.df=12))
    dev.off()
    }
    return()
}

plot.pd.gbm <-
function (object, saData, nVar = 4, res = 1500, rescale = F)
{
    require(gbm)
    rel.inf <- summary(object, plotit = F)[1:nVar, ]
    var <- as.character(rel.inf[, 1])
    plotData <- plot(object, var[1], continuous.resolution = res,
        return.grid = T)
    names(plotData)[1] <- "x"
    if (rescale == T) {
        plotData$x <- rms(saData[which(names(saData) == as.character(rel.inf$var[1]))]) *
            plotData$x
    }
    for (i in 2:nVar) {
        tmp <- plot(object, var[i], continuous.resolution = res,
            return.grid = T)
        names(tmp)[1] <- "x"
        if (rescale == T) {
            tmp$x <- rms(saData[which(names(saData) == as.character(rel.inf$var[i]))]) *
                tmp$x
            print(rms(saData[which(names(saData) == as.character(rel.inf$var[i]))]))
            print(summary(saData[which(names(saData) == as.character(rel.inf$var[i]))]))
        }
        plotData <- rbind(plotData, tmp)
    }
    plotData$var <- var
    plotData$SI <- rep(rel.inf[, 2], each = res)
    label <- paste(rel.inf[, 1], "\n SI = ", signif(rel.inf[,
        2], 3), sep = "")
    label <- factor(rep(label, each = res), levels = label, ordered = T)
    plotData$label <- label
    plotData
}

sa.pd.plot <-
function (plotModel, saData, nVars = 4, den.adj=2, spl.df=12, rescale = F, logem = F)
{
    if (logem == T) {
        plotModel$sa$y <- exp(plotModel$sa$y)
    }
    plotModel$sa$y[plotModel$sa$y < 0] <- 0
    saPlot <- xyplot(y ~ x | label, data = plotModel$sa, subset = SI >=
        rev(sort(unique(SI)))[nVars], groups = var, las = 1,
        panel = function(x, y, groups, subscripts, ...) {
            p.num <- panel.number()
            xData <- as.data.frame(matrix(plotModel$data$x, nrow = length(plotModel$data$y)))
            names(xData) <- plotModel$var.names
            xData <- xData[, groups[p.num]]
            if (rescale == T) {
                xData <- rms(saData[which(names(saData) == groups[p.num])]) *
                  xData
            }
            dens <- density(xData, from = min(xData), to = max(xData),
                adjust = den.adj)#length(xData)/25)
            dens$y <- (dens$y - min(dens$y))/diff(range(dens$y)) *
                diff(range(y)) + min(y)
            dens$y[dens$y <= 0] <- 0.01 * max(dens$y)
            panel.grid(v = -1, h = 0, lty = 1)
            panel.xyplot(dens$x, dens$y, type = "h", col = "palegreen",
                lw = 5, ...)
            panel.xyplot(x, smooth.spline(x, y, df = spl.df)$y, type = "l",
                col = "darkblue", lw = 3, ...)
        }, par.strip.text = list(lines = 2, cex = 0.75), xlab = "",
        ylab = "partial dependence", between = list(x = 1, y = 1),
        scale = list(x = list(relation = "free", cex = 0.8),
            y = list(relation = "free", rot = 0)), layout = c(2,
            nVars/2), as.table = T)
    return(saPlot)
}
