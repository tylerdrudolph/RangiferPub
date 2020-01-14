##########################################################################################################################################
## 1. Reduce your dataset to one observation per day per unique collared animal. Here is a function that will do that for you, though you 
##    may want to modify it to accommodate your data set. You will need to have the 'chron' package installed. The input data frame must have 
##    a variable called 'date' that can be coerced to POSIXct format, a variable called 'id' (unique collar ID) and a variable called 'day' of 'Date' format.

  ## Rarify data.frame to one observation per day (reference = 12 noon)
  rarifyGPS = function(df, dt=NULL, id=NULL, reftime = times(c("12:00:00"))) {
  
    if(is.null(dt)) stop("Must provide a date/time field 'dt'")
    if(!is.null(dt) & !inherits(dt, "POSIXct")) stop("'dt' must be in 'POSIXct' format")
  
    df$dt <- dt
  
    if(is.null(id)) stop("Must provide a unique 'id' column")
  
    df$id <- id
  
    ## Produce a 'time' field in 'times' format (package::chron)
    library(chron)
    df$time = times(format(df$dt, format="%H:%M:%S"))
  
    ## Subtract date-times from reftime (noon) to derive median values
    df$diff <- abs(reftime - df$time)
    df$diff <- as.numeric(df$diff)
  
    ## Order date times by id, then day, then diff
    df <- df[order(df$id, df$dt, df$diff), ]
  
    ## Extract the observation(s) situated closest to reftime (noon by default) for each unique id/day
    df.min <- aggregate(df$diff, list(id = df$id, date = as.Date(df$dt)), min, na.rm = TRUE)
      colnames(df.min) <- c("id", "date", "diff")  # rename columns
    df.subset <- subset(df, paste(df$id, as.Date(df$dt), df$diff) %in% paste(df.min$id, df.min$date, df.min$diff))
  
    ## Identify and eliminate exceptions consisting of > 1 observation per id/day
    df.subset.tab <- data.frame(table(df.subset$id, df.subset$date))
      colnames(df.subset.tab) <- c("id", "date", "Freq")
    df.subset.tab <- subset(df.subset.tab, df.subset.tab$Freq > 1)
  
    ## Eliminate multiple records where applicable
    if(nrow(df.subset.tab) > 0) df.subset <- df.subset[!duplicated(paste(df.subset$id, df.subset$date)), ]
  
    return(df.subset)
  
  }

#################################################################################################################################################
## 2. Add a variable corresponding to Julian day that begins anew every year (i.e. range = c(1, 365)). This kind of thing can be tricky when you 
##    have multiple years of data and potentially leap years, but the key is proper application of the 'origin' argument. Here's a generic example 
##    of how you might want to do that assuming you have a data frame ('df') consisting of multiple years of GPS telemetry monitoring:

  attach(df)
  
  df <- do.call(rbind, lapply(split(df, year, drop=T), function(x) {
    x$jday <- julian(x$date, origin=as.Date(paste(format(x$date, "%Y)[1]-1, "-12-31", sep="")))
    paste("01-01-", format(x$date, "%Y)[1], sep=""))
    return(x)
  }))

#################################################################
## 3. Smooth data to eliminate random noise (if desired)

  ## order by id and then julian day
  df = df[order(df$id, df$judy),]

  ## smooth with moving window (I chose 4 days, you may prefer to smooth more or less)
  library(zoo)
  gps = do.call(rbind, lapply(split(df, df$id, drop=TRUE), function(x) {
    x$movwin = rollapply(x$mov.rate, 4, mean, fill=NA)
    return(x)
  }))
  row.names(df) = 1:nrow(df)

  # log-transform movement rate (assuming this makes sense in your case)
  df$l_movwin = log(df$movwin)

#########################################################################################################################################
## 4. At this point, plot the average daily movement rate, and potentially net displacement, over the course of a year 
##    from your sample population. If you don't yet have those parameters in your data then you can generate them using this function I borrowed and modified:
  
  ###########################################################################################
  ## Calculate movement parameters as per Clement Calenge (adehabitatLT::as.ltraj)

  movparms = function(xy, date, id) {
    x <- xy[, 1]
    y <- xy[, 2]
    res <- split(data.frame(x = x, y = y, date = date), id)
    res <- lapply(res, function(y) y[order(y$date), , drop = FALSE])
    foo <- function(x) {
      x1 <- x[-1, ]
      x2 <- x[-nrow(x), ]
      dist <- c(sqrt((x1$x - x2$x)^2 + (x1$y - x2$y)^2), NA)
      R2n <- (x$x - x$x[1])^2 + (x$y - x$y[1])^2
      dt <- c(unclass(x1$date) - unclass(x2$date), NA)
      dx <- c(x1$x - x2$x, NA)
      dy <- c(x1$y - x2$y, NA)
      abs.angle <- ifelse(dist < 1e-07, NA, atan2(dy, dx))
      so <- cbind.data.frame(dx = dx, dy = dy, dist = dist, 
                  dt = dt, R2n = R2n, abs.angle = abs.angle)
      return(so)
    }
  
  speed <- lapply(res, foo)
  res <- lapply(1:length(res), function(i) cbind(res[[i]], speed[[i]]))
  ang.rel <- function(df, slspi = slsp) {
    ang1 <- df$abs.angle[-nrow(df)]
    ang2 <- df$abs.angle[-1]
    dist <- c(sqrt((df[-nrow(df), "x"] - df[-1, "x"])^2 + 
                     (df[-nrow(df), "y"] - df[-1, "y"])^2), NA)
    wh.na <- which(dist < 1e-07)
    if (length(wh.na) > 0) {
      no.na <- (1:length(ang1))[!(1:length(ang1)) %in% wh.na]
      for (i in wh.na) {
        indx <- no.na[no.na < i]
        ang1[i] <- ifelse(length(indx) == 0, NA, ang1[max(indx)])
      }
    }
    res <- ang2 - ang1
    res <- ifelse(res <= (-pi), 2 * pi + res, res)
    res <- ifelse(res > pi, res - 2 * pi, res)
    return(c(NA, res))
  }
  rel.angle <- lapply(res, ang.rel)
  res <- lapply(1:length(res), function(i) data.frame(id=levels(as.factor(id))[i], res[[i]], rel.angle = rel.angle[[i]]))
  return(do.call(rbind, res))
  
}

######################################################################################################################################################
## 5. Something like this should then produce what you need to plot the mean relationship of interest. Try also with net displacement, taking care to 
##    set your spatial 'anchor point' of reference to something biologically meaningful:
  
  # means
  mean.daily.movrates <- sapply(1:365, function(j) mean(df$mov.rate[df$jday==j], na.rm=T))
  # std deviations
  sd.daily.movrates <- sapply(1:365, function(j) sd(df$mov.rate[df$jday==j], na.rm=T))

  # now visualize (you can also generate sample confidence intervals uing sd.daily.movrates)
  plot(mean.daily.movrates~1:365, type="l")

#########################################################################################################################################################
## 6. By now you will have a good idea whether or not the calving period can be discerned from your data set. Try transforming the variables (e.g. log) 
##    for a better result. At this point you can try a population-level RE-EM tree model (random effects expectation maximization), and seek out the 
##    partitions that most closely approximate the transitions associated with the calving period that you are presumably able to observe in the previous plots.

  library(rpart)
  rmod = REEMtree(formula=l_movwin ~ julian.day, data=df[!is.na(gps$l_movwin),], random= ~ 1 | id/year, method="REML")

  ## Function to extract ablines from rpart model  
  rtreeInflection = function(rmod) {  
    micro<-rmod$frame[rmod$frame[,1]!="<leaf>",]
    micro<-cbind.data.frame(rmod$splits[,4],micro$complexity)
    micro <- micro[order(micro[,2], decreasing=T),]
    names(micro) = c("abline","complexity")
    return(micro)
  }

  # Derive inflection points and plot
  plot(mean.daily.movrates~1:365, type="l")
  abline(v=rtreeInflection(rmod))
  
#######################################################################################################################################################
## 7. You may also wish to fit the recursive partitioning model to vital rates from individual collared animals (individual-based delineation), which
##    can be accomplished using rpart (see example for id=="C001"). 

  rmod = rpart(l_movwin ~ julian.day, data=df[!is.na(gps$l_movwin) & df$id=="C001",])