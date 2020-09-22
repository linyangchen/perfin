#Lin Yangchen
#Laboratory of Computational Philately
#Coconut Academy of Sciences
#September 2020


perfin <- function(file)
{
    x <- readPNG(paste(file, ".png", sep = ""))
    
    
    #==================================================
    #OBJECT SEGMENTATION
    #==================================================
    
    rind <- which(x[,,1] < rthresh)
    gind <- which(x[,,2] < gthresh)
    bind <- which(x[,,3] < bthresh)
    ind <- intersect(rind, intersect(gind, bind))
    
    #sum pixel values of three colour channels
    x <- apply(x, c(1,2), sum)
    
    threshold <- x
    threshold[ind] <- 1
    threshold[-ind] <- 0

    y <- Image(threshold)

    #segment and label objects
    y <- EBImage::opening(y, makeBrush(size = 1, shape='Gaussian'))
    y <- fillHull(y)
    y <- bwlabel(y)

    #hole positions and sizes
    holepos <- computeFeatures.moment(y)
    holesize <- computeFeatures.shape(y)
    
    
    srmcol <- which(colnames(holesize) == "s.radius.mean")
    mcxcol <- which(colnames(holepos) == "m.cx")
    mcycol <- which(colnames(holepos) == "m.cy")
    
    #exclude stray specks
    ind <- which(holesize[,srmcol] < speckthresh)
    
    y <- rmObjects(y, ind)
    holesize <- holesize[-ind,]
    holepos <- holepos[-ind,]

    numholes <- nrow(holepos)

    #actual hole diameters
    holesize <- cbind(holesize, holesize[,srmcol]*2*mmperpixel)
    colnames(holesize)[ncol(holesize)] <- "actualdiam"

    adcol <- which(colnames(holesize) == "actualdiam")

    #distribution of hole sizes
    sizedist <- stats::density(holesize[,adcol])
    chardiam <- sizedist$x[which.max(sizedist$y)]

    
    
    
    #total hole area
    pixelarea <- mmperpixel^2
    holearea <- length(which(imageData(y)!=0))*pixelarea
    
    
    
    
    
    
    #pixel coordinates of centroid
    xcentroid <- mean(holepos[,mcxcol])
    ycentroid <- mean(holepos[,mcycol])
    
    #zero the centroid
    holepos <- cbind(holepos, holepos[,mcxcol] - xcentroid)
    holepos <- cbind(holepos, holepos[,mcycol] - ycentroid)
    colnames(holepos)[(ncol(holepos)-1):ncol(holepos)] <- c("zx", "zy")
    
    zxcol <- which(colnames(holepos) == "zx")
    zycol <- which(colnames(holepos) == "zy")
    
    
    
    #perfin bounding box
    top <- max(holepos[,zxcol])
    bottom <- min(holepos[,zxcol])
    left <- min(holepos[,zycol])
    right <- max(holepos[,zycol])
    
    #cap height
    capheight <- top - bottom
    capheight <- capheight*mmperpixel




    #==================================================
    #NEAREST NEIGHBOURS
    #==================================================

    nndistances <- nndist(holepos[,c(zxcol, zycol)], k = 1:(nrow(holepos)-1))
    ann <- apply(nndistances, 2, mean)*mmperpixel

    
    
    randompoints <- NULL
    for(nncounter in 1:randomsets)
    {
        #generate random points
        random1 <- runif(numholes, min = bottom, max = top)
        random2 <- runif(numholes, min = left, max = right)
    
        randomnn <- nndist(random1, random2, k = 1:(nrow(holepos)-1))
        randomann <- apply(randomnn, 2, mean)*mmperpixel
        
        randompoints <- rbind(randompoints, randomann)
    }
    
    randomavg <- apply(randompoints, 2, mean)
    rmse <- sqrt(mean((randomavg - ann)^2))
    
    
    
    
    strokeclar <- capheight/(chardiam*ann[1])
    names(strokeclar) <- NULL
    
    
        
    
    

    #==================================================
    #PLOTS
    #==================================================
    
    gfxfile <- paste("plots_", file, ".pdf", sep = "")
    pdf(gfxfile, height = 5, width = 5, family = font)
    
    par(pty = "s")
    
    allticks <- seq(0, 1, length = dim(raster(x))[1] + 1)
    ticks <- seq(0, dim(raster(x))[1], by = 500)
    

    #brightness map
    par(mar = c(0,3,0,1))
    plot(raster(x),
    col = timPalette(1000),
    axes = F, box = F,
    mgp = mgp, tck = tck,
    ylab = "Pixels"
    #,main = paste(file, "perfin\nRelative pixel brightness map")
    )

    axis(2, mgp = mgp, tck = tck, at = allticks[ticks + 1], labels = ticks)
    
    
    
    #thresholded image
    par(mar = c(0.1,3,2,0))
    plot(raster(threshold),
    axes = F, box = F,
    mgp = mgp, tck = tck,
    ylab = "Pixels",
    col = grey.colors(2, start = 0, end = 1), legend = F)
    
    axis(2, mgp = mgp, tck = tck, at = allticks[ticks + 1], labels = ticks)
    
    text(0.5, 0.025,
    paste("Channel thresholds", rval, gval, bval),
    col = "white")
    

    
    
    
    
    

    par(mar = c(2,2,2,2))

    #segmented image
    plot(EBImage::rotate(flip(y), 90))
    text(dim(y)[1]/2, dim(y)[2]*0.025, "Segmented objects", col = "white")

    #graduated scale 10 mm
    scalebar(y)
        
        
    
    
    
    

    #hole positions
    
    xdata <- holepos[,zycol]*mmperpixel
    ydata <- -holepos[,zxcol]*mmperpixel
    
    #same scaling for both axes
    wholerange <- range(c(xdata, ydata))

    par(mar = c(3,3,3,1))
    plot(xdata, ydata,
    xlab = "Distance from centroid (mm)", ylab = "Distance from centroid (mm)",
    xlim = wholerange, ylim = wholerange,
    axes = F, mgp = mgp, tck = tck,
    col = "transparent",
    main = "Hole centre coordinates")
    grid(50,50, lwd = 0.25, lty = 1, col = "black")
    points(0, 0, pch = 3, cex = 5, col = "red")
    points(xdata, ydata,
    pch = 13, cex = 1.5, col = "black")
    axis(1, mgp = mgp, tck = tck)
    axis(2, mgp = mgp, tck = tck, las = 2)
    box(bty = "L")
    
    

    #distribution of hole sizes
    plot(sizedist$x, sizedist$y, type = "l",
    axes = F, mgp = mgp, tck = intck,
    xaxs = "i",
    xlab = "Hole diameter (mm)",
    ylab = "Probability density",
    main = "Hole diameter distribution")
    
    abline(v = chardiam, lty = 3, lwd = 1)
    text(chardiam*0.99, max(sizedist$y)/3,
    paste(
    round(chardiam, 3)
    , "mm"), srt = 90)
    
    axis(1, mgp = mgp, tck = intck)
    axis(2, mgp = mgp, tck = intck, las = 2)
    box(bty = "L")



    #ANN signature
    par(xpd = T)
    matplot(t(randompoints), type = "l", lwd = 0.18, lty = 1, col = "black",
    xlab = expression(paste(italic(n), 'th nearest neighbour')),
    ylab = "Average distance to neighbour (mm)",
    axes = F, mgp = mgp, tck = intck,
    xaxs = "i",
    main = "Average nearest neighbours")
    
    points(ann,
    col = "red"
    )

    lines(ann, col = "red")
    
    axis(1, mgp = mgp, tck = intck)
    axis(2, mgp = mgp, tck = intck, las = 2)
    box(bty = "L")
    
    legend(ncol(randompoints)*0.16, max(randompoints)*0.9,
    legend = c("perfin holes", "random holes"),
    col = c("red", "black"),
    lwd = c(1, 1), pch = c(1, NA),
    bty = "n")
    
    text(ncol(randompoints)*0.64, max(randompoints)*0.16, paste("RMSE =",
    round(rmse, 3)
    ))
    
    

    

    dev.off()
    embedFonts(gfxfile)




    
    return(
    list(
    file = file,
    rthresh = rthresh,
    gthresh = gthresh,
    bthresh = bthresh,
    image = x,
    threshold = threshold,
    objects = y,
    holepos = holepos,
    centroid = c(xcentroid, ycentroid),
    holesize = holesize,
    numholes = numholes,
    holearea = holearea,
    chardiam = chardiam,
    capheight = capheight,
    strokeclar = strokeclar,
    nndistances = nndistances,
    ann = ann,
    randompoints = randompoints,
    rmse = rmse
    )
    )


}







#==================================================
#EUCLIDEAN TRANSFORM
#==================================================


eutrans <- function(file, obj, holepos, centroids)
{        
    #translate second image to match centroids
    trans1 <- centroids[[1]][1] - centroids[[2]][1]
    trans2 <- centroids[[1]][2] - centroids[[2]][2]
    trans <- c(trans1, trans2)
    obj[[2]] <- translate(obj[[2]], v = trans)

    

    zxcol <- which(colnames(holepos[[1]]) == "zx")
    zycol <- which(colnames(holepos[[1]]) == "zy")
    pos1 <- holepos[[1]][,c(zxcol, zycol)]
    pos2 <- holepos[[2]][,c(zxcol, zycol)]

    
    #Find corresponding pairs of hole centres using PCA

    pca1 <- prcomp(pos1)
    pca2 <- prcomp(pos2)


    #flip orientation of PC1 loading if necessary

    if(pca1$rotation[1,1] < 0)
    {
        pca1$rotation[,1] <- -pca1$rotation[,1]
        pos1 <- pos1 %*% pca1$rotation
        pos1 <- rotrad(pos1, pi)
    } else
    {pos1 <- pos1 %*% pca1$rotation}
    
    if(pca2$rotation[1,1] < 0)
    {
        pca2$rotation[,1] <- -pca2$rotation[,1]
        pos2 <- pos2 %*% pca2$rotation
        pos2 <- rotrad(pos2, pi)
    } else
    {pos2 <- pos2 %*% pca2$rotation}
    

    
    
    
    #corresponding pairs of objects
    orderind <- NULL
    for(corrcounter in 1:nrow(pos1))
    {
        xdist <- abs(pos1[corrcounter,1] - pos2[,1])
        ydist <- abs(pos1[corrcounter,2] - pos2[,2])
        xydist <- xdist + ydist
        orderind <- c(orderind, which.min(xydist))    
    }


    function(){
    #visually check correspondence
    gfxfile = paste("correspondence_", file[1], "-", file[2], ".pdf", sep = "")
    pdf(gfxfile, height = 5, width = 9, family = font)
    par(mfrow = c(1,2), pty = "s", mar = c(0.1,0.1,2,0.1))
    plot(holepos[[1]][,c(zxcol, zycol)], col = "transparent",
    main = file[1], axes = F)
    box(bty = "O")
    text(holepos[[1]][,zxcol], holepos[[1]][,zycol], 1:nrow(holepos[[1]]))
    plot(holepos[[2]][orderind,c(zxcol, zycol)], col = "transparent",
    main = file[2], axes = F)
    box(bty = "O")
    text(holepos[[2]][orderind,zxcol], holepos[[2]][orderind,zycol], 1:nrow(holepos[[2]]))
    dev.off()
    embedFonts(gfxfile)
    } #END OF COMMENT-OUT FUNCTION
    
    
    
    #covariance matrix of corresponding object centres
    covmat <- cov(holepos[[1]][,c(zxcol, zycol)], holepos[[2]][orderind,c(zxcol, zycol)])
    
    #singular value decomposition
    svded <- svd(covmat)
    R <- svded$v %*% t(svded$u)
    
    #recover angle from rotation matrix
    rot <- rad2deg(asin(R[1,2]))
    
    obj[[2]] <- EBImage::rotate(obj[[2]], angle = rot, output.origin = centroids[[1]])
    
    
    #misalignment map
    for(mscounter in 1:2)
    {
        obj[[mscounter]] <- imageData(obj[[mscounter]])
        obj[[mscounter]][which(obj[[mscounter]] != 0)] <- 1
    }
    obj3 <- obj[[1]] - obj[[2]]
    
    save(obj3, file =
    paste(
    "map_", file[1], "-", file[2], ".RData", sep = ""
    ))
    
    #percentage misalignment
    bwpixels <- length(which(obj3 == -1 | obj3 == 1))
    holepixels <- sum(length(which(obj[[1]] == 1)), length(which(obj[[2]] == 1)))
    misalign <- bwpixels/holepixels*100






    #plots
    
    
    gfxfile <- paste("matched_", file[1], "-", file[2], ".pdf", sep = "")
    pdf(gfxfile, family = font)
    par(mar = c(0,0,0,0))
    


    #misalignment map
    
    plot(raster(obj3), col = grey.colors(1000, start = 0, end = 1),
    legend = F, axes = F, box = F)
    
    text(0.5, 0.975, "Misalignment map", col = "black")
    
    text(0.5, 0.025,
    paste(
    round(misalign, 3)
    , "%", sep = ""))
    
    
    
    
    
    
    dev.off()
    embedFonts(gfxfile)
    
    
    
    return(
    list(
    file1 = file[1],
    file2 = file[2],
    correspondence = orderind,
    rotation = rot,
    misalign = misalign
    
    )
    )

    
}











#==================================================
#AUXILIARY FUNCTIONS
#==================================================


rad2deg <- function(rad) {(rad * 180) / (pi)}


rotrad <- function(x, rad)
{    
    rotmat <- matrix(
    nrow = 2,
    c(
    cos(rad),
    -sin(rad),
    sin(rad),
    cos(rad)
    ),
    byrow = T
    )
    
    return(x %*% rotmat)
}




scalebar <- function(y)
{

    #barheight <- dim(y)[2]*0.975
    barheight <- dim(y)[2]
    midpt <- dim(y)[1]/2
    startpt <- midpt - 5*pixelsonemm
    endpt <- midpt + 5*pixelsonemm
    majticks <- seq(startpt, endpt, length = 11)
    medticks <- seq(startpt, endpt, length = 21)
    minticks <- seq(startpt, endpt, length = 101)
    #majtickheight <- dim(y)[2]*0.965
    #mintickheight <- dim(y)[2]*0.970
    majtickheight <- dim(y)[2]*0.985
    medtickheight <- dim(y)[2]*0.990
    mintickheight <- dim(y)[2]*0.995
    #textheight <- dim(y)[2]*0.965
    textheight <- dim(y)[2]*0.995
    labs = 0:10


    #horizontal bar
    #segments(startpt, barheight, endpt, barheight, col = "white")

    #ticks
    segments(majticks, barheight, majticks, majtickheight, col = "white", lwd = 0.5)
    segments(medticks, barheight, medticks, medtickheight, col = "white", lwd = 0.5)
    segments(minticks, barheight, minticks, mintickheight, col = "white", lwd = 0.5)

    #labels
    text(majticks, textheight, labels = labs, col = "white", pos = 3)
    text(majticks[length(majticks)]*1.064, textheight, "mm", col = "white", pos = 3)

}


