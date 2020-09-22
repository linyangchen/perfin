#Lin Yangchen
#Laboratory of Computational Philately
#Coconut Academy of Sciences
#September 2020




rm(list = ls())

#change path to your own folder
setwd("/Users/yangchen/philately/perfin/analysis")

require(stringr)

require(png)
require(raster)

require(plot3D)
require(lattice)
require(fBasics)
require(colorspace)

require(EBImage)
require(spatstat)
#require(dbscan)


require(extrafont)
#font_import()
#loadfonts(device = 'pdf')
#fonts()
font <- 'Courier New'


source("functions.R")

#===============================
#SETTINGS

#seed the pseudorandom number generator (optional)
set.seed(12345)


#put in a subfolder the images of perfin samples to be matched
matchfolder <- "match"



#brightness thresholds for hole detection

rval <- 50
gval <- 50
bval <- 50

rgblevels <- 256
rthresh <- rval/rgblevels
gthresh <- gval/rgblevels
bthresh <- bval/rgblevels


#radius threshold for deleting stray objects
speckthresh <- 25


#calibrated resolution of imaging system
mmperpixel <- 35/5123
pixelsonemm <- 1/mmperpixel


#how many sets of randomly distributed holes for ANN comparison
randomsets <- 250



#plotting
mgp <- c(2,0.5,0)
tck <- -0.01
intck <- 0.01


#===============================

#retrieve file names to process
files <- list.files()[str_detect(list.files(), ".png")]
files <- str_remove(files, ".png")

perfindata <- list(NULL)
for(counter in 1:length(files))
{
    perfindata[[counter]] <- perfin(files[counter])
    
    print(paste("Analysis:", counter, "of", length(files), "perfins done"))
}

names(perfindata) <- files
save(perfindata, file = "perfindata.RData")

#load("perfindata.Rdata")



#table of everything
perfintable <- NULL
for(counter in 1:length(perfindata))
{
    perfintable <- rbind(perfintable,
    data.frame(
    Perfin = names(perfindata)[counter],
    Capheight = round(perfindata[[counter]]$capheight, 3),
    Holes = perfindata[[counter]]$numholes,
    Diameter = round(perfindata[[counter]]$chardiam, 3),
    Area = round(perfindata[[counter]]$holearea, 3),
    Interhole = round(perfindata[[counter]]$ann[1], 3),
    Congestion = round(perfindata[[counter]]$chardiam/perfindata[[counter]]$ann[1], 3),
    Stroke = round(perfindata[[counter]]$strokeclar, 3),
    Readability = round(perfindata[[counter]]$rmse, 3)
    )
    )
}




perfintable[which(perfintable$Perfin == "HSBC" | perfintable$Perfin == "YSB"),
which(names(perfintable) == "Capheight" | names(perfintable) == "Stroke")] <- NA

#export table
write.csv(perfintable, file = "perfintable.csv", row.names = F)




#panel plot of segmented perfins
pdf("perfins.pdf", family = font, height = 5, width = 5)
par(mfrow = c(3,3))

for(counter in 1:length(files))
{
    y <- perfindata[[counter]]$objects
    plot(EBImage::rotate(flip(y), 90))
    
    #if(names(perfindata)[counter] != "HSBC" &
    #names(perfindata)[counter] != "YSB")
    if(counter == 2 | counter == 6)
    {
        barheight <- dim(y)[2]*0.975
        midpt <- dim(y)[1]/2
        startpt <- midpt - 5*pixelsonemm
        endpt <- midpt + 5*pixelsonemm
        majticks <- seq(startpt, endpt, length = 11)
        medticks <- seq(startpt, endpt, length = 21)
        minticks <- seq(startpt, endpt, length = 101)
        majtickheight <- dim(y)[2]*0.960
        medtickheight <- dim(y)[2]*0.965
        mintickheight <- dim(y)[2]*0.970
        textheight <- dim(y)[2]*0.9

        #horizontal bar
        segments(startpt, barheight, endpt, barheight, col = "white", lwd = 0.25)
    
        #ticks
        segments(majticks, barheight, majticks, majtickheight, col = "white", lwd = 0.25)
        segments(medticks, barheight, medticks, medtickheight, col = "white", lwd = 0.25)
        segments(minticks, barheight, minticks, mintickheight, col = "white", lwd = 0.25)

        #labels
        text(midpt, textheight, "10 mm", col = "white")
    }    
}

#fill in 9th panel
plot(Image())

dev.off()
embedFonts("perfins.pdf")




#=============================
#Perfin matching

setwd(matchfolder)
#setwd("matchHSBC")

matchfiles <- list.files()[str_detect(list.files(), ".png")]
matchfiles <- str_remove(matchfiles, ".png")



matchperfindata <- list(NULL)
for(counter in 1:length(matchfiles))
{
    matchperfindata[[counter]] <- perfin(matchfiles[counter])
    
    print(paste("Analysis:", counter, "of", length(matchfiles), "perfins done"))
}

names(matchperfindata) <- matchfiles
save(matchperfindata, file = "matchperfindata.RData")

#load("matchperfindata.Rdata")





filepairs <- combn(matchfiles, 2)

matchdata <- list(NULL)
for(counter in 1:ncol(filepairs))
{
    ind1 <- which(names(matchperfindata) == filepairs[1,counter])
    ind2 <- which(names(matchperfindata) == filepairs[2,counter])
    
    matchdata[[counter]] <- eutrans(
    names(matchperfindata)[c(ind1, ind2)],
    list(matchperfindata[[ind1]]$objects, matchperfindata[[ind2]]$objects),
    list(matchperfindata[[ind1]]$holepos, matchperfindata[[ind2]]$holepos),
    list(matchperfindata[[ind1]]$centroid, matchperfindata[[ind2]]$centroid))
    
    
    if(ncol(filepairs)==1)
    {
        print(paste("Matching:", counter, "of", ncol(filepairs), "pair done"))
    } else
    {
        print(paste("Matching:", counter, "of", ncol(filepairs), "pairs done"))
    }
}

save(matchdata, file = "matchdata.RData")

#load("matchdata.RData")




misaligns <- NULL
for(counter in 1:length(matchdata))
{
    misaligns <- c(misaligns, matchdata[[counter]]$misalign)
}



#=======================================================
#sample-wise distribution of misalignment percentages

bws <- c(0.5, 1.5)

samplewise <- NULL

for(supercounter in 1:length(bws))
{
    for(counter in 1:length(matchfiles))
    {
        ind <- which(filepairs[1,] == matchfiles[counter] | filepairs[2,] == matchfiles[counter])
        sampledistrib <- stats::density(misaligns[ind], bw = bws[supercounter])
        samplewise <- rbind(samplewise,
        data.frame(
        file = rep(matchfiles[counter], length(sampledistrib$x)),
        bw = bws[supercounter],
        x = sampledistrib$x,
        y = sampledistrib$y
        )
        )
    }
}

samplewise$bw <- paste("bandwidth", samplewise$bw)







#=======================================================
#clustering


#dbscan()





#================
#plots

pdf("anomaly.pdf", family = font, height = 5, width = 5)

par(mar = c(3,3,3,1), pty = "s")
plot(rev(sort(misaligns)),
ylab = "Misalignment (%)", xlab = "Ranked perfin sample pairs",
axes = F, mgp = mgp, tck = intck,
main = "Misalignment distribution")

axis(1, mgp = mgp, tck = intck)
axis(2, mgp = mgp, tck = intck, las = 2)
box(bty = "L")



xyplot(y~x | bw, groups = file, data = samplewise, type = "l",
col = adjustcolor(timPalette(length(matchfiles)), alpha.f = 0.5),
xlab = "Misalignment (%)", ylab = "Probability density",
aspect = 1, as.table = T,
scales = list(x = list(at = c(5,10,15))),
par.settings = list(strip.background = list(col = "transparent")),
main = "Sample-wise distributions")



dev.off()
embedFonts("anomaly.pdf")





#back to root directory
#change path to your own folder
setwd("/Users/yangchen/philately/perfin/analysis")








