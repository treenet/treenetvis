# treenetvis (v0.1)
# Authors: Roman Zweifel & Micah Wilhelm
# Modified: 09.09.2021

# Parameters to be set by user #####################################################################################
path        <- "~/treenet/treenetvis"                               # path of treenetvis directory
y_sel       <- 2021                                                 # year selected for data download and processing
plot.format <- "pdf"                                                # output format of plots ("pdf" or "jpg")
####################################################################################################################


# PREAMBLE #########################################################################################################

# Install necessary packages if not present
cran.pkgs           <- c("dplyr", "ggplot2", "gridExtra", "lubridate", "raster", "remotes")
cran.pkgs.install   <- cran.pkgs[!(cran.pkgs %in% installed.packages()[,"Package"])]
if (length(cran.pkgs.install))   install.packages(cran.pkgs.install)

# Load packages
library(dplyr)
library(ggplot2)
library(lubridate)

# Set variables and function
setwd(path)
timestart  <- paste0("time >= '", y_sel-1, "-01-01 00:00:00'")
reso       <- 10
lvls.gro <- c("off season", "high GRO", "above avg. GRO", "below avg. GRO", "poor GRO")
lvls.twd <- c("high TWD", "above avg. TWD","below avg. TWD", "low TWD")


stat.smooth <- function(color) {
  ggplot2::stat_smooth(
    data = NULL,
    geom = "smooth",
    position = "identity",
    method = "loess",
    formula = y ~ x,
    se = F,
    n = 365,
    span = 0.25,
    fullrange = F,
    level = 0.95,
    method.args = list(),
    na.rm = T,
    orientation = NA,
    show.legend = F,
    inherit.aes = T,
    color=color) 
}

# Read TreeNet metadata
metadata <- readRDS("metadata.Rds")

# Intitialise dataframes for result aggregation
dfsignal <- setNames(data.frame(matrix(ncol = 22, nrow = 0)),
                     c("Series", "Site_ID", "Site_XCor","Site_YCor","Region","Species","RelGRO","RelGRO_1","RelGRO_2","RelGRO_3","RelGRO_4", "GRO","GRO_1","GRO_2","GRO_3","GRO_4", "TWD","TWD_1","TWD_2","TWD_3","TWD_4","last_doy"))
dfsignalactual <- dfsignal


# PROCESS DATA #####################################################################################################

# Plot GRO and TWD of selected year and historical range ----
for (ii in 1:dim(metadata)[1]) {
  
  # Download best available data set with LM data combined with new L2 data
  data_L2M <- readRDS(paste0(metadata$Seriesname[ii],".Rds"))
  
  # Check latest data point in data_L2M
  lastdata_L2Mpoint <- as.numeric(yday(data_L2M$ts[length(data_L2M$ts)]))
  doy_today <- as.numeric(yday(Sys.Date()))
  writeLines(paste0("Latest doy in data set ",metadata$Seriesname[ii],": ",lastdata_L2Mpoint,". Today is DOY: ", doy_today))
  
  data_L2M$series <- metadata$Seriesname[ii]
  data_L2M <- data_L2M %>%
    mutate(year=as.numeric(format.Date(ts, "%Y")),
           doy=as.numeric(yday(ts)),
           hour=as.numeric(format.Date(ts, "%H")),
           Time=format.Date(ts, "%H:%M"),
           dec.date = decimal_date(ts))
  
  data_L2M <- data_L2M %>%
    mutate(dec.time = as.numeric(dec.date)-as.numeric(year))  
  
  
  myplot <- ggplot(data_L2M)+
    geom_point(aes(ts,gro_yr, colour=year), na.rm=T)
  
  myplot %+% subset(data_L2M, year %in% c("1997","2008"))
  
  data_L2M <- data_L2M %>%
    mutate(gro_yr = case_when(year == 1997 ~ gro_yr+750, year > 1997 ~ gro_yr ))
  
  # Summarise data  
  
  df <- data_L2M %>%
    mutate(year=format.Date(ts, "%Y"),
           doy=yday(ts),
           hour=format.Date(ts, "%H"),
           Time=format.Date(ts, "%H:%M"),
           dec.date = decimal_date(ts))
  
  df <- df %>%
    mutate(dec.time = as.numeric(dec.date)-as.numeric(year))
  
  
  # Summarise historic data with all years before the selected year
  dfhist <- df %>%
    filter(year < y_sel)
  
  grostart.med <- median(dfhist$gro_start, na.rm=T)
  grostart.0 <- quantile(dfhist$gro_start, probs=0.05, na.rm=T)
  grostart.25 <- quantile(dfhist$gro_start, probs=0.25, na.rm=T)
  grostart.75 <- quantile(dfhist$gro_start, probs=0.75, na.rm=T)
  grostart.100 <- quantile(dfhist$gro_start, probs=0.95, na.rm=T)
  
  df.avg <- dfhist %>%
    group_by(doy) %>%
    summarise(GRO.med=median(gro_yr, na.rm=T),
              GRO.0=quantile(gro_yr, probs=0.05, na.rm=T),
              GRO.25=quantile(gro_yr, probs=0.25, na.rm=T),
              GRO.75=quantile(gro_yr, probs=0.75, na.rm=T),
              GRO.100=quantile(gro_yr, probs=0.95, na.rm=T),
              TWD.med=median(twd, na.rm=T),
              TWD.0=quantile(twd, probs=0.05, na.rm=T),
              TWD.25=quantile(twd, probs=0.25, na.rm=T),
              TWD.75=quantile(twd, probs=0.75, na.rm=T),
              TWD.100=quantile(twd, probs=0.95, na.rm=T)
    )
  
  df.avg <- df.avg %>%
    mutate(time = doy)
  
  # Align to 365 days
  df.avg <- df.avg[df.avg$doy<=365,]
  df.365 <- data.frame(doy = c(1:365)) 
  df.avg <- full_join(df.365,df.avg)
  
  # Calculate smoothed curves with ggplot and extract the data with ggplot_build
  options(warn=-1)
  p <- ggplot(data = df.avg, aes(x = doy, y = GRO.0)) +
    stat.smooth(color="darkgreen") 
  Smooth_GRO.0 <- ggplot_build(p)$data[[1]][["y"]]
  
  p <- ggplot(data = df.avg, aes(x = doy, y = GRO.25)) +
    stat.smooth(color="darkgreen") 
  Smooth_GRO.25 <- ggplot_build(p)$data[[1]][["y"]]
  
  p <- ggplot(data = df.avg, aes(x = doy, y = GRO.med)) +
    stat.smooth(color="darkgreen") 
  Smooth_GRO.med <- ggplot_build(p)$data[[1]][["y"]]
  
  p <- ggplot(data = df.avg, aes(x = doy, y = GRO.75)) +
    stat.smooth(color="lightgreen")
  Smooth_GRO.75 <- ggplot_build(p)$data[[1]][["y"]]
  
  p <- ggplot(data = df.avg, aes(x = doy, y = GRO.100)) +
    stat.smooth(color="darkgreen") 
  Smooth_GRO.100 <- ggplot_build(p)$data[[1]][["y"]]
  
  p <- ggplot(data = df.avg, aes(x = doy, y = TWD.0)) +
    stat.smooth(color="darkgreen") 
  Smooth_TWD.0 <- ggplot_build(p)$data[[1]][["y"]]
  
  p <- ggplot(data = df.avg, aes(x = doy, y = TWD.25)) +
    stat.smooth(color="darkgreen") 
  Smooth_TWD.25 <- ggplot_build(p)$data[[1]][["y"]]
  
  p <- ggplot(data = df.avg, aes(x = doy, y = TWD.med)) +
    stat.smooth(color="darkgreen") 
  Smooth_TWD.med <- ggplot_build(p)$data[[1]][["y"]]
  
  p <- ggplot(data = df.avg, aes(x = doy, y = TWD.75)) +
    stat.smooth(color="lightgreen")
  Smooth_TWD.75 <- ggplot_build(p)$data[[1]][["y"]]
  
  p <- ggplot(data = df.avg, aes(x = doy, y = TWD.100)) +
    stat.smooth(color="darkgreen") 
  Smooth_TWD.100 <- ggplot_build(p)$data[[1]][["y"]]
  
  options(warn=0)
  df.avg <- df.avg %>%
    mutate(Smooth_GRO.0=Smooth_GRO.0,
           Smooth_GRO.25=Smooth_GRO.25,
           Smooth_GRO.med=Smooth_GRO.med,
           Smooth_GRO.75=Smooth_GRO.75,
           Smooth_GRO.100=Smooth_GRO.100,
           Smooth_TWD.0=Smooth_TWD.0,
           Smooth_TWD.25=Smooth_TWD.25,
           Smooth_TWD.med=Smooth_TWD.med,
           Smooth_TWD.75=Smooth_TWD.75,
           Smooth_TWD.100=Smooth_TWD.100) 
  
  # Substitute values < zero with zero
  df.avg <- df.avg %>% 
    mutate(Smooth_GRO.0 = if_else(Smooth_GRO.0 < 0, 0, Smooth_GRO.0),
           Smooth_GRO.25 = if_else(Smooth_GRO.25 < 0, 0, Smooth_GRO.25),
           Smooth_GRO.med = if_else(Smooth_GRO.med < 0, 0, Smooth_GRO.med),
           Smooth_GRO.75 = if_else(Smooth_GRO.75 < 0, 0, Smooth_GRO.75),
           Smooth_GRO.100 = if_else(Smooth_GRO.100 < 0, 0, Smooth_GRO.100),
           Smooth_TWD.0 = if_else(Smooth_TWD.0 < 0, 0, Smooth_TWD.0),
           Smooth_TWD.25 = if_else(Smooth_TWD.25 < 0, 0, Smooth_TWD.25),
           Smooth_TWD.med = if_else(Smooth_TWD.med < 0, 0, Smooth_TWD.med),
           Smooth_TWD.75 = if_else(Smooth_TWD.75 < 0, 0, Smooth_TWD.75),
           Smooth_TWD.100 = if_else(Smooth_TWD.100 < 0, 0, Smooth_TWD.100),
    )
  
  # Prepare selected year's data
  dfsel <- df %>%
    mutate(year=format.Date(ts, "%Y"),
           doy=yday(ts),
           hour=format.Date(ts, "%H"),
           Time=format.Date(ts, "%H:%M"),
           dec.date = decimal_date(ts))
  
  dfsel <- dfsel%>%
    mutate(dec.time = as.numeric(dec.date)-as.numeric(year))
  
  
  dfsel <- dfsel %>%
    filter(year == y_sel)%>%
    group_by(doy)%>%
    summarise(year=y_sel,
              gro_sel=median(gro_yr, na.rm=T),
              rel_gro_sel=gro_sel/max(df.avg$Smooth_GRO.med, na.rm=T),
              twd_sel=min(twd, na.rm=T))
  
  Linenumber <- as.numeric(length(dfsel$doy))
  z1=24*60/reso
  z2=0
  for (iii in 1:length(dfsel$doy)-2) 
  {
    z1 <- z1+24*60/reso
    z2 <- z1-2*24*60/reso
    
    if(length(df$ts) > z1){
      dfsel$twd_sel[Linenumber] <- min(df$twd[(length(df$ts)-z1):(length(df$ts)-z2)], na.rm=T)
    }
    
    Linenumber <- Linenumber-1
  }
  
  doytoday <- as.numeric(tail(dfsel$doy, n=1))
  twdtoday <- tail(dfsel$twd_sel, n=1)
  grotoday <- tail(dfsel$gro_sel, n=1)
  groreltoday <- tail(dfsel$rel_gro_sel, n=1)
  
  
  if (length(doytoday)==0) {
    warning("Skipping - no current data found")
    next()
  }  
  
  dfcom <- full_join(df.avg,dfsel)
  dfcom <- dfcom %>%
    mutate(signalGRO = case_when(
      doy < grostart.med-20 ~ "off season",
      gro_sel > Smooth_GRO.75 ~ "high GRO",
      gro_sel > Smooth_GRO.med ~ "above avg. GRO",
      gro_sel > Smooth_GRO.25 ~ "below avg. GRO",
      gro_sel >= 0 ~ "poor GRO"
      
    ),
    signalTWD = case_when(
      twd_sel > Smooth_TWD.75 ~ "high TWD",
      twd_sel > Smooth_TWD.med ~ "above avg. TWD",
      twd_sel > Smooth_TWD.25 ~ "below avg. TWD",
      twd_sel >= 0 ~ "low TWD"
    )
    )
  
  dfcom$signalGRO <- as.factor(dfcom$signalGRO)
  dfcom$signalTWD <- as.factor(dfcom$signalTWD)
  
  title <- as.character(metadata$S_Series_ID_Species[ii])
  title <- gsub("/", "_", title)
  file.name <- paste0(metadata$Site[ii], "_", metadata$Seriesname[ii])
  
  
  # Plot GRO
  dfcom$signalGRO <- as.factor(dfcom$signalGRO)
  dfcom$signalGRO <- factor(dfcom$signalGRO, levels=lvls.gro)
  
  options(warn=-1)
  p1 <- ggplot(data = dfcom, aes(doy)) +
    xlab("DOY") +
    ylab("Annual growth [µm]")+
    ggtitle(paste0(title, " (",y_sel,")"))+
    geom_ribbon(aes(x = doy, ymin = Smooth_GRO.0, ymax = Smooth_GRO.100),fill="grey70", alpha=0.2) +
    geom_ribbon(aes(x = doy, ymin = Smooth_GRO.25, ymax = Smooth_GRO.75),fill="grey70", alpha=0.5) +
    geom_line(aes(x = doy, y = Smooth_GRO.med), size=1, color="grey50")+
    geom_line(aes(x = doy, y = gro_sel), size=0.1)+
    geom_point(aes(x = doy, y = gro_sel, color=signalGRO), size=1)+
    
    geom_point(aes(x = doytoday, y = grotoday, color=signalGRO), size=4)+
    scale_color_manual(values = c("high GRO" = "darkgreen", "above avg. GRO" = "green4","below avg. GRO"="green1", "poor GRO"= "yellow", "off season"= "grey80"),
                       drop = F,
                       na.translate = F)+
    annotate("rect", ymin = -Inf, ymax = Inf, xmin = grostart.25, xmax = grostart.75, fill="yellowgreen", alpha=0.2) +
    geom_vline(xintercept = grostart.med, size=0.25)+
    geom_text(aes(x=grostart.med, label=paste0("\n", grostart.med)), y=0.8*max(Smooth_GRO.100), colour="black", angle=90, text=element_text(size=10)) +
    geom_text(aes(x=doytoday, label=paste0(doytoday)), y=grotoday+0.08*max(Smooth_GRO.100), colour="black", angle=90, text=element_text(size=10)) +
    theme_bw() + 
    theme(plot.title = element_text(size=11)) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    
    theme(legend.title=element_blank())+
    theme(legend.position=c(0.15, 0.75))+
    theme(legend.background = element_rect(fill=F))+
    theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))
  
  # Plot TWD
  
  dfcom$signalTWD <- as.factor(dfcom$signalTWD)
  dfcom$signalTWD <- factor(dfcom$signalTWD, levels=lvls.twd)
  
  
  p2 <- ggplot(data = dfcom, aes(doy)) +
    xlab("DOY") +
    ylab("TWD [µm]")+
    geom_ribbon(aes(x = doy, ymin = Smooth_TWD.0, ymax = Smooth_TWD.100),fill="grey70", alpha=0.2) +
    geom_ribbon(aes(x = doy, ymin = Smooth_TWD.25, ymax = Smooth_TWD.75),fill="grey70", alpha=0.5) +
    geom_line(aes(x = doy, y = Smooth_TWD.med), size=1, color="grey50")+
    geom_line(aes(x = doy, y = twd_sel), size=0.1)+
    geom_point(aes(x = doy, y = twd_sel, color=signalTWD), size=1)+
    geom_point(aes(x = doytoday, y = twdtoday, color=signalTWD), size=4)+
    scale_color_manual(values = c("high TWD" = "red", "above avg. TWD" = "orange","below avg. TWD"="royalblue1", "low TWD"= "darkblue"),
                       drop = F,
                       na.translate = F)+
    annotate("rect", ymin = -Inf, ymax = Inf, xmin = grostart.25, xmax = grostart.75, fill="yellowgreen", alpha=0.2) +
    geom_vline(xintercept = grostart.med, size=0.25)+
    geom_text(aes(x=grostart.med, label=paste0("\n", grostart.med)), y=0.8*max(Smooth_TWD.100), colour="black", angle=90, text=element_text(size=10)) +
    geom_text(aes(x=doytoday, label=paste0(doytoday)), y=twdtoday+0.08*max(Smooth_TWD.100), colour="black", angle=90, text=element_text(size=10)) +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(legend.title=element_blank())+
    theme(legend.position=c(0.85, 0.85))+
    theme(legend.background = element_rect(fill=F))+
    theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))
  
  # Check if values are plausible
  if (grotoday> 2*Smooth_GRO.100[doytoday]&grotoday>0.025*max(Smooth_GRO.med)|twdtoday>2*Smooth_TWD.100[doytoday])  {
    GROTWD_filename <- paste0(path, "/NA_GROTWD_", file.name)
    warning("Implausible values found in ", GROTWD_filename)
    a <- 0
  } else {
    GROTWD_filename <- paste0(path, "/GROTWD_", file.name)
    a <- 1
  }  
  
  
  # Save plot in selected format
  if (plot.format == "pdf") {
    pdf(paste0(GROTWD_filename, ".pdf"), width = 8, height = 9)
    gridExtra::grid.arrange(p1, p2, ncol=1)
    dev.off()
  } else if (plot.format == "jpg") {
    ggsave(
      filename = paste0(GROTWD_filename, ".jpg"),
      plot = cowplot::plot_grid(p1, p2, ncol=1, align="v"),
      device = "jpeg",
      scale = 1,
      width = 8,
      height = 9,
      units = "in",
      dpi = 300,
      limitsize = T
    )
  }
  options(warn=0)
  
  # Save last 5 days of data for map
  if (a==1){
    dfsignalactual[1,1]=metadata$Seriesname[ii]
    dfsignalactual[1,2]=metadata$Site_ID[ii]
    dfsignalactual[1,3]=metadata$Site_XCor[ii]
    dfsignalactual[1,4]=metadata$Site_YCor[ii]
    dfsignalactual[1,5]=metadata$Region[ii]
    dfsignalactual[1,6]=paste0(metadata$Genus[ii]," ",metadata$Species[ii])
    dfsignalactual[1,7:11]=as.numeric(dfcom$rel_gro_sel[doytoday:(doytoday-4)])
    dfsignalactual[1,12:16]=as.character(dfcom$signalGRO[doytoday:(doytoday-4)])
    dfsignalactual[1,17:21]=as.character(dfcom$signalTWD[doytoday:(doytoday-4)])
    dfsignalactual[1,22]=tail(dfsel$doy, n=1)
    dfsignal <- dfsignal %>%
      rbind(dfsignalactual)
  }
}

# Plot map of current GRO and TWD values aggregated by site ----
download.file("https://github.com/lukasmartinelli/swissdem/releases/download/v1.0/switzerland_dem.tif", "switzerland_dem.tif")
DEM <- raster::raster(paste0(path,"/switzerland_dem.tif"))

# Find the most frequent value/factor per Site_ID
dfsignalplotsite <- dfsignal %>%
  filter(last_doy >= doytoday-1) %>%
  group_by(Site_ID) %>% 
  summarise(Site_XCor = names(which(table(Site_XCor) == max(table(Site_XCor)))[1]),
            Site_YCor = names(which(table(Site_YCor) == max(table(Site_YCor)))[1]),
            Species = names(which(table(Species) == max(table(Species)))[1]),
            RelGRO = mean(RelGRO,na.rm=T),
            RelGRO_1 = mean(RelGRO_1,na.rm=T),
            RelGRO_2 = mean(RelGRO_2,na.rm=T),
            RelGRO_3 = mean(RelGRO_3,na.rm=T),
            RelGRO_4 = mean(RelGRO_4,na.rm=T),
            GRO = names(which(table(GRO) == max(table(GRO)))[1]),
            GRO_1 = names(which(table(GRO_1) == max(table(GRO_1)))[1]),
            GRO_2 = names(which(table(GRO_2) == max(table(GRO_2)))[1]),
            GRO_3 = names(which(table(GRO_3) == max(table(GRO_3)))[1]),
            GRO_4 = names(which(table(GRO_4) == max(table(GRO_4)))[1]),
            TWD = names(which(table(TWD) == max(table(TWD)))[1]),
            TWD_1 = names(which(table(TWD_1) == max(table(TWD_1)))[1]),
            TWD_2 = names(which(table(TWD_2) == max(table(TWD_2)))[1]),
            TWD_3 = names(which(table(TWD_3) == max(table(TWD_3)))[1]),
            TWD_4 = names(which(table(TWD_4) == max(table(TWD_4)))[1]),
            last_doy = names(which(table(last_doy) == max(table(last_doy)))[1]))

dfsignalplotsite$GRO <- as.factor(dfsignalplotsite$GRO)
dfsignalplotsite$GRO <- factor(dfsignalplotsite$GRO, levels=lvls.gro)
dfsignalplotsite$GRO_1 <- as.factor(dfsignalplotsite$GRO_1)
dfsignalplotsite$GRO_1 <- factor(dfsignalplotsite$GRO_1, levels=lvls.gro)
dfsignalplotsite$GRO_2 <- as.factor(dfsignalplotsite$GRO_2)
dfsignalplotsite$GRO_2 <- factor(dfsignalplotsite$GRO_2, levels=lvls.gro)
dfsignalplotsite$GRO_3 <- as.factor(dfsignalplotsite$GRO_3)
dfsignalplotsite$GRO_3 <- factor(dfsignalplotsite$GRO_3, levels=lvls.gro)
dfsignalplotsite$GRO_4 <- as.factor(dfsignalplotsite$GRO_4)
dfsignalplotsite$GRO_4 <- factor(dfsignalplotsite$GRO_4, levels=lvls.gro)

dfsignalplotsite$TWD <- as.factor(dfsignalplotsite$TWD)
dfsignalplotsite$TWD <- factor(dfsignalplotsite$TWD, levels=lvls.twd)
dfsignalplotsite$TWD_1 <- as.factor(dfsignalplotsite$TWD_1)
dfsignalplotsite$TWD_1 <- factor(dfsignalplotsite$TWD_1, levels=lvls.twd)
dfsignalplotsite$TWD_2 <- as.factor(dfsignalplotsite$TWD_2)
dfsignalplotsite$TWD_2 <- factor(dfsignalplotsite$TWD_2, levels=lvls.twd)
dfsignalplotsite$TWD_3 <- as.factor(dfsignalplotsite$TWD_3)
dfsignalplotsite$TWD_3 <- factor(dfsignalplotsite$TWD_3, levels=lvls.twd)
dfsignalplotsite$TWD_4 <- as.factor(dfsignalplotsite$TWD_4)
dfsignalplotsite$TWD_4 <- factor(dfsignalplotsite$TWD_4, levels=lvls.twd)

datetoday  <- as.Date(doytoday, origin = paste0(y_sel,"-01-01"))
datefigs <- datetoday -1

# Save plots in selected format
if (plot.format == "jpg") {
  jpeg(paste0(path, "Map_Gro_TWD_all_today_", doytoday,".jpg"), width = 5, height = 7, units = "in", res = 300) 
  {
    par(mfrow=c(2, 1), mai = c(0, 0.1, 0.1, 0), mar = c(0, 2, 0.1, 0))
    
    raster::plot(DEM, xlim=c(5.9,10.6), ylim=c(45.75,47.85), col=adjustcolor(grey.colors(100), alpha.f=0.5), zlim=c(0.1,5000),legend =F, axes=F,box = F)
    mtext( paste0("Tree growth performance (all species)"), line= -1,cex=0.8, adj=0)
    palette(c("grey80","darkgreen", "green4","#99FF99", "#FFFFCC") )
    legend('topright', legend = levels(dfsignalplotsite$GRO), col = 1:6, cex = 0.6, pch = 16, bty = "n")
    points(dfsignalplotsite$Site_XCor~dfsignalplotsite$Site_YCor, pch=21, col="darkgrey", cex=1.5)
    points(dfsignalplotsite$Site_XCor~dfsignalplotsite$Site_YCor, pch=16, col=as.factor(dfsignalplotsite$GRO), cex=1.4)
    mtext( paste0(datefigs), cex=0.6, adj=0, line = -2)
    
    raster::plot(DEM, xlim=c(5.9,10.6), ylim=c(45.75,47.85), col=adjustcolor(grey.colors(100), alpha.f=0.5), zlim=c(0.1,5000),legend =F, axes=F,box = F)
    mtext( paste0("Tree water deficit (all species)"), line=-1,cex=0.8, adj=0)
    palette(c("red","orange","royalblue1","darkblue") )
    legend('topright', legend = levels(dfsignalplotsite$TWD), col = 1:5, cex = 0.6, pch = 16, bty = "n")
    points(dfsignalplotsite$Site_XCor~dfsignalplotsite$Site_YCor, pch=21, col="darkgrey", cex=1.5)
    points(dfsignalplotsite$Site_XCor~dfsignalplotsite$Site_YCor, pch=16, col=as.factor(dfsignalplotsite$TWD), cex=1.4)
    mtext( paste0(datefigs), cex=0.6, adj=0, line = -2)
    
    dev.off() 
  }
} else if (plot.format == "pdf") {
  pdf(paste0("Map_Gro_TWD_all_today_", doytoday, ".pdf"),width=5, height=7) 
  {
    par(mfrow=c(2, 1), mai = c(0, 0.1, 0.1, 0), mar = c(0, 2, 0.1, 0))
    
    raster::plot(DEM, xlim=c(5.9,10.6), ylim=c(45.75,47.85), col=adjustcolor(grey.colors(100), alpha.f=0.5), zlim=c(0.1,5000),legend =F, axes=F,box = F)
    mtext( paste0("Tree growth performance (all species)"), line= -1,cex=0.8, adj=0)
    palette(c("grey80","darkgreen", "green4","#99FF99", "#FFFFCC") )
    legend('topright', legend = levels(dfsignalplotsite$GRO), col = 1:6, cex = 0.6, pch = 16, bty = "n")
    points(dfsignalplotsite$Site_XCor~dfsignalplotsite$Site_YCor, pch=21, col="darkgrey", cex=1.5)
    points(dfsignalplotsite$Site_XCor~dfsignalplotsite$Site_YCor, pch=16, col=as.factor(dfsignalplotsite$GRO), cex=1.4)
    mtext( paste0(datefigs), cex=0.6, adj=0, line = -2)
    
    raster::plot(DEM, xlim=c(5.9,10.6), ylim=c(45.75,47.85), col=adjustcolor(grey.colors(100), alpha.f=0.5), zlim=c(0.1,5000),legend =F, axes=F,box = F)
    mtext( paste0("Tree water deficit (all species)"), line=-1,cex=0.8, adj=0)
    palette(c("red","orange","royalblue1","darkblue") )
    legend('topright', legend = levels(dfsignalplotsite$TWD), col = 1:5, cex = 0.6, pch = 16, bty = "n")
    points(dfsignalplotsite$Site_XCor~dfsignalplotsite$Site_YCor, pch=21, col="darkgrey", cex=1.5)
    points(dfsignalplotsite$Site_XCor~dfsignalplotsite$Site_YCor, pch=16, col=as.factor(dfsignalplotsite$TWD), cex=1.4)
    mtext( paste0(datefigs), cex=0.6, adj=0, line = -2)
    
    dev.off() 
  }
}

