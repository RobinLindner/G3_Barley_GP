## ---------------------------
##
## Script name: Procrustes-PCA-plot.R
##
## Purpose of script: 
##        Plot the procrustes rotation of a PCA plot to approximate geographic sampling distribution
##
## Author: M.Sc. Robin Lindner
##
## Date Created: 2026-09-29
##
## Copyright (c) Robin Lindner, 2026
## Email: robin.lindner@uni-potsdam.de
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## ---------------------------

## load up the packages we will need:  

library(dplyr)
library(readxl)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(viridisLite)
library(vegan)

## ---------------------------
source("0_utils.R")

## Define output path:
figure_out = paste0(figure_dir,"Procrustes_PCA.png")

## ---------------------------
  
## load up our data into memory:
# Used
sample_pheno = read.csv(phenotype_nonHSR_file)

# Used
K = read.csv(GRM_path,row.names = 1)

# Used
location_data = read_excel(sampling_coord_file,sheet = 2,skip=1)

## ---------------------------    

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

## ---------------------------


## 1. - Get the sampling sites with coordinates from the Prusty et al. dataset
location_data = location_data[,1:11]

# format the site names to match genotype identifiers
location_data$Site = sub(x=location_data$Site,pattern = "-",replacement = "")

# generate a genotype to location map using the phenotypng data
location_map = data.frame(Genotype = unique(sample_pheno$Genotype)) %>%
  mutate(LocationID = substr(Genotype,1,5))

# Update the genotype to location map with the location coordinates
location_data = merge(location_map,location_data,by.x = "LocationID",by.y="Site")


## 2. - Map geographic coordinates to a specific CRS grid.

# 2. Your geographic coordinates (same row order/individuals as pc_scores!)

loc_filtered = location_data %>%
  dplyr::select(Genotype,Location,X,Y) %>%       # Extract Genotype, Location and Coordinates
  filter(Genotype %in% rownames(K)) %>%   # Filter for Genotypes in the GRM
  arrange(Genotype,rownames(K))           # Arrange Genotypes according to the rownames in K

# Add the second location identifier of the PSI data into the data.frame
sec_loc_map = sample_pheno %>%
  dplyr::select(Genotype,Location) %>%
  distinct() %>%
  rename(Loc=Location)

loc_filtered = merge(loc_filtered,sec_loc_map,by="Genotype")

# Project the X and Y coordinates onto a specific CRS grid

points_sf <- st_as_sf(loc_filtered, coords = c("X", "Y"), crs = 32636)  # projects the coordinates onto the israeli grid
points_wgs84 <- st_transform(points_sf, crs = 4326)                     # use a general wgs grid 

# Extract back to lat/long columns
coords_latlong <- st_coordinates(points_wgs84)      # Extract longitude and lattitude
loc_filtered$long <- coords_latlong[, 1]
loc_filtered$lat  <- coords_latlong[, 2]

# Store geographical coordinates in seperate data.frame
geo_coords <- data.frame(long = loc_filtered$long, lat = loc_filtered$lat)

## 3. - Run a PCA on the GRM and perform procrustes analysis
K = K[order(rownames(K)),order(rownames(K))]


# 1. Run PCA on your genetic/marker data
pca <- prcomp(K, scale. = TRUE)
pc_scores <- pca$x[, 1:2]  # PC1 and PC2

# 3. Procrustes rotation: rotates/scales/reflects PCA to best match geography
proc <- procrustes(X = geo_coords, Y = pc_scores, scale = TRUE, truemean=TRUE)
summary(proc)

pt <- protest(X = geo_coords, Y = pc_scores, permutations = 999)
summary(pt)
# extract rotated pca result
rotated_pca = proc$Yrot

# Base map (adjust to the region, e.g. Levant/Fertile Crescent for B1K)
world <- ne_countries(scale = "medium", returnclass = "sf")

# Build a data frame combining true coords + rotated PCA coords
plot_df <- data.frame(
  id        = rownames(geo_coords),
  long_true = geo_coords$long,
  lat_true  = geo_coords$lat,
  long_pca  = proc$Yrot[, 1] / (proc$scale*10) + proc$xmean[1],
  lat_pca   = proc$Yrot[, 2] / (proc$scale*10) + proc$xmean[2],
  region = loc_filtered$Loc,
  location = loc_filtered$Location
)


# Add lines to indicate PC axis
axis1_dir <- proc$scale * proc$rotation[1,1:2]
axis2_dir <- proc$scale * proc$rotation[2,1:2]

slope_axis1 <- proc$rotation[1, 2] / proc$rotation[1, 1]
slope_axis2 <- proc$rotation[2, 2] / proc$rotation[2, 1]

origin_x <- proc$translation[1]
origin_y <- proc$translation[2]

intercept_axis1 <- origin_y - (slope_axis1 * origin_x)
intercept_axis2 <- origin_y - (slope_axis2 * origin_x)

# Coloring and marker of the 51 locations:
# - Color changes based on y position (deciles)
# - marker changes within the deciles
df <- plot_df %>%
  dplyr::select(long_true,lat_true,location) %>%
  distinct() %>%
  mutate(
    y_decile = cut(
      lat_true,
      breaks = quantile(lat_true, probs = seq(0, 1, 0.1), na.rm = TRUE),
      include.lowest = TRUE,
      labels = paste0("D", 1:10)
    )
  ) %>%
  arrange(long_true) %>%
  group_by(y_decile) %>%
  mutate(shape_id = row_number()) %>%
  ungroup() %>%
  mutate(
    point_id = factor(
      paste0(y_decile, "-", shape_id),
      levels = paste0("D", rep(1:10, each = max(shape_id)), "-", rep(1:max(shape_id), 10))
    )
  )



# Build named lookup vectors: one color/shape per point_id level
decile_colors <- setNames(viridis(10), paste0("D", 1:10))
shape_pool    <- c(16, 17, 15, 3, 8, 4, 1, 2, 5, 6, 7, 9, 10, 11, 12, 13)

# Assign each point a shape and color combination
lookup <- df %>%
  distinct(point_id, y_decile, shape_id) %>%
  arrange(point_id) %>%
  mutate(
    color_val = decile_colors[as.character(y_decile)],
    shape_val = shape_pool[shape_id]
  )

# Named vectors = Dictionaries for lookup
color_map <- setNames(lookup$color_val, df$location[match(lookup$point_id,df$point_id)])
shape_map <- setNames(lookup$shape_val, df$location[match(lookup$point_id,df$point_id)])


# Plot the original sampling locations
p1 <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70") +
  coord_sf(
    xlim = range(c(plot_df$long_true, plot_df$long_pca)) + c(-0.5, 0.5),
    ylim = range(c(plot_df$lat_true, plot_df$lat_pca)) + c(-0.5, 0.5)
  ) +
  geom_point(data = plot_df, aes(long_true, lat_true, color=factor(location,levels=names(color_map)),shape=factor(location,levels=names(color_map))), size = 2) +

  scale_color_manual(values = color_map, name = "Location") +
  scale_shape_manual(values = shape_map, name = "Location") +
  scale_x_continuous(breaks = seq(33.5,37,1)) +
  labs(x="Longitude",
       y = "Latitude",
       #subtitle = "Location of origin",
       tag="B")+
  theme_minimal() +
  theme(legend.title = element_text(size = 8),
        legend.text  = element_text(size = 7),
        legend.spacing.y = unit(0.1, 'cm')) +
  guides(fill=guide_legend(ncol=5,byrow=F),
         color=guide_legend(ncol=5,byrow=F))

# Plot the rotated PCA points 
p2 <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70") +
  coord_sf(
    xlim = range(c(plot_df$long_true, plot_df$long_pca)) + c(-0.5, 0.5),
    ylim = range(c(plot_df$lat_true, plot_df$lat_pca)) + c(-0.5, 0.5)
  )+
  geom_abline(intercept = intercept_axis1,slope = slope_axis1,color="blue") +
  geom_abline(intercept = intercept_axis2,slope = slope_axis2,color="red") +
  geom_point(data = plot_df, aes(long_pca, 
                                 lat_pca, 
                                 color=factor(location,levels=names(color_map)),
                                 shape=factor(location,levels=names(color_map))), 
             size = 2) +
  annotate("text",
           x = (30 - intercept_axis1)/slope_axis1,
           y = 30,
           label = "PC1",
           color="blue",
           hjust=1) +
  annotate("text",
           x = 33.7,
           y = 33.7*slope_axis2 + intercept_axis2,
           label = "PC2",
           color="red",
           vjust=1) +
  scale_color_manual(values = color_map, name = "Location") +
  scale_shape_manual(values = shape_map, name = "Location") +
  scale_x_continuous(breaks = seq(33.5,37,1)) +
  labs(x="Longitude",
       y = "Latitude",
       #subtitle = "Rotated PCA coordinates",
       tag="C")+
  theme_minimal() 




pca = prcomp(K)
df=data.frame(Genotype=rownames(K),PC1=pca$rotation[,1],PC2=pca$rotation[,2],PC3=pca$rotation[,3])

# Extract variance explained by each PC
variance_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100
cumulative_variance <- cumsum(variance_explained)

# Number of PCs to display in the scree plot
num_pcs <- 10  # Adjust this to the desired number of PCs
num_pcs <- min(num_pcs, length(variance_explained))  # Ensure it doesn't exceed total PCs

# Create a data frame for plotting
scree_df <- data.frame(
  PC = 1:num_pcs,
  Variance_Explained = variance_explained[1:num_pcs],
  Cumulative_Variance = cumulative_variance[1:num_pcs]
)
scree_df$Included = "No"
scree_df$Included[scree_df$Variance_Explained>5] = "Yes"

p3 <-ggplot(scree_df, aes(x = PC, y = Variance_Explained,fill=factor(Included,levels=c("Yes","No")))) +
  geom_bar(stat = "identity", alpha = 1) +
  scale_x_continuous(breaks = 1:num_pcs) +
  labs(
    title = "",
    x = "Principal Component",
    y = "Variance Explained (%)"
  ) +
  labs(tag="A") +
  scale_fill_manual(values = c("Yes"="lightgreen","No"="grey"),name="Considered in GWAS")+
  #geom_vline(xintercept = 4.5,linetype = "dotted",color = "black")+
  theme_linedraw()+
  theme(plot.tag = element_text(),
        legend.position = "bottom")



mylegend<-g_legend(p1)

#p4 <- grid.arrange(arrangeGrob(p3,
#                               p1 + theme(legend.position="none"),
#                               mylegend,
#                               p2 + theme(legend.position="none"), 
#                               ncol=2,nrow=2))

p4 <- grid.arrange(arrangeGrob(p3,
                               mylegend,
                               ncol=1,nrow=2),
                   arrangeGrob(p1 + theme(legend.position="none"),
                               p2 + theme(legend.position="none"),
                               ncol=2,nrow=1),
                   ncol=2,nrow=1)

ggsave(figure_out,plot = p4,width = 12,height=6,bg="white")


##### SOMETHING IS NOT WORKING

