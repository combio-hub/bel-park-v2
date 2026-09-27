source("scripts/00. setup.R")
# install.packages("sf")
library(sf)
sf_use_s2(FALSE)
# install.packages("terra")
library(terra)

## What does a point need? To be inside a park or/and a tree
## Interesting predictors:
##    01. Tree density per park / per point (Lucas)
##    02. Distance to the nearest tree (Lucas)
##    03. Area of the nearest tree (Lucas)
##    04. Probability of connectivity (Lucas)
##    05. Integral index of connectivity (Lucas)
##    06. Minimum Cumulative Resistance (Lucas)
##    07. Patch/park density
##    08. Patch/park size
##    09. Land cover proportion
##    10. Patch/park edge density
##    11. Distance from source habitat
##    12. Distance to the nearest patch/park

br_crop = rast("dataset/spatial/brazil_coverage-col4_10m_2025.tif")
plot(br_crop)
legend = read.csv("dataset/legend_code_mapbiomas_brazil_collection_11.csv")
legend$class_name_en = gsub(" ","_",legend$class_name_en)
legend$class_name_pt_br = gsub(" ","_",legend$class_name_pt_br)
legend$class_name_pt_br = gsub("\\,","",legend$class_name_pt_br)
legend$class_name_en = gsub("\\,","",legend$class_name_en)

# Download large tree file
googledrive::drive_download("https://drive.google.com/file/d/1fi3qKzZzW-u4gj49OmLBbk2sd_A_AxGW/view?usp=sharing","dataset/spatial/trees_gsat_bel.kml")

tree_large = vect("dataset/spatial/trees_gsat_bel.kml")
tree_simple = vect("dataset/spatial/simple_trees.shp")
nrow(tree_large) == nrow(tree_simple)
tree_large$treeID = sprintf("%06d",1:nrow(tree_large))
tree_simple$treeID = sprintf("%06d",1:nrow(tree_simple))

plot(tree_simple[1:10000,])
park_bel = read_sf("dataset/spatial/praças_belem.kml")
park_ana = read_sf("dataset/spatial/praças_ananindeua.kml")
plot(park_bel$geometry)
plot(park_ana$geometry,add=T)
park_bel = park_bel[st_geometry_type(park_bel) == "POLYGON",]
park_ana = park_ana[st_geometry_type(park_ana) == "POLYGON",]



parks = rbind(park_bel[,colnames(park_bel) %in% colnames(park_ana)],park_ana[,colnames(park_ana) %in% colnames(park_bel)[(colnames(park_bel) %in% colnames(park_ana))]])
#par(mar = c(0,0,0,0))
plot(parks$geometry)

# install.packages("landscapemetrics")
library(landscapemetrics)
# install.packages("ggplot2")
library(ggplot2)

# Criar um grid de pontos para calcular as métricas da paisagem -----
tree_simple = st_as_sf(tree_simple)
tree_t = st_transform(tree_simple, crs = 31982)

tree_c = st_centroid(tree_t)
tree_c$group = "tree"
tree_c

parks_t = st_transform(parks, crs = 31982)
parks_c = st_centroid(parks_t)
parks_c$group = "park"
parks_c = st_transform(parks_c, crs = 31982)
plot(parks_c$geometry)

bufs = c(100,500,1000)

bufz = st_buffer(parks_c, dist = bufs[3])
plot(bufz$geometry)

bufz_t = st_transform(bufz,4326)
ld_crop = crop(br_crop,bufz_t[1,],mask=T)
tree_crop = st_intersection(tree_simple,bufz_t[1,])

# Match the vector CRS to the raster
tree_vect <- vect(tree_crop)

# 3 = tree present; 0 = no tree
tree_layer <- terra::rasterize(
  tree_vect,
  ld_crop[[1]],
  field = 3,
  background = 0,
  touches = TRUE
)

# Keep NA outside the cropped study area
tree_layer <- mask(tree_layer, ld_crop[[1]])
names(tree_layer) <- "tree_presence"
tree_val = values(tree_layer)
ld_val = values(ld_crop)
ld_crop[!is.na(tree_val)&(tree_val==3|ld_val==3|ld_val==6|ld_val==4|ld_val==7|ld_val==5)] = 3
# Add it as a new raster layer
plot(ld_crop)

##  Probability of connectivity
# install.packages("remotes")
#remotes::install_github("connectscape/Makurhini")
library(Makurhini)
ld_new = ld_crop
ld_new[!is.na(values(ld_new))&values(ld_new)!=3] = NA
ld_new = as.polygons(ld_new)
ld_new = st_cast(st_as_sf(ld_new), "POLYGON", do_split = TRUE)
ld_new$ID = 1:nrow(ld_new)
PC = MK_dPCIIC(ld_new,  metric = c("PC"), area_unit = "m2",overall =TRUE,onlyoverall=TRUE,distance = list(type = "centroid"),distance_thresholds = c(10,50,100,500))
PC

PC_view = MK_dPCIIC(ld_new,  metric = c("PC"), area_unit = "m2",overall =F,onlyoverall=F,distance = list(type = "centroid"),distance_thresholds = c(10,50,100,500))
PC_view

library(ggplot2)
ggplot()+
  geom_sf(data = PC_view$d500,aes(fill=log(dPC+1)))+
  scale_fill_viridis_c()

# Park metrics ----
## Patch/park size
park_data = data.frame(matrix(NA,nrow(parks),1))
colnames(park_data) = "parkID"
parks$ID = sprintf("%03d",1:nrow(parks))
park_data$parkID = parks$ID
park_data$name = parks$Name
park_data[park_data$name=="",]$name = NA
park_coords = st_coordinates(st_transform(parks_c,4326))
park_data$long = park_coords[,1]
park_data$lat = park_coords[,2]

parks_t = st_transform(parks,crs = 31982)
park_data$area = as.numeric(st_area(parks_t))/10000 # in hectares
hist(park_data$area)
hist(log10(park_data$area))

tree_simple_t = st_transform(tree_simple,31982)

park_met = pbapply::pblapply(1:nrow(parks_t),function(x){
  park_inter = st_intersection(parks_t[x,],tree_c)
  sel_trees = tree_large[as.numeric(park_inter$treeID),]
  
  # plot(st_transform(parks_t[x,]$geometry,4326))
  # plot(sel_trees,add=TRUE,col="red")
  
  park_land = crop(br_crop,parks[x,],mask=T)
  park_land = project(park_land,"EPSG:31982",method="near")
  park_land = crop(park_land,parks_t[x,],mask=T)
  
  park_val = values(park_land)
  tried = try(values(park_land)[!is.na(park_val)&(park_val==3|park_val==6|park_val==4|park_val==7|park_val==5|park_val==11)],silent = T)
  if(is(tried)[1] != "try-error"){
    values(park_land)[!is.na(park_val)&(park_val==3|park_val==6|park_val==4|park_val==7|park_val==5|park_val==11)] = 3
    park_val = values(park_land)
    
    if(length(sel_trees)>0){
      tree_layer <- terra::rasterize(
        project(sel_trees,"EPSG:31982"),
        park_land[[1]],
        field = 3,
        background = 0,
        touches = TRUE
      )
      tree_layer = crop(tree_layer,parks_t[x,],mask=T)
      tree_values = values(tree_layer)
      values(park_land)[!is.na(park_val)&(tree_values==3&park_val!=3)] = 3
      #plot(park_land)
    }
    
  } 
  lsms = calculate_lsm(park_land,what = c("lsm_c_pd","lsm_c_pland","lsm_c_ed","lsm_c_clumpy","lsm_c_area_mn"))
  rbind(lsms[lsms$metric!="pland",][lsms[lsms$metric!="pland",]$class == 3,],
        lsms[lsms$metric=="pland",]) -> lsms
  lsms$id = parks$ID[x]
  
  return(data.frame(class = lsms$class,
                    id = parks$ID[x],
                    metric = lsms$metric,
                    value = lsms$value))
})
park_met = do.call("rbind",park_met)
cast_data = reshape2::dcast(data = park_met,formula =  id ~ metric + class,measure.var = "value",fun.aggregate = mean)
park_data = cbind(park_data,cast_data[,2:ncol(cast_data)])
park_data[is.na(park_data)] = 0

distz = as.matrix(as.dist(st_distance(parks_t)))
min_distances = unlist(lapply(1:nrow(distz),function(x){min(distz[x,][distz[x,]!=0])}))
park_data$nearest = min_distances
plot(br_crop)

br_val = values(br_crop)
values(br_crop)[!is.na(br_val)&(br_val==3|br_val==6|br_val==4|br_val==7|br_val==5|br_val==11)] = 3
br_val = values(br_crop)
values(br_crop)[!is.na(br_val)&(br_val!=3)] = NA
plot(br_crop)

vec_br = st_as_sf(stars::st_as_stars(br_crop),merge = TRUE)
vec_br_t = st_transform(vec_br,31982)
area_vec = st_area(vec_br_t)
hec_area = as.numeric(area_vec)/10000
sum(hec_area > 10000)
hec_50 = vec_br_t[hec_area > 50,]
hec_1000 = vec_br_t[hec_area > 1000,]
hec_100 = vec_br_t[hec_area > 100,]

dist_100 = st_distance(parks_t,hec_100)
dist_100 = apply(dist_100,1,min)

dist_1000 = st_distance(parks_t,hec_1000)
dist_1000 = apply(dist_1000,1,min)

dist_10000 = st_distance(parks_t,hec_10000)
dist_10000 = apply(dist_10000,1,min)

park_data$mainland_100ha = as.numeric(dist_100)
park_data$mainland_1000ha = as.numeric(dist_1000)
park_data$mainland_10000ha = as.numeric(dist_10000)

br_crop = rast("dataset/spatial/brazil_coverage-col4_10m_2025.tif")
park_scale = list()
for(y in bufs){
  print(y)
  park_met2 = pbapply::pblapply(1:nrow(parks_t),function(x){
    buf_temp = st_buffer(parks_t[x,],y)
    park_inter = st_intersection(buf_temp,tree_c)
    sel_trees = tree_large[as.numeric(park_inter$treeID),]
    
    # plot(st_transform(buf_temp$geometry,4326))
    # plot(sel_trees,add=TRUE,col="red")
    
    park_land = crop(br_crop,st_transform(buf_temp,4326),mask=T)
    park_land = project(park_land,"EPSG:31982",method="near")
    park_land = crop(park_land,buf_temp,mask=T)
    
    park_val = values(park_land)
    tried = try(values(park_land)[!is.na(park_val)&(park_val==3|park_val==6|park_val==4|park_val==7|park_val==5|park_val==11)],silent = T)
    if(is(tried)[1] != "try-error"){
      values(park_land)[!is.na(park_val)&(park_val==3|park_val==6|park_val==4|park_val==7|park_val==5|park_val==11)] = 3
      park_val = values(park_land)
      
      if(length(sel_trees)>0){
        tree_layer <- terra::rasterize(
          project(sel_trees,"EPSG:31982"),
          park_land[[1]],
          field = 3,
          background = 0,
          touches = FALSE
        )
        tree_layer = crop(tree_layer,buf_temp,mask=T)
        tree_values = values(tree_layer)
        values(park_land)[!is.na(park_val)&(tree_values==3&park_val!=3)] = 3
        #plot(park_land)
      }
      
    } 
    lsms = calculate_lsm(park_land,what = c("lsm_c_pd","lsm_c_pland","lsm_c_ed","lsm_c_clumpy","lsm_c_area_mn"))
    rbind(lsms[lsms$metric!="pland",][lsms[lsms$metric!="pland",]$class == 3,],
          lsms[lsms$metric=="pland",]) -> lsms
    lsms$id = parks$ID[x]
    
    return(data.frame(class = lsms$class,
                      id = parks$ID[x],
                      metric = lsms$metric,
                      value = lsms$value,
                      scale = y))
  })
  park_met2 = do.call("rbind",park_met2)
  park_scale[[length(park_scale)+1]] = park_met2
}
park_scales = do.call("rbind",park_scale)

cast_data2 = reshape2::dcast(data = park_scales,formula =  id ~ metric + class + scale,measure.var = "value",fun.aggregate = mean)
cast_data2[,2:ncol(cast_data2)][is.na(cast_data2[,2:ncol(cast_data2)])] = 0
all_data_park = cbind(park_data,cast_data2[,2:ncol(cast_data2)])
all_data_park =all_data_park[,!grepl("group",colnames(all_data_park))]

## Gradiente de paisagem para as praças -----
park_scale = all_data_park
park_scale$area = log10(park_scale$area)
pca = prcomp(scale(park_scale[,5:ncol(park_scale)]))
summ <- summary(pca)
sum(summ$importance[2,1:4])

pca_data = data.frame(pca$x)
pca_data$id = park_data$parkID
pca_data$lat = park_data$lat
pca_data$long = park_data$long
var_scores = data.frame(pca$rotation)

ggplot()+
  geom_point(data = pca_data, aes(x = PC1, y = PC2))+
  geom_segment(data = var_scores, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10), color = "black", linewidth = 0.5, arrow = arrow(length = unit(0.06, "inches"))) +
  ggrepel::geom_text_repel(data = var_scores, aes(x = PC1 * 10, y = PC2 * 10, label = rownames(var_scores)), size = 3.5, fontface = "italic", color = "black") 

library(vegan)  
all_data_park = all_data_park[-362,]
eu_dist = vegdist(scale(all_data_park[,5:ncol(all_data_park)]),method="euclidean")
clust = hclust(eu_dist)
plot(clust)
samp_groups = cutree(clust,k=30)
all_data_park$group = samp_groups
hist(all_data_park$group)
all_data_park[all_data_park$group==30,]

sel_parks = list()
sel_parks[[1]] = st_as_sf(all_data_park[356,],coords=c("long","lat"),crs=4326)

for(x in 1:29){
  print(x)
  temp_samp = all_data_park[all_data_park$group==x,][sample(1:nrow(all_data_park[all_data_park$group==x,]),1),]
  temp_samp = st_as_sf(temp_samp,coords=c("long","lat"),crs=4326)
  d = 0
  while(min(as.numeric(st_distance(temp_samp,do.call("rbind",sel_parks))))<2000&d<999){
    temp_samp = all_data_park[all_data_park$group==x,][sample(1:nrow(all_data_park[all_data_park$group==x,]),1),]
    temp_samp = st_as_sf(temp_samp,coords=c("long","lat"),crs=4326)
    d = d +1
    #print(d)
    if(d==999){
      print("999")
    }
  }
  sel_parks[[length(sel_parks)+1]] = temp_samp
}
sel_parks = do.call("rbind",sel_parks)
plot(sel_parks$geometry)

hist(sel_parks$pland_3)
hist(sel_parks$pland_3_100)
hist(sel_parks$pland_3_1000)
hist(log10(sel_parks$area))

plot(parks[as.numeric(sel_parks$parkID),]$geometry)

## Variáveis para as árvores ----
bufs = c(100,500,1000)

bufz = st_buffer(tree_c, dist = bufs[length(bufs)])
bufz_t = st_transform(bufz,4326)
bufz_uni = st_union(bufz_t)
ld_crop_all = crop(br_crop,bufz_uni,mask=T)
ld_crop_all = terra::mask(ld_crop_all,vect(bufz_uni))
plot(ld_crop_all)

tree_vect <- terra::vect(tree_simple)
tree_layer <- terra::rasterize(
  tree_vect,
  ld_crop_all[[1]],
  field = 3,
  background = 0,
  touches = TRUE
)
# Keep NA outside the cropped study area
tree_layer <- terra::mask(tree_layer, ld_crop_all[[1]])
names(tree_layer) <- "tree_presence"
tree_val = values(tree_layer)
ld_val = values(ld_crop_all)
ld_crop_all[!is.na(tree_val)&(tree_val==3|ld_val==3|ld_val==6|ld_val==4|ld_val==7|ld_val==5|ld_val==11)] = 3
ld_crop_all = terra::project(ld_crop_all,"EPSG:31982",method="near")
plot(ld_crop_all)

#writeRaster(ld_crop_all,"dataset/spatial/tree_bel_classified.tif")
library(stars)
library(Makurhini)
library(landscapemetrics)
ld_crop_all = read_stars("dataset/spatial/tree_bel_classified.tif")

#install.packages("osmdata")
library(osmdata)
available_features()

PA = read_sf("dataset/spatial/PA_Municipios_2025/PA_Municipios_2025.shp")
PA = PA[PA$NM_MUN=="Belém"|PA$NM_MUN=="Ananindeua",]

PA_highway <- st_bbox(PA) %>%
  opq() %>%
  add_osm_feature(key = "highway") %>%
  osmdata_sf()
PA_highway = PA_highway$osm_lines

plot(PA_highway$geometry)

PA_highway_t = st_transform(PA_highway,31982)
buf_high = st_buffer(PA_highway_t,5)
parks

inter_tree = st_intersects(tree_c,buf_high)
inter_park = st_intersects(tree_c,parks_t)

listed_tree = lapply(inter_tree,function(x){length(x)>0})
listed_park = lapply(inter_park,function(x){length(x)>0})
street_trees = tree_c[unlist(listed_tree)|unlist(listed_park),]
street_trees

#random_trees = street_trees[sample(1:nrow(street_trees),9999),]
#write_sf(random_trees,"dataset/spatial/random_trees.shp")
random_trees = read_sf("dataset/spatial/random_trees.shp")
random_trees = st_transform(random_trees,31982)
plot(random_trees$geometry)

#y=100
all_scales = {}
for(y in bufs){
  bufz = st_buffer(random_trees, dist = y)
  cl = makeCluster(8)
  clusterEvalQ(cl,library(sf))
  clusterEvalQ(cl,library(stars))
  clusterEvalQ(cl,library(landscapemetrics))
  clusterExport(cl,list("ld_crop_all","bufz"))
  tree_mets = pbapply::pblapply(cl=cl,1:nrow(random_trees),function(x){
    ld_new = st_crop(ld_crop_all,bufz[x,],mask=T)
    ld_pol = ld_new
    
    ld_pol[[1]][!is.na(ld_pol[[1]])&ld_pol[[1]] != 3] = NA
    ld_pol = st_as_sf(ld_pol,merge = TRUE)
    ld_pol$ID = 1:nrow(ld_pol)
    if(nrow(ld_pol)<2){
      PC = NA
      lsms_tre = calculate_lsm(ld_new,what = c("lsm_c_pd","lsm_c_pland","lsm_c_ed","lsm_c_clumpy","lsm_c_area_mn"))
      rbind(lsms_tre[lsms_tre$metric!="pland",][lsms_tre[lsms_tre$metric!="pland",]$class == 3,],
            lsms_tre[lsms_tre$metric=="pland",]) -> lsms_tre
      lsms_tre$id = bufz[x,]$treeID
      return(data.frame(class = lsms_tre$class,
                        id = bufz[x,]$treeID,
                        metric = lsms_tre$metric,
                        value = lsms_tre$value))
      
    } else {
      PC = suppressMessages(Makurhini::MK_dPCIIC(ld_pol,  metric = c("PC"), area_unit = "m2",overall =TRUE,onlyoverall=TRUE,distance = list(type = "centroid"),distance_thresholds = c(10,50,100)))
      
      PC = reshape2::melt(PC)
      PC$id = x
      PC$Index = paste0(PC$Index,"_",PC$L1)
      PC$Index = gsub("___|__|____|_____|______","_",gsub("\\(|\\)|%","_",PC$Index))
      lsms_tre = calculate_lsm(ld_new,what = c("lsm_c_pd","lsm_c_pland","lsm_c_ed","lsm_c_clumpy","lsm_c_area_mn"))
      rbind(lsms_tre[lsms_tre$metric!="pland",][lsms_tre[lsms_tre$metric!="pland",]$class == 3,],
            lsms_tre[lsms_tre$metric=="pland",]) -> lsms_tre
      lsms_tre$id = bufz[x,]$treeID
      return(data.frame(class = c(lsms_tre$class,rep(3,nrow(PC))),
                        id = bufz[x,]$treeID,
                        metric = c(lsms_tre$metric,PC$Index),
                        value = c(lsms_tre$value,PC$value)))
    }
    
  })
  all_mets = do.call("rbind",tree_mets)
  all_mets$scale = y
  all_scales[[length(all_scales)+1]] = all_mets
  stopCluster(cl)
  gc()
}

all_scales_list = all_scales
all_scales = do.call("rbind",all_scales)
write.csv(all_scales,"dataset/all_scales.csv")
all_scales = read.csv("dataset/all_scales.csv",row.names = 1)

all_scales$class = as.factor(all_scales$class)
all_scales$scale = as.factor(all_scales$scale)
all_scales$id = as.factor(all_scales$id)
cast_scales = reshape2::recast(all_scales,id ~ metric + class + scale,measure.var = "value")
cast_scales[is.na(cast_scales)] = 0

distz_tree = st_distance(random_trees,tree_c)

nearest_trees = pbapply::pblapply(1:nrow(distz_tree),function(x){
  return(as.numeric(min(distz_tree[x,][as.numeric(distz_tree[x,])!=0])))
})

nearest_trees = unlist(nearest_trees)
hist(nearest_trees)

dist50 = st_distance(random_trees,hec_50)
dist100 = st_distance(random_trees,hec_100)
dist1000 = st_distance(random_trees,hec_1000)

min_dist50 = apply(dist50,1,min)
min_dist100 = apply(dist100,1,min)
min_dist1000 = apply(dist1000,1,min)

cast_scales$nearest = nearest_trees
cast_scales$mainland50 = min_dist50
cast_scales$mainland100 = min_dist100
cast_scales$mainland1000 = min_dist1000

library(vegan)
cast_scales = cast_scales[,!grepl("group",colnames(cast_scales))]
eu_dist = vegdist(scale(cast_scales[,2:ncol(cast_scales)]),method="euclidean",na.rm=T)
clust = hclust(eu_dist)
plot(clust)
samp_groups = cutree(clust,k=40)
cast_scales$group = samp_groups
hist(cast_scales$group)
cast_scales[cast_scales$group==30,]
cast_scales$id = sprintf("%06d",as.numeric(cast_scales$id))

sel_trees = list()
#sel_trees[[1]] = st_as_sf(tree_c,coords=c("long","lat"),crs=4326)

for(x in 1:40){
  #print(x)
  d = 0
  temp_samp = cast_scales[cast_scales$group==x,][sample(1:nrow(cast_scales[cast_scales$group==x,]),1),]$id
  temp_samp = tree_c[tree_c$treeID == temp_samp,]
  if(length(sel_trees)!=0){
    while(min(as.numeric(st_distance(temp_samp,do.call("rbind",sel_trees))))<2000&d<99){
      temp_samp = cast_scales[cast_scales$group==x,][sample(1:nrow(cast_scales[cast_scales$group==x,]),1),]$id
      temp_samp = tree_c[tree_c$treeID == temp_samp,]
      d = d +1
      #print(d)
      if(d==999){
        print("999")
      }
      #print(d)
    }
  }
  #print(nrow(temp_samp))
  sel_trees[[length(sel_trees)+1]] = temp_samp
  print(length(sel_trees))
}
sel_trees = do.call("rbind",sel_trees)
plot(sel_trees$geometry)

sel_trees = cbind(sel_trees,cast_scales[match(sel_trees$treeID,cast_scales$id),])

length(sel_trees$treeID)

hist(sel_trees$pland_3_100)
hist(sel_trees$pland_3_1000)
hist(sel_trees$PCintra_d50_3_100)
hist(cast_scales$PCintra_d50_3_100)
hist(log10(sel_trees$area))

plot(sel_trees[as.numeric(sel_trees$parkID),]$geometry)

sel_trees_t = st_transform(sel_trees,4326)
sel_trees_df = cbind(st_coordinates(sel_trees_t),as.data.frame(sel_trees_t)[,1:(ncol(sel_trees_t)-2)])
write.csv(sel_trees_df,"dataset/sel_trees_df.csv",row.names = F)
write_sf(sel_trees_t,"dataset/spatial/sel_trees.kml",row.names = F)
