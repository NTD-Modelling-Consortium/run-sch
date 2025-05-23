
library("AMISforInfectiousDiseases")
library("vioplot")
library("useful")
library("dplyr")
library("viridis")
library("patchwork")
library("sf")
library("readxl")
library("tidyr")
library("ggplot2")
library("gridExtra")

plot_trajectories <- T

# uncomment for relevant species
#species <- "ascaris"
#failed_ids = c(1007,1127,1128,1333,252,269)
#failed_ids_sigma0.025 = c(1007,1127,1128,1333,252,269)

# species <- "hookworm"
# failed_ids = c(269,512)
# failed_ids_sigma0.025 = c(269,512)

species <- "trichuris"
failed_ids = c(1106,1168,1306,326,47,83)
failed_ids_sigma0.025 = c(1106,1168,1306,326,47,83)


# loading 'iu_task_lookup' (batches-IUs look up table for the fitting)
if(species == "trichuris"){
  load("../Maps-STH/iu_task_lookup_trichuris.rds")
}else if (species %in% c("ascaris","hookworm")){
  load("../Maps-STH/iu_task_lookup_sth.rds")
}

num_batches <- max(iu_task_lookup$TaskID)

cat(paste0("species: ",species, " \n"))
cat(paste0("num_batches: ",num_batches, " \n"))

ids_sample_pars = setdiff(1:num_batches,failed_ids_sigma0.025)


load(paste0('../Maps-STH/',species,'_maps.rds')) # load species_map_allyears

if(species=="ascaris"){
  species_map_allyears <- ascaris_map_allyears
}else if(species=="hookworm"){
  species_map_allyears <- hookworm_map_allyears
}else if(species=="trichuris"){
  species_map_allyears <- trichuris_map_allyears
}

length(species_map_allyears) # T: number of time points
str(species_map_allyears[[1]]) # L x #samples + 2
str(species_map_allyears[[2]])
str(species_map_allyears[[3]])
M_l <- 200 # just to simulate from normal to run vioplots

load(paste0("../trajectories/trajectories_",1,"_",species,".Rdata")) # load 'trajectories' of batch 1
total_num_years <- ncol(trajectories)

map_years <- c(2000, 2013, 2018)  # years in the map samples amis fitted to
all_years <- 2000:2018

# Look up table that links each IU to its index in the corresponding batch
if(species %in% c("ascaris","hookworm")){
  table_iu_idx <- read.csv("../Maps-STH/table_iu_idx_STH.csv") 
} else {
  table_iu_idx <- read.csv("../Maps-STH/table_iu_idx_trichuris.csv") 
}

dim(table_iu_idx)
head(table_iu_idx)

countries <- sort(unique(table_iu_idx$country))


# for(country in countries){

  cat(paste0(" \n"))
  cat(paste0("=============================================== \n"))
  cat(paste0(" \n"))
  cat(paste0("Producing results for all countries \n"))

    
# wh <- which(table_iu_idx$country==country) # specific country
# Selecting ALL countries
ius_for_failed_batches = table_iu_idx %>% filter(TaskID %in% failed_ids_sigma0.025)
ius_for_failed_batches_index = which(table_iu_idx$IU_CODE %in% ius_for_failed_batches$IU_CODE)
wh <- setdiff(1:nrow(table_iu_idx),ius_for_failed_batches_index)
table_country <- table_iu_idx[wh,]
head(table_country)
L <- nrow(table_country)  # number of locations in the country

# Folder where figures are saved
pathForPlots <- paste0("outputs_all_countries_",species,"/")
if (!dir.exists(pathForPlots)) {dir.create(pathForPlots)}


n_samples <- 1000 # number of samples per iteration
num_sub_samples <- 200

panel_ncols <- 4
panel_nrows <- ceiling(L/panel_ncols)   ####


# decide order in which IUs are plotted in trajectories plots
ess_all_iu <- data.frame()
ius_requiring_sigma0.025 = data.frame()
for (id in ids_sample_pars){

  #### Load AMIS output (this is just to get IU names and initial ESS)
  if(!id %in% failed_ids){
    load(paste0("../AMIS_output/",species,"_amis_output",id,".Rdata")) # loads amis_output
  } else {
    load(paste0("../AMIS_output/",species,"_amis_output",id,"_sigma0.025.Rdata")) # loads amis_output
  }

  ess = amis_output$ess
  iu_names <- rownames(amis_output$prevalence_map[[1]]$data)

  if(!id %in% failed_ids){
    iu_names_lt200 = iu_names[ess<200]
    iu_names_ge200 = iu_names[!ess<200]
  } else {
    iu_names_lt200 = iu_names # if failed when sigma=0.0025 then use sigma=0.025 for all IUs
  }

  ess_iu = data.frame(IU_CODE = iu_names, ess = rep(NA,length(iu_names)))

  ius_requiring_sigma0.025 = c(ius_requiring_sigma0.025,iu_names_lt200)

  # get ESS for ius in iu_names_ge200
  if (file.exists(paste0("../../AMIS_output/",species,"_amis_output",id,".Rdata"))) {
    load(paste0("../AMIS_output/",species,"_amis_output",id,".Rdata")) # loads amis_output
    ess_iu$ess[which(iu_names %in% iu_names_ge200)] = amis_output$ess[which(iu_names %in% iu_names_ge200)]
  }
  
  # get ESS for ius in iu_names_lt200
  if (file.exists(paste0("../AMIS_output/",species,"_amis_output",id,"_sigma0.025.Rdata"))) {
    load(paste0("../AMIS_output/",species,"_amis_output",id,"_sigma0.025.Rdata")) # loads amis_output
    ess_iu$ess[which(iu_names %in% iu_names_lt200)] = amis_output$ess[which(iu_names %in% iu_names_lt200)]
  }

  ess_all_iu <- rbind(ess_all_iu,ess_iu)
}
ess_all_iu$ess = as.numeric(ess_all_iu$ess)
ess_all_iu$IU_CODE = as.numeric(ess_all_iu$IU_CODE)
# reorder to align with table_country 
ess_all_iu = ess_all_iu[sapply(1:nrow(table_country), function(j) which(ess_all_iu$IU_CODE == table_country$IU_CODE[j])),]
ixd_ord_traj_plots <- order(ess_all_iu[,"ess"])


if (species %in% c("ascaris","hookworm","trichuris")){
  source(paste0("../run-sch/",species,"_prior.R"))
  prior = Prior
} 
prior_samples = prior$rprior(100000)

priorkmed <- quantile(prior_samples[,2],probs=0.5)
priorklwr <- quantile(prior_samples[,2],probs=0.025)
priorkupr <- quantile(prior_samples[,2],probs=0.975)

priorR0med <- quantile(prior_samples[,1],probs=0.5)
priorR0lwr <- quantile(prior_samples[,1],probs=0.025)
priorR0upr <- quantile(prior_samples[,1],probs=0.975)


species_plot=species


# graphical parameters
par(mar=c(4,4,2,1)+.5)
cex <- 3
lwd <- 3
cex.lab <- 1.5
cex.axis <- 1.5
cex.main <- 1.5
breaks <- 500
main <- NULL
param_names <- c('R0', 'k')
ylab <- NULL
xlab <- NULL

load(paste0("InputPars_MTP_",species,"/InputPars_MTP_allIUs.rds"))
colnames(sampled_params_all) = c("IU_ID","location","seed","R0","k","prev_t1","prev_t2","prev_t3")

dat <- sampled_params_all %>% 
  mutate_if(is.character, as.numeric) %>%
  left_join(ess_all_iu %>% mutate(IU_ID=as.numeric(IU_CODE))) %>%
  filter(ess >= 200)

mnR0 <- aggregate(R0~IU_ID, data=dat, FUN="median")
lwrR0 <- aggregate(R0~IU_ID, data=dat, FUN=quantile, probs=0.05)
uprR0 <- aggregate(R0~IU_ID, data=dat, FUN=quantile, probs=0.95)

R0 <- (cbind(mnR0, lwrR0[,2], uprR0[,2]))
colnames(R0) <- c("IU_ID", "R0", "R0lwr", "R0upr")

mnK <- aggregate(k~IU_ID, data=dat, FUN="median")
lwrK <- aggregate(k~IU_ID, data=dat, FUN=quantile, probs=0.05)
uprK <- aggregate(k~IU_ID, data=dat, FUN=quantile, probs=0.95)

k <- (cbind(mnK, lwrK[,2], uprK[,2]))
colnames(k) <- c("IU_ID", "k", "klwr", "kupr")

mnprev1 <- aggregate(prev_t1~IU_ID, data=dat, FUN="median")
lwrprev1 <- aggregate(prev_t1~IU_ID, data=dat, FUN=quantile, probs=0.05)
uprprev1 <- aggregate(prev_t1~IU_ID, data=dat, FUN=quantile, probs=0.95)

prev1 <- (cbind(mnprev1, lwrprev1[,2], uprprev1[,2]))
colnames(prev1) <- c("IU_ID", "prev1", "prev1lwr", "prev1upr")

mnprev2 <- aggregate(prev_t2~IU_ID, data=dat, FUN="median")
lwrprev2 <- aggregate(prev_t2~IU_ID, data=dat, FUN=quantile, probs=0.05)
uprprev2 <- aggregate(prev_t2~IU_ID, data=dat, FUN=quantile, probs=0.95)

prev2 <- (cbind(mnprev2, lwrprev2[,2], uprprev2[,2]))
colnames(prev2) <- c("IU_ID", "prev2", "prev2lwr", "prev2upr")

mnprev3 <- aggregate(prev_t3~IU_ID, data=dat, FUN="median")
lwrprev3 <- aggregate(prev_t3~IU_ID, data=dat, FUN=quantile, probs=0.05)
uprprev3 <- aggregate(prev_t3~IU_ID, data=dat, FUN=quantile, probs=0.95)

prev3 <- (cbind(mnprev3, lwrprev3[,2], uprprev3[,2]))
colnames(prev3) <- c("IU_ID", "prev3", "prev3lwr", "prev3upr")

df <- cbind(R0, k[,2:4], prev1[,2:4],prev2[,2:4],prev3[,2:4]) #%>% 
df <- df[order(df$prev1),]
df$IUN <- seq(1,nrow(df))


pR0 <- ggplot(data=df, aes(y=IUN, x=R0))+
  geom_segment(aes(x=R0lwr, xend=R0upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
  geom_point() +
  xlab("R0") +
  theme(axis.text.y = element_blank()) +
  theme_bw() +
  scale_y_continuous(name="IU")

#labels = scales::trans_format("log10", scales::math_format(10^.x)))

pk <- ggplot(data=df, aes(y=IUN, x=k))+
  geom_segment(aes(x=klwr, xend=kupr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
  geom_point() +
  xlab("k") +
  theme(axis.text.y = element_blank()) +
  theme_bw() +
  scale_y_continuous(name="IU")


pprev1 <- ggplot(data=df, aes(y=IUN, x=prev1))+
  geom_segment(aes(x=prev1lwr, xend=prev1upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
  geom_point() +
  xlab("Prevalence 2002") +
  theme(axis.text.y = element_blank()) +
  theme_bw() +
  scale_y_continuous(name="IU") +
  xlim(0,1)

pprev2 <- ggplot(data=df, aes(y=IUN, x=prev2))+
  geom_segment(aes(x=prev2lwr, xend=prev2upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
  geom_point() +
  xlab("Prevalence 2013") +
  theme(axis.text.y = element_blank()) +
  theme_bw() +
  scale_y_continuous(name="IU") +
  xlim(0,1)

pprev3 <- ggplot(data=df, aes(y=IUN, x=prev3))+
  geom_segment(aes(x=prev3lwr, xend=prev3upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
  geom_point() +
  xlab("Prevalence 2022") +
  theme(axis.text.y = element_blank()) +
  theme_bw() +
  scale_y_continuous(name="IU") +
  xlim(0,1)


posterior_plots = grid.arrange(
  pR0+ 
    annotate(geom="segment", x=priorR0lwr, xend=priorR0upr, y=mean(df$IUN), yend=mean(df$IUN), col="red") +
    annotate(geom="point", x=priorR0med, y=mean(df$IUN), col="red"),
  
  pk+   annotate(geom="segment", x=priorklwr, xend=priorkupr, y=mean(df$IUN), yend=mean(df$IUN), col="red") +
    annotate(geom="point", x=priorkmed, y=mean(df$IUN), col="red"),
  ggplot(),
  pprev1, pprev2, pprev3,
  
  nrow=2)
ggsave(posterior_plots,file=paste0(pathForPlots,"posteriors_multipletimepts_",species_plot,".png"))


# Map median posterior prevalence, R0 and k
shape <- read_sf(dsn = "ESPEN_IU_2021/", layer = "ESPEN_IU_2021")

#Prevalence
prev_t1_medians = mnprev1 %>%
  mutate_if(is.character, as.numeric)
colnames(prev_t1_medians) = c("IU_ID","prev_t1")
espen_prev_t1 = shape %>%
  dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
  left_join(prev_t1_medians, by="IU_ID")
map_prev_t1 = st_as_sf(espen_prev_t1)

prev_t2_medians = mnprev2 %>%
  mutate_if(is.character, as.numeric)
colnames(prev_t2_medians) = c("IU_ID","prev_t2")
espen_prev_t2 = shape %>%
  dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
  left_join(prev_t2_medians, by="IU_ID")
map_prev_t2 = st_as_sf(espen_prev_t2)

prev_t3_medians = mnprev3 %>%
  mutate_if(is.character, as.numeric)
colnames(prev_t3_medians) = c("IU_ID","prev_t3")
espen_prev_t3 = shape %>%
  dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
  left_join(prev_t3_medians, by="IU_ID")
map_prev_t3 = st_as_sf(espen_prev_t3)

#R0
R0_medians = mnR0 %>%
  mutate_if(is.character, as.numeric)
colnames(R0_medians) = c("IU_ID","R0")
espen_R0 = shape %>%
  dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
  left_join(R0_medians, by="IU_ID")
map_R0 = st_as_sf(espen_R0)

#k
K_medians = mnK %>%
  mutate_if(is.character, as.numeric)
colnames(K_medians) = c("IU_ID","k")
espen_K = shape %>%
  dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
  left_join(K_medians, by="IU_ID")
map_K = st_as_sf(espen_K)

# Prevalences in Evandro's style
th <- theme(plot.title = element_text(size=12, hjust=0.5))
lwd_borders <- 0
colour <- "D"

p_prev1 <- ggplot(map_prev_t1) +
  geom_sf(aes(fill=prev_t1), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
  theme_bw() + theme(legend.position="right") + 
  guides(fill=guide_colorbar(title="Median estimated prevalence")) +
  ggtitle(map_years[1]) + th
#p_prev1

p_prev2 <- ggplot(map_prev_t2) +
  geom_sf(aes(fill=prev_t2), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
  theme_bw() + theme(legend.position="right") + 
  guides(fill=guide_colorbar(title="Median estimated prevalence")) +
  ggtitle(map_years[2]) + th
#p_prev2

p_prev3 <- ggplot(map_prev_t3) +
  geom_sf(aes(fill=prev_t3), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
  theme_bw() + theme(legend.position="right") + 
  guides(fill=guide_colorbar(title="Median estimated prevalence")) +
  ggtitle(map_years[3]) + th
#p_prev3

plot_prevs <- p_prev1 + p_prev2 + p_prev3 & theme(legend.position = "right")
plot_prevs <- plot_prevs + plot_layout(guides = "collect")

fileName <- paste0(pathForPlots,"map_median_prevs_all_countries_",species_plot,".png")
ggsave(fileName, plot_prevs, width = 120*3, height = 120, units = "mm")

# Parameters estimates in Evandro's style
p_R0 <- ggplot(map_R0) +
  geom_sf(aes(fill=R0), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") + 
  guides(fill=guide_colorbar(title="Median R0")) +
  ggtitle("R0") + th
#p_R0

p_k <- ggplot(map_K) +
  geom_sf(aes(fill=k), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80") +
  theme_bw() + theme(legend.position="right") + 
  guides(fill=guide_colorbar(title="Median k")) +
  ggtitle("k") + th
#p_k

plot_pars <- p_R0 + p_k & theme(legend.position = "right")

fileName <- paste0(pathForPlots,"parameters_all_countries_",species_plot,".png")
ggsave(fileName, plot_pars, width = 120*2, height = 120, units = "mm")


# plot geostatistical map prevalences for comparison
logistic <- function(x){return(1/(1+exp(-x)))}
mapped_prev_t1_means = data.frame(species_map_allyears[[1]]$data[,c("IU_2021","mu")])%>%
  mutate(IU_ID = as.numeric(IU_2021),
         mapped_prev_t1 = logistic(mu)) %>%
  select(IU_ID,mapped_prev_t1)
espen_prev_t1 = shape %>%
  dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
  left_join(mapped_prev_t1_means, by="IU_ID")
geomap_prev_t1 = st_as_sf(espen_prev_t1)

mapped_prev_t2_means = data.frame(species_map_allyears[[2]]$data[,c("IU_2021","mu")])%>%
  mutate(IU_ID = as.numeric(IU_2021),
         mapped_prev_t1 = logistic(mu)) %>%
  select(IU_ID,mapped_prev_t1)
colnames(mapped_prev_t2_means) = c("IU_ID","mapped_prev_t2")
espen_prev_t2 = shape %>%
  dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
  left_join(mapped_prev_t2_means, by="IU_ID")
geomap_prev_t2 = st_as_sf(espen_prev_t2)

mapped_prev_t3_means = data.frame(species_map_allyears[[3]]$data[,c("IU_2021","mu")])%>%
  mutate(IU_ID = as.numeric(IU_2021),
         mapped_prev_t1 = logistic(mu)) %>%
  select(IU_ID,mapped_prev_t1)
colnames(mapped_prev_t3_means) = c("IU_ID","mapped_prev_t3")
espen_prev_t3 = shape %>%
  dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
  left_join(mapped_prev_t3_means, by="IU_ID")
geomap_prev_t3 = st_as_sf(espen_prev_t3)

p_mapped_prev1 <- ggplot(geomap_prev_t1) +
  geom_sf(aes(fill=mapped_prev_t1), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
  theme_bw() + theme(legend.position="right") + 
  guides(fill=guide_colorbar(title="Mean prevalence from \n geostatistical map")) +
  ggtitle(map_years[1]) + th
#p_mapped_prev1

p_mapped_prev2 <- ggplot(geomap_prev_t2) +
  geom_sf(aes(fill=mapped_prev_t2), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
  theme_bw() + theme(legend.position="right") + 
  guides(fill=guide_colorbar(title="Mean prevalence from \n geostatistical map")) +
  ggtitle(map_years[2]) + th
#p_mapped_prev2

p_mapped_prev3 <- ggplot(geomap_prev_t3) +
  geom_sf(aes(fill=mapped_prev_t3), lwd=lwd_borders) +
  scale_fill_viridis(option=colour, direction = 1, limits=c(0,1), na.value="gray80") +
  theme_bw() + theme(legend.position="right") + 
  guides(fill=guide_colorbar(title="Mean prevalence from \n geostatistical map")) +
  ggtitle(map_years[3]) + th
#p_mapped_prev3

plot_mapped_prevs <- p_mapped_prev1 + p_mapped_prev2 + p_mapped_prev3 & theme(legend.position = "right")
plot_mapped_prevs <- plot_mapped_prevs + plot_layout(guides = "collect")

fileName <- paste0(pathForPlots,"geomap_mean_prevs_all_countries_",species_plot,".png")
ggsave(fileName, plot_mapped_prevs, width = 120*3, height = 120, units = "mm")


if(plot_trajectories){
  
  #### Plot draws from posterior
  
  pdf(paste0(pathForPlots, "/traj_post_all_countries_",species,".pdf"), 
      width=10, height=panel_nrows*2.5)
  par(mfrow=c(panel_nrows,panel_ncols))
  for (l in ixd_ord_traj_plots){
    
    iu <- table_country$IU_CODE[l]
    id <- table_country$TaskID[l]
    # idx <- table_country$idx[l] # wrong
    
    sampled_params_iu <- sampled_params_all %>%
      filter(IU_ID == iu)
    
    
    if (iu %in% ius_requiring_sigma0.025){
      load(paste0("../trajectories/trajectories_",id,"_",species,"_sigma0.025.Rdata")) # load 'trajectories'
      load(paste0("../AMIS_output/",species,"_amis_output",id,"_sigma0.025.Rdata")) # loads amis_output
    } else {
      load(paste0("../trajectories/trajectories_",id,"_",species,".Rdata")) # load 'trajectories'
      load(paste0("../AMIS_output/",species,"_amis_output",id,".Rdata")) # loads amis_output
    }
    prevalence_map <- amis_output$prevalence_map
    
    plot(x=all_years, y=trajectories[sampled_params_iu$seed[1],],type="l", ylim=c(0,1), xlim=range(all_years) + 2*c(-1,1),
         main=paste0(iu, "\n ESS ",round(amis_output$ess[which(names(amis_output$ess)==iu)],digits=2)), xlab="year", ylab="prevalence",xaxt="n",yaxt="n")
    axis(1, at=all_years)
    axis(2, at=c(0,0.5,1))
    sub_samp = sampled_params_iu$seed[-1]
    #sub_samp <- sample(sampled_params_iu$seed[-1], num_sub_samples, replace=F)
    for(i in seq_along(sub_samp)){
      lines(x=all_years, y=trajectories[sub_samp[i],])
    }
    
    
    prevalence_map_l <- lapply(1:length(prevalence_map), function(t) {
      map_t <- prevalence_map[[t]]$data
      return(list(data=map_t[which(rownames(map_t)==iu),,drop=F]))
    })
    logistic <- function(x){return(1/(1+exp(-x)))}
    for (year_ind in 1:length(prevalence_map_l)){
      year <- map_years[year_ind]
      params <- prevalence_map_l[[year_ind]]$data
      lo <- logistic(params[1,"mu"] - 1.96*params[1,"sigma"])
      up <- logistic(params[1,"mu"] + 1.96*params[1,"sigma"])
      mu <- logistic(params[1,"mu"])
      points(y=c(mu,lo,up),x=rep(year,3),col="red", cex=2, type = "b")
    }
  }
  dev.off()
  
#### Plot draws from prior
pdf(paste0(pathForPlots, "/traj_prior_all_countries_",species,".pdf"), 
    width=10, height=panel_nrows*2.5)
sub_samp <- sample(2:n_samples, num_sub_samples, replace=T) # show same subsamples for all locations as it's the same prior for all locations
par(mfrow=c(panel_nrows,panel_ncols))

for (l in ixd_ord_traj_plots){
  
  iu <- table_country$IU_CODE[l]
  id <- table_country$TaskID[l]
  # idx <- table_country$idx[l] # wrong
  
  if (iu %in% ius_requiring_sigma0.025){
    load(paste0("../trajectories/trajectories_",id,"_",species,"_sigma0.025.Rdata")) # load 'trajectories'
    load(paste0("../AMIS_output/",species,"_amis_output",id,"_sigma0.025.Rdata")) # loads amis_output
  } else {
    load(paste0("../trajectories/trajectories_",id,"_",species,".Rdata")) # load 'trajectories'
    load(paste0("../AMIS_output/",species,"_amis_output",id,".Rdata")) # loads amis_output
  }
  
  plot(x=all_years, y=trajectories[1,],type="l", ylim=c(0,1), xlim=range(all_years) + 2*c(-1,1),
       main=iu, xlab="year", ylab="prevalence",xaxt="n",yaxt="n")
  axis(1, at=all_years)
  axis(2, at=c(0,0.5,1))
  for(i in seq_along(sub_samp)){
    lines(x=all_years, y=trajectories[sub_samp[i],])
  }
  
  prevalence_map <- amis_output$prevalence_map
  
  prevalence_map_l <- lapply(1:length(prevalence_map), function(t) {
    map_t <- prevalence_map[[t]]$data
    return(list(data=map_t[which(rownames(map_t)==iu),,drop=F]))
  })
  logistic <- function(x){return(1/(1+exp(-x)))}
  for (year_ind in 1:length(prevalence_map_l)){
    year <- map_years[year_ind]
    params <- prevalence_map_l[[year_ind]]$data
    lo <- logistic(params[1,"mu"] - 1.96*params[1,"sigma"])
    up <- logistic(params[1,"mu"] + 1.96*params[1,"sigma"])
    mu <- logistic(params[1,"mu"])
    points(y=c(mu,lo,up),x=rep(year,3),col="red", cex=2, type = "b")
  }
}
dev.off()
}


