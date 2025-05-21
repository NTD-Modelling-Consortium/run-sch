
#### This script depends on outputs from prepare_histories.R and prepare_histories_trichuris.R ####

#library(haven)
library(readxl)
library(dplyr)
#library(ggplot2)
#library(gridExtra)
library(pracma)
#setwd("~/Documents/STH-endgame/Maps")

# Functions to calculate mu and sigma 
# (obtained by taking upper and lower values of the distribution and plugging into the relevant CDF)
calculate_mu = function(pL,pU){
  (qlogis(pL)+qlogis(pU))/2
}

calculate_sigma = function(pL,pU){
  mu=(qlogis(pL)+qlogis(pU))/2
  sigma=(qlogis(pU)-mu)/(sqrt(2)*erfinv(0.95))
}

logitN_likelihood = function(data, sim_prev, log){
  logit_prev = qlogis(sim_prev)
  mu = data[1]
  sigma = data[2]
  lh = dnorm(logit_prev, mu, sigma, log=log)
  return(lh)
}

# Read in data and remove problematic IUs
sth_distributions = read.csv("sartorious_2021_withTaskID.csv") %>%
  mutate(logitp_hook_mu = calculate_mu(p1L,p1U),
         logitp_hook_sigma= calculate_sigma(p1L,p1U),
         logitp_ascaris_mu = calculate_mu(p2L,p2U),
         logitp_ascaris_sigma = calculate_sigma(p2L,p2U)) 

#### Hookworm 

# Get relevant parts of the data frame for each year
hook_distributions_2000 = sth_distributions %>%
  filter(Year==2000) %>%
  select(IU_2021,TaskID,logitp_hook_mu,logitp_hook_sigma) 
colnames(hook_distributions_2000) = c("IU_2021","TaskID","mu","sigma")
rownames(hook_distributions_2000) = hook_distributions_2000$IU_2021

hook_distributions_2013 = sth_distributions %>%
  filter(Year==2013) %>%
  select(IU_2021,TaskID,logitp_hook_mu,logitp_hook_sigma) 
colnames(hook_distributions_2013) = c("IU_2021","TaskID","mu","sigma")
rownames(hook_distributions_2013) = hook_distributions_2013$IU_2021

hook_distributions_2018 = sth_distributions %>%
  filter(Year==2018) %>%
  select(IU_2021,TaskID,logitp_hook_mu,logitp_hook_sigma) 
colnames(hook_distributions_2018) = c("IU_2021","TaskID","mu","sigma")
rownames(hook_distributions_2018) = hook_distributions_2018$IU_2021

# Generate samples for each year and transform back to prevalence scale
hook_samples_2000 = cbind(hook_distributions_2000[,1:2],
                      t(sapply(1:nrow(hook_distributions_2000), function(x){
                        plogis(rnorm(1000,mean=as.numeric(hook_distributions_2000[x,"mu"]),sd=as.numeric(hook_distributions_2000[x,"sigma"])))
                      })))

hook_samples_2013 = cbind(hook_distributions_2013[,1:2],
                      t(sapply(1:nrow(hook_distributions_2013), function(x){
                        plogis(rnorm(1000,mean=as.numeric(hook_distributions_2013[x,"mu"]),sd=as.numeric(hook_distributions_2013[x,"sigma"])))
                      })))

hook_samples_2018 = cbind(hook_distributions_2018[,1:2],
                      t(sapply(1:nrow(hook_distributions_2018), function(x){
                        plogis(rnorm(1000,mean=as.numeric(hook_distributions_2018[x,"mu"]),sd=as.numeric(hook_distributions_2018[x,"sigma"])))
                      })))

# Combine maps and save
hookworm_map_allyears = list(list(data=hook_distributions_2000, likelihood=logitN_likelihood),
                         list(data=hook_distributions_2013, likelihood=logitN_likelihood),
                         list(data=hook_distributions_2018, likelihood=logitN_likelihood))
save(hookworm_map_allyears,file="hookworm_maps.rds")


#### Ascaris 

# Get relevant parts of the data frame for each year
ascaris_distributions_2000 = sth_distributions %>%
  filter(Year==2000) %>%
  select(IU_2021,TaskID,logitp_ascaris_mu,logitp_ascaris_sigma) 
colnames(ascaris_distributions_2000) = c("IU_2021","TaskID","mu","sigma")
rownames(ascaris_distributions_2000) = ascaris_distributions_2000$IU_2021

ascaris_distributions_2013 = sth_distributions %>%
  filter(Year==2013) %>%
  select(IU_2021,TaskID,logitp_ascaris_mu,logitp_ascaris_sigma) 
colnames(ascaris_distributions_2013) = c("IU_2021","TaskID","mu","sigma")
rownames(ascaris_distributions_2013) = ascaris_distributions_2013$IU_2021

ascaris_distributions_2018 = sth_distributions %>%
  filter(Year==2018) %>%
  select(IU_2021,TaskID,logitp_ascaris_mu,logitp_ascaris_sigma) 
colnames(ascaris_distributions_2018) = c("IU_2021","TaskID","mu","sigma")
rownames(ascaris_distributions_2018) = ascaris_distributions_2018$IU_2021

# Generate samples for each year and transform back to prevalence scale
ascaris_samples_2000 = cbind(ascaris_distributions_2000[,1:2],
                      t(sapply(1:nrow(ascaris_distributions_2000), function(x){
                        plogis(rnorm(1000,mean=as.numeric(ascaris_distributions_2000[x,"mu"]),sd=as.numeric(ascaris_distributions_2000[x,"sigma"])))
                      })))

ascaris_samples_2013 = cbind(ascaris_distributions_2013[,1:2],
                      t(sapply(1:nrow(ascaris_distributions_2013), function(x){
                        plogis(rnorm(1000,mean=as.numeric(ascaris_distributions_2013[x,"mu"]),sd=as.numeric(ascaris_distributions_2013[x,"sigma"])))
                      })))

ascaris_samples_2018 = cbind(ascaris_distributions_2018[,1:2],
                      t(sapply(1:nrow(ascaris_distributions_2018), function(x){
                        plogis(rnorm(1000,mean=as.numeric(ascaris_distributions_2018[x,"mu"]),sd=as.numeric(ascaris_distributions_2018[x,"sigma"])))
                      })))

# removing map samples of the year 2013 for specific IUs
# IUs_2013_to_rm <- c(36874, 36869, 36870, 49710, 18686)
IUs_2013_to_rm <- c(36874, 36869, 36870, 49710, 18686, 36865, 37001, 37009, 36870, 49778, 18686, 50884)
wh <- which(ascaris_distributions_2013$IU_2021%in%IUs_2013_to_rm)
ascaris_distributions_2013[wh, c("mu", "sigma")] <- NA

# Combine maps and save
ascaris_map_allyears = list(list(data=ascaris_distributions_2000, likelihood=logitN_likelihood),
                            list(data=ascaris_distributions_2013, likelihood=logitN_likelihood),
                            list(data=ascaris_distributions_2018, likelihood=logitN_likelihood))
save(ascaris_map_allyears,file="ascaris_maps.rds")



#### Trichuris 

# Read in data and remove problematic IUs (for trichuris)
sth_distributions_tri = read.csv("sartorious_2021_withTaskID_trichuris.csv") %>%
  mutate(logitp_tri_mu =  calculate_mu(p3L,p3U),
         logitp_tri_sigma = calculate_sigma(p3L,p3U)) 

# Get relevant parts of the data frame for each year
tri_distributions_2000 = sth_distributions_tri %>%
  filter(Year==2000) %>%
  select(IU_2021,TaskID,logitp_tri_mu,logitp_tri_sigma) 
colnames(tri_distributions_2000) = c("IU_2021","TaskID","mu","sigma")
rownames(tri_distributions_2000) = tri_distributions_2000$IU_2021

tri_distributions_2013 = sth_distributions_tri %>%
  filter(Year==2013) %>%
  select(IU_2021,TaskID,logitp_tri_mu,logitp_tri_sigma) 
colnames(tri_distributions_2013) = c("IU_2021","TaskID","mu","sigma")
rownames(tri_distributions_2013) = tri_distributions_2013$IU_2021

tri_distributions_2018 = sth_distributions_tri %>%
  filter(Year==2018) %>%
  select(IU_2021,TaskID,logitp_tri_mu,logitp_tri_sigma) 
colnames(tri_distributions_2018) = c("IU_2021","TaskID","mu","sigma")
rownames(tri_distributions_2018) = tri_distributions_2018$IU_2021

# Generate samples for each year and transform back to prevalence scale
tri_samples_2000 = cbind(tri_distributions_2000[,1:2],
                         t(sapply(1:nrow(tri_distributions_2000), function(x){
                           plogis(rnorm(1000,mean=as.numeric(tri_distributions_2000[x,"mu"]),sd=as.numeric(tri_distributions_2000[x,"sigma"])))
                         })))

tri_samples_2013 = cbind(tri_distributions_2013[,1:2],
                         t(sapply(1:nrow(tri_distributions_2013), function(x){
                           plogis(rnorm(1000,mean=as.numeric(tri_distributions_2013[x,"mu"]),sd=as.numeric(tri_distributions_2013[x,"sigma"])))
                         })))

tri_samples_2018 = cbind(tri_distributions_2018[,1:2],
                         t(sapply(1:nrow(tri_distributions_2018), function(x){
                           plogis(rnorm(1000,mean=as.numeric(tri_distributions_2018[x,"mu"]),sd=as.numeric(tri_distributions_2018[x,"sigma"])))
                         })))

# removing map samples of the year 2000 for IUs 19448 and 19517
wh <- which(tri_distributions_2000$IU_2021==19448 | tri_distributions_2000$IU_2021==19517)
tri_distributions_2000[wh, c("mu", "sigma")] <- NA

# Combine maps and save
trichuris_map_allyears = list(list(data=tri_distributions_2000, likelihood=logitN_likelihood),
                              list(data=tri_distributions_2013, likelihood=logitN_likelihood),
                              list(data=tri_distributions_2018, likelihood=logitN_likelihood))
save(trichuris_map_allyears,file="trichuris_maps.rds")

cat("Finished running prepare_maps_allspecies.R \n")

# #### Plot distributions by year
# # hookworm
# ggplot(data=sth_distributions) + 
#   geom_segment(aes(x=p1L,xend=p1U,y=idlim,alpha=0.1)) + 
#   geom_point(aes(x=p1,y=idlim,alpha=0.05),colour="red") + 
#   facet_wrap(~Year)
# ggsave("plots/hookworm_maps.pdf")
# # ascaris
# ggplot(data=sth_distributions) + 
#   geom_segment(aes(x=p2L,xend=p2U,y=idlim,alpha=0.1)) + 
#   geom_point(aes(x=p2,y=idlim,alpha=0.05),colour="red") + 
#   facet_wrap(~Year)
# ggsave("plots/ascaris_maps.pdf")
# # trichuris
# ggplot(data=sth_distributions) + 
#   geom_segment(aes(x=p3L,xend=p3U,y=idlim,alpha=0.1)) + 
#   geom_point(aes(x=p3,y=idlim,alpha=0.05),colour="red") + 
#   facet_wrap(~Year)
# ggsave("plots/trichuris_maps.pdf")


# # Compare 2000 and 2013 (we don't have MDA history between these years) and 2018
# library(sf)
# library(magrittr)
# library(viridis)
# library(gridExtra)
# 
# shape <- read_sf(dsn = "../../ESPEN_IU_2021/", layer = "ESPEN_IU_2021")
# shape %<>% filter(ADMIN0ISO3 %in% unique(sth_distributions$ADMIN0ISO3))
# 
# #Hookworm
# hook_2000_means = sth_distributions %>% filter(Year == 2000) %>% select(IU_2021,logitp_hook_mu)
# colnames(hook_2000_means) = c("IU_ID","hook_2000")
# espen_hook_2000 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(hook_2000_means, by="IU_ID")
# map_hook_2000 = st_as_sf(espen_hook_2000)
# 
# hook_2000_map=ggplot(map_hook_2000) +
#   geom_sf(aes(fill=hook_2000), lwd = 0)+
#   xlab("hook_2000") +
#   scale_fill_viridis(option = "C", direction = -1, limits=c(-15,5)) +
#   scale_colour_manual(na.value="gray") + 
#   theme_bw() + theme(legend.position="bottom")
# 
# 
# hook_2013_means = sth_distributions %>% filter(Year == 2013) %>% select(IU_2021,logitp_hook_mu)
# colnames(hook_2013_means) = c("IU_ID","hook_2013")
# espen_hook_2013 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(hook_2013_means, by="IU_ID")
# map_hook_2013 = st_as_sf(espen_hook_2013)
# 
# hook_2013_map=ggplot(map_hook_2013) +
#   geom_sf(aes(fill=hook_2013), lwd = 0)+
#   xlab("hook_2013") +
#   scale_fill_viridis(option = "C", direction = -1, limits=c(-15,5)) +
#   scale_colour_manual(na.value="gray") + 
#   theme_bw() + theme(legend.position="bottom")
# 
# hook_2018_means = sth_distributions %>% filter(Year == 2018) %>% select(IU_2021,logitp_hook_mu)
# colnames(hook_2018_means) = c("IU_ID","hook_2018")
# espen_hook_2018 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(hook_2018_means, by="IU_ID")
# map_hook_2018 = st_as_sf(espen_hook_2018)
# 
# hook_2018_map=ggplot(map_hook_2018) +
#   geom_sf(aes(fill=hook_2018), lwd = 0)+
#   xlab("hook_2018") +
#   scale_fill_viridis(option = "C", direction = -1, limits=c(-15,5)) +
#   scale_colour_manual(na.value="gray") + 
#   theme_bw() + theme(legend.position="bottom")
# 
# png("plots/maps_hook.png" , width=21, height=10,units ="in",res=1500)
# grid.arrange(hook_2000_map,hook_2013_map,hook_2018_map,nrow=1)
# dev.off()
# 
# 
# # Ascaris
# ascaris_2000_means = sth_distributions %>% filter(Year == 2000) %>% select(IU_2021,logitp_ascaris_mu)
# colnames(ascaris_2000_means) = c("IU_ID","ascaris_2000")
# espen_ascaris_2000 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(ascaris_2000_means, by="IU_ID")
# map_ascaris_2000 = st_as_sf(espen_ascaris_2000)
# 
# ascaris_2000_map=ggplot(map_ascaris_2000) +
#   geom_sf(aes(fill=ascaris_2000), lwd = 0)+
#   xlab("ascaris_2000") +
#   scale_fill_viridis(option = "C", direction = -1, limits=c(-15,5)) +
#   scale_colour_manual(na.value="gray") + 
#   theme_bw() + theme(legend.position="bottom")
# 
# 
# ascaris_2013_means = sth_distributions %>% filter(Year == 2013) %>% select(IU_2021,logitp_ascaris_mu)
# colnames(ascaris_2013_means) = c("IU_ID","ascaris_2013")
# espen_ascaris_2013 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(ascaris_2013_means, by="IU_ID")
# map_ascaris_2013 = st_as_sf(espen_ascaris_2013)
# 
# ascaris_2013_map=ggplot(map_ascaris_2013) +
#   geom_sf(aes(fill=ascaris_2013), lwd = 0)+
#   xlab("ascaris_2013") +
#   scale_fill_viridis(option = "C", direction = -1, limits=c(-15,5)) +
#   scale_colour_manual(na.value="gray") + 
#   theme_bw() + theme(legend.position="bottom")
# 
# ascaris_2018_means = sth_distributions %>% filter(Year == 2018) %>% select(IU_2021,logitp_ascaris_mu)
# colnames(ascaris_2018_means) = c("IU_ID","ascaris_2018")
# espen_ascaris_2018 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(ascaris_2018_means, by="IU_ID")
# map_ascaris_2018 = st_as_sf(espen_ascaris_2018)
# 
# ascaris_2018_map=ggplot(map_ascaris_2018) +
#   geom_sf(aes(fill=ascaris_2018), lwd = 0)+
#   xlab("ascaris_2018") +
#   scale_fill_viridis(option = "C", direction = -1, limits=c(-15,5)) +
#   scale_colour_manual(na.value="gray") + 
#   theme_bw() + theme(legend.position="bottom")
# 
# png("plots/maps_ascaris.png" , width=21, height=10,units ="in",res=1500)
# grid.arrange(ascaris_2000_map,ascaris_2013_map,ascaris_2018_map,nrow=1)
# dev.off()
# 
# # Trichuris
# tri_2000_means = sth_distributions %>% filter(Year == 2000) %>% select(IU_2021,logitp_tri_mu)
# colnames(tri_2000_means) = c("IU_ID","tri_2000")
# espen_tri_2000 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(tri_2000_means, by="IU_ID")
# map_tri_2000 = st_as_sf(espen_tri_2000)
# 
# tri_2000_map=ggplot(map_tri_2000) +
#   geom_sf(aes(fill=tri_2000), lwd = 0)+
#   xlab("tri_2000") +
#   scale_fill_viridis(option = "C", direction = -1, limits=c(-15,5)) +
#   scale_colour_manual(na.value="gray") + 
#   theme_bw() + theme(legend.position="bottom")
# 
# 
# tri_2013_means = sth_distributions %>% filter(Year == 2013) %>% select(IU_2021,logitp_tri_mu)
# colnames(tri_2013_means) = c("IU_ID","tri_2013")
# espen_tri_2013 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(tri_2013_means, by="IU_ID")
# map_tri_2013 = st_as_sf(espen_tri_2013)
# 
# tri_2013_map=ggplot(map_tri_2013) +
#   geom_sf(aes(fill=tri_2013), lwd = 0)+
#   xlab("tri_2013") +
#   scale_fill_viridis(option = "C", direction = -1, limits=c(-15,5)) +
#   scale_colour_manual(na.value="gray") + 
#   theme_bw() + theme(legend.position="bottom")
# 
# 
# tri_2018_means = sth_distributions %>% filter(Year == 2018) %>% select(IU_2021,logitp_tri_mu)
# colnames(tri_2018_means) = c("IU_ID","tri_2018")
# espen_tri_2018 = shape %>%
#   dplyr::select(IU_ID,Shape_Leng,Shape_Area,geometry) %>%
#   left_join(tri_2018_means, by="IU_ID")
# map_tri_2018 = st_as_sf(espen_tri_2018)
# 
# tri_2018_map=ggplot(map_tri_2018) +
#   geom_sf(aes(fill=tri_2018), lwd = 0)+
#   xlab("tri_2018") +
#   scale_fill_viridis(option = "C", direction = -1, limits=c(-15,5)) +
#   scale_colour_manual(na.value="gray") + 
#   theme_bw() + theme(legend.position="bottom")
# 
# png("plots/maps_tri.png" , width=21, height=10,units ="in",res=1500)
# grid.arrange(tri_2000_map,tri_2013_map,tri_2018_map,nrow=1)
# dev.off()


# # Compare the simulated maps against the distribution provided in the data
# 
# test=cbind(hook_samples_2000[,1:2],
#            t(apply(hook_samples_2000[,3:1002],1,function(x) quantile(x,prob=c(0.025,0.5,0.975),na.rm=T))),
#            mean=apply(hook_samples_2000[,3:1002],1,function(x) mean(x,na.rm=T))) %>%
#   left_join(sth_distributions %>% select(IU_2021,Year,p1L,p1,p1U))%>%
#   select(IU_2021,Year,'2.5%',p1L,'mean',p1,'97.5%',p1U)
# colnames(test) = c("IU_2021","Year","quant_2.5","p1L","mean","p1","quant_97.5","p1U")
# 
# 
# hookworm_plot = test %>% 
#   arrange(p1) %>%
#   mutate(id = row_number()) %>%
#   ggplot() +
#   geom_ribbon(aes(x=id,ymin=quant_2.5, ymax=quant_97.5,alpha=0.5),fill="red")  +
#   geom_ribbon(aes(x=id,ymin=p1L, ymax=p1U,alpha=0.5),fill="blue")+
#   geom_point(aes(x=id,y=mean,alpha=0.5,colour="Transformed")) +
#   geom_point(aes(x=id,y=p1,alpha=0.5,colour="Original"))+
#   ggtitle("Hookworm") + 
#   ylab("Prevalence") +
#   scale_color_manual(values = c("Original"="blue","Transformed"="red"),
#                      labels = c("Original","Transformed"))
# 
# 
# test2=cbind(ascaris_samples_2000[,1:2],
#             t(apply(ascaris_samples_2000[,3:1002],1,function(x) quantile(x,prob=c(0.025,0.5,0.975),na.rm=T))),
#             mean=apply(ascaris_samples_2000[,3:1002],1,function(x) mean(x,na.rm=T))) %>%
#   left_join(sth_distributions %>% select(IU_2021,Year,p2L,p2,p2U))%>%
#   select(IU_2021,Year,'2.5%',p2L,'mean',p2,'97.5%',p2U) 
# colnames(test2) = c("IU_2021","Year","quant_2.5","p2L","mean","p2","quant_97.5","p2U")
# 
# ascaris_plot = test2 %>% 
#   arrange(p2) %>%
#   mutate(id = row_number()) %>%
#   ggplot() +
#   geom_ribbon(aes(x=id,ymin=quant_2.5, ymax=quant_97.5,alpha=0.5),fill="red")  +
#   geom_ribbon(aes(x=id,ymin=p2L, ymax=p2U,alpha=0.5),fill="blue")+
#   geom_point(aes(x=id,y=mean,alpha=0.5,colour="Transformed")) +
#   geom_point(aes(x=id,y=p2,alpha=0.5,colour="Original")) +
#   ggtitle("Ascaris") + 
#   ylab("Prevalence") +
#   scale_color_manual(values = c("Original"="blue","Transformed"="red"),
#                      labels = c("Original","Transformed"))
# 
# test3=cbind(tri_samples_2000[,1:2],
#             t(apply(tri_samples_2000[,3:1002],1,function(x) quantile(x,prob=c(0.025,0.975),na.rm=T))),
#             mean=apply(tri_samples_2000[,3:1002],1,function(x) mean(x,na.rm=T))) %>%
#   left_join(sth_distributions %>% select(IU_2021,Year,p3L,p3,p3U)) %>%
#   select(IU_2021,Year,'2.5%',p3L,'mean',p3,'97.5%',p3U)
# colnames(test3) = c("IU_2021","Year","quant_2.5","p3L","mean","p3","quant_97.5","p3U")
# 
# tri_plot = test3 %>% 
#   arrange(p3) %>%
#   mutate(id = row_number()) %>%
#   ggplot() +
#   geom_ribbon(aes(x=id,ymin=quant_2.5, ymax=quant_97.5,alpha=0.5),fill="red")  +
#   geom_ribbon(aes(x=id,ymin=p3L, ymax=p3U,alpha=0.5),fill="blue")+
#   geom_point(aes(x=id,y=mean,alpha=0.5,colour="Transformed")) +
#   geom_point(aes(x=id,y=p3,alpha=0.5,colour="Original")) +
#   ggtitle("Trichuris trichiura") + 
#   ylab("Prevalence") +
#   scale_color_manual(values = c("Original"="blue","Transformed"="red"),
#                      labels = c("Original","Transformed"))
# 
# pdf("~/Documents/STH-endgame/Maps/plots/transformed_maps_v2.pdf", width=7, height=15)
# grid.arrange(hookworm_plot,ascaris_plot,tri_plot,nrow=3)
# dev.off()
# 
# 
# ggplot(data=sth_distributions %>% filter(Year==2000)) +
#   geom_point(aes(x=qlogis(p1),y=qlogis(p1L),alpha=0.1),colour="red") +
#   geom_point(aes(x=qlogis(p1),y=qlogis(p1U),alpha=0.1),colour="blue")
# 
# ggplot(data=sth_distributions %>% filter(Year==2000)) +
#   geom_point(aes(x=qlogis(p1),y=(qlogis(p1U)-qlogis(p1L))/2, alpha=0.1),colour="blue") 
# 
# 

