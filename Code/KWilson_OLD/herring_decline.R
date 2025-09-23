# instructions
# - convert numbers at age vector to age comp vector by dividing numbers at age vector by sum of numbers at age. Then:
# - maturity vector * age comp vector = mature age comp vector
# - mature age comp vector * ssb = ssb by age vector
# - ssb by age vector / weight at age vector = mature individuals at age vector

# For generation length: Age at 50% maturity is ~2
# mean age is ~4 from eyeballing Figure 5 in the CSAS doc (this could be calculated from the .dat files for each SAR). Age structure has almost certainty been truncated from over a half a century of size selective harvest, but I am not sure there is data to inform pre-commercial fishery age structure. So all to say I think we could probably get away with assuming 5 years, with some caveated language to accompany assumption.  

setwd("~/My Drive/CCIRA - AQB/COSEWIC Marine Fishes SSC/Data")
columns <- c("Year","Gear","Area","Group","Sex",2:10)
CC <- read.csv("iscam_sbt_mcmc_CC.csv",header=TRUE)
HG <- read.csv("iscam_sbt_mcmc_HG.csv",header=TRUE)
PRD <- read.csv("iscam_sbt_mcmc_PRD.csv",header=TRUE)
SoG <- read.csv("iscam_sbt_mcmc_SoG.csv",header=TRUE)
WCVI <- read.csv("iscam_sbt_mcmc_WCVI.csv",header=TRUE)

gen_length <- 5 # from DFO Stock Status update 2022
# Central Coast

declines <- function(region,post,file,gear_line,year_line,wt_line,mat_line,year_line2,num_age_line)
{
  ngears <- scan(file=file,skip=gear_line,nlines=1,sep=" ",comment.char="#")
  gear_types <- factor(c(1,2,3),labels=c("Other","RoeSN","RoeGN"))
  nyears <- scan(file=file,skip=year_line,nlines=1,sep=" ",comment.char="#")
  wt_age <- scan(file=file,skip=wt_line,nlines=nyears,sep="\t",comment.char="#")
  wt_age <- matrix(wt_age,nrow=nyears,ncol=14,byrow=TRUE)
  wt_age <- as.data.frame(wt_age)
  colnames(wt_age) <- columns
  mat_vec <- scan(file=file,skip=mat_line,nlines=1,n=9,sep=",",comment.char="#")
  nyears <- scan(file=file,skip=year_line2,nlines=1,sep="\t",comment.char="#",n=ngears)
  num_age <- scan(file=file,skip=num_age_line,nlines=sum(nyears),sep="\t",comment.char="#")
  num_age <- matrix(num_age,nrow=sum(nyears),ncol=length(columns),byrow=TRUE)
  num_age <- as.data.frame(num_age)
  colnames(num_age) <- columns
  
  prop_age <- num_age[num_age$Gear!=3,]
  dup_years <- prop_age$Year[duplicated(prop_age$Year)]
  prop_age <- prop_age[!(prop_age$Year%in%dup_years & prop_age$Gear!=2),]
  prop_age[,-(1:5)] <- prop_age[,-(1:5)]/rowSums(prop_age[,-(1:5)])
  years <- c(min(prop_age$Year):max(prop_age$Year))
  miss_years <- years[!years%in%unique(prop_age$Year)]
  prop_age_miss <- data.frame("Year"=miss_years,"Gear"=NA,"Area"=NA,"Group"=NA,"Sex"=NA)
  miss_dat <- cbind(prop_age_miss,matrix(NA,nrow=length(miss_years),ncol=9))
  colnames(miss_dat)[-(1:5)] <- names(prop_age[,-(1:5)])
  miss_dat[,-(1:5)] <- t(sapply(miss_dat$Year,function(x){colMeans(prop_age[which.min(abs(x-prop_age$Year)),-(1:5)])}))
  prop_age <- merge(prop_age,miss_dat,all=TRUE)
  mat_age <- prop_age
  mat_age[,-(1:5)] <- mat_age[,-(1:5)]*mat_vec
  
  ssb_mn <- apply(post,2,median)
  names(ssb_mn) <- gsub("sbt1_","",names(ssb_mn))
  ssb_df <- data.frame("Region"=region,"Year"=as.numeric(names(ssb_mn)),"ssb"=ssb_mn)
  ssb_df <- ssb_df[ssb_df$Year!=2023,]
  ssb_age <- mat_age
  ssb_age[,-(1:5)] <- mat_age[,-(1:5)]*ssb_df$ssb
  adult_age <- ssb_age
  adult_age[,-(1:5)] <- ssb_age[,-(1:5)]/wt_age[,-(1:5)]
  adult_time <- data.frame("Region"=region,"Year"=adult_age$Year,"N"=rowSums(adult_age[,-(1:5)]),"ssb"=ssb_df$ssb)
  adult_time <- adult_time[order(adult_time$Year),]
  gen_3 <- adult_time[(nrow(adult_time)-3*gen_length+1):nrow(adult_time),]
  ten_yr <- adult_time[(nrow(adult_time)-10+1):nrow(adult_time),]
  long_term <- adult_time[1:nrow(adult_time),]
  gen_3$pred <- predict(glm(N~Year,data=gen_3,family=gaussian(link="log")),data=gen_3,type="r")
  ten_yr$pred <- predict(glm(N~Year,data=ten_yr,family=gaussian(link="log")),data=ten_yr,type="r")
  long_term$pred <- predict(glm(N~Year,data=long_term,family=gaussian(link="log")),data=long_term,type="r")
  
  decline <- 100*((gen_3$pred[nrow(gen_3)]-gen_3$pred[1])/(gen_3$pred[1]))
  
  ten_decline <- 100*((ten_yr$pred[nrow(ten_yr)]-ten_yr$pred[1])/(ten_yr$pred[1]))
  total_decline <- 100*((long_term$pred[nrow(long_term)]-long_term$pred[1])/(long_term$pred[1]))
  ts_decline <- 100*sapply((3*gen_length+1):nrow(adult_time),function(x){(adult_time$N[x]-adult_time$N[x-3*gen_length])/adult_time$N[x-3*gen_length]})

return(list("ts_data"=adult_time,"Recent three gen. % change"=decline,"10-year % change"=ten_decline,"Total % change"=total_decline,"Annual_change"=data.frame("Region"=region,"Year"=years[(3*gen_length+1):nrow(adult_time)],"3-gen. decline"=ts_decline),"Adults"=adult_time))
}
CC_decline <- declines(region="CC",post=CC,file="HerringCC2022.dat",gear_line=219,year_line=327,wt_line=329,mat_line=29,year_line2=220,num_age_line=226)

HG_decline <- declines(region="HG",post=HG,file="HerringHG2022.dat",gear_line=195,year_line=282,wt_line=284,mat_line=29,year_line2=196,num_age_line=202)

PRD_decline <- declines(region="PRD",post=PRD,file="HerringPRD2022.dat",gear_line=254,year_line=389,wt_line=391,mat_line=29,year_line2=255,num_age_line=261)

SoG_decline <- declines(region="SoG",post=SoG,file="HerringSoG2022.dat",gear_line=290,year_line=467,wt_line=469,mat_line=29,year_line2=291,num_age_line=297)

WCVI_decline <- declines(region="WCVI",post=WCVI,file="HerringWCVI2022.dat",gear_line=208,year_line=313,wt_line=315,mat_line=29,year_line2=209,num_age_line=215)

decline_metrics <- merge(WCVI_decline$Annual_change,merge(SoG_decline$Annual_change,merge(PRD_decline$Annual_change,merge(CC_decline$Annual_change,HG_decline$Annual_change,all=TRUE),all=TRUE),all=TRUE),all=TRUE)

decline_metrics$X3.gen..decline[decline_metrics$X3.gen..decline>=100] <- 100
library(ggplot2)
herring <- ggplot(data=decline_metrics,aes(x=Year,y=X3.gen..decline,col=Region,fill=Region)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Region,scales="free_x") +
  theme_minimal() +
  scale_color_brewer(type="qual",palette = 2) +
  scale_fill_brewer(type="qual",palette = 2) +
  labs(y="% change in mature individuals (3 gen.)") +
  theme(legend.position = "top") +
  ylim(-100,100) +
  ggtitle(label="Annual change in Pacific herring over 3-generations",subtitle="Note: changes >100% shown as 100% & 3 gen. = 15 years")

ggsave(filename="pacific_herring.jpeg",plot=herring,units="in",width=7,height=5,dpi=600)

adults_change <- merge(WCVI_decline$Adults,merge(SoG_decline$Adults,merge(PRD_decline$Adults,merge(CC_decline$Adults,HG_decline$Adults,all=TRUE),all=TRUE),all=TRUE),all=TRUE)

library(ggplot2)
herring_adults <- ggplot(data=adults_change,aes(x=Year,y=N,col=Region,fill=Region)) +
  geom_point() +
  geom_line() +
  geom_smooth(data=adults_change,aes(x=Year,y=N),col="dodgerblue",fill="dodgerblue",method="glm",method.args=list(family=gaussian(link="log")))+
  geom_smooth(data=adults_change[adults_change$Year>=(2023-15),],aes(x=Year,y=N),col="black",fill="black",method="glm",method.args=list(family=gaussian(link="log")))+
  facet_wrap(~Region,scales="free_x") +
  theme_minimal() +
  scale_color_brewer(type="qual",palette = 2) +
  scale_fill_brewer(type="qual",palette = 2) +
  labs(y="Mature individuals") +
  theme(legend.position = "top") +
  ggtitle(label="Pacific herring mature individuals over time")

ggsave(filename="pacific_herring_adults.jpeg",plot=herring_adults,units="in",width=7,height=5,dpi=600)

herring_adults <- ggplot(data=adults_change,aes(x=Year,y=ssb,col=Region,fill=Region)) +
  geom_point() +
  geom_line() +
  geom_smooth(data=adults_change,aes(x=Year,y=ssb),col="dodgerblue",fill="dodgerblue",method="glm",method.args=list(family=gaussian(link="log")))+
  geom_smooth(data=adults_change[adults_change$Year>=(2023-15),],aes(x=Year,y=ssb),col="black",fill="black",method="glm",method.args=list(family=gaussian(link="log")))+
  facet_wrap(~Region,scales="free_x") +
  theme_minimal() +
  scale_color_brewer(type="qual",palette = 2) +
  scale_fill_brewer(type="qual",palette = 2) +
  labs(y="Spawner biomass") +
  theme(legend.position = "top") +
  ggtitle(label="Pacific herring spawner biomass over time")

ggsave(filename="pacific_herring_ssb.jpeg",plot=herring_adults,units="in",width=7,height=5,dpi=600)

decline_criteria <- data.frame("Region"=c("CC","HG","PRD","SoG","WCVI"),
                               "Recent 3-gen change"=round(c(CC_decline$`Recent three gen. % change`,HG_decline$`Recent three gen. % change`,PRD_decline$`Recent three gen. % change`,SoG_decline$`Recent three gen. % change`,WCVI_decline$`Recent three gen. % change`)),
                               "Long-term change"=round(c(CC_decline$`Total % change`,HG_decline$`Total % change`,PRD_decline$`Total % change`,SoG_decline$`Total % change`,WCVI_decline$`Total % change`)))

openxlsx::write.xlsx(decline_criteria,file="herring_decline.xlsx")


ts_data <- merge(WCVI_decline$ts_data,merge(SoG_decline$ts_data,merge(PRD_decline$ts_data,merge(CC_decline$ts_data,HG_decline$ts_data,all=TRUE),all=TRUE),all=TRUE),all=TRUE)
openxlsx::write.xlsx(ts_data,file="herring_time_series.xlsx")
