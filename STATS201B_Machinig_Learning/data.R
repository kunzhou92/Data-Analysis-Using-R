ecological_zones = read.csv(file="E:\\Papers\\Statistics\\202B Matrix Algebra and Optimization\\project\\data\\hcapi-agro-ecological_zones-ssa.csv", header = T)
ecological_zones = ecological_zones[c("AEZ16_CLAS")]
#Agro-Ecological Zones 
climate = read.csv(file="E:\\Papers\\Statistics\\202B Matrix Algebra and Optimization\\project\\data\\hcapi-climate-ssa.csv", header = T)
climate = climate[c("LGP_AVG", "LGP_CV", "pre_cv", "pre_mean")]
#Growing Period (mean), Growing Period (c.v.), Annual rainfall (c.v.), Annual rainfall (mean)
ELEVATION = read.csv(file="E:\\Papers\\Statistics\\202B Matrix Algebra and Optimization\\project\\data\\hcapi-ELEVATION--ssa.csv", header = T)
ELEVATION = ELEVATION[c("ELEVATION")]
#Elevation
farming_systems = read.csv(file="E:\\Papers\\Statistics\\202B Matrix Algebra and Optimization\\project\\data\\hcapi-farming_systems-ssa.csv", header = T)
farming_systems = farming_systems[c("FS_2012_TX")]
#Farming System Name
health_and_nutrition = read.csv(file="E:\\Papers\\Statistics\\202B Matrix Algebra and Optimization\\project\\data\\hcapi-health_and_nutrition-ssa.csv", header = T)
health_and_nutrition = health_and_nutrition[c("bmi_rur", "child_mortality_rur", "diarrhea_rur", "infant_mortality_rur", "stunted_low_rur", "underweight_low_rur")]
#Rural BMI, Rural Child Mortality, Children with Diarrhea - rural, Rural Infant Mortality, Rural Stunting, low (z-score -6 to -1), Rural Underweight, low (z-score -6 to -1)
income_and_poverty = read.csv(file="E:\\Papers\\Statistics\\202B Matrix Algebra and Optimization\\project\\data\\hcapi-income_and_poverty-ssa.csv", header = T)
income_and_poverty =  income_and_poverty[c("RPOV_GINI", "RPOV_PCEXP")]
#Rural Gini '05, Rural Poverty Per Cap. Exp. '05
land_cover_and_use = read.csv(file="E:\\Papers\\Statistics\\202B Matrix Algebra and Optimization\\project\\data\\hcapi-land_cover_and_use-ssa.csv", header = T)
land_cover_and_use = land_cover_and_use[c("AREA_CR_LC")]
#Cropland '00, 
livestock = read.csv(file="E:\\Papers\\Statistics\\202B Matrix Algebra and Optimization\\project\\data\\hcapi-livestock-ssa.csv", header = T)
livestock = livestock[c("AD05_CATT", "AD05_CHIC", "AD05_GOAT", "AD05_PIG", "AD05_SHEE")]
#Cattle density '05, Poultry density '05, Goat density '05, Pig density '05, Sheep density '05
MSH_50K_TX = read.csv(file="E:\\Papers\\Statistics\\202B Matrix Algebra and Optimization\\project\\data\\hcapi-MSH_50K_TX--ssa.csv", header = T)
MSH_50K_TX = MSH_50K_TX[c("MSH_50K_TX")]
#Marketsheds 50+
population = read.csv(file="E:\\Papers\\Statistics\\202B Matrix Algebra and Optimization\\project\\data\\hcapi-population-ssa.csv", header = T)
population = population[c("PD05_RUR")]
#Rural Pop. density '05
portshed = read.csv(file="E:\\Papers\\Statistics\\202B Matrix Algebra and Optimization\\project\\data\\hcapi-portshed-ssa.csv", header = T)
portshed = portshed[c("PSH_TX")]
#Portshed
production = read.csv(file="E:\\Papers\\Statistics\\202B Matrix Algebra and Optimization\\project\\data\\hcapi-production-ssa.csv", header = T)
names = names(production)
index = which(!grepl("_i", names) & !grepl("_r", names))
production = production[index[-(1:7)]]
#production
soil_resource = read.csv(file="E:\\Papers\\Statistics\\202B Matrix Algebra and Optimization\\project\\data\\hcapi-soil_resources-ssa.csv", header = T)
soil_resource = soil_resource[8:26]
#soil resource
travel_time = read.csv(file="E:\\Papers\\Statistics\\202B Matrix Algebra and Optimization\\project\\data\\hcapi-travel_time-ssa.csv", header = T)
travel_time = travel_time[c("TT_20K", "TT_PORT")]
#Travel Time 20K+, Travel Time Nearest Port
value_of_production = read.csv(file="E:\\Papers\\Statistics\\202B Matrix Algebra and Optimization\\project\\data\\hcapi-value_of_production-ssa.csv", header = T)
value_of_production = value_of_production[c("vp_food_ar", "vp_food", "vp_nonf", "vp_nonf_ar")]
#Food val. prod. per hectare '05, Food val. prod. '05, Non-food val. prod. '05, Non-food val. prod. per hectare '05
location = read.csv(file="E:\\Papers\\Statistics\\202B Matrix Algebra and Optimization\\project\\data\\hcapi-agro-ecological_zones-ssa.csv", header = T)
location = location[c("X", "Y")]
#Longitude, Latitude)


bonus = cbind(health_and_nutrition, income_and_poverty)
agriculture = cbind(ecological_zones, livestock, production, value_of_production)
climate_geology = cbind(climate, ELEVATION, land_cover_and_use, soil_resource)
social_factor = cbind(MSH_50K_TX, portshed, travel_time)
index = complete.cases(bonus) & complete.cases(agriculture) & complete.cases(climate_geology) & complete.cases(social_factor)
bouns.c = bonus[index,]
agriculture.c = agriculture[index,]
climate_geology.c = climate_geology[index,]
social_factor.c = social_factor[index,]
location.c = location[index,]
#########################
set.seed(26)
s = sample(1:nrow(index), size=floor(0.3*nrow(index)))
index.train = index
index.train[s] = F
index.test = index
index.test[-s] = F
bonus.c.tra = bonus[index.train,]
bonus.c.tes = bonus[index.test,]
agriculture.c.tra = agriculture[index.train,]
agriculture.c.tes = agriculture[index.test,]
climate_geology.c.tra = climate_geology[index.train,]
climate_geology.c.tes = climate_geology[index.test,]
social_factor.c.tra = social_factor[index.train,]
social_factor.c.tes = social_factor[index.test,]
location.c.tra = location[index.train,]
location.c.tes = location[index.test,]
########################################################
#OLS
temp = livestock[index,]
livestock_cli = lm(cbind(AD05_CATT,  AD05_CHIC, AD05_GOAT, AD05_PIG, AD05_SHEE)~
LGP_AVG+LGP_CV+pre_cv+pre_mean+ELEVATION+AREA_CR_LC
, data = cbind(agriculture.c.tra, climate_geology.c.tra))
summary(livestock_cli)

data = cbind(agriculture.c.tra, climate_geology.c.tra)
data = data[data$AD05_CATT!=0,]
catt.lm = lm(log(AD05_CATT)~LGP_AVG+LGP_CV+pre_cv+pre_mean+ELEVATION+AREA_CR_LC, 
             data = data)
summary(catt.lm)

catt.pred = exp(predict.lm(catt.lm, newdata=cbind(agriculture.c.tes, climate_geology.c.tes)))
wilcox.test(agriculture.c.tes$AD05_CATT, catt.pred, paired=TRUE) 
quantile(agriculture.c.tes$AD05_CATT, 0.7)
quantile(catt.pred, 0.7)


data = cbind(agriculture.c.tra, climate_geology.c.tra)
data = data[data$AD05_GOAT!=0,]
goat.lm = lm(log(AD05_GOAT)~LGP_AVG+LGP_CV+pre_cv+pre_mean+ELEVATION+AREA_CR_LC, 
             data = data)
summary(goat.lm)
goat.pred = exp(predict.lm(goat.lm, newdata=cbind(agriculture.c.tes, climate_geology.c.tes)))
wilcox.test(agriculture.c.tes$AD05_GOAT, goat.pred, paired=TRUE) 
quantile(agriculture.c.tes$AD05_GOAT, 0.7)
quantile(goat.pred, 0.7)

###########tls


