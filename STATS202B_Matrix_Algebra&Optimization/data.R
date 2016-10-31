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
s = sample(1:length(index), size=floor(0.3*length(index)))
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
temp = cbind(climate, ELEVATION, land_cover_and_use)
cor(temp[index.train,])

soil_resource.tra = soil_resource[index.train,]
soil_resource.pca = prcomp(soil_resource.tra)
sum((soil_resource.pca$sdev[1:5])^2)/sum((soil_resource.pca$sdev)^2)
soil_resource.pca$rotation
a = scale(soil_resource.tra, center=T, scale=F)%*% soil_resource.pca$rotation[,1:5] 
climate_geology.p.tra = cbind(climate, ELEVATION, land_cover_and_use)
climate_geology.p.tra = climate_geology.p.tra[index.train,]
climate_geology.p.tra =cbind(climate_geology.p.tra, a)

soil_resource.tes = soil_resource[index.test,]
c = scale(soil_resource.tes, center=T, scale=F)%*% soil_resource.pca$rotation[,1:5]
climate_geology.p.tes = cbind(climate, ELEVATION, land_cover_and_use)
climate_geology.p.tes = climate_geology.p.tra[index.test,]
climate_geology.p.tes =cbind(climate_geology.p.tes, c)


livestock_cli = lm(cbind(AD05_CATT,  AD05_CHIC, AD05_GOAT, AD05_PIG, AD05_SHEE)~
                     LGP_AVG+LGP_CV+pre_cv+pre_mean+ELEVATION+AREA_CR_LC+PC1+PC2+PC3+PC4+PC5,
                   data = cbind(agriculture.c.tra, climate_geology.p.tra))
summary(livestock_cli)
b = lm(AD05_SHEE~LGP_AVG+LGP_CV+pre_cv+pre_mean+ELEVATION+AREA_CR_LC+PC1+PC2+PC3+PC4+PC5,
      data = cbind(agriculture.c.tra, climate_geology.p.tra))
summary(b)

plot(agriculture.c.tra$AD05_SHEE[1:300], b$fitted.values[1:300],xlab="real sheep density", ylab="fitted values")
plot(agriculture.c.tra$AD05_SHEE[1:300], b$residuals[1:300],xlab="real sheep density", ylab="residuals")
cor(agriculture.c.tra$AD05_SHEE[1:300], b$fitted.values[1:300])
sum(b$fitted.values<0) / length(b$fitted.values)

sheep.disc = agriculture.c.tra$AD05_SHEE
sheep.disc[sheep.disc>0] = 1
data=cbind(sheep.disc, climate_geology.p.tra)
sheep.glm = glm(sheep.disc~LGP_AVG*LGP_CV*pre_cv*pre_mean+ELEVATION+AREA_CR_LC+PC1+PC2+PC3+PC4+PC5, data=data)
summary(sheep.glm)

v = sheep.glm$fitted.values
v[v>0.5] = 1
v[v<=0.5] = 0
sum(v==sheep.disc)/length(sheep.disc)

sheep.disc = agriculture.c.tra$AD05_SHEE
data = cbind(sheep.disc, climate_geology.p.tra)
data = data[sheep.disc>0,]
g = lm(log(sheep.disc)~LGP_AVG+LGP_CV+pre_cv+pre_mean+ELEVATION+AREA_CR_LC+PC1+PC2+PC3+PC4+PC5,
       data = data)
summary(g)
h = sample(nrow(data), 300)
plot(log(data[h, "sheep.disc"]), g$fitted.values[h], xlab="real sheep density", ylab= "fitted values")
points(c(-20,20), c(-20,20), type="l", col = 2)
plot(log(data[h, "sheep.disc"]), g$residuals[h])


g = lm(log(AD05_SHEE)~LGP_AVG+LGP_CV+pre_cv+pre_mean+ELEVATION+AREA_CR_LC,
       data = data)
summary(g)


######################################




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
livestock.tra.mat = log(as.matrix(livestock[index.train,])+1)
livestock.tes.mat = log(as.matrix(livestock[index.test,])+1)
climate_geology.tra.mat = as.matrix(climate_geology.c.tra)
climate_geology.tes.mat = as.matrix(climate_geology.c.tes)
live_cg.tls = tls(climate_geology.tra.mat,livestock.tra.mat,kappa=10,lambda=1)

live_cg.tes.pred = climate_geology.tes.mat %*% live_cg.tls$b
plot(livestock.tes.mat[,1], live_cg.tes.pred[,1])

live_cg.tra.pred = climate_geology.tra.mat %*% live_cg.tls$b
plot(livestock.tra.mat[,1], live_cg.tra.pred[,1])

index = sample(nrow(livestock.tes.mat), 300)
test.true = apply(livestock.tes.mat[index,], 1, sum)
test.fit = apply(live_cg.tes.pred[index,], 1, sum)
plot(test.true, test.fit,  xlab="actual values", ylab="fitted values")
points(c(-20,20), c(-20,20), type="l", col = 2)
cor(test.true, test.fit)
plot(test.true, test.true-test.fit, xlab="actural values", yalb="residuals")
points(c(-20,20), c(0,0), type="l", col = 2)
plot(density(test.true-test.fit), xlab="Residual", main="Density Function")
qqnorm(test.true-test.fit)
qqline(test.true-test.fit, col = 1)
###########################pca
library(MASS)
climate_geology.tra.pca = prcomp(climate_geology.tra.mat, cor = F)
a=prcomp(soil_resource[index.train,])


#########################

TT_20K = social_factor.c.tra[,3]
TT_20K.factor = TT_20K %/% 5
TT_20K.factor[TT_20K.factor>6] = 7

vp_ar.tra = as.matrix(value_of_production[index.train,"vp_food_ar"]+
                             value_of_production[index.train,"vp_nonf_ar"])
vp_ar.factor = vp_ar.tra %/% 300
vp_ar.factor[vp_ar.factor>9] = 10 

sort(unique(vp_ar.factor))
plot(density(vp_ar.tra))
quantile(vp_ar.factor, 0.99)

a = table(vp_ar.factor, TT_20K.factor)
a = matrix(a, nrow = 11)
corresp(a, nf = 2)
biplot(co)

cor(climate[index.train,])
######################################################
length(unique(MSH_50K_TX$MSH_50K_TX[index.train]))
length(unique(portshed$PSH_TX[index.train]))

market.tra = MSH_50K_TX$MSH_50K_TX[index.train]
port.tra = portshed$PSH_TX[index.train]
TT_20K.factor
vp_ar.factor

table(transportation_time = TT_20K.factor, production_value_per_area= vp_ar.factor)
library(MASS)
a = table(vp_ar.factor, TT_20K.factor)
a = matrix(a, nrow = 11)
chisq.test(a) 
a.cor = corresp(a, nf=2)
biplot(a.cor)
abline(v = 0, h = 0, lty = 3)


b = table(vp_ar.factor, market.tra)
temp = b
b = matrix(b, nrow = 11)
colnames(b) = colnames(temp)
rownames(b) = rownames(temp)
temp = colSums(b)
b = b[, temp != 0]
temp =apply(b ,1, sum)
chisq.test(b) 
b.cor = corresp(b, nf=2)
biplot(b.cor)


c = table(vp_ar.factor, port.tra)
temp = c
c = matrix(c, nrow = 11)
colnames(c) = colnames(temp)
rownames(c) = rownames(temp)
temp =apply(c ,2, sum)
c = c[, temp != 0]
chisq.test(c) 
c.cor = corresp(c, nf=2)
biplot(c.cor)
abline(v = 0, h = 0, lty = 3)


################aov
market = factor(market.tra)
port = factor(port.tra)
tran.time = factor(TT_20K.factor)
prod.value = as.matrix(value_of_production[index.train,"vp_food_ar"]+
                         value_of_production[index.train,"vp_nonf_ar"])
q = aov(prod.value~market+port*tran.time)
summary(q)
plot(density(q$residuals), xlab="residual", xlim=c(-1000, 1000), main="")
qqnorm(q$residuals)
qqline(q$residuals)


#########################
library(randomForest)
names = unique(ecological_zones[index.train,])
ecol = rep(0, length(names))
for(i in 1:length(names))
{
  ecol[ecological_zones[index.train,] == names[i]] = i
}

ecol = factor(ecol)
cli = as.data.frame(climate_geology.p.tra) 
data = cbind(ecol, cli)

temp = data[1:5000,]
r = randomForest(ecol~., data=temp, mtry=6, ntree = 2000, na.action = na.omit)

climate_geology.p.tes

#######################
names.test = unique(ecological_zones[index.test,])
ecol.test = rep(0, length(names.test))
for(i in 1:length(names.test))
{
  ecol.test[ecological_zones[index.test,] == names.test[i]] = i
}

sum(r$y ==round(r$predicted,0)) / length(r$y)



