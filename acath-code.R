library(Hmisc)
library(tidyverse)
library(rigr)
library(naniar)
library(gridExtra)
library(knitr)
library(Amelia)
library(mice)

getHdata(acath)
acath1 <- subset(acath, age>19)
acath1$high_choleste <- ifelse(acath1$choleste>238, 1, 0)
acath_sigdz <- subset(acath1, sigdz==1)

tab1 <- acath1 %>% descrip()
tab1 <- tab1[, 1:9]
tab1 <- data.frame(tab1[, c("N", "Msng", "Mean", "Std Dev", " Min", " Max")])
names(tab1) <- c("N", "Msng", "Mean/Prop", "Std Dev", "Min", "Max")
tab1 <- tab1 %>% mutate(`Mean/Prop` = as.character(round(`Mean/Prop`, 2)),
                        `Std Dev` = as.character(round(`Std Dev`, 2)),
                        Max = as.character(round(Max, 1)),
                        Min = as.character(round(Min, 1)))
rownames(tab1) <- c("Sex (Female)", "Age (years)", 
                    "Duration of Chest Pain (days)", "Cholesterol Level (mg/dl)",
                    "Has Significant Coronary Disease",
                    "Has Severe Coronary Disease", "Has High Cholesterol (>238 mg/dl)")
tab1$`Std Dev`[c(1,5:7)] <- "--"
tab1$Min[c(1,5:7)] <- "--"
tab1$Max[c(1,5:7)] <- "--"

tab2 <- acath_sigdz %>% select(-c(sigdz)) %>% descrip() 
tab2 <- tab2[, 1:9]
tab2 <- data.frame(tab2[, c("N", "Msng", "Mean", "Std Dev", " Min", " Max")])
names(tab2) <- c("N", "Msng", "Mean/Prop", "Std Dev", "Min", "Max")
tab2 <- tab2 %>% mutate(`Mean/Prop` = as.character(round(`Mean/Prop`, 2)),
                        `Std Dev` = as.character(round(`Std Dev`, 2)),
                        Max = as.character(round(Max, 1)),
                        Min = as.character(round(Min, 1)))
rownames(tab2) <- c("Sex (Female)", "Age (years)", 
                    "Duration of Chest Pain (days)", "Cholesterol Level (mg/dl)",
                    "Has Severe Coronary Disease", "Has High Cholesterol (>238 mg/dl)")
tab2$`Std Dev`[c(1,5:6)] <- "--"
tab2$Min[c(1,5:6)] <- "--"
tab2$Max[c(1,5:6)] <- "--"

# Covariate of interest transformation
ggplot(data=acath1, 
       aes(x=cad.dur)) +
  geom_histogram() +
  xlab("Duration of Chest Pain (Days)")

acath1$logcad.dur <- log(acath1$cad.dur+1)
acath_sigdz <- subset(acath1, sigdz==1)

ggplot(data=acath1, 
       aes(x=logcad.dur)) +
  geom_histogram() +
  xlab("Duration of Chest Pain (log-Days)")

# Boxplots
p1 <- ggplot(data=acath1, 
       aes(x=logcad.dur, y=factor(sigdz))) +
  geom_boxplot() +
  xlab("Duration of Chest Pain (log-Days)") + 
  ylab("Significant Coronary Disease")
p2 <- ggplot(acath_sigdz %>% filter(!is.na(tvdlm)), 
       aes(x=logcad.dur, y=as.factor(tvdlm))) +
  geom_boxplot() + 
  xlab("Duration of Chest Pain (log-Days)") + 
  ylab("Severe Coronary Disease")
grid.arrange(p1,p2)

p3 <- ggplot(data=na.omit(acath1), 
             aes(x=logcad.dur, y=as.factor(sigdz), fill=as.factor(high_choleste))) +
  geom_boxplot() + 
  facet_wrap(~high_choleste) +
  xlab("Duration of Chest Pain (log-Days)") + 
  ylab("Significant Coronary Disease")  
p4 <- ggplot(data=na.omit(acath_sigdz), 
             aes(x=logcad.dur, y=as.factor(tvdlm), fill=as.factor(high_choleste))) +
  geom_boxplot() + 
  facet_wrap(~high_choleste) +
  xlab("Duration of Chest Pain (log-Days)") + 
  ylab("Severe Coronary Disease") 
grid.arrange(p3,p4)

missmap(acath1 %>% select(-c(high_choleste, logcad.dur)), main = "Missing vs Observed Values")
md.pattern(acath1)

plot1 <- ggplot(
  data = acath1, aes(x=choleste, y=age)) +
  geom_miss_point() +
  xlab("Cholesterol Levels (mg/dl)") +
  ylab("Age (years)")
plot2 <- ggplot(
  data = acath1, aes(x=choleste, y=logcad.dur)) +
  geom_miss_point() +
  xlab("Cholesterol Levels (mg/dl)") +
  ylab("Duration of Chest Pain (log-Days)")
grid.arrange(plot1,plot2)

acath_na <- subset(acath1, is.na(acath1$choleste))
table(acath_na$sex)[]/nrow(acath_na)

# Stochastic regression imputation via mice
acath_pred <- subset(acath, age>19)
acath_imp <- mice(acath_pred, meth = 
                    c("norm.nob", "norm.nob", "norm.nob", "norm.nob", "norm.nob", "logreg"), m = 5)
# Store data
acath_imp <- complete(acath_imp)

acath_imp$age20 <- acath_imp$age-20
acath_imp$high_choleste <- ifelse(acath_imp$choleste>238, 1, 0)
acath_imp$logcad.dur <- log(acath_imp$cad.dur+1)
acath_sigdz_imp <- subset(acath_imp, sigdz==1)

r1 <- regress("odds", sigdz~logcad.dur+choleste+sex+age20, acath_imp)
r2 <- regress("odds", tvdlm~logcad.dur+choleste+sex+age20, acath_sigdz_imp)

log1 <- glm(sigdz~logcad.dur+choleste+sex+age20, family="binomial", acath_imp)
summary(log1)
log2 <- glm(tvdlm~logcad.dur+choleste+sex+age20, family="binomial", acath_sigdz_imp)
summary(log2)

log1_null <- glm(sigdz~choleste+sex+age20, family="binomial", acath_imp)
test1 <- anova(log1_null, log1, test="LRT")
log2_null <- glm(tvdlm~choleste+sex+age20, family="binomial", acath_sigdz_imp)
test2 <- anova(log2_null, log2, test="LRT")

r3 <- regress("odds", sigdz~logcad.dur*high_choleste+sex+age20, acath_imp)
r4 <- regress("odds", tvdlm~logcad.dur*high_choleste+sex+age20, acath_sigdz_imp)

log3 <- glm(sigdz~logcad.dur*high_choleste+sex+age20, family="binomial", acath_imp)
summary(log3)
log4 <- glm(tvdlm~logcad.dur*high_choleste+sex+age20, family="binomial", acath_sigdz_imp)
summary(log4)

test3 <- anova(log3, test="LRT")
test4 <- anova(log4, test="LRT")

tab3 <- r1$transformed
tab3 <- cbind(tab3[-1, 1:3], test1[2, 5])
tab3 <- as.data.frame(tab3)
names(tab3) <- c("Estimate", "Lower 95% CI", "Upper 95% CI", "LRT P-Value")
rownames(tab3) <- c("Duration of Chest Pain (log-days)", "Cholesterol Level (mg/dl)", 
                    "Sex (Female)", "Age (years over 20)")
tab3 <- tab3 %>% round(3)
tab3$`LRT P-Value`[2:4] <- "--"

tab4 <- r2$transformed
tab4 <- cbind(tab4[-1, 1:3], test2[2, 5])
tab4 <- as.data.frame(tab4)
names(tab4) <- c("Estimate", "Lower 95% CI", "Upper 95% CI", "LRT P-Value")
rownames(tab4) <- c("Duration of Chest Pain (log-days)", "Cholesterol Level (mg/dl)", 
                    "Sex (Female)", "Age (years over 20)")
tab4 <- tab4 %>% mutate(`Estimate` = round(`Estimate`, 3),
                        `Lower 95% CI` = round(`Lower 95% CI`, 3),
                        `Upper 95% CI` = round(`Upper 95% CI`, 3),
                        `LRT P-Value` = round(`LRT P-Value`, 21))
tab4$`LRT P-Value`[2:4] <- "--"

tab5 <- r3$transformed
tab5 <- cbind(tab5[-1, 1:3], test3[6, 5])
tab5 <- as.data.frame(tab5)
names(tab5) <- c("Estimate", "Lower 95% CI", "Upper 95% CI", "LRT P-Value")
rownames(tab5) <- c("Duration of Chest Pain (log-days)", "Cholesterol Level (High)", 
                    "Sex (Female)", "Age (years over 20)",
                    "Duration of Chest Pain (log-days):Cholesterol Level (High)")
tab5 <- tab5 %>% round(3)
tab5$`LRT P-Value`[1:4] <- "--"

tab6 <- r4$transformed
tab6 <- cbind(tab6[-1, 1:3], test4[6, 5])
tab6 <- as.data.frame(tab6)
names(tab6) <- c("Estimate", "Lower 95% CI", "Upper 95% CI", "LRT P-Value")
rownames(tab6) <- c("Duration of Chest Pain (log-days)", "Cholesterol Level (High)", 
                    "Sex (Female)", "Age (years over 20)",
                    "Duration of Chest Pain (log-days):Cholesterol Level (High)")
tab6 <- tab6 %>% round(3)
tab6$`LRT P-Value`[1:4] <- "--"

t <- data.frame(`Model 1` = c("*sigdz*", "*cad.dur*", "*sex=Female*", "*age-20*", "*choleste*", "--"),
                `Model 2` = c("$tvdlm|sigdz=1$", "*cad.dur*", "*sex=Female*", "*age-20*", "*choleste*", "--"),
                `Model 3` = c("*sigdz*", "*cad.dur*", "*sex=Female*", "*age-20*", 
                              "*choleste=high*", "$cad.dur \times choleste=high$"),
                `Model 4` = c("$tvdlm|sigdz=1$", "*cad.dur*", "*sex=Female*", "*age-20*", 
                              "*choleste=high*", "$cad.dur \times choleste=high$"))
names(t) <- c("**Model 1**", "**Model 2**", "**Model 3**", "**Model 4**")
rownames(t) <- c("**Response**", "**Covariates**", "1", "2", "3", "4")
knitr::kable(t, caption = "t")

library(table1)
label(acath1$sex) <- "Sex (Female)"
label(acath1$age) <- "Age (years)"
label(acath1$cad.dur) <- "Duration of Chest Pain (days)"
label(acath1$choleste) <- "Cholesterol Level (mg/dl)"
label(acath1$high_choleste) <- "Has High Cholesterol (>238 mg/dl)"
t1 <- table1(~sex + age + cad.dur + choleste + high_choleste| sigdz, data = acath1)

label(acath_sigdz$sex) <- "Sex (Female)"
label(acath_sigdz$age) <- "Age (years)"
label(acath_sigdz$cad.dur) <- "Duration of Chest Pain (days)"
label(acath_sigdz$choleste) <- "Cholesterol Level (mg/dl)"
label(acath_sigdz$high_choleste) <- "Has High Cholesterol (>238 mg/dl)"
t2 <- table1(~sex + age + cad.dur + choleste + high_choleste| tvdlm, data = acath_sigdz)
knitr::kable(t2, "pipe")

table1(~cad.dur+sigdz | high_choleste, data = acath1)
table1(~cad.dur+tvdlm | high_choleste, data = acath_sigdz)
