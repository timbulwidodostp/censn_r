# Olah Data Semarang
# WhatsApp : +6285227746673
# IG : @olahdatasemarang_
# Censored skew-normal regression approach with delayed entry Use censn With (In) R Software
install.packages("knitr")
install.packages("haven")
install.packages("sn")
install.packages("curl")
library(knitr)
library(haven)
library("sn")
library("curl")
censn_ = read.csv("https://raw.githubusercontent.com/timbulwidodostp/censn_r/main/censn/censn.csv",sep = ";")
source(paste("https://raw.githubusercontent.com/timbulwidodostp/censn_r/main/censn/censn_r.r",sep = ";"))
# Estimation Censored skew-normal regression approach with delayed entry Use censn With (In) R Software
### No censoring, no delayed-entry
censn_1 <- censn(age1~1, ltrun=NULL, data=censn_, weights = rep(5, nrow(censn_)))
summary.censn(censn_1)
### Censoring, but no delayed-entry
censn_2 <- censn(age1~1, failure=died, ltrun=NULL, data=censn_, opt.method="BFGS")
summary.censn(censn_2)
### Censoring and delayed entry
censn_3 <- censn(age1~1, failure=died, ltrun=age, data=censn_)
summary.censn(censn_3)
### Censoring and delayed entry, covariate drug
censn_4 <- censn(age1~factor(drug), failure=died, ltrun=age, data=censn_)
summary.censn(censn_4)
# Censored skew-normal regression approach with delayed entry Use censn With (In) R Software
# Olah Data Semarang
# WhatsApp : +6285227746673
# IG : @olahdatasemarang_
# Finished