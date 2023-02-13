
# -------------------------------------------------------------------------#
## ------------ BIAS ANALYSIS MULTINOMIAL ------------------------------- ##
## --------- Model and calculate sensitivity and specificity ------------ ##
# -------------------------------------------------------------------------#

# create some fake quintiles 
set.seed(1234)
df$neighb <- sample(1:5, 1846, replace=TRUE, prob=c(0.2, 0.2, 0.2, 0.2, 0.2))
df$neighb <- factor(df$neighb)

set.seed(123)
df$ind <- sample(1:5, 1846, replace=TRUE, prob=c(0.19, 0.2, 0.2, 0.2, 0.21))
df$ind <- factor(df$ind)

# make dummies, REF = 5
df$ind_1 <- ifelse(df$ind==1, 1, 0)
df$ind_2 <- ifelse(df$ind==2, 1, 0)
df$ind_3 <- ifelse(df$ind==3, 1, 0)
df$ind_4 <- ifelse(df$ind==4, 1, 0)

df$neighb_1 <- ifelse(df$neighb==1, 1, 0)
df$neighb_2 <- ifelse(df$neighb==2, 1, 0)
df$neighb_3 <- ifelse(df$neighb==3, 1, 0)
df$neighb_4 <- ifelse(df$neighb==4, 1, 0)
df$neighb_5 <- ifelse(df$neighb==5, 1, 0)


# -------------------------------------------------------------------------#
# First calculate crude sens and spec by hand for each quintile
# Se = true positives / true positives + false negatives
# Sp = true negatives / true negatives + false negatives 
# -------------------------------------------------------------------------#
# get 5*5 table 
tab <- table(df$neighb,df$ind) # ind is on the top

# values needed for sensitivity 
tq1 <- tab[1,1]
tq2 <- tab[2,2]
tq3 <- tab[3,3]
tq4 <- tab[4,4]
tq5 <- tab[5,5]
fnq1 <- tab[2,1] + tab[3,1] + tab[4,1] + tab[5,1]
fnq2 <- tab[1,2] + tab[3,2] + tab[4,2] + tab[5,2]
fnq3 <- tab[1,3] + tab[2,3] + tab[4,3] + tab[5,3]
fnq4 <- tab[1,4] + tab[2,4] + tab[3,4] + tab[5,4]
fnq5 <- tab[1,5] + tab[2,5] + tab[3,5] + tab[4,5]

# calculate sensitivity 
# se Q1
seq1 <- tq1 / (tq1+fnq1)
seq2 <- tq2 / (tq2+fnq2)
seq3 <- tq3 / (tq3+fnq3)
seq4 <- tq4 / (tq4+fnq4)
seq5 <- tq5 / (tq5+fnq5)

# values needed for specificity 
tab
tnq1 <- tab[2,2]+tab[2,3]+tab[2,4]+tab[2,5]+
        tab[3,2]+tab[3,3]+tab[3,4]+tab[3,5]+
        tab[4,2]+tab[4,3]+tab[4,4]+tab[4,5]+
        tab[5,2]+tab[5,3]+tab[5,4]+tab[5,5] 

tnq2 <- tab[1,1]+tab[1,3]+tab[1,4]+tab[1,5]+
        tab[3,1]+tab[3,3]+tab[3,4]+tab[3,5]+
        tab[4,1]+tab[4,3]+tab[4,4]+tab[4,5]+
        tab[5,1]+tab[5,3]+tab[5,4]+tab[5,5] 

tnq3 <- tq1+tq2+tq4+tq5+tab[1,2]+tab[1,4]+tab[1,5]+tab[2,1]+tab[2,4]+tab[2,5]+tab[4,1]+tab[4,2]+tab[4,5]+tab[5,1]+tab[5,2]+tab[5,4]
tnq4 <- tq1+tq2+tq3+tq5+tab[1,2]+tab[1,3]+tab[1,5]+tab[2,1]+tab[2,3]+tab[2,5]+tab[3,1]+tab[3,2]+tab[3,5]+tab[5,1]+tab[5,2]+tab[5,3]
tnq5 <- tq1+tq2+tq3+tq4+tab[1,2]+tab[1,3]+tab[1,4]+tab[2,1]+tab[2,3]+tab[2,4]+tab[3,1]+tab[3,2]+tab[3,4]+tab[4,1]+tab[4,2]+tab[4,3]

fpq1 <- tab[1,2]+tab[1,3]+tab[1,4]+tab[1,5]
fpq2 <- tab[2,1]+tab[2,3]+tab[2,4]+tab[2,5]
fpq3 <- tab[3,1]+tab[3,2]+tab[3,4]+tab[3,5]
fpq4 <- tab[4,1]+tab[4,2]+tab[4,3]+tab[4,5]
fpq5 <- tab[5,1]+tab[5,2]+tab[5,3]+tab[5,4]

# calculate specificity 
spq1 <- tnq1/(tnq1+fpq1)
spq2 <- tnq2/(tnq2+fpq2)

spq3 <- tnq3/(tnq3+fpq3)
spq4 <- tnq4/(tnq4+fpq4)
spq5 <- tnq5/(tnq5+fpq5)

# -------------------------------------------------------------------------#
# Use multinomial regression to calculate sensitivity and specificity 
# See banack 2018 for stata code
# outcome = measure exposure (neighbourhood income)
# exposure = true exposure (individual income)
# confounders = rural residence, possibly age and sex 
# -------------------------------------------------------------------------#

#### CALCULATE SENSITIVITY USING MODELS #### 
### use outcome = quintile of interest vs rest Q1 (i.e. use logistic regression)

### QUINTILE 1 ###
mq1 <- glm(neighb_1 ~ ind_1 + ind_2 + ind_3 + ind_4, data=df, family='binomial')
summary(mq1)
coef_mq1 <- coef(mq1)
coef_mq1
# sensivity for Q1 -- matches crude 
m_seq1 <- (1) / (1+(exp(-(coef_mq1[1]+1*(coef_mq1[2]))))) 
m_seq1
seq1
# specificity for Q1 -- doesn't match crude?
m_spq1 <- (1 - ((1) / (1+(exp(-(coef_mq1[1]+0*(coef_mq1[2])+1*(coef_mq1[3])+1*(coef_mq1[4])+1*(coef_mq1[5])))))  ) )
m_spq1
spq1

### QUINTILE 2 ###
mq2 <- glm(neighb_2 ~ ind_1 + ind_2 + ind_3 + ind_4, data=df, family='binomial')
summary(mq2)
coef_mq2 <- coef(mq2)
# sensitivity q2 -- matches crude 
m_seq2 <- (1) / (1+(exp(-(coef_mq2[1]+1*(coef_mq2[3]))))) 
m_seq2
seq2
# specificity q2
m_spq2 <- 1 - ((1) / (1+(exp(-(coef_mq2[1]+0*(coef_mq2[3]))))) )
m_spq2
spq2

### QUINTILE 3 ###
mq3 <- glm(neighb_3 ~ ind_1 + ind_2 + ind_3 + ind_4, data=df, family='binomial')
summary(mq3)
coef_mq3 <- coef(mq3)
# sensitivity q3 -- matches crude 
m_seq3 <- (1) / (1+(exp(-(coef_mq3[1]+1*(coef_mq3[4]))))) 
m_seq3
seq3
# specificity q3
m_spq3 <- 1 -( (1) / (1+(exp(-(coef_mq3[1]+1*(coef_mq3[4]))))) )
m_spq3
spq3

### QUINTILE 4 ###
mq4 <- glm(neighb_4 ~ ind_1 + ind_2 + ind_3 + ind_4, data=df, family='binomial')
summary(mq4)
coef_mq4 <- coef(mq4)
# sensitivity q4 -- matches crude 
m_seq4 <- (1) / (1+(exp(-(coef_mq4[1]+1*(coef_mq4[5]))))) 
m_seq4
seq4
# specificity q4

### QUINTILE 5 ###
mq5 <- glm(neighb_5 ~ ind_1 + ind_2 + ind_3 + ind_4, data=df, family='binomial')
summary(mq5)
coef_mq5 <- coef(mq5)
# snesitivity q5 -- matches crude 
m_seq5 <- (1) / (1+(exp(-(coef_mq5[1])))) 
m_seq5
seq5
# specificity q5










