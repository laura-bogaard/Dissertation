## Read in data

## Format data

## Check for variable collinearity in linear model

## Adjust model spec

## 

#########
df <- padat
df$pa <- as.numeric(df$pa)
df$treat <- as.factor(df$treatment)
df$obs <- as.factor(df$obs)
df$x.pos <- df$x
df$y.pos <- df$y

# add numeric ID for each experimental unit (sesh id)
seshlevels <-  levels(as.factor(df$session_id)) 
numericsesh <- data.frame(session_id = levels(as.factor(df$session_id)),
                          block_as_session = as.factor(seq(1, length(seshlevels), 1)))

df <- left_join(df, numericsesh, by = "session_id" )
arrange(df, session_id)

str(df)
#######
#EDA
qplot(df$treat, df$pa)



# start with full model 
m0 <- gam(pa ~ treatment + s(x, y) + s(rl) + s(fish) + obs + s(dist) + day + hour, method = "REML", family = binomial)
m0
summary(m0)
gam.check(m0)
# does not converge

# Try smooths for rl
m1 <- gam(pa ~ treatment + s(x, y) + s(rl) + fish + obs + dist + day + hour, method = "REML", family = binomial)
m1
summary(m1)
gam.check(m1)
# not bad, kick out observer  
m1b <- gam(pa ~ treatment + s(x, y) + s(rl) + fish + dist + day + hour, method = "REML", family = binomial)
m1b
summary(m1b)
gam.check(m1b)
# hm worse

# not bad, kick out observer + dist 
m1c <- gam(pa ~ treatment + s(x, y) + s(rl) + fish + day + hour, method = "REML", family = binomial)
m1c
summary(m1c)
gam.check(m1c)
# hm worse

# not bad, kick out observer + dist 
m1d <- gam(pa ~ treatment + s(x, y) + s(rl) + fish + hour, method = "REML", family = binomial)
m1d
summary(m1d)
gam.check(m1d)
# hm worse


AIC(m1b, m1c, m1d)