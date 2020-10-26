# examples
# lm()
y <- cars$speed
x <- cars$dist
fit <- lm(y~x)
ib(fit)
fit_ib <- ib(fit,var=TRUE)
summary(fit_ib)

# glm(poisson())
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
data.frame(treatment, outcome, counts) # showing data
glm_fit <- glm(counts ~ outcome + treatment, family = poisson())
fit_ib <- ib(glm_fit, control=list(H=100,verbose=F))
summary(fit_ib)

glm_fit <- MASS::glm.nb(counts ~ outcome + treatment)
fit_ib <- ib(glm_fit, control=list(H=100,verbose=T))


# glm(Gamma())
clotting <- data.frame(
  u = c(5,10,15,20,30,40,60,80,100),
  lot1 = c(118,58,42,35,27,25,21,19,18),
  lot2 = c(69,35,26,21,18,16,13,12,12))
fit_gamma <- glm(lot2 ~ log(u), data = clotting, family = Gamma(link = "inverse"))
summary(fit_gamma)
fit_ib <- ib(fit_gamma,control=list(H=1e3,verbose=T,tol=1e-8))
summary(fit_ib)

fit_ib2 <- ib(fit_gamma,control=list(H=1e3,verbose=T,tol=1e-8),shape=TRUE)
summary(fit_ib2)

# negative binomial
require(MASS)
quine.nb1 <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)
fit_ib <- ib(quine.nb1,control=list(H=1e3,verbose=T))
fit_ib <- ib(quine.nb1,control=list(H=1e2,verbose=T),overdispersion = T)
summary(fit_ib)

# lmer
require(lme4)
# data("sleepstudy")
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
summary(fm1)# (with its own print method; see class?merMod % ./merMod-class.Rd
fit_ib1 <- ib(fm1, control = list(H=5,verbose=T))
summary(fit_ib)

fm2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)
fit_ib2 <- ib(fm2, control = list(H=5,verbose=T), Sigma=TRUE)

fm3 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML=FALSE)
fit_ib3 <- ib(fm3, control = list(H=5,verbose=T))

fm4 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML=FALSE)
fit_ib4 <- ib(fm4, control = list(H=100, verbose=T, maxit=100), Sigma=TRUE)


cnms <- getME(fm2,"cnms")
nc <- lengths(cnms)
theta <- getME(fm2,"theta")
sigma <- sigma(fm2)
fl <- getME(fm2,"flist")
nms <- names(fl)[attr(fl, "assign")]
vc <- mkVarCorr(sigma, cnms , nc , theta, nms)
