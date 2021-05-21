################################################################################
# This code calculates the halflife of neutralising antibodies observed        #
# between May and November, considering only subjects who were positive in May #
# and not increasing between May and November.                                 #
#                                                                              #
# Cite as:                                                                     #
# Dorigatti I et al, SARS-CoV-2 antibody dynamics, within-household            #
# transmission and the impact of contact tracing from community-wide           #
# serological testing in the Italian municipality of Vo'.                      #
# Pre-print available at: http://dx.doi.org/10.2139/ssrn.3792892               #
#                                                                              #
# Description:                                                                 #
# See README.md                                                                #
################################################################################

cat("\n\n### Running 'scripts/Neutralization_decay_script.R'\n\n")

# source scripts
source("R/functions_association_antibody_slopes.R")


################################################################################
# functions to transform categorical neutralization values into continuous
################################################################################

median_FORW.fun = function(d, i, rep, burn)
{
  dd = d[i,]
  hl.sim = vector("numeric", length=rep)
  for(idx in 1:rep)
  {
    d3 = dd
    d3[d3 == 1] = sample(1:9, sum(d3==1), replace=TRUE, prob=dgeom(1:9, 0.5))
    d3[d3 == 10] = sample(10:19, sum(d3==10), replace=TRUE, prob=dgeom(10:19-10, 0.5))
    d3[d3 == 20] = sample(20:39, sum(d3==20), replace=TRUE, prob=dgeom(20:39-20, 0.25))
    d3[d3 == 40] = sample(40:79, sum(d3==40), replace=TRUE, prob=dgeom(40:79-40, 0.125))
    d3[d3 == 80] = sample(80:159, sum(d3==80), replace=TRUE, prob=dgeom(80:159-80, 0.0625))
    d3[d3 == 160] = sample(160:319, sum(d3==160), replace=TRUE, prob=dgeom(160:319-160, 0.03125))
    d3[d3 == 320] = sample(320:639, sum(d3==320), replace=TRUE, prob=dgeom(320:639-320, 0.015625))
    d3[d3 == 640] = sample(640:1279, sum(d3==640), replace=TRUE, prob=dgeom(640:1279-640, 0.0078125))
    d3[d3 == 1280] = sample(1280:2559, sum(d3==1280), replace=TRUE, prob=dgeom(1280:2559-1280, 0.00390625))
    hl.sim[idx] = median(212/log(d3[,1]/d3[,2],2))
  }
  median(hl.sim[(burn+1):rep])
}

sample_FORW.fun = function(x)
{
  test = x/10
  val = 0

  if(test != 0)
    if(test < 1) val = sample(x:9, 1, prob=dgeom(x:9, 0.5)) else val = sample(x:(2*x-1), 1, prob=dgeom(x:(2*x-1)-x, 0.5/(2^log(test,2))))
  val
}

################################################################################
# end of function
################################################################################

# # read data
ds <- read.csv("data/Vo_serology_data.csv", header = TRUE)


# select baseline groud truth
ds <- ds[which(ds$Groundtruth == 1), ]

# keep only subjects positive in May and not increasing between May and November

Neu_1 <- as.character(ds$Neutralization)
Neu_1[Neu_1=="NA"] = NA
Neu_1 = as.factor(Neu_1)
Neu_1 = factor(Neu_1, levels = levels(Neu_1), labels = paste(c(1, 10, 1280, 160, 20, 320, 40, 640, 80, 0)))


Neu_2 <- as.character(ds$Neutralization_november_2020)
Neu_2[Neu_2=="NA"] = NA
Neu_2 = as.factor(Neu_2)
Neu_2 = factor(Neu_2, levels = levels(Neu_2), labels = paste(c(1, 320, 10, 160, 20, 320, 40, 80)))

Neu <- data.frame(Neu_1, Neu_2)
Neu = Neu[!apply(Neu, 1, function(x) any(is.na(x))),]

Neu[,1] = as.numeric(as.character(Neu[,1]))
Neu[,2] = as.numeric(as.character(Neu[,2]))

Neu_pos_May_notIncreasing = subset(Neu, Neu$Neu_1 >40 & Neu$Neu_2 <= Neu$Neu_1)

# number of days between serosurveys
days_surveys <- as.numeric(as.Date("11/29/2020", "%m/%d/%Y") -
                             as.Date("05/01/2020", "%m/%d/%Y"))

# calculate halflife and CI

dt = Neu_pos_May_notIncreasing
hl_db <- vector()

HL_N.sim = vector("numeric", length=4999)
for(idx in 1:4999)
{
  dt = data.frame(matrix(unlist(lapply(as.matrix(Neu_pos_May_notIncreasing), sample_FORW.fun)), ncol=ncol(Neu_pos_May_notIncreasing), nrow=nrow(Neu_pos_May_notIncreasing)))
  HL_N.sim[idx] = median(days_surveys/log(dt[,1]/dt[,2],2))
}

cat("Median = ", median(HL_N.sim), "\n")

cat("95% CI : ", quantile(HL_N.sim, c(0.025, 0.975)))

