library(readxl)
library(metafor)

# Read coding sheet
data <- as.data.frame(
  read_excel("Coding sheet.xlsx"),
  check.names = FALSE
)

names(data) <- trimws(gsub("[\r\n]+", "", names(data)))

# Assumed pre-post correlation
corr <- 0.50


# ----------------------------------------------------------
# Post = 1: IGRM
# ----------------------------------------------------------

# Pooled pretest SD across treatment and control groups
data$s2Pre.Pooled <- with(
  data,
  ((pre.sd.tx^2) * (pre.n.tx - 1) +
     (pre.sd.ctl^2) * (pre.n.ctl - 1)) /
    (pre.n.tx + pre.n.ctl - 2)
)

data$sPre.Pooled <- sqrt(data$s2Pre.Pooled)


# Control-group standardized change
data <- escalc(
  measure = "SMCR",
  m1i = post.m.ctl,
  m2i = pre.m.ctl,
  sd1i = sPre.Pooled,
  ni = pre.n.ctl,
  ri = corr,
  data = data,
  correct = FALSE,
  replace = FALSE,
  var.names = c("delta_IG_ctl", "v_IG_ctl")
)


# Treatment-group standardized change
data <- escalc(
  measure = "SMCR",
  m1i = post.m.tx,
  m2i = pre.m.tx,
  sd1i = sPre.Pooled,
  ni = pre.n.tx,
  ri = corr,
  data = data,
  correct = FALSE,
  replace = FALSE,
  var.names = c("delta_IG_tx", "v_IG_tx")
)


# Difference in improvement between treatment and control
data$delta_IG <- data$delta_IG_tx - data$delta_IG_ctl
data$v_IG <- data$v_IG_tx + data$v_IG_ctl


# ----------------------------------------------------------
# Post = 0: Posttest SMD
# ----------------------------------------------------------

post_ES <- escalc(
  measure = "SMD",
  m1i = post.m.tx,
  sd1i = post.sd.tx,
  n1i = post.n.tx,
  m2i = post.m.ctl,
  sd2i = post.sd.ctl,
  n2i = post.n.ctl,
  data = data,
  correct = FALSE
)


# ----------------------------------------------------------
# Select effect-size metric
# Post = 1: IGRM
# Post = 0: Posttest SMD
# ----------------------------------------------------------

data$ES_d <- NA_real_
data$v_d <- NA_real_

data$ES_d[data$Post == 1] <- data$delta_IG[data$Post == 1]
data$v_d[data$Post == 1]  <- data$v_IG[data$Post == 1]

data$ES_d[data$Post == 0] <- post_ES$yi[data$Post == 0]
data$v_d[data$Post == 0]  <- post_ES$vi[data$Post == 0]


# ----------------------------------------------------------
# Convert to Hedges' g
# ----------------------------------------------------------

data$df <- data$post.n.tx + data$post.n.ctl - 2
data$J <- 1 - 3 / (4 * data$df - 1)

data$ES_g <- data$J * data$ES_d
data$v_g <- data$J^2 * data$v_d


# Study identifier for RVE
data$StudyID_clean <- trimws(as.character(data$StudyID))