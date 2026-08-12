library(robumeta)

# ----------------------------------------------------------
# Moderator variables
# ----------------------------------------------------------

data$Intervention_format <- factor(
  data$MP,
  levels = c(0, 1),
  labels = c("Standalone MV", "MV-embedded")
)

data$MD_status <- factor(
  data$MD,
  levels = c(0, 1),
  labels = c("TD", "MD")
)

data$Proximality <- factor(
  data$Proximal,
  levels = c(0, 1),
  labels = c("Distal", "Proximal")
)

# Fundamental is the reference group
data$Content_type <- factor(
  data$Content,
  levels = c(1, 2),
  labels = c("Fundamental", "Higher-order")
)

data$Study_design <- factor(
  ifelse(data$Design == 3, 1, 0),
  levels = c(0, 1),
  labels = c("Group design", "Single-case design")
)

data$Publication_type <- factor(
  data$Publication,
  levels = c(0, 1),
  labels = c("Non-peer-reviewed", "Peer-reviewed")
)


# ----------------------------------------------------------
# Prepare data
# Main analyses: Ctl = 1, 3, 4
# Ctl = 2: component-added comparisons
# ----------------------------------------------------------

main_data <- subset(data, Ctl %in% c(1, 3, 4))
ctl2_data <- subset(data, Ctl == 2)

MV <- subset(main_data, Outcome == 0)
MP <- subset(main_data, Outcome == 1)


# ----------------------------------------------------------
# RVE function
# ----------------------------------------------------------

RVE <- function(dat, moderator = NULL, rho = .80) {
  
  formula <- if (is.null(moderator)) {
    ES_g ~ 1
  } else {
    as.formula(paste("ES_g ~", moderator))
  }
  
  robu(
    formula,
    data = dat,
    studynum = StudyID,
    var.eff.size = v_g,
    modelweights = "CORR",
    rho = rho,
    small = TRUE
  )
}


# ----------------------------------------------------------
# Pooled effects
# ----------------------------------------------------------

MV_overall <- RVE(MV)
MP_overall <- RVE(MP)

MV_overall
MP_overall


# ----------------------------------------------------------
# Moderator analyses
# One moderator at a time
# ----------------------------------------------------------

MV_mod <- list(
  Intervention_format = RVE(MV, "Intervention_format"),
  Proximality         = RVE(MV, "Proximality"),
  MD_status           = RVE(MV, "MD_status"),
  Age                 = RVE(MV, "Age"),
  Content_type        = RVE(MV, "Content_type"),
  MV_dosage           = RVE(MV, "MVDosage"),
  MV_percentage       = RVE(MV, "`MV perccentage`"),
  Study_design        = RVE(MV, "Study_design"),
  Publication_type    = RVE(MV, "Publication_type")
)

MP_mod <- list(
  Intervention_format = RVE(MP, "Intervention_format"),
  Proximality         = RVE(MP, "Proximality"),
  MD_status           = RVE(MP, "MD_status"),
  Age                 = RVE(MP, "Age"),
  Content_type        = RVE(MP, "Content_type"),
  MV_dosage           = RVE(MP, "MVDosage"),
  MV_percentage       = RVE(MP, "`MV perccentage`"),
  MP_dosage           = RVE(MP, "MPdose"),
  Study_design        = RVE(MP, "Study_design"),
  Publication_type    = RVE(MP, "Publication_type")
)


# ----------------------------------------------------------
# View key moderator results
# ----------------------------------------------------------

MV_mod$MD_status
MV_mod$Age
MV_mod$Content_type

MP_mod$MD_status
MP_mod$Age
MP_mod$Content_type


# ----------------------------------------------------------
# Subgroup pooled effects
# ----------------------------------------------------------

# MD status
RVE(subset(MV, MD_status == "TD"))
RVE(subset(MV, MD_status == "MD"))

RVE(subset(MP, MD_status == "TD"))
RVE(subset(MP, MD_status == "MD"))

# Mathematics content
RVE(subset(MV, Content_type == "Fundamental"))
RVE(subset(MV, Content_type == "Higher-order"))

RVE(subset(MP, Content_type == "Fundamental"))
RVE(subset(MP, Content_type == "Higher-order"))


# ----------------------------------------------------------
# Ctl = 2: Component-added effects
# ----------------------------------------------------------

MV_added <- RVE(subset(ctl2_data, Outcome == 0))
MP_added <- RVE(subset(ctl2_data, Outcome == 1))

MV_added
MP_added


# ----------------------------------------------------------
# Sensitivity analyses
# ----------------------------------------------------------

# Remove single-case designs
RVE(subset(MV, Design != 3))
RVE(subset(MP, Design != 3))

# Sensitivity to rho
for (rho in c(.60, .70, .80, .90)) {
  
  cat("\nMV: rho =", rho, "\n")
  print(RVE(MV, rho = rho)$reg_table)
  
  cat("\nMP: rho =", rho, "\n")
  print(RVE(MP, rho = rho)$reg_table)
}