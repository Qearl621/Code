library(readxl)
library(metaSEM)

# ----------------------------------------------------------
# Read correlation matrices
# ----------------------------------------------------------

raw <- as.data.frame(
  read_excel("Martix.xlsx", col_names = FALSE),
  check.names = FALSE
)

var_names <- c("treat", "PostMV", "PostMP")


# Identify rows containing matrix headers
header_rows <- which(
  apply(raw, 1, function(x) {
    x <- gsub("\\s+", "", as.character(x))
    all(var_names %in% x)
  })
)


# Extract correlation matrices, StudyID, and N
mats <- lapply(header_rows, function(i) {
  
  mat <- apply(
    raw[(i + 1):(i + 3), 5:7],
    2,
    as.numeric
  )
  
  mat <- matrix(mat, nrow = 3)
  dimnames(mat) <- list(var_names, var_names)
  
  mat
})

StudyID <- as.numeric(raw[header_rows, 1][[1]])
MatrixID <- as.numeric(raw[header_rows, 2][[1]])
N <- as.numeric(raw[header_rows, 3][[1]])


# ----------------------------------------------------------
# Combine dependent matrices within studies
# StudyID 2, 5, and 6
# ----------------------------------------------------------

fisher_average <- function(mats, n) {
  
  out <- diag(1, 3)
  dimnames(out) <- list(var_names, var_names)
  
  for (i in 1:2) {
    for (j in (i + 1):3) {
      
      r <- sapply(mats, function(x) x[i, j])
      
      z <- weighted.mean(
        atanh(r),
        w = n - 3
      )
      
      out[i, j] <- out[j, i] <- tanh(z)
    }
  }
  
  out
}


study_mats <- list()
study_n <- numeric()

for (id in unique(StudyID)) {
  
  index <- which(StudyID == id)
  
  if (length(index) > 1) {
    
    study_mats[[length(study_mats) + 1]] <-
      fisher_average(mats[index], N[index])
    
    # Retain one study-level sample size
    study_n <- c(study_n, max(N[index]))
    
  } else {
    
    study_mats[[length(study_mats) + 1]] <- mats[[index]]
    study_n <- c(study_n, N[index])
  }
}

names(study_mats) <- unique(StudyID)
study_ids <- unique(StudyID)


# ----------------------------------------------------------
# Two-stage TSSEM
# ----------------------------------------------------------

run_TSSEM <- function(exclude_study1 = FALSE) {
  
  keep <- if (exclude_study1) {
    study_ids != 1
  } else {
    rep(TRUE, length(study_ids))
  }
  
  # Stage 1: pool correlation matrices
  stage1 <- tssem1(
    Cov = study_mats[keep],
    n = study_n[keep],
    method = "REM",
    RE.type = "Diag"
  )
  
  # Stage 2: mediation/path model
  model <- "
    PostMV ~ a*treat
    PostMP ~ cprime*treat + b*PostMV

    treat ~~ 1*treat
    PostMV ~~ eMV*PostMV
    PostMP ~~ eMP*PostMP
  "
  
  RAM <- lavaan2RAM(model)
  
  stage2 <- tssem2(
    stage1,
    RAM = RAM,
    intervals.type = "z"
  )
  
  list(
    stage1 = stage1,
    stage2 = stage2
  )
}


# ----------------------------------------------------------
# Main model
# ----------------------------------------------------------

model_main <- run_TSSEM()

summary(model_main$stage1)
summary(model_main$stage2)


# ----------------------------------------------------------
# Sensitivity analysis: exclude StudyID 1
# ----------------------------------------------------------

model_no1 <- run_TSSEM(exclude_study1 = TRUE)

summary(model_no1$stage1)
summary(model_no1$stage2)


# ----------------------------------------------------------
# Indirect and total effects
# ----------------------------------------------------------

indirect_effect <- function(model) {
  
  coef_table <- as.data.frame(
    summary(model$stage2)$coefficients
  )
  
  coef_table$parameter <- rownames(coef_table)
  
  a <- coef_table[coef_table$parameter == "a", ]
  b <- coef_table[coef_table$parameter == "b", ]
  cprime <- coef_table[coef_table$parameter == "cprime", ]
  
  indirect <- a$Estimate * b$Estimate
  
  indirect_SE <- sqrt(
    b$Estimate^2 * a$Std.Error^2 +
      a$Estimate^2 * b$Std.Error^2
  )
  
  z <- indirect / indirect_SE
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  
  total <- cprime$Estimate + indirect
  
  data.frame(
    effect = c(
      "Intervention -> MV",
      "MV -> MP",
      "Intervention -> MP",
      "Indirect effect",
      "Total effect"
    ),
    beta = c(
      a$Estimate,
      b$Estimate,
      cprime$Estimate,
      indirect,
      total
    ),
    SE = c(
      a$Std.Error,
      b$Std.Error,
      cprime$Std.Error,
      indirect_SE,
      NA
    ),
    p = c(
      a[["Pr(>|z|)"]],
      b[["Pr(>|z|)"]],
      cprime[["Pr(>|z|)"]],
      p,
      NA
    )
  )
}


main_paths <- indirect_effect(model_main)
sensitivity_paths <- indirect_effect(model_no1)

main_paths
sensitivity_paths