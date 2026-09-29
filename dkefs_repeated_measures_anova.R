# =====================================================================
# D-KEFS Color-Word Interference Test -- one-way repeated-measures ANOVA
# PSY 410 Experiment I
#
# REBUILT analysis script. This is NOT the original Sept 17 file --
# it was reconstructed to reproduce the verified results below.
# Run it yourself (base R only, no extra packages needed) and confirm
# the printed values match the verified results at the bottom.
#
# Usage: set your working directory to the folder containing
#        Fall_26_Experiment_1_Data_Set_Sheet1.csv, then source() this
#        file or run it line by line.
# =====================================================================

# ---- 1. Load data -----------------------------------------------------
data_file <- "Fall_26_Experiment_1_Data_Set_Sheet1.csv"

raw <- read.csv(data_file, check.names = FALSE)
names(raw) <- c("pid", "color_naming", "word_reading",
                "inhibition", "inhib_switch")

# Keep only participants with complete data on all four conditions
df <- raw[complete.cases(raw[, c("color_naming", "word_reading",
                                 "inhibition", "inhib_switch")]), ]
n <- nrow(df)
cat("Complete participants:", n, "\n\n")

conds      <- c("word_reading", "color_naming", "inhibition", "inhib_switch")
cond_names <- c("Word Reading", "Color Naming",
                "Inhibition", "Inhibition/Switching")
k <- length(conds)

# ---- 2. Descriptive statistics ----------------------------------------
cat("--- Condition means (seconds) ---\n")
means <- sapply(df[conds], mean)
names(means) <- cond_names
print(round(means, 2))
cat("\n--- Condition SDs ---\n")
print(round(sapply(df[conds], sd), 2))
cat("\n")

# ---- 3. Reshape to long format ----------------------------------------
long <- data.frame(
  pid       = rep(df$pid, times = k),
  condition = factor(rep(conds, each = n), levels = conds,
                     labels = cond_names),
  score     = unlist(df[conds], use.names = FALSE)
)

# ---- 4. Repeated-measures ANOVA ---------------------------------------
fit <- aov(score ~ condition + Error(pid/condition), data = long)
s   <- summary(fit)
print(s)

within  <- s[["Error: pid:condition"]][[1]]
ss_cond <- within["condition", "Sum Sq"]
df_cond <- within["condition", "Df"]
ss_err  <- within["Residuals", "Sum Sq"]
df_err  <- within["Residuals", "Df"]

F_val <- (ss_cond / df_cond) / (ss_err / df_err)
p_val <- pf(F_val, df_cond, df_err, lower.tail = FALSE)
eta2p <- ss_cond / (ss_cond + ss_err)

cat(sprintf(paste0("\nRepeated-measures ANOVA: F(%d,%d) = %.2f, ",
                   "p = %.2e, partial eta-squared = %.2f\n"),
            df_cond, df_err, F_val, p_val, eta2p))

# ---- 5. Sphericity ----------------------------------------------------
# Covariance matrix of the four conditions
S <- cov(as.matrix(df[conds]))

# Mauchly's W, computed on orthonormal contrasts of the conditions
C <- contr.poly(k)                       # k x (k-1), orthonormal columns
Sc <- t(C) %*% S %*% C
mauchly_W <- det(Sc) / (sum(diag(Sc)) / (k - 1))^(k - 1)
rho <- (n - 1) - (2 * (k - 1)^2 + (k - 1) + 2) / (6 * (k - 1))
chisq_m <- -rho * log(mauchly_W)
df_m <- k * (k - 1) / 2 - 1
p_mauchly <- pchisq(chisq_m, df_m, lower.tail = FALSE)
cat(sprintf(paste0("Mauchly's test of sphericity: W = %.4f, ",
                   "chi2(%d) = %.2f, p = %.3g\n"),
            mauchly_W, df_m, chisq_m, p_mauchly))

# Greenhouse-Geisser epsilon: computed on the covariance matrix of the
# orthonormal contrasts (this matches SPSS / afex / ezANOVA)
gg_eps <- sum(diag(Sc))^2 / ((k - 1) * sum(Sc^2))
df1_gg <- gg_eps * df_cond
df2_gg <- gg_eps * df_err
p_gg <- pf(F_val, df1_gg, df2_gg, lower.tail = FALSE)
cat(sprintf(paste0("Greenhouse-Geisser epsilon = %.4f; ",
                   "corrected F(%.2f,%.2f) = %.2f, p = %.2e\n"),
            gg_eps, df1_gg, df2_gg, F_val, p_gg))

# ---- 6. Post-hoc pairwise comparisons (Bonferroni) ---------------------
cat("\n--- Bonferroni post-hoc pairwise t-tests (paired) ---\n")
print(pairwise.t.test(long$score, long$condition,
                      paired = TRUE, p.adjust.method = "bonferroni"))

# ---- 7. Self-check against the verified results ------------------------
# Verified 2026-09-23: F(3,42)=46.98, p=1.77e-13, partial eta^2=.77,
# Mauchly's W=.1685 (p=.000393), GG epsilon=.5255 (corrected p<.001),
# means 21.64 / 31.27 / 52.80 / 59.12 s.
checks <- c(
  n == 15,
  abs(F_val - 46.98) < 0.05,
  abs(eta2p - 0.77) < 0.01,
  abs(mauchly_W - 0.1685) < 0.002,
  abs(gg_eps - 0.5255) < 0.002,
  all(abs(unname(means) - c(21.64, 31.27, 52.80, 59.12)) < 0.02)
)
if (all(checks)) {
  cat("\nAll checks passed: output matches the verified results.\n")
} else {
  warning("Some values differ from the verified results -- inspect the output above.")
}
