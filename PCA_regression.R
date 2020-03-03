#  Principal component regression analysis plots
lmmatrix<- function (y, cov){
  p <- matrix(NA, ncol = ncol(cov), nrow = ncol(y))
  for (j in 1:ncol(cov)) {
    x <- cov[, j]
    for (i in 1:ncol(y)) {
      fit <- summary(lm(y[, i] ~ x, na.action = na.omit))
      f <- fit$fstatistic
      p[i, j] <- pf(f["value"], f["numdf"], f["dendf"], lower.tail = FALSE)
    }
  }
  colnames(p) <- names(cov)
  return(p)
}

lmmatrix2<- function (y, cov){
  p <- matrix(NA, ncol = ncol(cov), nrow = ncol(y))
  #for (j in 1:ncol(cov)) {
  #x <- cov[, j]
  for (i in 1:ncol(y)) {
    fit <- anova(lm(y[, i] ~ ., cov, na.action = na.omit))
    p[i, ] <- fit$`Pr(>F)`[1:ncol(cov)]
  }
  #}
  colnames(p) <- names(cov)
  return(p)
}

plotp2 <- function (p, yaxis, xmax, title){
  plot(1, xlim = c(0, xmax), ylim = c(0, length(yaxis) + 1), type = "n", bty = "n", axes = FALSE, xlab = "Principal Component", ylab = "", main = title)
  axis(1, at = c(1:xmax), pos = 0.5, las = 1, lwd = 3)
  for (i in 1:length(yaxis)) {
    text(0.3, i, yaxis[i], xpd = TRUE, adj = 1)
  }
  for (i in 1:ncol(p)) {
    for (j in 1:nrow(p)) {
      pp <- p[j, i]
      colcode <- "white"
      if (pp <= 1e-09) {
        colcode = "darkred"
      }
      else if (pp <= 1e-04) {
        colcode = "red"
      }
      else if (pp <= 0.01) {
        colcode = "orange"
      }
      else if (pp <= 0.05) {
        colcode = "pink"
      }
      polygon(c(j - 0.5, j - 0.5, j + 0.5, j + 0.5), c(i - 0.5, i + 0.5, i + 0.5, i - 0.5), col = colcode, border = NA)
    }
  }
  legend("topright", c("<0.05", "<0.01", "<10E-5", "<10E-10"), col = c("pink", "orange", "red", "darkred"), pch = 15, pt.cex = 2, bty = "o", horiz = TRUE, xpd = TRUE)
}
