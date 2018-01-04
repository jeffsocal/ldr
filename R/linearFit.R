################################################################################
# Copyright (2016) SoCal Bioinformatics Inc. All rights reserved.
# This script is the confidential and proprietary product of SoCal
# Bioinformatics Inc. Any unauthorized reproduction or transfer of
# the contents herein is strictly prohibited.
#
################################################################################
# AUTH:     Jeff Jones | SoCal Bioinofrmatics Inc
# DATE:     2016.11.17
# OWNER:    Jeff Jones | SoCal Bioinofrmatics Inc
# PROJECT:  C. Hughey | JMU | NegativeElectrospray
# DESC:     linearFit
##################################################################################


#' linearFit
#'
#' @param data data.frame of dilution~response values
#' @param col_dil column containing the variable for dilution concentration
#' @param col_res column containing the variable for dilution response
#' @param col_uid column containing the variable for unique id
#' @param rjct_pre rejection criteria for min precision, default 0.25,
#' @param rjct_acc rejection criteria for min accuracy, default 0.25,
#' @param rjct_rsq rejection criteria for min linearity Pearson corrleation, default 0.90,
#' @param rjct_pts rejection criteria for dilution points retained, default 0.75,
#' @param filter filter dataset to top solution passing reject criteria, default TRUE
#' @keywords linearFit
#' @export
#' @examples
#' linearFit


linearFit <- function(data,
                      col_dil = 'dilution_conc',
                      col_res = 'dilution_resp',
                      col_uid = 'uid',
                      rjct_pre = 0.25,
                      rjct_acc = 0.25,
                      rjct_rsq = 0.90,
                      rjct_pts = 0.75,
                      filter = T){

  # sanitize data
  fat <- linearDataClean(data, col_dil, col_res)

  l_fat <- linearPrecision(fat)

  # collect unique elements in DF for iteration
  v_uid <- unique(l_fat$dat[,col_uid])
  l_uid <- length(v_uid)

  # est final DF
  dat_dil <- c()

  cat("Calculate [", l_uid, "] linear fits ... \n")
  pb <- txtProgressBar(max=l_uid, style=3, min=0)

  for ( i in 1:l_uid ){
    setTxtProgressBar(pb, i)

    dat_this <- c()

    # set raw DF & summary DF
    this_fat_a <- l_fat$dat[l_fat$dat$uid == v_uid[i],]

    # all possible dilution spans of 3 or greater dilution points
    ap <- t(combn(min(this_fat_a$dilution_n):max(this_fat_a$dilution_n),2))
    ap <- ap[which((ap[,2] - ap[,1]) > 2),]

    for( dr in 1:dim(ap)[1] ){

      fit_fail_cat <- ""
      this_fat <- this_fat_a[this_fat_a$dilution_n >= ap[dr,1] & this_fat_a$dilution_n <= ap[dr,2],]

      this_lm <- lm(log10(this_fat$dilution_resp) ~ log10(this_fat$dilution_conc))
      this_int <- as.numeric(this_lm$coeff[1])
      this_slope <- as.numeric(this_lm$coeff[2])
      this_rsq <- as.numeric(summary(this_lm)$r.sq)

      # residual calcs based on ordinal value
      this_fat$residuals_norm <- abs(this_fat$dilution_resp - 10^this_lm$fitted.values) / this_fat$dilution_resp

      this_cvs <- linearPrecision(this_fat)[['dat_cv']]

      this_res <- merge(this_fat, this_cvs, by=c(col_uid, col_dil))

      this_n_total      <- dim(this_res)[1]
      this_n_precision  <- dim(this_res[this_res$resp_cv <= rjct_pre,])[1]
      this_n_accuracy   <- dim(this_res[this_res$residuals_norm <= rjct_acc,])[1]

      fit_fail_cat <- ''
      if(is.na(this_rsq) | this_rsq < rjct_rsq)
        fit_fail_cat <- paste0('linearity: rqs < ', rjct_rsq, " ")

      if(this_n_accuracy / this_n_total < rjct_pts)
        fit_fail_cat <- paste0('accuracy: ', round((this_n_accuracy / this_n_total) * 100,1), '% ')

      if(this_n_precision / this_n_total < rjct_pts)
        fit_fail_cat <- paste0('precision: ', round((this_n_precision / this_n_total) * 100,1), '% ')

      this_res <- this_res[this_res$resp_cv <= rjct_pre & this_res$residuals_norm <= rjct_acc,]

      if(dim(this_res)[1] / dim(this_fat)[1] <= rjct_pts)
        fit_fail_cat <- paste0('prec and accu: ', round((this_n_precision / this_n_total) * 100,1), '% ')

      if(length(unique(this_res$dilution_n)) < 3 | dim(this_cvs)[1] == 0) {
        fit_fail_cat <- '3 or less conc'
      } else {
        # if(min(this_res$dilution_conc) != min(this_cvs$dilution_conc) |
        #    max(this_res$dilution_conc) != max(this_cvs$dilution_conc))
        #   fit_fail_cat <- "out of range"
      }

      dat_this <- rbind(dat_this, data.frame(
        unique_id = v_uid[i],
        fit_remain = dim(this_res)[1] / dim(this_fat)[1],
        fit_points = length(unique(this_res$dilution_n)),
        fit_span = ap[dr,2] - ap[dr,1] + 1,
        fit_slope = this_slope,
        fit_intercept = this_int,
        fit_rsquared = this_rsq,
        fit_cv_mean = mean(this_res$resp_cv, rm.na=T),
        fit_res_mean = mean(abs(this_fat$residuals_norm)),
        fit_dil_min = min(this_res$dilution_conc),
        fit_dil_max = max(this_res$dilution_conc),
        fit_lloq = round(min(this_res$resp_mean)),
        fit_uloq = round(max(this_res$resp_mean)),
        fit_fail = fit_fail_cat != '',
        fit_fail_cat = fit_fail_cat))
    }

    dat_this <- dat_this[order(dat_this$fit_fail, -round(dat_this$fit_rsquared), -dat_this$fit_points),]

    # Filter Data

    if( filter == T){
      dat_this_temp <- dat_this

      dat_this <- dat_this[dat_this$fit_rsquared >= rjct_rsq,]
      dat_this <- dat_this[dat_this$fit_res_mean <= rjct_acc,]
      dat_this <- dat_this[dat_this$fit_cv_mean <= rjct_pre,]
      dat_this <- dat_this[dat_this$fit_points >= dat_this$fit_span-1 ,]
      dat_this <- dat_this[order(dat_this$fit_fail, -dat_this$fit_points),]

      if(dim(dat_this)[1] == 0)
        dat_this <- dat_this_temp
    }

    dat_dil <- rbind(dat_dil, dat_this)
  }

  return(dat_dil)
}
