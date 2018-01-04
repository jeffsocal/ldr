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
# DESC:     linearPrecision
##################################################################################


# linearPrecision
#
# @param data
# @param col_dil
# @param col_res
# @param col_uid
# @keywords linearPrecision
# @export
# @examples
# linearPrecision

linearPrecision <- function(data,
                           col_dil = 'dilution_conc',
                           col_res = 'dilution_resp',
                           col_uid = 'uid'){

  dat <- data
  dat_cv <- ddply(dat, c(col_uid, col_dil, col_dil), summarise,
                  resp_count     = length(dilution_resp),
                  resp_mean  = mean(dilution_resp),
                  resp_sd    = sd(dilution_resp),
                  resp_cv    = sd(dilution_resp) / mean(dilution_resp))

  d_out <- list()
  d_out$dat <- dat
  d_out$dat_cv <- dat_cv

  return(d_out)
}
