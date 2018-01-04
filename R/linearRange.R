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
# DESC:     linearRange
##################################################################################


# linearRange
#
# @param data
# @param col_dil
# @param col_res
# @param col_uid
# @keywords linearRange
# @export
# @examples
# linearRange

linearRange <- function(data,
                        col_dil = 'dilution_conc',
                        col_res = 'dilution_resp',
                        col_uid = 'uid') {

  # sanitize data
  data <- linearDataClean(data, col_dil, col_res)


  dats <- ddply(data, c('uid','dilution_conc'), summarize,
                res_mean = mean(dilution_resp))

  dats <- dats[order(dats$uid, dats$dilution_conc),]

  uid_cols <- unique(dats$uid)
  for( uid_col in uid_cols){
    this_dat <- dats[dats$uid == uid_col,]

    # remove upper dilution roll off
    w_max <- which(this_dat$res_mean == max(this_dat$res_mean))
    w_min <- which(this_dat$res_mean == min(this_dat$res_mean))

    w_max_rm <- which(this_dat$res_mean < this_dat[w_max,]$res_mean &
                        this_dat$dilution_conc > this_dat[w_max,]$dilution_conc)
    w_min_rm <- which(this_dat$res_mean > this_dat[w_min,]$res_mean &
                        this_dat$dilution_conc < this_dat[w_min,]$dilution_conc)

    if(length(unique(w_max_rm, w_min_rm)) == 0) next()
    dil_rm <- this_dat[unique(w_max_rm, w_min_rm),]$dilution_conc

    w_rm <- which(data$uid == uid_col & data$dilution_conc %in% dil_rm)

    data <- data[-w_rm,]

  }

  return(data)
}


