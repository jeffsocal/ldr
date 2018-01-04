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
# DESC:     linearTTest
##################################################################################


# linearTTest
#
# @param data
# @param col_dil
# @param col_res
# @param col_uid
# @param rjct_tt
# @keywords linearTTest
# @export
# @examples
# linearTTest

linearTTest <- function(data,
                        col_dil = 'dilution_conc',
                        col_res = 'dilution_resp',
                        col_uid = 'uid',
                        rjct_tt = 0.05){

  # sanitize data
  data <- linearDataClean(data, col_dil, col_res)

  uid_cols <- unique(data[,col_uid])
  for( uid_col in uid_cols){
    this_dat <- data[data[,col_uid] == uid_col,]
    uid_dils <- unique(this_dat$dilution_n)
    uid_dils <- uid_dils[order(uid_dils)]

    uid_dil_last <- 0
    for( uid_dil in uid_dils ){
      if(uid_dil_last == 0){
        uid_dil_last <- uid_dil
        next()
      }

      tt_out <- t.test(this_dat[this_dat$dilution_n == uid_dil_last,col_res], this_dat[this_dat$dilution_n == uid_dil,col_res])

      cat(uid_dil_last, " -> ", uid_dil, "\n")
      print(this_dat[this_dat$dilution_n %in% c(uid_dil_last, uid_dil),])
      print(tt_out)

      uid_dil_last <- uid_dil
    }
  }

}
