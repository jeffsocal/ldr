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
# DESC:     linearDataClean
##################################################################################


# linearDataClean
#
# @param data
# @param col_dil
# @param col_res
# @keywords linearDataClean
# @export
# @examples
# linearDataClean

linearDataClean <- function(data,
                            col_dil,
                            col_res){
  dat <- data
  u_dil <- data.frame(dilution_conc = sort(unique(dat[,col_dil])),
                      dilution_n = 1:length(unique(dat[,col_dil]))
  )

  w <- which(colnames(dat) == col_dil)
  colnames(dat)[w] <- 'dilution_conc'
  w <- which(colnames(dat) == col_res)
  colnames(dat)[w] <- 'dilution_resp'

  dat <- dat[!is.na(dat$dilution_resp),]

  if(!'dilution_n' %in% colnames(dat))
    dat <- merge(u_dil, dat)

  return(dat)

}
