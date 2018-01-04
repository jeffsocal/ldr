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
# DESC:     calculate linear dynamic range for the dilution profiles
################################################################################


forceColName <- function(data, colnow, colset){
  coln <- colnames(data)
  w <- which(coln == colnow)
  colnames(data)[w] <- colset
  return(data)
}
