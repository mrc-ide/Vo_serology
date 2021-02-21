################################################################################
# Function reurning case counts for Figure 1 and Table 4                       #
# Frequency of assay result combinations (e.g. A+D+R+, etc.)                   #
################################################################################

output_case_counts <- function(ds){

  # ------------------------------ May 2020 -----------------------------------#
  message("subjects with at least one serological result (May)")
  
  Dt <- which(!is.na(ds$Diasorin_IgG_qualitative))
  Rt <- which(!is.na(ds$Roche_Total_qualitative))
  At <- which(!is.na(ds$Abbot_qualitative))   
  idx <- unique(c(Dt, Rt, At))
  
  print(length(idx))

  # ---------------------------------------------------------------------------#
  message("subjects with both PCR and at least one sero result (May)")
  
  tab <- table(ds$positive_swab[idx], useNA = "no")
  
  print(c(sum(tab["FALSE"]+tab["TRUE"]), 
      round((sum(tab["FALSE"]+tab["TRUE"])/length(idx))*100, 1)))
  
  # ---------------------------------------------------------------------------#
  message("subjects negative to all three assays (May)") 
  
  tab <- table(ds$Diasorin_IgG_qualitative,
               ds$Abbot_qualitative,
               ds$Roche_Total_qualitative, 
               useNA = "always")
  
    print(c(tab["Negative", "Negative", "Negative"], 
          round((tab["Negative", "Negative", "Negative"]/length(idx))*100, 1)))
   
  # ---------------------------------------------------------------------------#         
  message("subjects negative to at least one assay, NA to others (May)")
  
  x <- tab["Negative",which(is.na(rownames(tab["Negative",,]))), 
           which(is.na(colnames(tab["Negative",,])))] + 
    tab[which(is.na(rownames(tab[,,"Negative"]))),"Negative","Negative"]+
    tab["Negative",which(is.na(colnames(tab[,,"Negative"]))),"Negative"]
  
  print(c(x, round((x/length(idx))*100,1)))
  
  # ---------------------------------------------------------------------------#
  message("subjects positive to at least one assay (May)") 
  
  y <- length(unique(ds$INTERNAL_ID[
      which(ds$Diasorin_IgG_qualitative == "Positive" | 
              ds$Roche_Total_qualitative == "Positive" |
              ds$Abbot_qualitative == "Positive")]))
  
  print(c(y, round((y/length(idx))*100,1)))
  
  # ---------------------------------------------------------------------------#
  message("subjects with neutralising titres (> 1:40) (May)")
  
  tab <- table(ds$Neutralization)
  
  z <- sum(tab[which(rownames(tab) %in% c("1:80", "1:160", "1:320", "1:640", 
                                          "1:1280"))])
  
  print(c(z, round((z/y)*100,1)))
  
  # ---------------------------------------------------------------------------#
  message("subjects with known results (excluding equivocal) for all assays (May)")
  
  tab <- table(ds$Abbot_qualitative, 
               ds$Diasorin_IgG_qualitative,
               ds$Roche_Total_qualitative, 
               useNA = "no")
  
  res <- data.frame(res = c("A+D+R+",
                            "A+D-R+",
                            "A+D+R-",
                            "A-D+R+",
                            "A+D-R-",
                            "A-D+R-",
                            "A-D-R+",
                            "A-D-R-"),
                    
                    freq = c(tab["Positive", "Positive", "Positive"],
                             tab["Positive", "Negative", "Positive"], 
                             tab["Positive", "Positive", "Negative"],
                             tab["Negative", "Positive", "Positive"],
                             tab["Positive", "Negative", "Negative"],
                             tab["Negative", "Positive", "Negative"],
                             tab["Negative", "Negative", "Positive"],
                             tab["Negative", "Negative", "Negative"]))
  
  print(c(sum(res$freq), round((sum(res$freq)/length(idx))*100,1)))
  
  # ---------------------------------------------------------------------------#
  message("assay combinations (May)") 
  
  res$prob <- round(res$freq/sum(res$freq)*100,1)
  
  print(res)
  
  # ---------------------------------------------------------------------------#
  message("PCR negative (May)") 
  
  tab2 <- table(ds$positive_swab,
                ds$Abbot_qualitative, 
                ds$Diasorin_IgG_qualitative, 
                ds$Roche_Total_qualitative,
                useNA = "no")
  
  res_pcr_neg <- data.frame(res = c("A+D+R+",
                            "A+D-R+",
                            "A+D+R-",
                            "A-D+R+",
                            "A+D-R-",
                            "A-D+R-",
                            "A-D-R+",
                            "A-D-R-"),
                    
                    freq = c(tab2["FALSE", "Positive", "Positive", "Positive"],
                             tab2["FALSE", "Positive", "Negative", "Positive"],
                            tab2["FALSE", "Positive", "Positive", "Negative"],
                            tab2["FALSE", "Negative", "Positive", "Positive"],
                            tab2["FALSE", "Positive", "Negative", "Negative"],
                            tab2["FALSE", "Negative", "Positive", "Negative"],
                            tab2["FALSE", "Negative", "Negative", "Positive"],
                            tab2["FALSE", "Negative", "Negative", "Negative"]))
 
  res_pcr_neg$prob <- round(res_pcr_neg$freq/sum(res_pcr_neg$freq)*100, 1)
  
  print(res_pcr_neg)
  
  # ---------------------------------------------------------------------------#
  message("PCR positve (May)") 
  
  res_pcr_pos <- 
    data.frame(res = c("A+D+R+",
                       "A+D-R+",
                       "A+D+R-",
                       "A-D+R+",
                       "A+D-R-",
                       "A-D+R-",
                       "A-D-R+",
                       "A-D-R-"),
               
               freq = c(tab2["TRUE", "Positive", "Positive", "Positive"],
                        tab2["TRUE", "Positive", "Negative", "Positive"],
                        tab2["TRUE", "Positive", "Positive", "Negative"],
                        tab2["TRUE", "Negative", "Positive", "Positive"],
                        tab2["TRUE", "Positive", "Negative", "Negative"],
                        tab2["TRUE", "Negative", "Positive", "Negative"],
                        tab2["TRUE", "Negative", "Negative", "Positive"],
                        tab2["TRUE", "Negative", "Negative", "Negative"]))
  
  res_pcr_pos$prob <- round(res_pcr_pos$freq/sum(res_pcr_pos$freq)*100, 1)
  
  print(res_pcr_pos)
  
  # ---------------------------------------------------------------------------#
  message("swab testing results (May)") 
  
  print(table(ds$Third_sampling_swab))
  
  # --------------------------- November 2020 ---------------------------------#
  
  message("subjects with at least one serological result (Nov)")
  Dt2 <- which(ds$Diasorin_qualiitative_november_2020 != "#N/A")
  Rt2 <- which(ds$Roche_qualitative_november_2020 != "#N/A")
  At2 <- which(ds$Abbott_qualitative_november_2020 != "#N/A")   
  idx2 <- unique(c(Dt2, Rt2, At2))
  
  print(length(idx2))
  
  # ---------------------------------------------------------------------------#
  message("number of subjects with known results for all three assays (Nov)")
  
  tab3 <- table(ds$Abbott_qualitative_november_2020, 
                ds$Diasorin_qualiitative_november_2020, 
                ds$Roche_qualitative_november_2020,
                useNA = "no")
  
  res_nov <- data.frame(res = c("A+D+R+",
                                    "A+D-R+",
                                    "A+D+R-",
                                    "A-D+R+",
                                    "A+D-R-",
                                    "A-D+R-",
                                    "A-D-R+",
                                    "A-D-R-"),
                            
                            freq = c(tab3["Positive", "Positive", "Positive"],
                                     tab3["Positive", "Negative", "Positive"], 
                                     tab3["Positive", "Positive", "Negative"],
                                     tab3["Negative", "Positive", "Positive"],
                                     tab3["Positive", "Negative", "Negative"],
                                     tab3["Negative", "Positive", "Negative"],
                                     tab3["Negative", "Negative", "Positive"],
                                     tab3["Negative", "Negative", "Negative"]))
  
  res_nov$prob <- round(res_nov$freq/sum(res_nov$freq)*100, 1)
  
  print(c(sum(res_nov$freq), round((sum(res_nov$freq)/length(idx2))*100,1)))
  
  # ---------------------------------------------------------------------------#
  message("subjects positive to at least one assay (Nov)") 
  
  y2 <- length(unique(ds$INTERNAL_ID[
    which(ds$Diasorin_qualiitative_november_2020 == "Positive" | 
            ds$Roche_qualitative_november_2020 == "Positive" |
            ds$Abbott_qualitative_november_2020 == "Positive")]))
  
  print(c(y2, round((y2/length(idx2))*100,1)))
  
  # ---------------------------------------------------------------------------#
  message("subjects negative to all three assays (Nov)") 
  print(c(tab3["Negative", "Negative", "Negative"], 
          round((tab3["Negative", "Negative", "Negative"]/length(idx2))*100, 1)))
  
  # ---------------------------------------------------------------------------#
  message("assay combinations (Nov)") 
  print(res_nov)
  
  # ---------------------------------------------------------------------------#
  message("swab testing results (Nov)") 
  print(table(ds$Swab_november_2020))
  
}
################################################################################
# retrieve standard color codes                                                #
################################################################################

gg_colour_hue <- function(n){
  
  hues= seq(15, 375, length = n + 1)
  
  hcl(h = hues, l = 65, c = 100)[1:n]
}

################################################################################
################################################################################