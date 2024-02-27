extractData <- function(data, use, immuneIDX){
  if(use == "s"){
    final <- data[1:54675,]
  } else if (use == "i"){
    raw <- data[54676:nrow(data),]
    subset <- raw[immuneIDX,]
    final <- subset$Expression.Value %>% as.matrix()
    rownames(final) <- subset[,1]
    colnames(final) <- "ex"
    final
  }
}

# immuneGood <- extractData(test_good, "i")
# head(immuneGood)
# dim(immuneGood)
# class(immuneGood)
# predict(finalKNN, t(immuneGood))
# 
# immuneBad <- extractData(test_bad, "i")
# head(immuneBad)
# dim(immuneBad)
# class(immuneBad)
# predict(finalKNN, t(immuneBad))
# # 
# # i <- match("NA.1", rownames(immuneBad))
# # 
# # ##hmmm
# # immuneGood[c("NA."),]
# # poopy[c("NA."),]
# # immuneBad[c("NA."),]
# # 
# # poopy <- immune_data[5,-1] %>% t()
# # head(poopy)
# # dim(poopy)
# # class(poopy)
# # predict(immuneKNN, t(poopy))
# # predict(immuneKNN, t(immuneBad))
