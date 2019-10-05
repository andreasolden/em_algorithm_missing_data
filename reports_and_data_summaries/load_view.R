#View
#Load

tot_time
apply(monteboot, 2, function(x) sum(is.na(x)))
colMeans(monteboot, na.rm = TRUE)
