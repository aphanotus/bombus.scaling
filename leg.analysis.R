# Explore the leg data

names(m)
unique(m$species)

legs <- read.csv("leg.measurements.csv")
dim(legs)
head(legs)
sum(duplicated(legs$specimen_id))
sum(grepl("duplicate",legs$specimen_id))/2

names(legs)
measurement.names <- names(legs)[-c(1:5,42)]

# Inspect distributions for outliers
{
  par(mfrow=c(6,6), mai = c(0.25, 0.25, 0.25, 0.25))
  lapply(measurement.names, function (name.i) {
    x <- legs[,name.i]
    hist(x, main = name.i)
  })
  par(mfrow=c(1,1), mai = (c(5, 4, 4, 2) + 0.1))
}

# Check for zeros
invisible(lapply(measurement.names, function (name.i) {
  x <- legs[,name.i]
  x <- x[!is.na(x)]
  cat(name.i,"\t",sum(x==0),"\n")
}))


