getwd()
#setwd("/Users/cnhirsch/Desktop/")

#########################################################################
#################### Siumulating Data ###################################
#########################################################################

data1<-as.data.frame(matrix(data=c(1,3,5,7,9,2,4,6,8,10), ncol=2))
data1

data2<-rnorm(1000, mean=50, sd=5)


#########################################################################
#################### Making Basic Plots #################################
#########################################################################

#Making a scatter plot with vanilla R
plot(data1$V1, data1$V2)

#Modifying the x and y labels and changing the symbol type
plot(data1$V1, data1$V2, xlab="Type 1", ylab="Type 2", pch=16)

#Modifying the x and y axis limits
plot(data1$V1, data1$V2, xlab="Type 1", ylab="Type 2", pch=16, xlim=c(0,10), ylim=c(0,10))

#Add a line with intercept 0 and slope of 1
abline(a=0, b=1)

#Change the color to be red
abline(a=0, b=1, col="red")

#Makeing a histogram with vanilla R
hist(data2)

#Modifying x-axis label, color, and number of breaks
hist(data2, xlab="Height", col="gray", breaks=20)

#Remove the main header
hist(data2, xlab="Height", col="gray", breaks=20, main="")

#Remove the space below the histogram
hist(data2, xlab="Height", col="gray", breaks=20, main="", xaxt="n")
axis(1, at=c(35,40,45,50,55,60,65,70), lab=c(35,40,45,50,55,60,65,70), pos=0)

#Opening a PDF file to write this figure to file
pdf.options(family="Helvetica")
pdf("hist.pdf", pointsize=10)
hist(data2, xlab="Height", col="gray", breaks=20, main="", xaxt="n")
axis(1, at=c(35,40,45,50,55,60,65,70), lab=c(35,40,45,50,55,60,65,70), pos=0)
dev.off()

#########################################################################
##################### Loading  ggplot2 ##################################
#########################################################################

install.packages("ggplot2")
library("ggplot2")  

#########################################################################
########################## Basic Plots ##################################
#########################################################################

#Making a scatter plot
qplot(data1$V1, data1$V2)

#Making a histogram
qplot(data2, geom="histogram")

#Another way to make a scatter plot
ggplot(data1, aes(V1, V2)) + geom_point()

#########################################################################
######################## Advanced Plots #################################
#########################################################################

#Creating a more complex data set
data3<-as.data.frame(matrix(data=c(1,2,3,4,5,1,3,5,7,9,2,4,6,8,10), ncol=3))
colnames(data3) <- c("Day", "Treatment1", "Treatment2")
data3

#Install and load reshape to convert data to long format
install.packages("reshape2")
library("reshape2")

#Melt data to long format
data3_long<-melt(data3, id="Day")
data3_long

#Lets make some more complex plots
ggplot(data3_long, aes(Day, value)) + geom_point()
ggplot(data3_long, aes(Day, value, color=variable)) + geom_point()
ggplot(data3_long, aes(Day, value, color=variable)) + geom_point() + geom_line()
ggplot(data3_long, aes(Day, value)) + geom_point() + geom_line(aes(color=variable))

