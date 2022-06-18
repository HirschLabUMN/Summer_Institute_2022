#######################################################################
####################  Setting up the file system  #####################
#######################################################################

#Determine current working directory
getwd()

#Set current working direcgtory
setwd("/Users/cnhirsch/Desktop")

#Check that current working directory is where you changed it to
getwd()

#See the list of files in the directory
list.files()

#######################################################################
####################  Built in Functions in R  ########################
#######################################################################

#Rounding of numbers
round(3.1414)

#Mathematical funcitons like factorial
factorial(3)

#Length of an object
name<-"alex"
length(name)
name<-c("Alex", "Tom")
length(name)

#R can also do nested functions going inside to outside
mean(1:6)
round(3.5)
mean_1_6<-mean(1:6)
round(mean_1_6)
round(mean(1:6))


#######################################################################
#######################  Getting help in R  ###########################
#######################################################################

#See a detailed manual for a function
?sample
help(sample)

#Get usage information for a function
args(sample)


#######################################################################
#################  Making your own function in R  #####################
#######################################################################

#This example function takes in a user defined input file, filters for a specific chromsome and then outputs to a user defined file
select_chromosome<-function(input_file_path,chrom,output_file_path){
  data=read.csv(input_file_path)
  filtered_data=subset(data, Chromosome==chrom)
  write.csv(filtered_data,file=output_file_path)
}

#Lets try out our function now
select_chromosome(chrom="chr3", input_file_path="snp.csv", output_file_path="snp_chr3.csv")


#######################################################################
#####################  Properties of a vector  ########################
#######################################################################

#create a vector
die <- c(1,2,3,4,5,6)

#see if it is a vector
is.vector(die)

#check what type of data is stored
typeof(die)

#check the length
length(die)

#parse a subset of the vector
die[2:4]
die[c(2,4,6)]

#doing math on a vector
die+1


#######################################################################
################  In class exercise - roll 2 dice  ####################
#######################################################################

roll <- function() {
  die <- 1:6
  dice <- sample(die, size = 2, replace = TRUE)
  sum(dice)
}

roll()


#######################################################################
#####################  Properties of a matrix  ########################
#######################################################################

#creating a 2x3 matrix from our die vector
m<-matrix(die, nrow=2)
m

#creating a 3x2 matrix from our die vector filling in by row
m<-matrix(die, ncol=2, byrow=TRUE)
m


#######################################################################
###################  Properties of a data frame  ######################
#######################################################################

#creating a data frame by hand
my_df<-data.frame(face = c("ace", "two", "six"), suit=c("clubs", "diamonds", "spades"), value=c(1,2,6))
my_df

#creating a data from by loading in data
mammals <- read.table("mammals_dataset.txt", sep="\t", stringsAsFactors=FALSE)

#Look at the structure of the data
str(mammals)

#See the first part of the data
head(mammals)

#Look at the dimensions of the dataset
dim(mammals)

#Look at the last few lines of the dataset
tail(mammals)

#Check for missing values in the adult_body_mass_g column
sum(is.na(mammals$adult_body_mass_g))

#Mean of the adult_body_mass_g column
mean(mammals$adult_body_mass_g, na.rm=TRUE)

#Check how many times each order appears in the dataset
table(mammals$order)
