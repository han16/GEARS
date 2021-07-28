## code to prepare `DATASET` dataset goes here
# Load raw data from .txt file
exampleData1 <- read.table("data-raw/Beta10nodes20edges200samp.txt")
exampleData2 <- read.table("data-raw/gams10nodes20edges200samp.txt")
exampleData3 <- read.table("data-raw/parents10nodes20edges200samp.txt")
exampleData4 <- read.table("data-raw/simuData10nodes20edges200samp.txt")
# Apply preprocessing...
# Save the cleaned data in the required R package location
usethis::use_data(exampleData1)
usethis::use_data(exampleData2)
usethis::use_data(exampleData3)
usethis::use_data(exampleData4)

