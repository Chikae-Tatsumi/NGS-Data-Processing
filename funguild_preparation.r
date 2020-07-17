setwd("~/R/Analysis/1_Test/ITS")
ASV.table <- read.table(file="rarefied_ASV_table.txt")
taxonomy <- paste(ASV.table$Kingdom, ASV.table$Phylum, ASV.table$Class, ASV.table$Order, ASV.table$Family, ASV.table$Genus, ASV.table$Species, sep="; ")
ASV.table <- ASV.table [,1:(ncol(ASV.table)-7)] 
ASV.table <- cbind (ASV.table, taxonomy)
ASV.table <- cbind (rownames(ASV.table),ASV.table)
colnames(ASV.table)[1] <- "OTU ID"
write.table (ASV.table, file="rarefied_ASV_table_funguild.txt",sep="\t",row.names=F, quote=F)

# Go to http://www.stbates.org/guilds/app.php
