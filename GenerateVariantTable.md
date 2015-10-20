## Get the Variant files

http://watson.nci.nih.gov/projects/nci60/wes/VCF/

	 wget -r -l1 --no-parent -A.vcf http://watson.nci.nih.gov/projects/nci60/wes/VCF/
	 wget -r -l1 --no-parent -A.idx http://watson.nci.nih.gov/projects/nci60/wes/VCF/
	 
	 mkdir VCF
	 mv watson.nci.nih.gov/projects/nci60/wes/VCF/* ./VCF
	 rm -r watson.nci.nih.gov/
	 
## annotate the variants 
	mkdir VCF_snpEff
	for i in ./VCF/*.vcf;\
	do \
	j=`echo $i | sed s/\.vcf/.snpEff.vcf/`;\
	java -Xmx8G -jar $SNPEFF_HOME/snpEff.jar eff -verbose -canon -no-downstream -no-upstream -config $SNPEFF_HOME/snpEff.config hg19 $i > $j;\
	done


## generate table in R

	setwd(".")
	library(VariantAnnotation)

	# now all of them
	listVCF <- list.files("./Data/VCF_snpEff", "*.vcf")
	CELLNAME <- sub(".snpEff.vcf", "", listVCF, fixed = TRUE)

	GeneLists <- vector("list", length(listVCF))
	names(GeneLists) <- CELLNAME

	for(ii in seq_along(listVCF)) {
		temp <- readVcf(file.path("./Data/VCF_snpEff", listVCF[ii]), genome = "hg19")
		print(listVCF[ii])
		temp_T2 <- temp[which(rowRanges(temp)$FILTER == "Type2")]
		GeneSymbol <- sapply(info(temp_T2)$EFF, function(xx) strsplit(xx[1], "|", fixed = TRUE)[[1]][6])
		info(temp_T2)$GeneSymbol <- GeneSymbol

		EFFECT <- sapply(info(temp_T2)$EFF, function(xx) strsplit(xx[1], "(", fixed = TRUE)[[1]][1])
		info(temp_T2)$EFFECT <- EFFECT
		print(table(info(temp_T2)$EFFECT))

		# keep only non synonymous
		# note this is based on previous version of snpEff, they have updated the EFFECT codes.
		# http://snpeff.sourceforge.net/SnpEff_manual.html#eff
		Eff_Keep <- c("CODON_CHANGE_PLUS_CODON_DELETION", "CODON_CHANGE_PLUS_CODON_INSERTION", "CODON_DELETION", "CODON_INSERTION", "FRAME_SHIFT", "NON_SYNONYMOUS_CODING", "STOP_GAINED", "START_LOST", "STOP_LOST", "CODON_CHANGE")

		temp_Type2_NS <- temp_T2[which(info(temp_T2)$EFFECT %in% Eff_Keep)]

		GeneLists[[ii]] <- sort(unique(info(temp_Type2_NS)$GeneSymbol))
		
		# alternative with all type 2 variants
		# test <- readVcf(file.path("./Data/VCF_snpEff", listVCF[ii]), genome = "hg19")
		# test_T2 <- test[which(rowData(test)$FILTER == "Type2")]
		# GeneSymbol <- sapply(info(test_T2)$EFF, function(xx) strsplit(xx[1], "|", fixed = TRUE)[[1]][6])
		# GeneLists[[ii]] <- sort(unique(GeneSymbol))	
	}

	# number of genes per cell line
	sort(sapply(GeneLists, length))

	AllGenes <- sort(unique(unlist(GeneLists)))
	# remove the NA gene names
	AllGenes <- AllGenes[(which(!(AllGenes %in% c(""))))]

	VariantTable <- matrix(0, nrow = length(CELLNAME), ncol = length(AllGenes))
	rownames(VariantTable) <- CELLNAME
	colnames(VariantTable) <- AllGenes
	for(ii in seq_along(CELLNAME)) {
		VariantTable[ii, which(AllGenes %in% GeneLists[[ii]])] <- 1
	}

	# save(VariantTable, file = "VariantTableByGene.RData")
	write.csv(VariantTable, file = "VariantTableByGene.csv")
