#!/applications/R/R-3.5.0/bin/Rscript

# Extract accession IDs of hexaploid wheat, and wild and dommesticated emmer wheat accessions
# used in exome-derived SNP analysis published in  He et al. (2019) Nat. Genet. 51: https://doi.org/10.1038/s41588-019-0382-2 ,
# grouped by geographic region

tblS1 <- read.table("Table_S1_41588_2019_382_MOESM3_ESM.tsv",
                    quote = "\"", sep = "\t", header = T,
                    check.names = T, stringsAsFactors = F)
# Extract Triticum aestivum subsp. aestivum accessions only
tblS1_DWemmer <- rbind(tblS1[which(tblS1$X21.populations.used.in.Admixture.analyses == "domesticatedE"),],
                       tblS1[which(tblS1$X21.populations.used.in.Admixture.analyses == "wildE"),])
tblS1 <- tblS1[which(tblS1$Species == "Triticum aestivum subsp. aestivum"),]
# OR
#tblS1 <- tblS1[which(tblS1$improvement.Status %in% c("cul", "lr", "unknown")),]

# NOTE: this script removes 5 accessions, as was done in the above paper:
# "On the basis of the patterns of clustering inconsistent with the clustering of wheat accessions
# having the same improvement status (wild, domesticated, landrace, cultivar), we have removed
# from analyses three wheat landraces (PI 345355, PI 131592, PI 534284) and
# two accessions of wild and domesticated emmer (PI 467000, PI 415152)." (He et al. 2019 Nat. Genet. 51)
tblS1 <- tblS1[-c(which(grepl("345355", tblS1$Accession_ID)),
                  which(grepl("131592", tblS1$Accession_ID)),
                  which(grepl("534284", tblS1$Accession_ID))),]
tblS1_DWemmer <- tblS1_DWemmer[-c(which(grepl("467000", tblS1_DWemmer$Accession_ID)),
                                  which(grepl("415152", tblS1_DWemmer$Accession_ID))),]

pop <- sort(unique(tblS1$X21.populations.used.in.Admixture.analyses))
pop <- pop[which(pop != "-")]

region <- sort(unique(tblS1$Continent))
region <- region[which(region != "-")]

impStat <- sort(unique(tblS1$improvement.Status))
impStat <- impStat[which(impStat != "-")] 

NorthAfrica <- tblS1[which(grepl(pattern = "AfricaNorth",
                                 x = tblS1$X21.populations.used.in.Admixture.analyses,
                                 ignore.case = F, fixed = T)),]
write.table(NorthAfrica$internal_ID,
            file = "NorthAfrica_accessions.txt",
            quote = F, row.names = F, col.names = F)

SubSaharanAfrica <- tblS1[which(grepl(pattern = "AfricaSouth",
                                      x = tblS1$X21.populations.used.in.Admixture.analyses,
                                      ignore.case = F, fixed = T)),]
write.table(SubSaharanAfrica$internal_ID,
            file = "SubSaharanAfrica_accessions.txt",
            quote = F, row.names = F, col.names = F)

Africa <- tblS1[which(grepl(pattern = "Africa",
                            x = tblS1$X21.populations.used.in.Admixture.analyses,
                            ignore.case = F, fixed = T)),]
write.table(Africa$internal_ID,
            file = "Africa_accessions.txt",
            quote = F, row.names = F, col.names = F)

CentralAmerica <- tblS1[which(tblS1$X21.populations.used.in.Admixture.analyses %in%
                              c("AmeCentral", "AmeMexico")),]
write.table(CentralAmerica$internal_ID,
            file = "CentralAmerica_accessions.txt",
            quote = F, row.names = F, col.names = F)

NorthAmerica <- tblS1[which(grepl(pattern = "AmeN",
                                  x = tblS1$X21.populations.used.in.Admixture.analyses,
                                  ignore.case = F, fixed = T)),]
write.table(NorthAmerica$internal_ID,
            file = "NorthAmerica_accessions.txt",
            quote = F, row.names = F, col.names = F)

SouthAmerica <- tblS1[which(grepl(pattern = "AmeS",
                                  x = tblS1$X21.populations.used.in.Admixture.analyses,
                                  ignore.case = F, fixed = T)),]
write.table(SouthAmerica$internal_ID,
            file = "SouthAmerica_accessions.txt",
            quote = F, row.names = F, col.names = F)

SouthAsia <- tblS1[which(grepl(pattern = "AsiaIndia",
                               x = tblS1$X21.populations.used.in.Admixture.analyses,
                               ignore.case = F, fixed = T)),]
write.table(SouthAsia$internal_ID,
            file = "SouthAsia_accessions.txt",
            quote = F, row.names = F, col.names = F)

EastAsia <- tblS1[which(grepl(pattern = "AsiaEast",
                              x = tblS1$X21.populations.used.in.Admixture.analyses,
                              ignore.case = F, fixed = T)),]
write.table(EastAsia$internal_ID,
            file = "EastAsia_accessions.txt",
            quote = F, row.names = F, col.names = F)

CentralAsia <- tblS1[which(grepl(pattern = "AsiaCentral",
                                 x = tblS1$X21.populations.used.in.Admixture.analyses,
                                 ignore.case = F, fixed = T)),]
write.table(CentralAsia$internal_ID,
            file = "CentralAsia_accessions.txt",
            quote = F, row.names = F, col.names = F)

MiddleEast <- tblS1[which(grepl(pattern = "MiddleEast",
                                x = tblS1$X21.populations.used.in.Admixture.analyses,
                                ignore.case = F, fixed = T)),]
write.table(MiddleEast$internal_ID,
            file = "MiddleEast_accessions.txt",
            quote = F, row.names = F, col.names = F)

Asia <- rbind(CentralAsia, SouthAsia, EastAsia)
write.table(Asia$internal_ID,
            file = "Asia_accessions.txt",
            quote = F, row.names = F, col.names = F)

FormerSU <- tblS1[which(grepl(pattern = "SOVIETUNION",
                              x = tblS1$X21.populations.used.in.Admixture.analyses,
                              ignore.case = F, fixed = T)),]
write.table(FormerSU$internal_ID,
            file = "FormerSU_accessions.txt",
            quote = F, row.names = F, col.names = F)

EasternEurope <- tblS1[which(tblS1$Continent == "EuropeEast"),]
write.table(EasternEurope$internal_ID,
            file = "EasternEurope_accessions.txt",
            quote = F, row.names = F, col.names = F)

WesternEurope <- tblS1[which(tblS1$Continent == "EuropeWest"),]
write.table(WesternEurope$internal_ID,
            file = "WesternEurope_accessions.txt",
            quote = F, row.names = F, col.names = F)

Oceania <- tblS1[which(tblS1$Continent == "OCEANIA"),]
write.table(Oceania$internal_ID,
            file = "Oceania_accessions.txt",
            quote = F, row.names = F, col.names = F)

DomesticatedEmmer <- tblS1_DWemmer[which(tblS1_DWemmer$X21.populations.used.in.Admixture.analyses == "domesticatedE"),]
write.table(DomesticatedEmmer$internal_ID,
            file = "DomesticatedEmmer_accessions.txt",
            quote = F, row.names = F, col.names = F)

WildEmmer <- tblS1_DWemmer[which(tblS1_DWemmer$X21.populations.used.in.Admixture.analyses == "wildE"),]
write.table(WildEmmer$internal_ID,
            file = "WildEmmer_accessions.txt",
            quote = F, row.names = F, col.names = F)
