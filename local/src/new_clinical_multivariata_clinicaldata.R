## clinical data multivariata

library(tidyverse)
library(sjPlot)
library(openxlsx)

cl_f <- snakemake@input[["clinical_data"]]
buoni_f <- snakemake@input[["pdo"]]
fra_f <- snakemake@input[["mut"]]
lmh_f <- snakemake@input[["umani"]]
msi_f <- snakemake@input[["msi"]]
res_circos <- snakemake@output[["df_circos"]]
plot_fit <- snakemake@output[["fit_plot"]]
res_fit <- snakemake@output[["results_fit"]]

## why new? modificati i buoni dopo essere ripartiti dalla table 2 di simo 

print("se clinical data non viene caricata controllare tramite visual studio code che 
alcune righe non contengano un invio non visibile su rstudio e cancellare l'invio
controllare che il tab sia corretto per mantere le colonne corrette")
## se clinical data non viene caricata controllare tramite visual studio code che 
## alcune righe non contengano un invio non visibile su rstudio e cancellare l'invio
## controllare che il tab sia corretto per mantere le colonne corrette

## carico i dati clinici, filtro per i pdo buoni e per le mutazioni annotate di fra
#cl_f <- "/scratch/trcanmed/biobanca/local/share/data/clinical_data_done_revision020424.tsv"
cl <- read.table(cl_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#buoni_f <- "/scratch/trcanmed/biobanca/local/share/data/biobanca_pdo_buoni.tsv"
#buoni_f <- "/scratch/trcanmed/biobanca/local/share/data/whoiswho_validation_xen_revision_derivation.tsv"
buoni <- read.table(buoni_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#names(buoni)[names(buoni) == "smodel"] <- "CASE"
names(buoni)[names(buoni) == "type"] <- "buoni"

merged <- merge(cl, buoni, by = "CASE")

#fra_f <- "/scratch/trcanmed/biobanca/dataset/V1/enrichment/fra_mutational_annotation.tsv"
fra <- read.table(fra_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(fra) <- fra$genes
fra$genes <- NULL
tfra <- as.data.frame(t(fra))
tfra$CASE <- rownames(tfra)
rownames(tfra) <- NULL

res <- merge(merged, tfra, by = "CASE", all.x = TRUE)
## remove hcc
#rownames(res) <- res$CASE
#row_names_df_to_remove<-c("CRC1379", "CRC1479", "CRC1737")
#res <- res[!(row.names(res) %in% row_names_df_to_remove),]

## creo un df con solo i dati necessari per la multivariata
res <- res[c("CASE", "AGE.AT.COLLECTION..years.", "SEX", "SITE.OF.PRIMARY", "STGE", "THERAPY.BEFORE.COLLECTION..Y.N.",
             "KRAS", "BRAF", "NRAS", "MSI_MSS", "buoni", "derivation_type")]
#names(res)[names(res) == "T"] <- "Classification_T"
#names(res)[names(res) == "N"] <- "Classification_N"

## sistemo il df per poter accorpare i siti del primario
rownames(res) <- res$CASE

res["CRC1241", "SITE.OF.PRIMARY"] <- "NONE"
res["CRC1575", "SITE.OF.PRIMARY"] <- "NONE"
res["CRC0578", "SITE.OF.PRIMARY"] <- "NONE"
res["CRC0277", "SITE.OF.PRIMARY"] <- "NONE"
res["CRC1376", "SITE.OF.PRIMARY"] <- "NONE"
res["CRC1474", "SITE.OF.PRIMARY"] <- "NONE"

for (i in seq(res$CASE)) {
  if (res[i, "SITE.OF.PRIMARY"] == "SIGMOID COLON" | res[i, "SITE.OF.PRIMARY"] == "ANO" | 
      res[i, "SITE.OF.PRIMARY"] == "LEFT COLON" | res[i, "SITE.OF.PRIMARY"] == "RECTUM" |
      res[i, "SITE.OF.PRIMARY"] == "SPLENIC FLEXURE"){ 
    res[i, "SITE.OF.PRIMARY"] <- "LEFT COLON"
  } else if (res[i, "SITE.OF.PRIMARY"] == "CAECUM" | res[i, "SITE.OF.PRIMARY"] == "TRANSVERSE COLON" | 
             res[i, "SITE.OF.PRIMARY"] == "HEPATIC FLEXURE" | res[i, "SITE.OF.PRIMARY"] == "RIGHT COLON" ) {
    res[i, "SITE.OF.PRIMARY"] <- "RIGHT COLON"
  } else {
    res[i, "SITE.OF.PRIMARY"] <- NA
  }
}

# ## sistemo T ed N per poterli rendere numerici senza perdere quelli con la lettera
# 
# res["CRC1241", "Classification_T"] <- "N"
# res["CRC1575", "Classification_T"] <- "N"
# res["CRC0578", "Classification_T"] <- "N"
# res["CRC1390", "Classification_T"] <- "N"
# 
# res["CRC1241", "Classification_N"] <- "N"
# res["CRC1575", "Classification_N"] <- "N"
# res["CRC0578", "Classification_N"] <- "N"
# res["CRC1390", "Classification_N"] <- "N"
# 
# for (i in seq(res$CASE)) {
#   if (res[i, "Classification_T"] == "4A") {
#     res[i, "Classification_T"] <- "4"
#   } else if (res[i, "Classification_T"] == "4B") {
#     res[i, "Classification_T"] <- "4"
#   } else {
#     res[i, "Classification_T"] <- res[i, "Classification_T"]
#   }
# }
# 
# for (i in seq(res$CASE)) {
#   if (res[i, "Classification_N"] == "1A") {
#     res[i, "Classification_N"] <- "1"
#   } else if (res[i, "Classification_N"] == "1B") {
#     res[i, "Classification_N"] <- "1"
#   } else if (res[i, "Classification_N"] == "2A") {
#     res[i, "Classification_N"] <- "2"
#   } else if (res[i, "Classification_N"] == "2B") {
#     res[i, "Classification_N"] <- "2"
#   } else {
#     res[i, "Classification_N"] <- res[i, "Classification_T"]
#   }
# }

## sistemo stage per non perdere i casi con la lettera
names(res)[names(res)== "STGE"] <- "STAGE"

res["CRC1241", "STAGE"] <- "N"
res["CRC1575", "STAGE"] <- "N"
res["CRC0578", "STAGE"] <- "N"
res["CRC1390", "STAGE"] <- "N"
res["CRC0277", "STAGE"] <- "N"
res["CRC1376", "STAGE"] <- "N"
res["CRC1474", "STAGE"] <- "N"
res["CRC1336", "STAGE"] <- "N"
res["CRC1501", "STAGE"] <- "N"

for (i in seq(res$CASE)) {
  if (res[i, "STAGE"] == "2A") {
    res[i, "STAGE"] <- "2"
  } else if (res[i, "STAGE"] == "2B") {
    res[i, "STAGE"] <- "2"
  } else if (res[i, "STAGE"] == "3A") {
    res[i, "STAGE"] <- "3"
  } else if (res[i, "STAGE"] == "3B") {
    res[i, "STAGE"] <- "3"
  } else if (res[i, "STAGE"] == "4A") {
    res[i, "STAGE"] <- "4"
  } else {
    res[i, "STAGE"] <- res[i, "STAGE"]
  }
}


## rendo tutto numerico o factor a seconda delle necessità

res$AGE.AT.COLLECTION <-  gsub('N', NA, res$AGE.AT.COLLECTION)
res$AGE.AT.COLLECTION <- as.numeric(res$AGE.AT.COLLECTION)

res$SEX <-  gsub('N', NA, res$SEX)
res$SEX <- as.factor(res$SEX)

# res$Classification_T <- gsub("N", NA, res$Classification_T)
# res$Classification_T <- as.numeric(res$Classification_T)
# 
# res$Classification_N <- gsub("N", NA, res$Classification_N)
# res$Classification_N <- as.numeric(res$Classification_N)

res$SITE.OF.PRIMARY <- gsub("NONE", NA, res$SITE.OF.PRIMARY)
res$SITE.OF.PRIMARY <- as.factor(res$SITE.OF.PRIMARY)

res$STAGE <-  gsub('N', NA, res$STAGE)
res$STAGE <- as.numeric(res$STAGE)

res$THERAPY.BEFORE.COLLECTION..Y.N. <- as.factor(res$THERAPY.BEFORE.COLLECTION..Y.N.)

#msi_fra <- "/scratch/trcanmed/biobanca/local/share/data/MSIstatus.xlsx"
msi <- read.xlsx(msi_f)
msi$N <- NULL
msi$model <- substr(msi$Genealogy.ID, 1, 7)
msi$status <- gsub("\\*", "", msi$status)
msi$Genealogy.ID <- NULL
names(msi)[names(msi)=="model"] <- "CASE"
rownames(msi) <- msi$CASE
msi$CASE <- NULL

for (i in rownames(msi)) {
  res[i, "MSI_MSS"] <- msi[i, "status"]
}

is.na(res$MSI_MSS) <- NA
res$MSI_MSS <-  gsub('NT', NA, res$MSI_MSS)
res$MSI_MSS <- as.factor(res$MSI_MSS)
res$MSI_MSS <- gsub("-H", "", res$MSI_MSS)

res$derivation_type <- as.factor(res$derivation_type)

## scrivo la tabella per i circos che comprenda anche i casi
write.table(res, file = res_circos, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

## rimuovere i 14 lmh
#lmh_f <- "/scratch/trcanmed/biobanca/local/share/data/lmh_detele_multivariate.tsv"
lmh <- read.table(lmh_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
lmh <- lmh$CASE

res <- res %>% filter(!CASE %in% lmh)

res <- res %>% filter(!buoni == "Validation not performed")
res$buoni <-  gsub('Validation successful', "True", res$buoni)
res$buoni <- gsub("Validation failed","False", res$buoni)
res$buoni <- gsub("Not established","False", res$buoni)
res$buoni <- as.factor(res$buoni)


res$CASE <- NULL
#res_prova <- res 
#res_prova <- res_prova %>% filter(!KRAS == "True")

# faccio il fit

sink(snakemake@log[['log']])
print(table(res$buoni))
sink()
#fit.full <- glm(buoni ~ SEX + AGE.AT.COLLECTION, data=res, family=binomial())
#fit.full <- glm(buoni ~ SEX + AGE.AT.COLLECTION + STAGE + THERAPY.BEFORE..Y.N. , data=res, family=binomial())
fit.full <- glm(buoni ~ SITE.OF.PRIMARY + STAGE + THERAPY.BEFORE.COLLECTION..Y.N. + KRAS + BRAF + NRAS + SEX + AGE.AT.COLLECTION, data = res, family = binomial())
#fit.full <- glm(buoni ~ SEX + AGE.AT.COLLECTION + THERAPY.BEFORE..Y.N. + SITE.OF.PRIMARY + STAGE + KRAS + NRAS + BRAF, data=res, family=binomial())
pdf(plot_fit, useDingbats=FALSE)
#plot_model(fit.full)
#name <- c("Sex", "Age", "Therapy", "Site of primary", "Stage", "KRAS", "NRAS", "BRAF")
#name <- c("Site of primary", "Stage", "Therapy", "KRAS", "BRAF", "NRAS", "Sex", "Age")
plot_model(fit.full, axis.lim = c(0.1, 10), title = "Validation")
dev.off()
#fit.full <- glm(buoni ~ SEX + AGE.AT.COLLECTION + THERAPY.BEFORE..Y.N. + SITE.OF.PRIMARY + BRAF + NRAS, data=res_prova, family=binomial())
fit <- as.data.frame(summary.glm(fit.full)$coefficients)
## intervallo di confidenza
conf_intervals <- confint(fit.full)
fit2 <- cbind(fit, exp(conf_intervals))
## odds ratio
# Obtain coefficients
coefficients <- coef(fit.full)
odds_ratios <- as.data.frame(exp(coefficients))
colnames(odds_ratios) <- "odds ratio"
fit3 <- cbind(fit2, odds_ratios)
rownames(fit3) <- c("Intercept", "Site of primary (right colon)", "Stage", "Therapy Before (yes)", 
                    "KRAS (mutant)", "BRAF (mutant)", "NRAS (mutant)", "Sex (male)", "Age at collection")

write.table(fit3, file = res_fit, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)