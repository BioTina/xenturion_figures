# Compute delta of gistic scores for genes
library(ggplot2)
library(ggrastr)

gistic_us_f <- snakemake@input[['us']]
gistic_msk_f <- snakemake@input[['msk']]
corrplot_f <- snakemake@output[['corrplot']]
corrs_f <- snakemake@output[['corrs']]
kind <- snakemake@wildcards[['kind']]
log <- snakemake@log[['log']]

osnakemake <- snakemake
load(snakemake@input[['Rimage']])
#eval(parse(text=myriad))
snakemake <- osnakemake

df_us <- read.table(gistic_us_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df_msk <- read.table(gistic_msk_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)
rownames(df_us) <- df_us$Gene.Symbol
rownames(df_msk) <- df_msk$Gene.Symbol
df_us$Locus.ID <- NULL
df_us$Gene.Symbol <- NULL
df_us$Cytoband <- NULL
df_msk$Locus.ID <- NULL
df_msk$Gene.Symbol <- NULL
df_msk$Cytoband <- NULL

#save.image('pippo.Rdata')


get_freq <- function(us, msk, direction, kind) {
  if (direction == 1) {
    freq_us <- as.data.frame(apply(us, 1, function(x) { sum(x >= direction)/length(x) }))
    freq_msk <- as.data.frame(apply(msk, 1, function(x) { sum(x >= direction)/length(x) }))
  } else {
    freq_us <- as.data.frame(apply(us, 1, function(x) { sum(x <= direction)/length(x) }))
    freq_msk <- as.data.frame(apply(msk, 1, function(x) { sum(x <= direction)/length(x) }))
  }
  colnames(freq_us) <- kind
  colnames(freq_msk) <- 'MSK'
  res <- merge(freq_us, freq_msk, by='row.names')
  sink(log)
  print(nrow(res))
  print(nrow(freq_us))
  print(nrow(freq_msk))
  sink()
  rownames(res) <- res$Row.names
  res$Row.names <- NULL
  return(res)
}

freqs_amp <- get_freq(df_us, df_msk, 1, kind)
freqs_del <- get_freq(df_us, df_msk, -1, kind)
freqs_amp$event <- 'amp'
freqs_del$event <- 'del'
freqs <- rbind(freqs_amp, freqs_del)

plot <- ggplot(data=freqs, aes_string(x='MSK', y=kind, color='event'))+rasterise(geom_point(alpha=0.5, size=0.1), dpi=300)+geom_smooth(method='lm', size=0.2)+
        scale_color_manual(values=c('#c84440','#185492'))

plotbis <- function(plot, theme_unmute, theme_mute, name, h=2.5, w=2.5, units='in', dpi=300) {
  unmute <- plot + theme_unmute  + theme(legend.position="none")
  ggsave(filename=name, plot=unmute, height=h, width=w, units=units, dpi=300)
  ext <- substr(name, nchar(name)-3, nchar(name)) 
  # this works only with 3 char extensions, TODO FIXME  https://stackoverflow.com/questions/29113973/get-filename-without-extension-in-r
  #name_mute <- paste0(substr(name, 0, nchar(name)-3), 'mute', ext)
  #mute <- plot + theme_mute  + theme(legend.position="none")
  #ggsave(filename=name_mute, plot=mute, height=h, width=w, units=units, dpi=300)
}

#plotbis(plot=plot, theme_unmute=unmute_theme, theme_mute=mute_theme, name=corrplot_f)
plotbis(plot=plot, theme_unmute=unmute_theme, theme_mute=mute_theme, name=corrplot_f, h=8, w=8)

cor_amp <- cor.test(freqs_amp[,kind], freqs_amp$MSK)
cor_del <- cor.test(freqs_del[,kind], freqs_del$MSK)

res <- data.frame(event=c('amp','del'), pearson=c(cor_amp$estimate, cor_del$estimate), pval=c(cor_amp$p.value, cor_del$p.value))
write.table(res, file=corrs_f, sep='\t', row.names=FALSE, quote=FALSE)