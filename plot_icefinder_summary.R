#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: plot_icefinder_summary.R COUNTS CONTEXTS OUTDIR")
counts <- read.delim(args[1], check.names = FALSE)
contexts <- read.delim(args[2], check.names = FALSE)
outdir <- args[3]

palette <- c(ICE="#E6864A", IME="#F2B179", rICE="#C46A76", AICE="#8D6A9F",
             Other="#7B8794", `non-ICE`="#9EAFBF")
types <- intersect(c("ICE", "IME", "rICE", "AICE", "Other"), names(counts))
counts$rank <- rank(-counts$total, ties.method="first")
long <- reshape(counts[c("genome_id", "rank", types)], varying=types,
                v.names="count", timevar="element_type", times=types,
                direction="long")
long$element_type <- factor(long$element_type, levels=types)

p1 <- ggplot(long, aes(rank, count, fill=element_type)) +
  geom_col(width=1, linewidth=0) +
  scale_fill_manual(values=palette, drop=FALSE) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=expansion(mult=c(0,0.04))) +
  labs(title="ICEfinder-opt predictions across E. faecium genomes",
       subtitle=paste0(sum(counts$total > 0), " of ", nrow(counts),
                       " genomes carry at least one predicted element"),
       x="Genomes ranked by total element count", y="Predicted elements", fill=NULL) +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold", size=16),
        plot.subtitle=element_text(size=11, colour="#4F4F4F"),
        legend.position="top", axis.line=element_line(colour="#303030"))

distribution <- as.data.frame(table(total=counts$total), stringsAsFactors=FALSE)
distribution$total <- as.integer(as.character(distribution$total))
p2 <- ggplot(distribution, aes(total, Freq)) +
  geom_col(fill="#E6864A", width=0.78) +
  scale_x_continuous(breaks=sort(unique(distribution$total))) +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  labs(x="Predicted elements per genome", y="Number of genomes") +
  theme_classic(base_size=12) +
  theme(axis.line=element_line(colour="#303030"))

type_totals <- data.frame(
  element_type=factor(types, levels=types),
  count=vapply(types, function(x) sum(counts[[x]]), numeric(1))
)
p_type <- ggplot(type_totals, aes(element_type, count, fill=element_type)) +
  geom_col(width=0.68, show.legend=FALSE) +
  geom_text(aes(label=format(count, big.mark=",")), vjust=-0.45, size=4) +
  scale_fill_manual(values=palette, drop=FALSE) +
  scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
  labs(x="Element type", y="Total predictions") +
  theme_classic(base_size=12) +
  theme(axis.line=element_line(colour="#303030"))

combined <- p1 / (p_type | p2) + plot_layout(heights=c(2.1,1))
ggsave(file.path(outdir, "icefinder_elements_per_genome.pdf"), combined,
       width=11.2, height=8.2, units="in", device="pdf", useDingbats=FALSE)
ggsave(file.path(outdir, "icefinder_elements_per_genome.png"), combined,
       width=11.2, height=8.2, units="in", dpi=320, bg="white")

context_types <- intersect(c("ICE", "IME", "rICE", "AICE", "Other", "non-ICE"), names(contexts))
ctx <- reshape(contexts[c("hotspot_id", context_types)], varying=context_types,
               v.names="count", timevar="context", times=context_types,
               direction="long")
ctx$context <- factor(ctx$context, levels=context_types)
ctx$hotspot_id <- factor(ctx$hotspot_id, levels=rev(contexts$hotspot_id))
p3 <- ggplot(ctx, aes(count, hotspot_id, fill=context)) +
  geom_col(position="fill", width=0.68) +
  scale_fill_manual(values=palette, drop=FALSE) +
  scale_x_continuous(labels=function(x) paste0(round(100*x), "%"), expand=c(0,0)) +
  labs(title="ICE/IME context of existing Phorific hotspots",
       x="Percentage of hotspot instances", y="Hotspot", fill=NULL) +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold", size=16), legend.position="top",
        axis.line=element_line(colour="#303030"))
ggsave(file.path(outdir, "existing_hotspot_ice_context.pdf"), p3,
       width=10.5, height=6.2, units="in", device="pdf", useDingbats=FALSE)
ggsave(file.path(outdir, "existing_hotspot_ice_context.png"), p3,
       width=10.5, height=6.2, units="in", dpi=320, bg="white")
