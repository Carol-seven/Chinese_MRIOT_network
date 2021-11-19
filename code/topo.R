#########################
#### Load R packages ####
#########################

library(circlize) ## generate chord graphs
library(ggplot2)  ## generate graphs (for paper)
library(ggpubr) ## arrange multiple ggplots 
library(igraph)  ## generate networks (for analysis)
library(RColorBrewer) ## generate a long color palette
library(tidyr) ## gather columns in the data frame into key-value pairs
## install the R package for weighted directed network
## library(devtools) ## install gitlab packages
## devtools::install_git(url = "https://gitlab.com/wdnetwork/wdnet.git")
library(wdnet) ## generate results of weighted directed network
library(xtable)  ## generate latex tables

#####################
#### Load .Rdata ####
#####################

file_names <- as.list(list.files(file.path("..", "data")))

## number of sectors per province
sec <- 30

## number of province
prov <- 30

## total number of nodes
tot <- sec*prov

## year list
year_list <- c(2007L, 2012L)
years <- length(year_list)

## create an array for network data
mriot_data <- array(NA, dim = c(tot, tot, years))

## create an array for total value added
tva <- array(NA, dim = c(tot, years))

## create an array for total final use
tfu <- array(NA, dim = c(tot, years))

## create an array for gross output/input
go <- array(NA, dim = c(tot, years))

## inflation rate from the World Bank
GDPdef <- c(1, 1.2715)

for (i in 1:years) {
  load(paste0(file.path("..", "data", file_names[[i]])))
  ## the intermediate flow matrix
  mriot_data[,,i] <- as.matrix(mriot[1:tot, 1:tot]/GDPdef[i])
  ## total value added
  tva[,i] <- as.numeric(mriot[which(rownames(mriot) == "Total value added"), 1:tot]/GDPdef[i])
  ## total final use
  tfu[,i] <- as.numeric(mriot[1:tot, which(colnames(mriot) == "Total.final.use")]/GDPdef[i])
  ## total output
  go[,i] <- as.numeric(mriot[1:tot, which(colnames(mriot) == "Total.output")]/GDPdef[i])
}

## indicator matrix: a block diagonal matrix, whose main diagonal blocks are ones matrices 
## and off-diagonal blocks are zeros matrices.
indicator_matrix <- array(0, dim = dim(mriot_data))
for (i in 1:years) {
  for (j in 1:prov) {
    indicator_matrix[((j - 1)*sec + 1):(j*sec), ((j - 1)*sec + 1):(j*sec), i] <- 1
  }
}

## the data at the intra-province level
mriot_data_intra <- mriot_data * indicator_matrix

## the data at the inter-province level
mriot_data_inter <- mriot_data * (1 - indicator_matrix)

province_name <- unique(gsub("[^A-z]", "", colnames(mriot)[1:tot]))

industry_list <- paste(gsub("[^A-z]", "", colnames(mriot)[1:tot]), 
                       formatC(as.numeric(gsub("\\D", "", colnames(mriot)[1:tot])), 
                               width = nchar(sec), flag = "0"))

##----------------------------------------------------------------------------------------

source("utility.R")
## the script contains the following functions: 
## backbone of a network (netbackbone), degree and strength distribution (gen_hist), 
## tail distribution (tail_dist), weighted PageRank (wpr).

if (! dir.exists("../image")) dir.create("../image")

##----------------------------------------------------------------------------------------

######################
#### Chord graphs ####
######################

## function matrix_split: split a matrix into many row-by-column submatrices.
## M: the original matrx
## r: number of rows in sub-matrix format
## c: number of columns in sub-matrix format

matrix_split <- function(M, r, c, byrow = FALSE) {
  r_tag <- (row(M)-1) %/% r + 1
  c_tag <- (col(M)-1) %/% c + 1
  if (byrow == TRUE) {
    m_tag <- (r_tag-1) * max(c_tag) + c_tag
  } else {
    m_tag <- (c_tag-1) * max(r_tag) + r_tag
  }
  dim3 <- prod(dim(M))/r/c
  result <- unlist(lapply(1:dim3, function(x) {M[m_tag == x]}))
  dim(result) <- c(r, c, dim3)
  return(result)
}

## aggregated data for sectors

mriot_sec <- array(NA, dim = c(sec, sec, years))

for (i in 1:years) {
  mriot_sec[,,i] <- as.matrix(apply(matrix_split(mriot_data[,,i], sec, sec), c(1,2), sum))
  mriot_sec[,,i] <- mriot_sec[,,i] - diag(diag(mriot_sec[,,i])) ## remove self-loops
}

colnames(mriot_sec) <- rownames(mriot_sec) <- formatC(1:sec, width = nchar(sec), flag = "0")

## aggregated data for provinces

mriot_prov <- array(NA, dim = c(prov, prov, years))

for (z in 1:years) {
  for (i in c(1:prov)) {
    for (j in c(1:prov)) {
      mriot_prov[i, j, z] <- sum(mriot_data[((i-1)*sec + 1):(i*sec), ((j-1)*sec + 1):(j*sec), z])
    }
  }
  mriot_prov[,,z] <- mriot_prov[,,z] - diag(diag(mriot_prov[,,z])) ## remove self-loops
}

colnames(mriot_prov) <- rownames(mriot_prov) <- formatC(1:prov, width = nchar(prov), flag = "0")

## function gen_chord: generate the chord graph

gen_chord <- function(mymat, mycolor) {
  par(mar = c(0, 0, 1, 0))
  circos.par(gap.after = rep(3, nrow(mymat)))
  chordDiagram(mymat, grid.col = mycolor, directional = 1, link.zindex = rank(mymat), 
               diffHeight = mm_h(3.5), target.prop.height = mm_h(2))
}

## the unit of economic flow is set at 1e7 * 10,000 = 100 billion Chinese Yuan (CNY).

mriot_sec <- mriot_sec / 1e7
mriot_prov <- mriot_prov / 1e7

## my_sec_color <- rand_color(sec, luminosity = "bright") ## the code can generate random colors
## my_prov_color <- rand_color(prov, luminosity = "bright")

## mycolor for inter-sectoral flow in the paper

my_sec_color <- c("#F830E7FF", "#F305CFFF", "#13DEAEFF", "#FDFA50FF", "#6DBF37FF", 
                  "#EE4699FF", "#A7E019FF", "#0C8037FF", "#4444E4FF", "#460676FF", 
                  "#62CF52FF", "#4B74ACFF", "#1C0B6EFF", "#FC028EFF", "#F79F30FF", 
                  "#FFD567FF", "#EACC5FFF", "#29CEADFF", "#EA5EBFFF", "#9D4DFCFF", 
                  "#070182FF", "#D65EF5FF", "#36BE51FF", "#79E909FF", "#0D8B38FF", 
                  "#FB6675FF", "#E63E8BFF", "#E3A352FF", "#C47C32FF", "#38819BFF")

## mycolor for inter-provincial flow in the paper

my_prov_color <- c("#B9D13DFF", "#1AA679FF", "#28D8BFFF", "#E55F06FF", "#962F04FF", 
                   "#E1844EFF", "#A4E54EFF", "#C05B3EFF", "#619514FF", "#104BD8FF", 
                   "#B02B42FF", "#EB7F44FF", "#DF31A5FF", "#C11915FF", "#F4B447FF", 
                   "#EC3E11FF", "#BF6CFCFF", "#AE2233FF", "#117B98FF", "#1C5099FF", 
                   "#7358C8FF", "#F72463FF", "#E046EDFF", "#C960E2FF", "#1D2979FF", 
                   "#44DD07FF", "#54E078FF", "#F747E6FF", "#4F2BDAFF", "#EF5A22FF")

## chord visualization

pdf(file.path("../image", "chord.pdf"), width = 10, height = 10)

par(mfrow = c(2, 2))

gen_chord(mriot_sec[,,1], my_sec_color)
title("inter-sectoral flow in 2007")
circos.clear()

gen_chord(mriot_sec[,,2], my_sec_color)
title("inter-sectoral flow in 2012")
circos.clear()

gen_chord(mriot_prov[,,1], my_prov_color)
title("inter-provincial flow in 2007")
circos.clear()

gen_chord(mriot_prov[,,2], my_prov_color)
title("inter-provincial flow in 2012")
circos.clear()

dev.off()

##----------------------------------------------------------------------------------------

########################################
#### Backbone network visualization ####
########################################

## regional GDP: total value added based on provinces
tva_prov <- array(NA, dim = c(prov, years))
for (i in 1:years) {
  for (j in 1:prov) {
    tva_prov[j,i] <- sum(tva[((j-1)*sec+1):(j*sec),i])
  }
}

## target provinces: the provinces with top 5 regional GDP in 2012
prov_target <- province_name[order(tva_prov[,2], decreasing = TRUE)[1:5]]

## the indexes for the target provinces
province.index <- which(province_name %in% prov_target)

## the indexes for the data of the target provinces
select.index <- as.vector(sapply((province.index - 1)*sec + 1, function(x) {seq(x, x + sec - 1)}))

## function gen_bn: the backbone network visualization from the data of target provinces, 
## at the significance level of alpha

gen_bn <- function(alpha) {
  par(mar = c(0, 0, 1.1, 0))
  ## direction is not considered; weight is preserved
  mriot_data_2012_ex <- netbackbone(mriot_data[,,2], alpha = alpha, weighted = TRUE, directed = FALSE)
  ## data for the target provinces
  mriot_adj2012_ex <- mriot_data_2012_ex$AdjacencyMatrix[select.index, select.index]
  ## Remove self-loops
  mriot_adj2012 <- mriot_adj2012_ex - diag(diag(mriot_adj2012_ex), sec*length(province.index))
  mriot_graph2012 <- graph_from_adjacency_matrix(mriot_adj2012, weighted = TRUE)
  E(mriot_graph2012)$arrow.mode <- 0
  V(mriot_graph2012)$size <- 1.3 * log(degree(mriot_graph2012, mode = "all") + 1.3)
  V(mriot_graph2012)$label <- NA
  for(i in 1:length(province.index)) {
    V(mriot_graph2012)[((i - 1)*sec + 1):(i*sec)]$color <- colors()[10*i]
  }
  ## Remove disconnected nodes in the graph
  mriot2012_clust <- clusters(mriot_graph2012, mode = "weak")
  mriot2012_cc <- induced.subgraph(mriot_graph2012, 
                                   V(mriot_graph2012)[which(mriot2012_clust$membership 
                                                            == which.max(mriot2012_clust$csize))])
  plot(mriot2012_cc, layout = layout_with_graphopt(mriot2012_cc))
}

## an example plot of the backbone MRIOT in 2012

pdf(file.path("../image", "backbone2012.pdf"), width = 7, height = 4)

layout(mat = matrix(c(1,2,3,3), nrow = 2, ncol = 2, byrow = TRUE), heights = c(5, 1))

gen_bn(alpha = 10^-3)
title(expression(delta == 10^{-3}))

gen_bn(alpha = 10^-4)
title(expression(delta == 10^{-4}))

## common legend
par(mar = c(1, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("bottom", inset = -0.01, 
       title = "Provinces", prov_target, horiz = TRUE, 
       pch = 19, col = c(colors()[10*(1:5)]))

dev.off()

##----------------------------------------------------------------------------------------

######################################
#### Degree/Strength distribution ####
######################################

###################
#### Histogram ####
###################

pdf(file.path("../image", "degdist.pdf"), height = 5)

gen_hist(mriot_data, c("2007", "2012"), mode = "degree")

dev.off()

pdf(file.path("../image", "strdist.pdf"), height = 5)

gen_hist(mriot_data, c("2007", "2012"), mode = "strength")

dev.off()

############################
#### Tail distributions ####
############################

## generate strength distributions

out_str2007 <- tail_dist(mriot_data[,,1], type = "out-strength")
out_str2012 <- tail_dist(mriot_data[,,2], type = "out-strength")
in_str2007 <- tail_dist(mriot_data[,,1], type = "in-strength")
in_str2012 <- tail_dist(mriot_data[,,2], type = "in-strength")
tot_str2007 <- tail_dist(mriot_data[,,1], type = "total-strength")
tot_str2012 <- tail_dist(mriot_data[,,2], type = "total-strength")

## function gen_df: generate the estimated parameters for power-law tails.
## input: the strength distribution
## thres: the threshold

gen_df <- function(input, thres = 0, type, year) {
  myinput <- input[input > thres]
  mystr <- sort(myinput)
  mytail <- 1 - ecdf(myinput)(sort(myinput))
  n <- length(myinput)
  mytype <- rep(type, n)
  myyear <- rep(year, n)
  if (thres == 0){
    datatype <- "full data"
  } else {
    datatype <- "tail"
  }
  mydatatype <- rep(datatype, n)
  res <- data.frame(mystr, mytail, mytype, myyear, mydatatype)
  return(res)
}

mystr_df <- rbind(gen_df(out_str2007$res, type = "out-strength", year = "2007"), 
                  gen_df(out_str2012$res, type = "out-strength", year = "2012"), 
                  gen_df(in_str2007$res, type = "in-strength", year = "2007"),
                  gen_df(in_str2012$res, type = "in-strength", year = "2012"),
                  gen_df(tot_str2007$res, type = "total-strength", year = "2007"),
                  gen_df(tot_str2012$res, type = "total-strength", year = "2012"),
                  gen_df(out_str2007$res, thres = out_str2007$xmin, type = "out-strength", year = "2007"), 
                  gen_df(out_str2012$res, thres = out_str2012$xmin, type = "out-strength", year = "2012"), 
                  gen_df(in_str2007$res, thres = in_str2007$xmin, type = "in-strength", year = "2007"),
                  gen_df(in_str2012$res, thres = in_str2012$xmin, type = "in-strength", year = "2012"),
                  gen_df(tot_str2007$res, thres = tot_str2007$xmin, type = "total-strength", year = "2007"),
                  gen_df(tot_str2012$res, thres = tot_str2012$xmin, type = "total-strength", year = "2012")
                  )
## xmin is the optimal lower cutoff obtained from a goodness-of-fit based approach

## remove zeros
rmzero <- intersect(which(mystr_df$mystr > 0), which(mystr_df$mytail > 0))

## plot tail distributions
## arrange two ggplot2 (facet_grid) graphs

fig1 <- ggplot(mystr_df[intersect(rmzero, which(mystr_df$mydatatype == "full data")),], 
               aes(x = mystr, y = mytail, color = myyear)) +
  geom_point(size = 0.7) +
  facet_grid(mydatatype~mytype, scales = "free_x") +
  scale_x_continuous(breaks = scales::log_breaks(), trans = scales::log_trans()) +
  scale_y_continuous(breaks = scales::log_breaks(), trans = scales::log_trans(),
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_colour_discrete("year") +
  xlab(NULL) + ylab(NULL) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        strip.background = element_rect(fill = "gray90", color = "black"))

fig2 <- ggplot(mystr_df[intersect(rmzero, which(mystr_df$mydatatype == "tail")),], 
               aes(x = mystr, y = mytail, color = myyear)) +
  geom_point(size = 0.7) +
  facet_grid(mydatatype~mytype, scales = "free_x") +
  scale_x_continuous(breaks = scales::log_breaks(), trans = scales::log_trans()) +
  scale_y_continuous(breaks = scales::log_breaks(), trans = scales::log_trans(),
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_colour_discrete("year") +
  xlab("strength") + ylab(NULL) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_rect(fill = "gray90", color = "black"))

figure <- ggarrange(fig1, fig2, ncol = 1, nrow = 2, heights = c(1.01, 1), 
                    common.legend = TRUE, legend = "bottom")

pdf(file.path("../image", "tail_str.pdf"), height = 5.5)

annotate_figure(figure, 
                left = text_grob("tail probability", rot = 90, size = 11))

dev.off()

## table of the the estimated parameters for power-law tails
## xmin: the threshold based on a goodness-of-fit test
## alpha: the exponent based on a goodness-of-fit test
## pvalue: the p-value of the power law tail beyond the estimated threshold obtained from bootstrapping

pl_df <- as.data.frame(
  rbind(cbind(in_str2007$xmin, out_str2007$xmin, tot_str2007$xmin, 
              in_str2012$xmin, out_str2012$xmin, tot_str2012$xmin) / 1e6,
        cbind(in_str2007$alpha, out_str2007$alpha, tot_str2007$alpha, 
              in_str2012$alpha, out_str2012$alpha, tot_str2012$alpha),
        cbind(in_str2007$pvalue, out_str2007$pvalue, tot_str2007$pvalue, 
              in_str2012$pvalue, out_str2012$pvalue, tot_str2012$pvalue)
))

rownames(pl_df) <- c("threshold", "exponent", "pvalue")

colnames(pl_df) <- c("in_str2007", "out_str2007", "tot_str2007", 
                     "in_str2012", "out_str2012", "tot_str2012")

print(xtable(pl_df, digits = c(0, rep(2, 6))))

##----------------------------------------------------------------------------------------

#######################
#### Assortativity ####
#######################

gen_assort <- function(data = mriot_data, alpha = 1) {
  
  ## generate weighted and undirected networks
  data_ud <- array(NA, dim = dim(data))
  data_ud[] <- apply(data, 3, function(M) {
    temp <- netbackbone(M, alpha = alpha, weighted = TRUE, directed = FALSE)
    return(temp$AdjacencyMatrix)})
  
  ## generate unweighted and directed networks
  data_uw <- array(NA, dim = dim(data))
  data_uw[] <- apply(data, 3, function(M) {
    temp <- netbackbone(M, alpha = alpha, weighted = FALSE, directed = TRUE)
    return(temp$AdjacencyMatrix)})
  
  ## generate unweighted and undirected networks
  data_uwd <- array(NA, dim = dim(data))
  data_uwd[] <- apply(data, 3, function(M) {
    temp <- netbackbone(M, alpha = alpha, weighted = FALSE, directed = FALSE)
    return(temp$AdjacencyMatrix)})
  
  ## unweighted assortativity
  assort_total_uw <- apply(data_uwd, 3, wdnet::dw_assort, type = "out-in")
  assort_outin_uw <- apply(data_uw, 3, wdnet::dw_assort, type = "out-in")
  assort_inin_uw <- apply(data_uw, 3, wdnet::dw_assort, type = "in-in")
  assort_outout_uw <- apply(data_uw, 3, wdnet::dw_assort, type = "out-out")
  assort_inout_uw <- apply(data_uw, 3, wdnet::dw_assort, type = "in-out")
  
  ## weighted assortativity
  assort_total <- apply(data_ud, 3, wdnet::dw_assort, type = "out-in")
  assort_outin <- apply(data, 3, wdnet::dw_assort, type = "out-in")
  assort_inin <- apply(data, 3, wdnet::dw_assort, type = "in-in")
  assort_outout <- apply(data, 3, wdnet::dw_assort, type = "out-out")
  assort_inout <- apply(data, 3, wdnet::dw_assort, type = "in-out")
  
  result <- rbind(assort_inin_uw, assort_inout_uw, assort_outin_uw, assort_outout_uw, assort_total_uw,
                  assort_inin, assort_inout, assort_outin, assort_outout, assort_total)
  colnames(result) <- year_list
  return(result)
}

## --------------------------------
##              2007
## --------------------------------
##  Overall     Intra       Inter    
## ---------- ---------- ---------- 
##   Uw  W      Uw  W       Uw  W
## --------------------------------

gen_assortdf <- function(year) {
  y <- which(year_list %in% year == TRUE)
  myvalue <- rbind(gen_assort(mriot_data)[,y],
                   gen_assort(mriot_data_intra)[,y],
                   gen_assort(mriot_data_inter)[,y])
  dim(myvalue) <- c(5, prod(dim(myvalue))/5)
  mydf <- data.frame(myvalue)
  rownames(mydf) <- c("in-in", "in-out", "out-in", "out-out", "total")
  return(mydf)
}

xtable(gen_assortdf(year_list), digits = c(0, rep(3, length(year_list)*3*2)))

## assortativity coefficients at a sequence of significance levels

gen_assort_seq <- function(x) {
  data <- mriot_data
  
  ## generate weighted and directed networks at significance level x
  data_x <- array(NA, dim = dim(data))
  data_x[] <- apply(data, 3, function(M) {
    temp <- netbackbone(M, alpha = x, weighted = TRUE, directed = TRUE)
    return(temp$AdjacencyMatrix)})
  
  ## generate weighted and undirected networks at significance level x
  data_ud_x <- array(NA, dim = dim(data))
  data_ud_x[] <- apply(data, 3, function(M) {
    temp <- netbackbone(M, alpha = x, weighted = TRUE, directed = FALSE)
    return(temp$AdjacencyMatrix)})
  
  ## weighted assortativity
  assort_total <- apply(data_ud_x, 3, wdnet::dw_assort, type = "out-in")
  assort_outin <- apply(data_x, 3, wdnet::dw_assort, type = "out-in")
  assort_inin <- apply(data_x, 3, wdnet::dw_assort, type = "in-in")
  assort_outout <- apply(data_x, 3, wdnet::dw_assort, type = "out-out")
  assort_inout <- apply(data_x, 3, wdnet::dw_assort, type = "in-out")
  
  result <- rbind(assort_inin, assort_inout, assort_outin, assort_outout, assort_total)
  return(result)
}

assort_seq <- mapply(gen_assort_seq, seq(0.001, 0.999, 0.001))
## Considering the long running time, the data has been saved.
## load("assort_seq.RData")

assort_seq <- data.frame(assort_seq)

colnames(assort_seq) <- seq(0.001, 0.999, 0.001)

assort_seq <- data.frame(type = rep(c("in-in", "in-out", "out-in", "out-out", "total")), 
                         gather(assort_seq, "alpha", "assort"), 
                         year = rep(year_list, each = 5))

assort_seq$alpha <- as.numeric(assort_seq$alpha)

pdf(file.path("../image", "assort_delta.pdf"), height = 4)

ggplot(assort_seq, aes(x = alpha, y = assort, color = type, group = type)) + 
  geom_line() + 
  xlab("delta") + 
  ylab("assortativity") + 
  facet_wrap(.~year) + 
  theme(panel.spacing = unit(1, "lines"), 
        legend.position = "bottom", 
        legend.spacing.x = unit(0.5, "cm"), 
        strip.background = element_rect(fill = "grey60", color = "grey60"))

dev.off()

##----------------------------------------------------------------------------------------

#########################
#### Clustering coef ####
#########################

cc_total <- apply(mriot_data, 3, function(M) {
  temp <- wdnet::dw_clustcoeff(M, method = "Clemente", isolates = "nan")$total
  return(temp)
})

cc_total_intra <- apply(mriot_data_intra, 3, function(M) {
  temp <- wdnet::dw_clustcoeff(M, method = "Clemente", isolates = "nan")$total
  return(temp)
})

cc_total_inter <- apply(mriot_data_inter, 3, function(M) {
  temp <- wdnet::dw_clustcoeff(M, method = "Clemente", isolates = "nan")$total
  return(temp)
})

cc_in <- apply(mriot_data, 3, function(M) {
  temp <- wdnet::dw_clustcoeff(M, method = "Clemente", isolates = "nan")$"in"
  return(temp)
})

cc_in_intra <- apply(mriot_data_intra, 3, function(M) {
  temp <- wdnet::dw_clustcoeff(M, method = "Clemente", isolates = "nan")$"in"
  return(temp)
})

cc_in_inter <- apply(mriot_data_inter, 3, function(M) {
  temp <- wdnet::dw_clustcoeff(M, method = "Clemente", isolates = "nan")$"in"
  return(temp)
})

cc_out <- apply(mriot_data, 3, function(M) {
  temp <- wdnet::dw_clustcoeff(M, method = "Clemente", isolates = "nan")$out
  return(temp)
})

cc_out_intra <- apply(mriot_data_intra, 3, function(M) {
  temp <- wdnet::dw_clustcoeff(M, method = "Clemente", isolates = "nan")$out
  return(temp)
})

cc_out_inter <- apply(mriot_data_inter, 3, function(M) {
  temp <- wdnet::dw_clustcoeff(M, method = "Clemente", isolates = "nan")$out
  return(temp)
})

cc_cyc <- apply(mriot_data, 3, function(M) {
  temp <- wdnet::dw_clustcoeff(M, method = "Clemente", isolates = "nan")$cycle
  return(temp)
})

cc_cyc_intra <- apply(mriot_data_intra, 3, function(M) {
  temp <- wdnet::dw_clustcoeff(M, method = "Clemente", isolates = "nan")$cycle
  return(temp)
})

cc_cyc_inter <- apply(mriot_data_inter, 3, function(M) {
  temp <- wdnet::dw_clustcoeff(M, method = "Clemente", isolates = "nan")$cycle
  return(temp)
})

cc_mid <- apply(mriot_data, 3, function(M) {
  temp <- wdnet::dw_clustcoeff(M, method = "Clemente", isolates = "nan")$middle
  return(temp)
})

cc_mid_intra <- apply(mriot_data_intra, 3, function(M) {
  temp <- wdnet::dw_clustcoeff(M, method = "Clemente", isolates = "nan")$middle
  return(temp)
})

cc_mid_inter <- apply(mriot_data_inter, 3, function(M) {
  temp <- wdnet::dw_clustcoeff(M, method = "Clemente", isolates = "nan")$middle
  return(temp)
})

## global clustering coefficients
## generate data frame for 2007 and 2012

gen_globalcc <- function(year) {
  myind <- which(year_list %in% year == TRUE)
  mydf <- cbind(sapply(myind, function(x){return(cc_total[[x]]$globalcc)}), 
                sapply(myind, function(x){return(cc_total_intra[[x]]$globalcc)}), 
                sapply(myind, function(x){return(cc_total_inter[[x]]$globalcc)}),
                sapply(myind, function(x){return(cc_cyc[[x]]$globalcc)}),
                sapply(myind, function(x){return(cc_cyc_intra[[x]]$globalcc)}),
                sapply(myind, function(x){return(cc_cyc_inter[[x]]$globalcc)}),
                sapply(myind, function(x){return(cc_mid[[x]]$globalcc)}),
                sapply(myind, function(x){return(cc_mid_intra[[x]]$globalcc)}),
                sapply(myind, function(x){return(cc_mid_inter[[x]]$globalcc)}),
                sapply(myind, function(x){return(cc_in[[x]]$globalcc)}),
                sapply(myind, function(x){return(cc_in_intra[[x]]$globalcc)}),
                sapply(myind, function(x){return(cc_in_inter[[x]]$globalcc)}),
                sapply(myind, function(x){return(cc_out[[x]]$globalcc)}),
                sapply(myind, function(x){return(cc_out_intra[[x]]$globalcc)}),
                sapply(myind, function(x){return(cc_out_inter[[x]]$globalcc)}))
  return(t(mydf))
}

globalcc_mat <- matrix(NA, nrow = 5, ncol = (3*length(year_list)))

for (i in 1:5){
  temp <- gen_globalcc(year_list)
  globalcc_mat[i,] <- t(as.vector(temp[c((3*(i - 1) + 1):(3*i)),]))
}

globalcc_df <- as.data.frame(globalcc_mat)

rownames(globalcc_df) <- c("total", "cyc", "mid", "in", "out")

print(xtable(globalcc_df, digits = c(0, rep(3, 6))))

## local clustering coefficients

gen_localcc <- function(year, top.set = 5, bottom.set = 5) {
  myind <- which(year_list %in% year == TRUE)
  local_tot <- rbind(
    sapply(myind, function(x){return(order(cc_total[[x]]$localcc, decreasing = TRUE)[1:top.set])}), 
    sapply(myind, function(x){return(order(cc_total[[x]]$localcc)[1:bottom.set])}))
  local_cyc <- rbind(
    sapply(myind, function(x){return(order(cc_cyc[[x]]$localcc, decreasing = TRUE)[1:top.set])}), 
    sapply(myind, function(x){return(order(cc_cyc[[x]]$localcc)[1:bottom.set])}))
  local_mid <- rbind(
    sapply(myind, function(x){return(order(cc_mid[[x]]$localcc, decreasing = TRUE)[1:top.set])}), 
    sapply(myind, function(x){return(order(cc_mid[[x]]$localcc)[1:bottom.set])}))
  local_in <- rbind(
    sapply(myind, function(x){return(order(cc_in[[x]]$localcc, decreasing = TRUE)[1:top.set])}), 
    sapply(myind, function(x){return(order(cc_in[[x]]$localcc)[1:bottom.set])}))
  local_out <- rbind(
    sapply(myind, function(x){return(order(cc_out[[x]]$localcc, decreasing = TRUE)[1:top.set])}), 
    sapply(myind, function(x){return(order(cc_out[[x]]$localcc)[1:bottom.set])}))
  mydf <- cbind(local_tot, local_cyc, local_mid, local_in, local_out)
  return(mydf)
}

localcc_df <- gen_localcc(year_list)

## nrow = top.set + bottom.set
localccname_df <- as.data.frame(matrix(industry_list[localcc_df], nrow = 10))

## only present top 5 in the manuscript
localccname_df_toponly <- localccname_df[c(1:5),]

localccname_df_topresent <- as.data.frame(
  rbind(as.matrix(localccname_df_toponly[,seq(1,9,2)]), 
        as.matrix(localccname_df_toponly[,seq(2,10,2)])))

print(xtable(localccname_df_topresent), include.rownames = FALSE)

##----------------------------------------------------------------------------------------

#############################
#### Community detection ####
#############################

## function gen_code: the community detection based on modularity maximization, and 
## province-sectors with the same code belong to the same community.

gen_code <- function(data = mriot_data, year) {
  data_ud <- array(NA, dim = c(tot, tot, years))
  data_ud[] <- apply(data, 3, function(M) {
    temp <- netbackbone(M, alpha = 1, weighted = TRUE, directed = FALSE)
    return(temp$AdjacencyMatrix)})
  y <- which(year_list == year)
  data_graph <- igraph::graph_from_adjacency_matrix(data_ud[,,y], mode = "undirected", weighted = TRUE)
  clust_res <- igraph::fastgreedy.community(data_graph)
  row.temp <- formatC(1:sec, width = nchar(sec), flag = "0")
  col.temp <- province_name
  data.temp <- expand.grid(sector = row.temp, province = col.temp)
  data.temp$Z <- as.factor(clust_res$membership)
  data.temp$year <- year
  return(data.temp)
}

code_df <- rbind(gen_code(mriot_data, 2007L), gen_code(mriot_data, 2012L))

## a long color palette

getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
colorCount <- length(unique(code_df$Z))
myPalette <- getPalette(colorCount)

pdf(file.path("../image", "cluster.pdf"), width = 13)

ggplot(code_df, aes(x = sector, y = province, fill = Z)) + 
  geom_tile(colour = "white", size = 0.5) + 
  facet_grid(.~year) + 
  scale_fill_manual(values = myPalette) + 
  theme(strip.background = element_rect(fill = "gray90", color = "black"), 
        panel.spacing = unit(0.7, "lines"), 
        axis.text.x = element_text(size = 7), legend.position = "NULL")

dev.off()

##----------------------------------------------------------------------------------------

#############################
#### Centrality analysis ####
#############################

## function gettop: get the top n province-sectors by sorting the scores.

gettop <- function(score, n) {
  score_ordered <- order(score, decreasing = TRUE)[1:n]
  return(industry_list[score_ordered])
}

## function pr_top: province-sectors with top "n" PageRank scores.
## prior: "tva" = using total value added as prior information
## "no" = no prior information; uniformly over all the nodes

pr_top <- function(data = mriot_data, theta = 1, prior = c("no", "tva"), 
                   year = year_list, n = 10) {
  if (prior == "no") {
    pr <- apply(data, 3, function(M) {
      temp <- wpr(M, theta = theta)[,2]
      return(temp)
    })
  } else if (prior == "tva") {
    pr <- array(NA, dim = c(tot, years))
    for (i in 1:years) {
      pr[,i] <- wpr(data[,,i], theta = theta, prior.info = tva[,i])[,2]
    }
  }
  y <- which(year_list %in% year == TRUE)
  mydf <- sapply(y, function(x) {gettop(pr[,x], n)})
  colnames(mydf) <- year
  return(mydf)
}

## --------------------------------------------------------------
##       theta = 0            theta = 1           theta = 1
##   (no prior info.)     (no prior info.)           (TVA)
## -------------------- -------------------- --------------------
##  year2007  year2012   year2007  year2012   year2007  year2012
## -------------------- -------------------- --------------------

pr_df <- cbind(pr_top(theta = 0, prior = "no"), 
               pr_top(prior = "no"), 
               pr_top(prior = "tva"))

xtable(pr_df)

## Generate network at the significance levels 0.001 and 0.0001.

mriot_data0.001 <- mriot_data0.0001 <-array(NA, dim = dim(mriot_data))

mriot_data0.001[] <- apply(mriot_data, 3, function(M) {
  temp <- netbackbone(M, alpha = 0.001, weighted = TRUE, directed = TRUE)
  return(temp$AdjacencyMatrix)})

mriot_data0.0001[] <- apply(mriot_data, 3, function(M) {
  temp <- netbackbone(M, alpha = 0.0001, weighted = TRUE, directed = TRUE)
  return(temp$AdjacencyMatrix)})

## ---------------------
##         alpha
## ---------------------
##  year2007   year2012
## ---------------------
## different significant levels: alpha = 1, 0.001, 0.0001.
## where province-sectors with * do not appear in pr_tva (raw weighted network).

pr_tva <- pr_top(prior = "tva")

pr0.001 <- pr_top(data = mriot_data0.001, prior = "tva")
for (i in 1:dim(pr_tva)[2]) {
  index <- which(is.na(match(pr0.001[,i], pr_tva[,i])))
  pr0.001[index, i] <- paste0(pr0.001[index, i], "*")
}

pr0.0001 <- pr_top(data = mriot_data0.0001, prior = "tva")
for (i in 1:dim(pr_tva)[2]) {
  index <- which(is.na(match(pr0.0001[,i], pr_tva[,i])))
  pr0.0001[index, i] <- paste0(pr0.0001[index, i], "*")
}

xtable(cbind(pr_tva, pr0.001, pr0.0001))
