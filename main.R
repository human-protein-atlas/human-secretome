
require(tidyverse)
require(magrittr)

setwd("C:/Users/max.karlsson/Documents/Scilifelab/Projects/Secretome/")

#Create result directory if non-existent:
result_folder <- paste0("./results/", Sys.Date())
dir.create(result_folder, showWarnings = FALSE)

# Function for creating savepath
save_path = function(savename) paste(result_folder, savename, sep = "/")

source("scripts/functions.R")

# ----- Color Palettes -----

location_colors <- 
  c("Blood" = "#B40001",
    "Brain locally" = "#FFDA00",
    "Female locally" = "#F4C1D4",
    "Gastric" = "#127FCA",
    "Intracellular/Membrane" = "#F8A265",
    "Male locally" = "#95D4F5",
    "Matrix" = "#7E6A9B",
    "Other locally" = "#FED480",
    "Secreted - no data" = "#A0A8AA") %>%
    {tibble(atlasloc = names(.),
            color = .)}
function_palette <- 
  c('Immunity,Defense' = '#F564E3', 
    'Enzymes,Enzyme inhibitors' = '#619CFF', 
    'Other' = 'darkgrey', 
    'Coagulation,Complement,Acute phase,APO' = '#F8766D', 
    'Hormones,Growth factors' = '#00BA38', 
    'Cytokines' = '#B79F00', 
    'Receptors,Transport' = '#00BFC4', 
    'No annotated function' = 'black')

function_palette_2 <- setNames(c(ggsci::pal_nejm()(6),
                                 'darkgrey', 
                                 'black'), 
                               c('Immunity,Defense',
                                 'Enzymes,Enzyme inhibitors',
                                 'Coagulation,Complement,Acute phase,APO',
                                 'Hormones,Growth factors',
                                 'Cytokines',
                                 'Receptors,Transport',
                                 'Other',
                                 'No annotated function'))

function_palette_3 <- setNames(c(ggsci::pal_jama()(6),
                                 'darkgrey', 
                                 'black'), 
                               c('Immunity,Defense',
                                 'Enzymes,Enzyme inhibitors',
                                 'Coagulation,Complement,Acute phase,APO',
                                 'Hormones,Growth factors',
                                 'Cytokines',
                                 'Receptors,Transport',
                                 'Other',
                                 'No annotated function'))

elevated.cat.cols <- rev(c("Not detected" = "grey",
                           "Low tissue specificity" = "grey40",
                           "Tissue enhanced" = "#984ea3",
                           "Group enriched" = "#FF9D00",
                           "Tissue enriched" = "#e41a1c"))

# ----- load data ----- 
secretome_genes <- 
  read_delim("data/2641_ensg_category.txt", delim = "\t") %>%
  left_join(location_colors, by = "atlasloc")


tissue_hierarchy <- 
  read_delim("data/tissue_atlas_hierarchy.txt", delim = "\t")

all.atlas.max <- 
  read_delim("data/consensus_max_of_all_groups_92.tsv", delim = "\t") %>%
  rename(consensus_content_name = 2, 
         max_tpm = tpm,
         max_ptpm = ptpm,
         max_nx = nx)

endothelial_atlas <- 
  read_delim("data/endothelium lims.txt", delim = "\t") %>%
  group_by(endo_ID) %>%
  mutate(ptpm = tpm / (sum(tpm) / 1e6)) %>%
  group_by(ensg_id) %>%
  summarise(ptpm = mean(ptpm), 
            log_ptpm = log10(ptpm + 1))

tissue_summation <- 
  c("liver" = "liver",
    "brain" = "brain",
    "blood" = "blood",
    "prostate" = "male tissues",
    "testis" = "male tissues",
    "ductus deferens" = "male tissues",
    "epididymis" = "male tissues",
    "seminal vesicle" = "male tissues",
    "vagina" = "female tissues", 
    "ovary" = "female tissues", 
    "fallopian tube" = "female tissues",
    "endometrium" = "female tissues", 
    "cervix, uterine" = "female tissues",
    "placenta" = "female tissues", 
    "breast" = "female tissues",
    "salivary gland" = "digestive system", 
    "esophagus" = "digestive system", 
    "tongue" = "digestive system",
    "pancreas" = "digestive system",
    "stomach" = "digestive system", 
    "intestine" = "digestive system"
    #"duodenum" = "digestive system",
    #"small intestine" = "digestive system",
    #"colon" = "digestive system", 
    #"rectum" = "digestive system"
  ) %>%
  {tibble(consensus_content_name = names(.),
          sum_name = .)}

loc_levels <- c('Blood', 
                'Female locally',
                'Male locally', 
                'Gastric', 
                'Brain locally', 
                'Matrix', 
                'Other locally', 
                'Secreted - no data', 
                'Intracellular/Membrane')


protein_names <- 
  read_delim("data/new_proteinclass_all_19670.txt", delim = "\t")




atlas_categories <- 
  read_delim("data/consensus_all_category_92.tsv", delim = "\t")

blood_functions <- 
  read_delim("data/730_blood_func.txt", delim = "\t")

blood_isoforms <- 
  read_delim("data/isoforms_for_730_blood.txt", delim = "\t")


# ----- data wrangling -----

# Calculate maximum expression in tissue groups
secretome_atlas <- 
  all.atlas.max %>%
  filter(consensus_content_name %in% tissue_summation$consensus_content_name) %>%
  filter(ensg_id %in% secretome_genes$ensg_id) %>%
  left_join(tissue_summation, by = "consensus_content_name") %>%
  group_by(ensg_id, sum_name) %>%
  summarise(max_ptpm = max(max_ptpm), 
            max_nx = max(max_nx)) %>%
  ungroup() %>%
  left_join(secretome_genes, by = "ensg_id") %>%
  mutate(log_ptpm = log10(max_ptpm + 1),
         log_nx = log10(max_nx + 1)) 

# Get labels for tissues with the number of genes included
loc_levels_n <- 
  secretome_atlas %>%
  group_by(atlasloc) %>%
  summarise(n_genes = paste0("n=", length(unique(ensg_id))),
            atlasloc_n = paste0(unique(atlasloc), " (", n_genes, ")")) %>%
  arrange(match(atlasloc, loc_levels)) %$%
  atlasloc_n

secretome_atlas <- 
  secretome_atlas %>%
  group_by(atlasloc) %>%
  mutate(n_genes = paste0("n=", length(unique(ensg_id))), 
         atlasloc_n = paste0(unique(atlasloc), " (", n_genes, ")")) %>%
  ungroup() %>%
  mutate(atlasloc_n = factor(atlasloc_n, levels = loc_levels_n), 
         atlasloc = factor(atlasloc, levels = loc_levels))




# ----- Line plot and facet heat map - expression trends -----
secretome_atlas %>% 
  mutate(sum_name = factor(sum_name, levels = c("blood", "brain", "digestive system", "female tissues", "male tissues", "liver"))) %>% {
    ggplot(., aes(sum_name, log_nx)) + 
      geom_line(aes(group = ensg_id),
                alpha = 0.1, 
                size = 0.7) +
      geom_line(data = group_by(., sum_name, atlasloc_n, n_genes) %>%
                  summarise(median_log_ptpm = median(log_nx)), 
                aes(sum_name, median_log_ptpm, group = atlasloc_n), 
                color = "firebrick2", size = 2) + 
      facet_wrap(~atlasloc_n)+
      scale_x_discrete(expand = c(0.05, 0.05))+
      scale_fill_manual(values = with(location_colors, setNames(color, atlasloc_n)))+
      theme(panel.background = element_rect(fill = NA, colour = NA), 
            #panel.border = element_rect(color = "black", fill = NA),
            plot.background = element_rect(fill = NA, color = NA),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "lightgray"), 
            strip.text = element_text(size = 10),
            legend.key = element_rect(colour = NA),
            #legend.position = "bottom",
            #legend.direction = "horizontal",
            legend.key.size= unit(0.3, "cm"),
            legend.title = element_text(face="italic"),
            axis.line = element_line(colour="black",size=0.5), 
            axis.text.x = element_text(angle = 60, hjust = 1)) + 
      xlab("Tissues") +
      ylab("Log10(NX + 1)")}
ggsave(save_path("line annotation trends.svg"), height = 8, width = 6)
ggsave(save_path("line annotation trends.pdf"), height = 8, width = 6)


secretome_atlas %>%
  Aasa_facet_heatmap(xcol = "ensg_id", 
                     ycol = "sum_name", 
                     value_col = "log_nx", 
                     facet_col = "atlasloc",
                     xband_color_cols = c("color"), 
                     x_clustering_method = "ward.D2",
                     y_clustering_method = "ward.D2",
                     x_margin = 0.1,
                     y_margin = 0.3,
                     x_dendrogram_nudge = 0,
                     y_dendrogram_nudge = 0,
                     value_title = "log10(NX + 1)") +
  xlab("Genes") +
  ylab("Tissues")
ggsave(save_path("heatmap facet annotations.svg"), height = 8, width = 10)
ggsave(save_path("heatmap facet annotations.pdf"), height = 8, width = 10)



# ----- Blood annotated genes plots -----

secretome_atlas_blood_all <- 
  all.atlas.max %>%
  filter(ensg_id %in% filter(secretome_atlas, atlasloc == "Blood")$ensg_id) %>%
  mutate(log_nx = log10(max_nx + 1)) 
 


# Create initial clusters
pheat_data <- 
  secretome_atlas_blood_all %>%
  select(ensg_id, consensus_content_name, log_nx) %>%
  spread(key = ensg_id, value = log_nx) 


n_clusters <- 9
cluster_colors = setNames(viridis::inferno(n_clusters), 1:n_clusters)
                   
cluster_colors_l <- setNames(viridis::viridis(n_clusters), LETTERS[1:n_clusters])

all_tissues_cluster <- 
  pheat_data %>%
  column_to_rownames("consensus_content_name") %>% 
  as.matrix() %>%
  pheatmap::pheatmap(cluster_cols = T,
                     cluster_rows = T,
                     fontsize_row = 10, 
                     show_colnames = F,
                     clustering_method = 'ward.D2',
                     cutree_col = n_clusters,
                     treeheight_col = 30,
                     color = viridis::inferno(n = 20, direction = -1))

pheat_data %>%
  column_to_rownames("consensus_content_name") %>% 
  as.matrix() %>%
  pheatmap::pheatmap(cluster_cols = T,
                     cluster_rows = T,
                     fontsize_row = 10, 
                     show_colnames = F,
                     clustering_method = 'ward.D2', 
                     annotation_col = cutree(all_tissues_cluster$tree_col, n_clusters) %>% 
                       as.data.frame() %>%
                       rownames_to_column("ensg_id") %>%
                       rename(cluster = ".") %>%
                       mutate(cluster = as.factor(cluster)) %>%
                       column_to_rownames("ensg_id"), 
                     annotation_colors = list(cluster = cluster_colors),
                     cutree_col = n_clusters,
                     treeheight_col = 30,
                     color = viridis::inferno(n = 20, direction = -1))



secretome_atlas_blood_all2 <- 
  secretome_atlas_blood_all %>%
  left_join(cutree(all_tissues_cluster$tree_col, n_clusters) %>% {
    tibble(ensg_id = names(.),
           all_cluster = .)
  },
  by = "ensg_id") %>%
  mutate(all_cluster_color = cluster_colors[match(all_cluster, names(cluster_colors))],
         all_cluster_l = case_when(all_cluster == 1 ~ "I",
                                   all_cluster == 2 ~ "F",
                                   all_cluster == 3 ~ "A",
                                   all_cluster == 4 ~ "E",
                                   all_cluster == 5 ~ "C",
                                   all_cluster == 6 ~ "G",
                                   all_cluster == 7 ~ "H",
                                   all_cluster == 8 ~ "B",
                                   all_cluster == 9 ~ "D")) 




# Create heatmap with cluster
exp_heatmap <- function() {
  pheat_data <- 
    secretome_atlas_blood_all2 %>%
    select(ensg_id, consensus_content_name, log_nx) %>%
    spread(key = ensg_id, value = log_nx) 
  
  pheat_meta <-
    secretome_atlas_blood_all2 %>%
    
    mutate(Cluster = as.factor(all_cluster_l)) %>%
    select(ensg_id, Cluster) %>%
    unique() %>%
    inner_join(atlas_categories %>%
                 select(ensg_id,
                        `Specificity category` = specificity_category), 
               by = "ensg_id") %>%
    inner_join(blood_functions %>%
                 select(ensg_id,
                        Function = plotfunc), 
               by = "ensg_id") %>%
    inner_join(blood_isoforms %>%
                 select(ensg_id,
                        Isoforms = type), 
               by = "ensg_id") %>%
    inner_join(endothelial_atlas %>%
                 select(`Endothelial log10(pTPM + 1)` = log_ptpm, 
                        ensg_id), 
               by = "ensg_id") %>%
    select(ensg_id, `Endothelial log10(pTPM + 1)`, Function, Isoforms, `Specificity category`, Cluster) %>%
    column_to_rownames("ensg_id") %>%
    as.data.frame()
  
  tiss_h <- 
    tissue_hierarchy %>%
    mutate(l = ifelse(l == "lymphoid system", 
                      "lymphoid tissue", 
                      tolower(l)))
  
  pheat_meta2 <-
    tiss_h %>%
    select(l, `Tissue group` = l1) %>%
    unique() %>%
    column_to_rownames("l") %>%
    as.data.frame()
  
  pheat_meta2_colors <- 
    tiss_h %>% 
    select(l1, colors) %>%
    unique() %$%
    setNames(colors, l1) 
  
  
                                                
  
  pheat_data %>%
    column_to_rownames("consensus_content_name") %>% 
    as.matrix() %>%
    pheatmap::pheatmap(cluster_cols = T,
                       cluster_rows = T,
                       fontsize_row = 10, 
                       show_colnames = F,
                       clustering_method = 'ward.D2',
                       annotation_legend = T,
                       annotation_colors = list(Cluster = cluster_colors_l, 
                                                `Tissue group` = pheat_meta2_colors, 
                                                `Specificity category` = elevated.cat.cols, 
                                                Isoforms = c("non-secreted isoforms" = "red", 
                                                             "secreted" = "white"), 
                                                Function = function_palette,
                                                `Endothelial log10(pTPM + 1)` = viridis::inferno(20, direction = -1)),
                       cutree_col = 9,
                       annotation_col = pheat_meta,
                       annotation_row = pheat_meta2,
                       treeheight_col = 50,
                       color = viridis::inferno(n = 20, direction = -1))
}

svg(save_path("Heatmap blood ann expression with allclusters.svg"), height = 8, width = 15)
exp_heatmap()
dev.off()
pdf(save_path("Heatmap blood ann expression with allclusters.pdf"), height = 8, width = 15)
exp_heatmap()
dev.off()

# Lineplot showing trends within clusters
secretome_atlas_blood_all2 %>%
  ggplot(aes(consensus_content_name, log_nx)) + 
  geom_line(aes(group = ensg_id),
            alpha = 0.1, 
            size = 0.7) +
  geom_line(data = secretome_atlas_blood_all2 %>% 
              group_by(consensus_content_name, all_cluster_l) %>%
              summarise(median_log_nx = median(log_nx)), 
            aes(consensus_content_name, median_log_nx, group = all_cluster_l), 
            color = "firebrick2", size = 2) + 
  ggrepel::geom_text_repel(data = secretome_atlas_blood_all2 %>% 
              group_by(consensus_content_name, all_cluster_l) %>%
              summarise(median_log_nx = median(log_nx)) %>%
              ungroup() %>%
              group_by(all_cluster_l) %>%
              mutate(rank = rank(median_log_nx)) %>%
              filter(rank > max(rank) - 2 | rank < min(rank) + 2), 
            aes(consensus_content_name, median_log_nx, group = all_cluster_l, label= consensus_content_name),
            nudge_x = 0, 
            nudge_y = 0.5) + 
  facet_wrap(.~all_cluster_l, ncol = 3)+
  scale_fill_manual(values = cluster_colors)+
  theme(panel.background = element_rect(fill = NA, colour = NA),
        plot.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.key = element_rect(colour = NA),
        #legend.position = "bottom",
        #legend.direction = "horizontal",
        legend.key.size= unit(0.3, "cm"),
        legend.title = element_text(face="italic"),
        axis.line = element_line(colour="black",size=0.5), 
        axis.text.x = element_text(angle = 60, hjust = 1))+
  xlab("Tissues")+
  ylab("Log10(NX + 1)")
ggsave(save_path("line trends in blood all clusters.pdf"), height = 10, width = 12)
ggsave(save_path("line trends in blood all clusters.svg"), height = 10, width = 12)





# ----- expression chord diagram -----
secretome_atlas_blood_all2 %>%
  select(ensg_id, all_cluster_l, consensus_content_name,  max_nx) %>%
  group_by(ensg_id, all_cluster_l) %>%
  summarise(median_NX = median(max_nx), 
            log10_median_NX = log10(median_NX + 1)) %>%
  write_csv(save_path("genes blood cluster.csv"))

exp_chord <- function() {
  chord_data <- 
    secretome_atlas_blood_all2 %>% 
    left_join(tissue_hierarchy %>%
                mutate(l = ifelse(l == "lymphoid system", 
                                  "lymphoid tissue", 
                                  tolower(l))), 
              by = c("consensus_content_name" = "l")) %>% 
    group_by(all_cluster_l) %>%
    mutate(n_cluster_genes = length(unique(ensg_id))) %>%
    group_by(l1, ensg_id) %>%
    mutate(max_nx = max(max_nx, na.rm = T)) %>%
    group_by(l1, all_cluster_l, n_cluster_genes) %>%
    summarise(mean_nx = mean(max_nx)) %>%
    ungroup() 
  
  to_subtitle <- 
    chord_data %>% 
    select(to = all_cluster_l, 
           subtitle = n_cluster_genes) %>% 
    unique() %>% 
    mutate(subtitle = paste0("n=", subtitle))
  
  from = chord_data$l1
  to = chord_data$all_cluster_l
  sizes = chord_data$mean_nx
  
  
  grid.col <- 
    c(tissue_hierarchy %>% 
        select(l1, colors) %>% 
        unique() %$% 
        setNames(colors, l1),
      cluster_colors_l)
  
  
  
  require(circlize) 
  
  factors.from <- unique(from)
  factors.to <- unique(to)
  factors <- c(factors.from, factors.to)
  
  groups = c(rep(1, length(factors.from)),
             rep(2, length(factors.to)))
  plot.order <- c(factors.from, factors.to)
  
  tb <- 
    tibble(from, to, sizes)
  
  #groups <- groups[plot.order]
  gap.after.par <- c()
  for(i in 1:(length(groups)-1)) {
    if(groups[i] == groups[i+1]) {
      gap.after.par <- c(gap.after.par, 2)
    } else {
      gap.after.par <- c(gap.after.par, 15)
    }
  }
  
  if(groups[length(groups)] == groups[1]) {
    gap.after.par <- c(gap.after.par, 2)
  } else {
    gap.after.par <- c(gap.after.par, 15)
  }
  
  circos.par(gap.after = gap.after.par)
  
  chord <-
    tb %>% 
    chordDiagram(grid.col = grid.col,
                 directional = 0,
                 annotationTrack="grid",
                 annotationTrackHeight = 0.05, 
                 preAllocateTracks = 1, 
                 order = plot.order, transparency = 0.2)
  
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    sector.name <- get.cell.meta.data("sector.index")
    sector.index <- get.cell.meta.data("sector.numeric.index")
    
    
    
    if(sector.name %in% factors.to) {
      circos.text(mean(xlim), 0.4, 
                  sector.name,
                  cex = 1.5, 
                  niceFacing = TRUE, 
                  facing = "bending")
      circos.text(mean(xlim), 0.1, 
                  to_subtitle$subtitle[match(sector.name, to_subtitle$to)],
                  cex = 1, 
                  niceFacing = TRUE, 
                  facing = "bending")
    } else {
      circos.text(mean(xlim), 0, 
                  sector.name, 
                  adj = degree(0), 
                  niceFacing = TRUE, 
                  facing = "clockwise")
    }
  }, bg.border = NA)
  
  circos.clear()
}



svg(save_path("chord mean NX tissues in clusters.svg"), width = 10, height = 10)
exp_chord()
dev.off()
pdf(save_path("chord mean NX tissues in clusters.pdf"), width = 10, height = 10)
exp_chord()
dev.off()


# ----- function chord diagram -----


func_chord <- function() {
  chord_data <- 
    secretome_atlas_blood_all2 %>% 
    inner_join(blood_functions %>%
                 select(ensg_id,
                        Function = plotfunc), 
               by = "ensg_id") %>%
    select(ensg_id, 
           to = all_cluster_l, 
           from = Function) %>%
    unique() %>%
    group_by(from, to) %>%
    summarise(sizes = length(ensg_id))
  
  to_subtitle <- 
    secretome_atlas_blood_all2 %>% 
    select(to = all_cluster_l, 
           subtitle = ensg_id) %>% 
    unique() %>% 
    group_by(to) %>%
    summarise(subtitle = paste0("n=", length(subtitle)))
  
  from = chord_data$from
  to = chord_data$to
  sizes = chord_data$sizes
  
  
  grid.col <- 
    c(function_palette,
      cluster_colors_l)
  
  
  
  require(circlize) 
  
  factors.from <- unique(from)
  factors.to <- unique(to)
  factors <- c(factors.from, factors.to)
  
  groups = c(rep(1, length(factors.from)),
             rep(2, length(factors.to)))
  plot.order <- c('No annotated function',
                  'Other',
                  'Immunity,Defense',
                  'Receptors,Transport', 
                  'Enzymes,Enzyme inhibitors',
                  'Hormones,Growth factors',
                  'Cytokines',
                  'Coagulation,Complement,Acute phase,APO',
                  LETTERS[1:9])
  
  tb <- 
    tibble(from, to, sizes)
  
  #groups <- groups[plot.order]
  gap.after.par <- c()
  for(i in 1:(length(groups)-1)) {
    if(groups[i] == groups[i+1]) {
      gap.after.par <- c(gap.after.par, 2)
    } else {
      gap.after.par <- c(gap.after.par, 15)
    }
  }
  
  if(groups[length(groups)] == groups[1]) {
    gap.after.par <- c(gap.after.par, 2)
  } else {
    gap.after.par <- c(gap.after.par, 15)
  }
  
  circos.par(gap.after = gap.after.par)
  
  chord <-
    tb %>% 
    chordDiagram(grid.col = grid.col,
                 directional = 0,
                 annotationTrack="grid",
                 annotationTrackHeight = 0.05, 
                 preAllocateTracks = 1, 
                 order = plot.order, transparency = 0.2)
  
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    sector.name <- get.cell.meta.data("sector.index")
    sector.index <- get.cell.meta.data("sector.numeric.index")
    
    
    
    if(sector.name %in% factors.to) {
      circos.text(mean(xlim), 0.4, 
                  sector.name,
                  cex = 1.5, 
                  niceFacing = TRUE, 
                  facing = "bending")
      circos.text(mean(xlim), 0.1, 
                  to_subtitle$subtitle[match(sector.name, to_subtitle$to)],
                  cex = 1, 
                  niceFacing = TRUE, 
                  facing = "bending")
    } else {
      circos.text(mean(xlim), 0, 
                  sector.name, 
                  adj = degree(0), 
                  niceFacing = TRUE, 
                  facing = "clockwise")
    }
  }, bg.border = NA)
  
  circos.clear()
}



svg(save_path("chord gene function clusters.svg"), width = 10, height = 10)
func_chord()
dev.off()
pdf(save_path("chord gene function clusters.pdf"), width = 10, height = 10)
func_chord()
dev.off()


# ----- endothelial expression analysis + plots incl. heatmap ------

endothelial_atlas_blood_cluster <- 
  left_join(endothelial_atlas, 
            secretome_atlas_blood_all2 %>%
              select(ensg_id, all_cluster_l) %>%
              unique(),
            by = "ensg_id") %>%
  filter(!is.na(all_cluster_l))

endo_exp_chord <- function() {
  chord_data <- 
    endothelial_atlas_blood_cluster %>% 
    mutate(content = "endothelial cells") %>%
    group_by(content, all_cluster_l) %>%
    summarise(mean_ptpm = mean(ptpm)) %>%
    filter(mean_ptpm > 0)
  
  
  
  from = chord_data$content
  to = chord_data$all_cluster_l
  sizes = chord_data$mean_ptpm
  
  
  grid.col <- 
    c("endothelial cells" = "chocolate",
      cluster_colors_l)
  
  
  
  require(circlize) 
  
  factors.from <- unique(from)
  factors.to <- unique(to)
  factors <- c(factors.from, factors.to)
  
  groups = c(rep(1, length(factors.from)),
             rep(2, length(factors.to)))
  plot.order <- c(factors.from, factors.to)
  
  tb <- 
    tibble(from, to, sizes)
  
  #groups <- groups[plot.order]
  gap.after.par <- c()
  for(i in 1:(length(groups)-1)) {
    if(groups[i] == groups[i+1]) {
      gap.after.par <- c(gap.after.par, 2)
    } else {
      gap.after.par <- c(gap.after.par, 15)
    }
  }
  
  if(groups[length(groups)] == groups[1]) {
    gap.after.par <- c(gap.after.par, 2)
  } else {
    gap.after.par <- c(gap.after.par, 15)
  }
  
  circos.par(gap.after = gap.after.par, start.degree = -20)
  
  chord <-
    tb %>% 
    chordDiagram(grid.col = grid.col,
                 directional = 0,
                 annotationTrack="grid",
                 annotationTrackHeight = 0.05, 
                 preAllocateTracks = 1, 
                 order = plot.order, transparency = 0.2)
  
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    sector.name <- get.cell.meta.data("sector.index")
    sector.index <- get.cell.meta.data("sector.numeric.index")
    
    
    # circos.segments(x0 = mean(xlim), x1 = mean(xlim), 
    #                 y0 = min(ylim), y1 = mean(ylim)-0.2 + adjustment, 
    #                 sector.name)
    
    if(sector.name %in% factors.to) {
      circos.text(mean(xlim), 0.2, 
                  sector.name,
                  cex = 1.5, 
                  niceFacing = TRUE, 
                  facing = "bending")
    } else {
      circos.text(mean(xlim), 0, 
                  sector.name, 
                  adj = degree(0), 
                  niceFacing = TRUE, 
                  facing = "clockwise")
    }
  }, bg.border = NA)
  
  circos.clear()
} 
svg(save_path("chord mean NX endothelial in clusters.svg"), width = 10, height = 10)
endo_exp_chord()
dev.off()

pdf(save_path("chord mean NX endothelial in clusters.pdf"), width = 10, height = 10)
endo_exp_chord()
dev.off()



# ---- endothelial ptpm cluster analysis + plots ----

# Cluster with endothelial and tissues on ptpm

pheat_data <- 
  secretome_atlas_blood_all2 %>%
  mutate(log_ptpm = log10(max_ptpm + 1)) %>%
  select(ensg_id, consensus_content_name, log_ptpm) %>%
  spread(key = consensus_content_name, value = log_ptpm) %>%
  inner_join(endothelial_atlas %>% 
               select(ensg_id, 
                      `endothelial cells` = log_ptpm), 
             by = "ensg_id") %>%
  gather(key = "consensus_content_name", value = "value", -1) %>%
  spread(key = ensg_id, value = value)

pheat_meta <-
  secretome_atlas_blood_all2 %>%
  
  mutate(Cluster = as.factor(all_cluster_l)) %>%
  select(ensg_id, Cluster) %>%
  unique() %>%
  column_to_rownames("ensg_id") %>%
  as.data.frame() 

pdf(save_path("Heatmap blood ann PTPM with endo.pdf"), height = 8, width = 10)
pheat_data %>%
  column_to_rownames("consensus_content_name") %>% 
  as.matrix() %>%
  pheatmap::pheatmap(cluster_cols = T,
                     cluster_rows = T,
                     fontsize_row = 10, 
                     show_colnames = F,
                     clustering_method = 'ward.D2',
                     annotation_legend = F,
                     annotation_colors = list(Cluster = cluster_colors_l),
                     cutree_col = 9,
                     annotation_col = pheat_meta,
                     treeheight_col = 50,
                     color = viridis::inferno(n = 20, direction = -1))
dev.off()


# ----- experimental enrichment analysis word clouds (not included in paper) -----

library(enrichR)

enrichR_data <- 
  secretome_atlas_blood_all2 %>%
  select(ensg_id, all_cluster_l) %>%
  unique() %>%
  inner_join(protein_names %>% 
               select(1:2), 
             by = "ensg_id") 
  

enr_levels <- sort(unique(enrichR_data$all_cluster_l))  
  
enr_results <- 
  lapply(enr_levels, 
         FUN = function(enr) {
           filter(enrichR_data,
                  all_cluster_l == enr) %$%
             enrichr(display_name, databases = c("GO_Biological_Process_2018",
                                                 "GO_Molecular_Function_2018")) %>%
             plyr::ldply(.id = "database")
         }) %>%
  set_names(enr_levels) %>%
  plyr::ldply(.id = "all_cluster_l")

enr_results_filtered <- 
  enr_results %>%
  as.tibble() %>%
  filter(Adjusted.P.value < 0.01) %>%
  mutate(log_score = round(-log10(Adjusted.P.value), 0),
         weighted_term = mapply(Term, log_score, FUN = function(t, r) rep(t, r) %>% paste(collapse = " ")))
  



require("tm")
require("SnowballC")
require("wordcloud")

enr_text <- 
  enr_results_filtered %>% 
  mutate(Term = gsub(" \\(.*\\)", "", Term)) %>% 
  group_by(all_cluster_l) %>%
  summarise(text = paste(weighted_term, collapse = " ")) %$% {
    lapply(text, 
           FUN = function(txt) {
             Corpus(VectorSource(txt)) %>%
               tm_map(content_transformer(tolower)) %>%
               tm_map(removeNumbers) %>%
               tm_map(removeWords, stopwords("english")) %>%
               tm_map(removeWords, c("regulation", "positive", "negative", "cell", "process", "activity", "response", "protein")) %>%
               tm_map(removePunctuation) %>%
               tm_map(stripWhitespace) %>%
               TermDocumentMatrix() %>%
               as.matrix() %>%
               {sort(rowSums(.),decreasing=TRUE)} %>%
               {tibble(word = names(.),
                       freq = .)}
           }) %>%
      set_names(all_cluster_l) %>%
      plyr::ldply(.id = "all_cluster_l")
  }

for(lvl in enr_levels) {
  pdf(save_path(paste0("worldcloud cluster ", lvl, ".pdf")), height = 5, width = 5)
  set.seed(42)
  enr_text %>%
    filter(all_cluster_l == lvl) %$%
    wordcloud(words = word, 
              freq = freq, 
              min.freq = 1,
              max.words=15, 
              random.order=FALSE, 
              rot.per=0, 
              colors=viridis::inferno(7)[-(6:7)])
  dev.off()
}



sapply(enr_levels, 
       FUN = function(a) {
         sapply(enr_levels, function(b) {
           inner_join(filter(enr_results_filtered, all_cluster_l == a),
                      filter(enr_results_filtered, all_cluster_l == b),
                      by = "Term") %>%
             nrow()
         })
       }) %>%
  as.tibble(rownames = "from") %>%
  gather(key = "to", value = "count", -from) %>%
  left_join(filter(., from == to) %>%
              select(from, from_total_count = count), 
            by = "from") %>%
  left_join(filter(., from == to) %>%
              select(to, to_total_count = count), 
            by = "to") %>%
  filter(match(from, LETTERS) < match(to, LETTERS)) %$%
  chord_classification(from = from, to = to, sizes = count, 
                       grid.col = cluster_colors_l, 
                       groups = rep(1, n_clusters), 
                       plot.order = enr_levels, 
                       line_expansion = 0)


