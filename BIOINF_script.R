install.packages("tidyverse")
install.packages("igraph")
library('tidyverse')
library('igraph')
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)
install.packages('expm')
library(expm)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(igraph)
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
# STRINGNetworkToolbox# STRING DATABASE ---------------------------------------------------------
string_raw <- read_delim("9606.protein.links.v12.0 STRING.txt", delim = " ")

# --- 2. PULIZIA E PREPARAZIONE ID ---
string_ensp <- string_raw %>%
  # Rimuoviamo il prefisso "9606." per isolare l'Ensembl Protein ID (ENSP)
  mutate(from_ensp = str_remove(protein1, "9606."),
         to_ensp = str_remove(protein2, "9606.")) %>%
  # Rimuoviamo i self-loops (Punto 1.1 del progetto)
  filter(from_ensp != to_ensp) %>%
  dplyr::select(from_ensp, to_ensp, combined_score)

# --- 3. MAPPATURA ENSP -> GENE SYMBOL (HGNC) ---
# Otteniamo la lista unica di tutti i codici ENSP nel dataset
all_ensp <- unique(c(string_ensp$from_ensp, string_ensp$to_ensp))

# Connessione ad Ensembl (biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Scarichiamo la tabella di conversione
string_mapping <- getBM(attributes = c('ensembl_peptide_id', 'hgnc_symbol'), 
                        filters = 'ensembl_peptide_id', 
                        values = all_ensp, 
                        mart = mart)

# Applichiamo il mapping al dataset originale
string_mapped <- string_ensp %>%
  left_join(string_mapping, by = c("from_ensp" = "ensembl_peptide_id")) %>%
  rename(from_symbol = hgnc_symbol) %>%
  left_join(string_mapping, by = c("to_ensp" = "ensembl_peptide_id")) %>%
  rename(to_symbol = hgnc_symbol)

# Pulizia finale: teniamo solo le coppie con simboli validi e non vuoti
string_final <- string_mapped %>%
  filter(!is.na(from_symbol) & from_symbol != "",
         !is.na(to_symbol) & to_symbol != "") %>%
  # Rimuoviamo ridondanze AB-BA (Punto 1.1)
  mutate(node_min = pmin(from_symbol, to_symbol), 
         node_max = pmax(from_symbol, to_symbol)) %>%
  distinct(node_min, node_max, .keep_all = TRUE) %>%
  dplyr::select(from = from_symbol, to = to_symbol)

# --- 4. CREAZIONE GRAFO E ISOLAMENTO LCC ---
g_string <- graph_from_data_frame(string_final, directed = FALSE)

# Identifichiamo i componenti connessi
comp_str <- components(g_string)

# Selezioniamo i nodi appartenenti alla componente gigante (LCC)
nodes_str_lcc <- V(g_string)[comp_str$membership == which.max(comp_str$csize)]

# Creiamo il sottografo LCC (Interactome LCC)
g_string_lcc <- induced_subgraph(g_string, nodes_str_lcc)

# --- 5. VERIFICA DATI PER TABELLA 1 ---
print(paste("Nodi LCC STRING:", vcount(g_string_lcc)))
print(paste("Archi LCC STRING:", ecount(g_string_lcc)))


#REACTOME DATABASE--------------------------------------------------------
reactome_raw <- read_tsv("FIsInGene_04142025_with_annotations_REACTOME.txt")

# Pulizia specifica Reactome (Step 1.1 del PDF)
reactome_clean <- reactome_raw %>%
  select(from = Gene1, to = Gene2) %>% # Usa le colonne Gene1 e Gene2 
  filter(from != to) %>% # Rimuovi self-loops [cite: 947]
  distinct(pmin(from, to), pmax(from, to), .keep_all = TRUE) # Rimuovi ridondanze [cite: 947]
# 1. Rimuovi le colonne extra e i self-loops
reactome_final <- reactome_clean %>%
  filter(from != to) %>%
  select(from, to) %>%
  distinct()

# 2. Crea il grafo e isola la LCC (Step 1.1 progetto)
g_reactome <- graph_from_data_frame(reactome_final, directed = FALSE)
comp <- components(g_reactome)
# Seleziona i nodi della componente più grande
nodes_lcc <- V(g_reactome)[comp$membership == which.max(comp$csize)]
g_reactome_lcc <- induced_subgraph(g_reactome, nodes_lcc)

# 3. SEGNA QUESTI DATI PER LA TABELLA 1
vcount(g_reactome_lcc) # Numero nodi LCC
ecount(g_reactome_lcc) # Numero archi LCC

#BIOGRID
biogrid_raw <- read_tsv("BIOGRID-ORGANISM-Homo_sapiens-5.0.252.tab3.txt")

# 2. Filtro e Pulizia (Step 1.1 del progetto)
biogrid_clean <- biogrid_raw %>%
  # Entrambi gli organismi devono essere 9606 (Homo sapiens) [cite: 932]
  filter(`Organism ID Interactor A` == 9606, `Organism ID Interactor B` == 9606) %>%
  # Solo interazioni fisiche [cite: 933]
  filter(`Experimental System Type` == "physical") %>%
  # Seleziona i Gene Symbols ufficiali [cite: 952]
  select(from = `Official Symbol Interactor A`, to = `Official Symbol Interactor B`) %>%
  filter(from != to) %>% # Rimuovi self-loops [cite: 947]
  distinct() # Rimuovi ridondanze

# 3. Estrazione LCC [cite: 947]
g_biogrid <- graph_from_data_frame(biogrid_clean, directed = FALSE)
comp_bio <- components(g_biogrid)
nodes_lcc_bio <- V(g_biogrid)[comp_bio$membership == which.max(comp_bio$csize)]
g_biogrid_lcc <- induced_subgraph(g_biogrid, nodes_lcc_bio)

# Segna i dati per la Tabella 1 
vcount(g_biogrid_lcc) # Numero nodi BioGRID LCC
ecount(g_biogrid_lcc) # Numero archi BioGRID LCC

#HURI
# 1. Caricamento HuRI (formato a due colonne) [cite: 936, 937]
huri_raw <- read_tsv("HuRI.tsv", col_names = c("from", "to"))

# 2. Pulizia (Step 1.1) [cite: 947]
huri_clean <- huri_raw %>%
  filter(from != to) %>%
  distinct()

# 3. Estrazione LCC [cite: 947]
g_huri <- graph_from_data_frame(huri_clean, directed = FALSE)
comp_huri <- components(g_huri)
nodes_lcc_huri <- V(g_huri)[comp_huri$membership == which.max(comp_huri$csize)]
g_huri_lcc <- induced_subgraph(g_huri, nodes_lcc_huri)

vcount(g_huri_lcc) # Numero nodi BioGRID LCC
ecount(g_huri_lcc)

# x questo dataset devo anche cercare gli esamble
# 1. Scarica la tabella di conversione (Mapping)
mapping_table <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                       filters = 'ensembl_gene_id', 
                       values = all_ensembl_ids, 
                       mart = mart)

# 2. Unisci i simboli e correggi l'errore select
huri_symbols <- huri_clean %>%
  left_join(mapping_table, by = c("from" = "ensembl_gene_id")) %>%
  rename(from_symbol = hgnc_symbol) %>%
  left_join(mapping_table, by = c("to" = "ensembl_gene_id")) %>%
  rename(to_symbol = hgnc_symbol)

# 3. Pulizia finale usando dplyr::select per evitare conflitti
huri_final <- huri_symbols %>%
  filter(!is.na(from_symbol) & from_symbol != "",
         !is.na(to_symbol) & to_symbol != "") %>%
  dplyr::select(from = from_symbol, to = to_symbol) %>% # <--- CORREZIONE QUI
  distinct()

# Ora ricrea la LCC con i nomi corretti
g_huri <- graph_from_data_frame(huri_final, directed = FALSE)
comp_huri <- components(g_huri)
nodes_lcc_huri <- V(g_huri)[comp_huri$membership == which.max(comp_huri$csize)]
g_huri_lcc <- induced_subgraph(g_huri, nodes_lcc_huri)




# TAKE DISEASE GENES ------------------------------------------------------
ms_genes_raw <- read_delim("diseases_", delim = "\t")


# Estrai i simboli unici dei geni
# Assicurati che la colonna si chiami "Associated genes" o adatta il nome
seed_genes <- ms_genes_raw %>%
  pull(`Associated genes`) %>%
  unique()

# Conta quanti sono inizialmente (per il report)
n_initial_genes <- length(seed_genes)
print(paste("Numero iniziale di geni GDA per Sclerosi Multipla:", n_initial_genes))


# DATA REACTOME FOR TABLE 1 ---------------------------------------------

# 1. Identifica i geni MS(multiple sclerosis) presenti nella LCC di Reactome
ms_in_reactome <- intersect(seed_genes, V(g_reactome_lcc)$name)
n_ms_in_reactome <- length(ms_in_reactome)
perc_ms_in_reactome <- (n_ms_in_reactome / n_initial_genes) * 100

# 2. Estrai il Disease Interactome (interazioni SOLO tra geni MS)
# (Punto 1.2 del PDF: identify the disease interactomes by getting the interactions among disease genes)
g_ms_reactome <- induced_subgraph(g_reactome_lcc, ms_in_reactome)

# 3. Trova la LCC del Disease Interactome
ms_comp <- components(g_ms_reactome)
g_ms_reactome_lcc <- induced_subgraph(g_ms_reactome, 
                                      V(g_ms_reactome)[ms_comp$membership == which.max(ms_comp$csize)])

# Interactome LCC size: vcount(g_reactome_lcc), ecount(g_reactome_lcc)
# Number of disease genes present: n_ms_in_reactome
# Percentage: perc_ms_in_reactome
# LCC size of disease interactome: vcount(g_ms_reactome_lcc)

# --- 1. FUNZIONE AUTOMATIZZATA PER LA TABELLA 1 ---
# Creiamo una funzione per non ripetere il codice 4 volte
calculate_tab1_metrics <- function(interactome_lcc, seed_genes, name) {
  
  # 1. Geni della malattia presenti nell'interattoma (Punto 1.3)
  ms_present <- intersect(seed_genes, V(interactome_lcc)$name)
  n_present <- length(ms_present)
  perc_present <- (n_present / length(seed_genes)) * 100
  
  # 2. Disease Interactome (solo interazioni tra geni MS) [cite: 955]
  g_ms_sub <- induced_subgraph(interactome_lcc, ms_present)
  
  # 3. LCC del Disease Interactome [cite: 957]
  ms_comp <- components(g_ms_sub)
  if(length(ms_present) > 0) {
    g_ms_lcc <- induced_subgraph(g_ms_sub, V(g_ms_sub)[ms_comp$membership == which.max(ms_comp$csize)])
    ms_lcc_size <- vcount(g_ms_lcc)
  } else {
    ms_lcc_size <- 0
  }
  
  # Restituisce una riga della tabella [cite: 969]
  return(data.frame(
    Interactome = name,
    Nodes_LCC = vcount(interactome_lcc),
    Links_LCC = ecount(interactome_lcc),
    MS_Genes_Present = n_present,
    Percentage = round(perc_present, 2),
    MS_LCC_Size = ms_lcc_size
  ))
}

# --- 2. ESECUZIONE PER TUTTI I DATASET ---

# Nota: Assicurati di aver convertito STRING e HuRI in Gene Symbols prima!
tab1_biogrid <- calculate_tab1_metrics(g_biogrid_lcc, seed_genes, "BioGRID")
tab1_string  <- calculate_tab1_metrics(g_string_lcc, seed_genes, "STRING")
tab1_huri    <- calculate_tab1_metrics(g_huri_lcc, seed_genes, "HuRI")
tab1_reactome <- calculate_tab1_metrics(g_reactome_lcc, seed_genes, "Reactome")

# Uniamo tutto nella Tabella 1 finale
tabella_1_finale <- rbind(tab1_biogrid, tab1_string, tab1_huri, tab1_reactome)
print(tabella_1_finale)


# --- TABLE 2 WITH STRING (Disease LCC) ---

# 1. Estraiamo specificamente la LCC della malattia per STRING
ms_present_string <- intersect(seed_genes, V(g_string_lcc)$name)
g_ms_string_sub <- induced_subgraph(g_string_lcc, ms_present_string)
ms_comp_str <- components(g_ms_string_sub)

# Questo è l'input corretto: la componente gigante della SCLEROSI MULTIPLA (235 nodi)
g_ms_string_lcc <- induced_subgraph(g_ms_string_sub, 
                                    V(g_ms_string_sub)[ms_comp_str$membership == which.max(ms_comp_str$csize)])

# 2. Ora usiamo g_ms_string_lcc come target
g_target <- g_ms_string_lcc

# Verifica dimensioni: dovresti vedere 235 nodi
print(paste("Nodi della Disease LCC:", vcount(g_target))) 

# 3. Calcolo metriche (ora sarà quasi istantaneo)
metrics <- data.frame(
  Gene_name = V(g_target)$name,
  Degree = degree(g_target),
  Betweenness = betweenness(g_target, directed = FALSE, normalized = TRUE), # [cite: 303]
  Eigenvector = eigen_centrality(g_target)$vector, # [cite: 327]
  Closeness = closeness(g_target, normalized = TRUE) # [cite: 315]
) %>%
  mutate(ratio_Betw_Degree = Betweenness / Degree) %>% # [cite: 964]
  arrange(desc(Degree)) 

# 4. Tabella 2 con i primi 20 geni [cite: 970]
tabella_2_centrality <- head(metrics, 20)
print(tabella_2_centrality)

# 5. Scatterplot richiesto [cite: 967]
library(ggplot2)
ggplot(metrics, aes(x = Degree, y = Betweenness)) +
  geom_point(color = "steelblue", alpha = 0.7) +
  geom_text(data = head(metrics, 5), aes(label = Gene_name), vjust = -1) + 
  theme_minimal() +
  labs(title = "Scatterplot: Node Degree vs Betweenness Centrality (MS LCC)",
       x = "Degree", y = "Betweenness Centrality")


# SAVE
save.image("progetto_bioinf_MS_completo.RData")
load("progetto_bioinf_MS_completo.RData")

# PART 2 OF THE PROJECT ---------------------------------------------------
# ==============================================================================
# PARTE 2: CROSS-VALIDATION E CONFRONTO ALGORITMI
# ==============================================================================

# 2. FUNZIONE DI PERFORMANCE (Punto 2.3)
calculate_performance <- function(predicted_genes, probe_set, top_k) {
  top_pred <- head(predicted_genes, top_k)
  hits <- length(intersect(top_pred, probe_set))
  
  precision <- hits / top_k
  recall <- hits / length(probe_set)
  f1 <- ifelse(precision + recall == 0, 0, 2 * (precision * recall) / (precision + recall))
  
  return(data.frame(Precision = precision, Recall = recall, F1 = f1))
}

# 3. CONFIGURAZIONE AMBIENTE
# Assicurati che il percorso sia corretto per il tuo PC
python_exe <- "C:/Users/Massimo/AppData/Local/Programs/Python/Python312/python.exe"
# INSTERSECTION LCC BIOGRID AND SEED GENES
seed_genes_biogrid <- intersect(seed_genes, V(g_reactome_lcc)$name)
# 4. PREPARAZIONE DATI (Esegui solo se non hai già questi oggetti)
set.seed(123)
shuffled_genes <- sample(seed_genes_biogrid)
folds <- cut(seq(1, length(shuffled_genes)), breaks = 5, labels = FALSE)


# SOLO DYMOND
# --- SOLO DIAMOnD ---
for (f in 1:5) {
  message("Lancio DIAMOnD - Fold ", f)
  
  training_set <- shuffled_genes[which(folds != f)]
  outfile <- paste0("diamond_reactome_fold_", f, ".txt")
  write.table(training_set, "current_seeds.txt", row.names=F, col.names=F, quote=F)
  
  # Lo esegue solo se il file non esiste già
  if (!file.exists(outfile)) {
    cmd <- paste(shQuote(python_exe), "DIAMOnD.py reactome_network.txt current_seeds.txt 200 1", outfile)
    system(cmd)
  }
}
message("Fase DIAMOnD completata!")

# --- SOLO DIFFUSIONE E ANALISI ---
tabella_performance <- data.frame()
times <- c(0.1, 1.0, 2.0, 5.0, 10.0)

run_diffusion_ultra_light <- function(graph, seeds, t, iterations = 20) {
  # 1. Matrice di Adiacenza sparsa
  W <- as_adjacency_matrix(graph, sparse = TRUE)
  
  # 2. Normalizzazione sparsa senza usare diag()
  # Dividiamo ogni riga per il grado del nodo corrispondente
  d <- degree(graph)
  d_inv <- ifelse(d > 0, 1/d, 0)
  W_norm <- Diagonal(x = d_inv) %*% W 
  
  # 3. Vettore iniziale (semi)
  p0 <- rep(0, vcount(graph))
  names(p0) <- V(graph)$name
  p0[intersect(seeds, names(p0))] <- 1
  if(sum(p0) > 0) p0 <- p0 / sum(p0) 
  
  # 4. Iterazione (Random Walk with Restart)
  # t controlla quanto il segnale si allontana (alpha alto = segnale lontano)
  alpha <- t / (1 + t) 
  pt <- as.matrix(p0) # Teniamo pt come vettore colonna
  
  for (i in 1:iterations) {
    pt <- alpha * (W_norm %*% pt) + (1 - alpha) * p0
  }
  
  # 5. Raggruppamento risultati
  pt_vec <- as.vector(pt)
  names(pt_vec) <- V(graph)$name
  
  # Rimuovi i semi dai risultati
  pt_vec[intersect(seeds, names(pt_vec))] <- -1
  
  return(sort(pt_vec, decreasing = TRUE))
}


for (f in 1:5) {
  message("\n Analisi Fold ", f)
  probe_set <- shuffled_genes[which(folds == f)]
  training_set <- shuffled_genes[which(folds != f)]
  
  # 1. DIAMOnD (Leggiamo i file che hai già creato)
  diamond_res <- read_tsv(paste0("diamond_reactome_fold_", f, ".txt"), comment="#", 
                          col_names=c("rank","gene","p"), show_col_types=F)$gene
  perf_diamond <- calculate_performance(diamond_res, probe_set, 50) %>%
    mutate(Algo="DIAMOnD", Parameter="alpha=1", Fold=f)
  tabella_performance <- rbind(tabella_performance, perf_diamond)
  
  # 2. Diffusione Iterativa (Veloce e leggera)
  for (t_val in times) {
    message("Diffusione ultra-light t = ", t_val)
    
    # Chiamata alla nuova funzione sparsa
    diff_predicted <- names(run_diffusion_ultra_light(g_reactome_lcc, training_set, t_val))
    
    # Validazione e salvataggio (Uguale a prima)
    perf_diff <- calculate_performance(diff_predicted, probe_set, 50) %>%
      mutate(Algo="Diffusion", Parameter=as.character(t_val), Fold=f)
    tabella_performance <- rbind(tabella_performance, perf_diff)
    
    gc() }}

# 3. Risultato Finale
print(tabella_performance %>% group_by(Algo, Parameter) %>% summarise(across(everything(), mean)))
# performance comparison HEAT DIFFUSION and DIAMOND on STRING
write_csv(tabella_performance, "performance_REACTOME.csv")


rm(list = setdiff(ls(), c("seed_genes", "python_exe", "calculate_performance", "run_diffusion_ultra_light")))
gc()

# TASK 2.3
# Carica i risultati (assicurati che i nomi dei file siano corretti)
res_string <- read_csv("performance_STRING.csv") %>% mutate(Network = "STRING")
res_biogrid <- read_csv("performance_BIOGRID.csv") %>% mutate(Network = "BIOGRID")
res_huri <- read_csv("performance_HURI.csv") %>% mutate(Network = "HURI")
res_reactome <- read_csv("performance_REACTOME.csv") %>% mutate(Network = "REACTOME")

# Unisci tutto e calcola le medie
tabella_comparativa <- bind_rows(res_string, res_biogrid, res_huri, res_reactome) %>%
  group_by(Network, Algo, Parameter) %>%
  summarise(Mean_Precision = mean(Precision), .groups = 'drop') %>%
  arrange(desc(Mean_Precision))

print(tabella_comparativa)


# COMPARISON PLOT
ggplot(tabella_comparativa, aes(x = Network, y = Mean_Precision, fill = Algo)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Comparison Performance",
       subtitle = "AVG Precision on first 50 genes predicted",
       y = "Mean Precision", x = "Rete") +
  scale_fill_brewer(palette = "Set1")



# PART 3 ------------------------------------------------------------------
# THE BEST COMBO NETWORK-ALGORITHM IS DIAMOND-STRING
# 1. Prepariamo tutti i semi (senza split)
# Filtriamo i geni MS presenti in STRING
final_seeds <- intersect(seed_genes, V(g_string_lcc)$name)
write.table(final_seeds, "final_ms_seeds.txt", row.names=F, col.names=F, quote=F)

# 2. Lanciamo DIAMOnD per l'ultima volta (200 iterazioni)
# Usiamo la rete STRING che è la più performante
message("Generazione dei 200 geni candidati definitivi...")
cmd_final <- paste(shQuote(python_exe), "DIAMOnD.py string_network.txt final_ms_seeds.txt 200 1 diamond_final_results.txt")
system(cmd_final)

# 3. Carichiamo i risultati in R
diamond_final <- read_tsv("diamond_final_results.txt", comment="#", 
                          col_names=c("rank","gene","p"), show_col_types=F)

# Visualizza i primi 10 geni candidati "scoperti"
print(head(diamond_final, 10))

# ENRICHMENT ANALYSIS
# 2. Conversione dei nomi dei geni candidati (DIAMOnD) in ENTREZ ID
genes_to_test <- diamond_final$gene

conversion <- bitr(genes_to_test, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

# 3. Esecuzione dell'analisi GO per i Biological Processes (BP)
go_enrich <- enrichGO(gene          = conversion$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

# 4. Visualizzazione dei risultati
dotplot(go_enrich, showCategory=15) + 
  ggtitle("Processi Biologici dei 200 Candidati DIAMOnD")


# PART 3.2
# 1. Conversione dei geni seme originali (usa la lista iniziale 'seed_genes')
conversion_seeds <- bitr(seed_genes, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)

# 2. Arricchimento GO per i semi
go_seeds <- enrichGO(gene          = conversion_seeds$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     readable      = TRUE)

# 3. Visualizzazione
dotplot(go_seeds, showCategory=15) + 
  ggtitle("Processi Biologici dei Geni Seme (MS Originali)")


# LAST PART OF POINT 3
# 1. Selezioniamo i primi 10 geni candidati (i "vincitori" di DIAMOnD)
top_10_candidates <- head(diamond_final$gene, 10)

# 2. Identifichiamo i semi originali presenti in STRING
seeds_in_network <- intersect(seed_genes, V(g_string_lcc)$name)

# 3. Estraiamo il sottografo
# Prendiamo i top 10 candidati e i semi a loro collegati
nodes_to_show <- c(top_10_candidates, seeds_in_network)
subgraph <- induced_subgraph(g_string_lcc, V(g_string_lcc)$name %in% nodes_to_show)

# Per rendere il grafico leggibile, teniamo solo i nodi che hanno almeno un legame
# con i nostri top 10 candidati
neighbors_of_top <- neighbors(g_string_lcc, V(g_string_lcc)$name %in% top_10_candidates)
final_nodes <- c(top_10_candidates, intersect(seeds_in_network, names(neighbors_of_top)))
subgraph_final <- induced_subgraph(g_string_lcc, V(g_string_lcc)$name %in% final_nodes)

# 4. Coloriamo i nodi: Blu per i Semi, Rosso per i Nuovi Candidati
V(subgraph_final)$color <- ifelse(V(subgraph_final)$name %in% top_10_candidates, "red", "skyblue")
V(subgraph_final)$size <- ifelse(V(subgraph_final)$name %in% top_10_candidates, 8, 4)

# 5. Disegniamo la rete
plot(subgraph_final, 
     vertex.label = ifelse(V(subgraph_final)$name %in% top_10_candidates, V(subgraph_final)$name, NA),
     vertex.label.color = "black",
     vertex.label.cex = 0.8,
     edge.arrow.size = 0.2,
     main = "Modulo Sclerosi Multipla: Semi (Blu) e Top 10 Candidati (Rossi)")
