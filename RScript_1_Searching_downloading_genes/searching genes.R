# Definir o diretório de trabalho
setwd("C:/Users/Emamnuell/Desktop/panoptosis/genes/genes_24_optosis")

# Instalar e carregar o pacote 'rentrez' se ainda não estiver instalado
if (!requireNamespace("rentrez", quietly = TRUE)) {
  install.packages("rentrez")
}

# Instalar e carregar o pacote 'rio' se ainda não estiver instalado
if (!requireNamespace("rio", quietly = TRUE)) {
  install.packages("rio")
}

library(rentrez)
library(rio)

# vias <- import("Integrated_Cell_Death_Definitions_cut.xlsx")
# 
# # Termos de pesquisa
# search_terms <- vias$`Cell Death Type`[1:19]
search_terms <- c("Apoptosis", "Autophagy", "Ferroptosis", "Cellular senescence", 
                  "Mitotic catastrophe", "Necroptosis", "Pyroptosis", "Anoikis", 
                  "Parthanatos", "Lysosome-dependent cell death", "Necrosis", 
                  "NETosis", "Entosis", "Mitoptosis", "Methuosis", "Paraptosis", 
                  "Alkaliptosis", "Oxeiptosis", "Cuproptosis", "Autosis", "Erebosis",
                  "Efferocytosis", "Immunogenic cell death", 
                  "Mitochondrial permeability transition", "Disulfidptosis")

# Definir as colunas desejadas
cols_to_keep <- c("uid", "name", "description", "chromosome", "summary")

# Inicializar uma lista para armazenar os termos que não retornaram resultados
termos_sem_resultados <- list()

# Inicializar uma lista para armazenar os dataframes com duplicatas
df_list_duplicatas <- list()

# Loop sobre os termos de pesquisa
for (term in search_terms) {
  # Mensagem indicando o termo atual
  cat("Pesquisando termo:", term, "\n")
  
  # Inicializar uma lista para armazenar os data frames de cada categoria de pesquisa
  gene_info_list <- list()
  
  # Loop sobre os tipos de pesquisa (cancer e tumor)
  for (category in c("cancer", "tumor")) {
    # Utilizar tryCatch para lidar com erros caso nenhum gene seja encontrado
    tryCatch({
      # Realizar uma pesquisa por palavra-chave para encontrar o gene de interesse
      search_results <- entrez_search(db = "gene", term = paste(term, "[ALL] AND Homo sapiens[Organism] AND", category, "NOT Neanderthal NOT Denisova"), retmax = 20000)
      
      # Verificar se a pesquisa retornou IDs de genes
      if (length(search_results$ids) > 0) {
        # Dividir a lista de IDs de genes em partes menores
        ids_chunks <- split(search_results$ids, ceiling(seq_along(search_results$ids) / 100))
        
        # Iterar sobre cada parte da lista de IDs de genes
        for (ids_chunk in ids_chunks) {
          # Obter os resumos dos genes para a parte atual
          gene_summaries <- entrez_summary(db = "gene", id = ids_chunk)
          
          # Adicionar os resumos dos genes à lista
          gene_info_list <- c(gene_info_list, gene_summaries)
          
          # Aguardar 1 segundo antes da próxima solicitação para evitar exceder o limite de taxa
          Sys.sleep(1)
        }
      }
    }, error = function(e) {
      cat("Erro ao buscar genes para o termo:", term, "\n")
    })
  }
  
  # Verificar se há resultados para pelo menos uma categoria
  if (length(gene_info_list) > 0) {
    # Combinar os resumos dos genes em um único data frame
    gene_info <- do.call(rbind, gene_info_list)
    
    # Subconjunto das colunas desejadas, mantendo apenas as colunas existentes
    cols_to_keep <- intersect(cols_to_keep, colnames(gene_info))
    gene_info <- gene_info[, cols_to_keep]
    
    # Adicionar coluna "Category_Search" com valores correspondentes
    num_rows <- nrow(gene_info)
    Category_Search <- rep(c("cancer", "tumor"), length.out = num_rows)
    
    # Adicionar a coluna "Category_Search" ao data frame
    gene_info <- cbind(gene_info, Category_Search)
    
    # Atribuir o objeto gene_info ao ambiente de trabalho
    assign(paste("gene_info_", gsub(" ", "_", term), sep = ""), gene_info, envir = .GlobalEnv)
    
    # Converter gene_info_Apoptosis para um data frame
    assign(paste("gene_info_", gsub(" ", "_", term), "_df", sep = ""), as.data.frame(get(paste("gene_info_", gsub(" ", "_", term), sep = ""))), envir = .GlobalEnv)
    
    # Adicionar a coluna "Pathway_Categorization" com os valores correspondentes
    Pathway_Categorization <- rep(term, nrow(get(paste("gene_info_", gsub(" ", "_", term), "_df", sep = ""))))
    
    # Adicionar a coluna "Pathway_Categorization" ao data frame
    assign(paste("gene_info_", gsub(" ", "_", term), "_df", sep = ""), cbind(get(paste("gene_info_", gsub(" ", "_", term), "_df", sep = "")), Pathway_Categorization), envir = .GlobalEnv)
    
    # Adicionar a coluna "Category_Search" ao data frame
    assign(paste("gene_info_", gsub(" ", "_", term), "_df", sep = ""), cbind(get(paste("gene_info_", gsub(" ", "_", term), "_df", sep = "")), Category_Search), envir = .GlobalEnv)
    
    # Adicionar o data frame à lista de dataframes com duplicatas
    df_list_duplicatas <- c(df_list_duplicatas, list(get(paste("gene_info_", gsub(" ", "_", term), "_df", sep = ""))))
  } else {
    # Se nenhum resultado for encontrado, adicionar o termo à lista de termos sem resultados
    termos_sem_resultados <- c(termos_sem_resultados, term)
    cat("Nenhum resultado encontrado para o termo:", term, "\n")
  }
}


# Salvar a lista de termos sem resultados em um arquivo TSV
termos_sem_resultados_df <- data.frame(Termos_Sem_Resultados = termos_sem_resultados)
rio::export(termos_sem_resultados_df, "termos_sem_resultados.tsv", sep = "\t")

# Verificar se todos os data frames na lista têm as mesmas colunas
colunas <- lapply(df_list_duplicatas, colnames)
todas_colunas <- unique(unlist(colunas))

# Ajustar os data frames para terem as mesmas colunas
df_list_duplicatas <- lapply(df_list_duplicatas, function(df) {
  missing_cols <- setdiff(todas_colunas, colnames(df))
  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      df[[col]] <- NA
    }
  }
  df <- df[todas_colunas]
  return(df)
})

# Combinar os data frames em um único dataframe com duplicatas
all_genes_duplicatas <- do.call(rbind, df_list_duplicatas)

# Exportar o dataframe com duplicatas como um arquivo TSV
rio::export(all_genes_duplicatas, "all_genes_duplicatas.tsv", sep = "\t")

# Remover as duplicatas na coluna "name"
all_genes_sem_duplicatas <- all_genes_duplicatas[!duplicated(all_genes_duplicatas$name), ]

# Exportar o dataframe sem duplicatas em "name" como um arquivo TSV
rio::export(all_genes_sem_duplicatas, "all_genes_sem_duplicatas.tsv", sep = "\t")
