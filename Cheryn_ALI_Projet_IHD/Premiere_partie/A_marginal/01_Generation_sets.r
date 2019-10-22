#!/usr/bin/env Rscript
    library(STRINGdb)
    library(tidyverse)
    library(igraph)
    library(neo4r)
    library(optparse)
#Arguments
 
    option_list = list(
        make_option(c("-u", "--neo4j_url"), type="character", default="http://localhost:7474", 
                help="neo4j URL", metavar="character"),
        
        make_option(c("-r", "--neo4j_usr"), type="character", default="neo4j", 
                help="neo4j user ", metavar="character"),
        
        make_option(c("-w", "--neo4j_psw"), type="character", default="bioinfo", 
                help="neo4j password ", metavar="character"),
        
        make_option(c("-i", "--sp_id"), type="integer", default=4113, 
                help="Stringdb spices id ", metavar="integer"),
        
        make_option(c("-p", "--output_path"), type="character", default="/home/guest/Desktop/Cheryn_ALI_Projet_IDH/Premiere_partie/", 
                help="output pathway", metavar="character"),
        
        make_option(c("-s", "--Compute_Stat"), type="logical", default=TRUE, 
                help="Generate tables and graphs", metavar="logical")              
    ); 
    
    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);

#Paths
    output_path = opt$output_path
    setwd(opt$output_path)

#Recuperation des données de l'espece d'interet (ex: Solanum tuberosum id=4113)
    sp_id = opt$sp_id
    go_molecular_function = "GO:0003674"
    go_biological_process = "GO:0008150"
    go_cellular_component = "GO:0005575"
    #STRINGdb
        species = get_STRING_species(version = '10', species_name=NULL) %>% filter(species_id == sp_id)
        species_compact_name = species$compact_name
        species_compact_name = strsplit(species_compact_name," ")
        species_compact_name = paste0(species_compact_name[[1]][1],"_",species_compact_name[[1]][2])

        sdb = STRINGdb$new(version='10', species=sp_id, score_threshold = 400, input_directory = paste0(output_path,"DATA/"))
        
        annot = sdb$get_annotations()

        GO_annot = annot[which((annot$category == 'Function') | (annot$category == 'Process') | (annot$category == 'Component' )),]
        GO_terms_all_unique = unique(GO_annot$term_id) #1 478 GO terms

        sp = sdb$get_proteins() #39 021 proteines

    #Neo4J
        con = neo4j_api$new(
        url=opt$neo4j_url,
        user=opt$neo4j_usr,
        password = opt$neo4j_psw
        )    

#Generation du fichier .sets pour les relation direct et implicite
    print(paste0('Generation of GOterm_prot_dir_imp_',species_compact_name,'.sets'))
        GOterm_prot_dir_imp_list = list(
            go_term = GO_terms_all_unique,
            go_info = seq(1,length(GO_terms_all_unique)), #type + description
            prot = seq(1,length(GO_terms_all_unique)) #String of all prot associated
        )
        nbprot_per_goterm_dir_imp_list = list(
            go_term = GO_terms_all_unique,
            go_type = seq(1,length(GO_terms_all_unique)),
            prot = seq(1,length(GO_terms_all_unique)) #Nbr of prot associated
        )
        nbgoterm_per_prot_dir_imp_dataFrame = data.frame(
            prot = sp$protein_external_id,
            nbGO = rep(0, length(sp$protein_external_id)) #Nbr of GO associated
        )

        for(i in seq(1:length(GO_terms_all_unique))){
            goterm = GOterm_prot_dir_imp_list$go_term[i]
                print(goterm)
            query_dir_imp = paste0( "MATCH (g:GOTerm {acc: '",goterm,"'})-[*] -> (p:Protein) return p")
            proteins_dir_imp = query_dir_imp %>% call_neo4j(con)

            query_dir_imp_info = paste0( "MATCH (g:GOTerm {acc: '",goterm,"'}) return g.term_type, g.name")
            GO_dir_imp_info = query_dir_imp_info %>% call_neo4j(con)
            GO_dir_imp_type = GO_dir_imp_info$g.term_type$value
            GO_dir_imp_name = GO_dir_imp_info$g.name$value

            if(length(GO_dir_imp_type)==0){
                GO_dir_imp_type = "Unknown"
            }
                print(GO_dir_imp_type)
            GO_dir_imp_info = paste0(GO_dir_imp_type,": ",GO_dir_imp_name)

            GOterm_prot_dir_imp_list$go_info[i] = GO_dir_imp_info
            protein_names_dir_imp= unique(proteins_dir_imp$p$protein_external_id)
            pt_n=''
            if(length(protein_names_dir_imp) != 0){
                for(j in seq(1:length(protein_names_dir_imp))){  
                    nbgoterm_per_prot_dir_imp_dataFrame[which(nbgoterm_per_prot_dir_imp_dataFrame[,1] == protein_names_dir_imp[j]),2] = nbgoterm_per_prot_dir_imp_dataFrame[which(nbgoterm_per_prot_dir_imp_dataFrame[,1] == protein_names_dir_imp[j]),2] + 1
                    pt = strsplit(protein_names_dir_imp[j], paste0(sp_id,".") ) 
                    pt = pt[[1]][2]
                    pt_n = paste0(pt_n,pt,"\t")
                }
            }
            GOterm_prot_dir_imp_list$prot[i] = pt_n

            nbprot_per_goterm_dir_imp_list$go_type[i] = GO_dir_imp_type
            nbprot_per_goterm_dir_imp_list$prot[i] = length(protein_names_dir_imp)

                cat("Number of proteines: ",nbprot_per_goterm_dir_imp_list$prot[i],"\n")
        }



        GOterm_prot_dir_imp_dataFrame = as.data.frame(GOterm_prot_dir_imp_list)
        nbprot_per_goterm_dir_imp_dataFrame = as.data.frame(nbprot_per_goterm_dir_imp_list)
        
        Unknowns_dir_imp = c() #Deleting the GOterms without data
        for(i in seq(1:nrow(GOterm_prot_dir_imp_dataFrame))){
            if(GOterm_prot_dir_imp_dataFrame[i,2] == "Unknown: "){
                Unknowns_dir_imp = c(Unknowns_dir_imp, i)
            }
        }


        nbr_empty_go_dir_imp = length(Unknowns_dir_imp)
        GOterm_prot_dir_imp_dataFrame = GOterm_prot_dir_imp_dataFrame[-Unknowns_dir_imp,]
        

        sink(paste0(output_path,"RESULTATS/GOterm_prot_dir_imp_",species_compact_name,".sets"))
            cat("# format: sets \n")
            cat("# version: 1.2 \n")
            cat(paste0("# strain: ",species_compact_name,"\n"))
            cat(paste0("# date: ",Sys.Date(),"\n"))
            cat("# comment: Gene Ontology terms with directly and implicitly associated proteins (External IDs)")
                for(i in seq(1:nrow(GOterm_prot_dir_imp_dataFrame))){
                cat("\n")
                cat(paste0(GOterm_prot_dir_imp_dataFrame[i,1],"\t",GOterm_prot_dir_imp_dataFrame[i,2],"\t",GOterm_prot_dir_imp_dataFrame[i,3]))
                }
        sink()

        nbprot_nbgoterm_per_prot_dir_imp_dataFrame = data.frame(
            nbgoterm_per_prot = unique(nbgoterm_per_prot_dir_imp_dataFrame$nbGO),
            nbprot = unique(nbgoterm_per_prot_dir_imp_dataFrame$nbGO)
        )

        for (nb in unique(nbgoterm_per_prot_dir_imp_dataFrame$nbGO)){
            nbprot_nbgoterm_per_prot_dir_imp_dataFrame[which(nbprot_nbgoterm_per_prot_dir_imp_dataFrame[,1]==nb),2] = nrow(nbgoterm_per_prot_dir_imp_dataFrame[which(nbgoterm_per_prot_dir_imp_dataFrame[,2]==nb),])
        }

        if(opt$Compute_Stat == TRUE){
            saveRDS(file= paste0(output_path, "DATA/GOterm_prot_dir_imp_dataFrame.RDS" ),  GOterm_prot_dir_imp_dataFrame)
            saveRDS(file= paste0(output_path, "DATA/nbprot_per_goterm_dir_imp_dataFrame.RDS" ),  nbprot_per_goterm_dir_imp_dataFrame)
            saveRDS(file= paste0(output_path, "DATA/nbgoterm_per_prot_dir_imp_dataFrame.RDS" ),  nbgoterm_per_prot_dir_imp_dataFrame)
            saveRDS(file= paste0(output_path, "DATA/nbprot_nbgoterm_per_prot_dir_imp_dataFrame.RDS" ),  nbprot_nbgoterm_per_prot_dir_imp_dataFrame)
        }

    print('DONE')


#Generation du fichier .sets pour les relation direct
    print(paste0('Generation of GOterm_prot_dir_',species_compact_name,'.sets'))
        GOterm_prot_dir_list = list(
            go_term = GO_terms_all_unique,
            go_info = seq(1,length(GO_terms_all_unique)),   #type + description
            prot = seq(1,length(GO_terms_all_unique))   #String of all prot associated
        )
        nbprot_per_goterm_dir_list = list(
            go_term = GO_terms_all_unique,
            go_type = seq(1,length(GO_terms_all_unique)),
            prot = seq(1,length(GO_terms_all_unique))   #Nbr of prot associated
        )
        nbgoterm_per_prot_dir_dataFrame = data.frame(
            prot = sp$protein_external_id,
            nbGO = rep(0, length(sp$protein_external_id))   #Nbr of GO associated
        )
        
        for(i in seq(1:length(GO_terms_all_unique))){
            goterm = GOterm_prot_dir_list$go_term[i]
                print(goterm)

            query_dir = paste0( "MATCH (g:GOTerm {acc: '",goterm,"'})--> (p:Protein) return p")
            proteins_dir = query_dir %>% call_neo4j(con)

            query_dir_info = paste0( "MATCH (g:GOTerm {acc: '",goterm,"'}) return g.term_type, g.name")
            GO_dir_info = query_dir_info %>% call_neo4j(con)
            GO_dir_type = GO_dir_info$g.term_type$value
            GO_dir_name = GO_dir_info$g.name$value

            if(length(GO_dir_type)==0){
                GO_dir_type = "Unknown"
            }
            GO_dir_info = paste0(GO_dir_type,": ",GO_dir_name)

            GOterm_prot_dir_list$go_info[i] = GO_dir_info
            
            protein_names_dir= unique(proteins_dir$p$protein_external_id)
            pt_n=''
            if(length(protein_names_dir) != 0){
                for(j in seq(1:length(protein_names_dir))){
                    nbgoterm_per_prot_dir_dataFrame[which(nbgoterm_per_prot_dir_dataFrame[,1] == protein_names_dir[j]),2] = nbgoterm_per_prot_dir_dataFrame[which(nbgoterm_per_prot_dir_dataFrame[,1] == protein_names_dir[j]),2] + 1
                    pt = strsplit(protein_names_dir[j], paste0(sp_id,".") )
                    pt = pt[[1]][2] 
                    pt_n = paste0(pt_n,pt,"\t")
                }
            }
            GOterm_prot_dir_list$prot[i] = pt_n

            nbprot_per_goterm_dir_list$go_type[i] = GO_dir_type
            nbprot_per_goterm_dir_list$prot[i] = length(protein_names_dir)
                
                print(GO_dir_type)
                cat("Number of proteines: ",nbprot_per_goterm_dir_list$prot[i],"\n")
        }

        GOterm_prot_dir_dataFrame = as.data.frame(GOterm_prot_dir_list)
        nbprot_per_goterm_dir_dataFrame = as.data.frame(nbprot_per_goterm_dir_list)

        Unknowns_dir = c()
        for(i in seq(1:nrow(GOterm_prot_dir_dataFrame))){
            if(GOterm_prot_dir_dataFrame[i,2] == "Unknown: "){
                Unknowns_dir = c(Unknowns_dir, i)
            }
        }
        nbr_empty_go_dir = length(Unknowns_dir)
        GOterm_prot_dir_dataFrame = GOterm_prot_dir_dataFrame[-Unknowns_dir,]


        sink(paste0(output_path,"RESULTATS/GOterm_prot_dir_",species_compact_name,".sets"))
            cat("# format: sets \n")
            cat("# version: 1.2 \n")
            cat(paste0("# strain: ",species_compact_name,"\n"))
            cat(paste0("# date: ",Sys.Date(),"\n"))
            cat("# comment: Genes Ontology terms with directly associated proteins (External IDs)")
                for(i in seq(1:nrow(GOterm_prot_dir_dataFrame))){
                cat("\n")
                cat(paste0(GOterm_prot_dir_dataFrame[i,1],"\t",GOterm_prot_dir_dataFrame[i,2],"\t",GOterm_prot_dir_dataFrame[i,3]))
                }
        sink()



        nbprot_nbgoterm_per_prot_dir_dataFrame=data.frame(
            nbgoterm_per_prot = unique(nbgoterm_per_prot_dir_dataFrame$nbGO),
            nbprot = unique(nbgoterm_per_prot_dir_dataFrame$nbGO)  
        )

        for (nb in unique(nbgoterm_per_prot_dir_dataFrame$nbGO)){
            nbprot_nbgoterm_per_prot_dir_dataFrame[which(nbprot_nbgoterm_per_prot_dir_dataFrame[,1]==nb),2] = nrow(nbgoterm_per_prot_dir_dataFrame[which(nbgoterm_per_prot_dir_dataFrame[,2]==nb),])
        }

        total_prot_annot_GO_dir = sum(nbprot_nbgoterm_per_prot_dir_dataFrame[-(which(nbprot_nbgoterm_per_prot_dir_dataFrame[,1]==0)),2])

        if(opt$Compute_Stat == TRUE){
            saveRDS(file= paste0(output_path, "DATA/GOterm_prot_dir_dataFrame.RDS" ),  GOterm_prot_dir_dataFrame)
            saveRDS(file= paste0(output_path, "DATA/nbprot_per_goterm_dir_dataFrame.RDS" ),  nbprot_per_goterm_dir_dataFrame)
            saveRDS(file= paste0(output_path, "DATA/nbgoterm_per_prot_dir_dataFrame.RDS" ),  nbgoterm_per_prot_dir_dataFrame)
            saveRDS(file= paste0(output_path, "DATA/nbprot_nbgoterm_per_prot_dir_dataFrame.RDS" ),  nbprot_nbgoterm_per_prot_dir_dataFrame)
        }

    print('DONE')

# Info.txt
    molecular_function_dataframe=data.frame(
        molecular_function_annot = c('dir_imp','dir'),
        nb_prot=c(nbprot_per_goterm_dir_imp_dataFrame[which(nbprot_per_goterm_dir_imp_dataFrame[,1]==go_molecular_function),3],
                    nbprot_per_goterm_dir_dataFrame[which(nbprot_per_goterm_dir_dataFrame[,1]==go_molecular_function),3])
    )
    biological_process_dataframe=data.frame(
        biological_process_annot = c('dir_imp','dir'),
        nb_prot=c(nbprot_per_goterm_dir_imp_dataFrame[which(nbprot_per_goterm_dir_imp_dataFrame[,1]==go_biological_process),3],
                    nbprot_per_goterm_dir_dataFrame[which(nbprot_per_goterm_dir_dataFrame[,1]==go_biological_process),3])
    )
    cellular_component_dataframe=data.frame(
        cellular_component_annot = c('dir_imp','dir'),
        nb_prot=c(nbprot_per_goterm_dir_imp_dataFrame[which(nbprot_per_goterm_dir_imp_dataFrame[,1]==go_cellular_component),3],
                    nbprot_per_goterm_dir_dataFrame[which(nbprot_per_goterm_dir_dataFrame[,1]==go_cellular_component),3])
    )

    sink(paste0(output_path,"RESULTATS/General_Infos_",species_compact_name,".txt"))
        cat(paste0("#Strain: ",species_compact_name, "\n"))
        cat("#Annotations: \n")
        print(
            annot %>%
                select(category) %>%
                group_by(category) %>%
                summarise(count=n()) %>%
                arrange(desc(count)))
        cat('\n')
        cat(paste0('#Total GO annotations = ',length(GO_annot$term_id),'\n'))
        cat(paste0('#Total unique GO annotations = ',length(GO_terms_all_unique),'\n'))
        cat(paste0('#Number of terms without data = ', nbr_empty_go_dir_imp,'\n'))
        cat(paste0('#Total proteins = ',nrow(sp),'\n'))
        cat(paste0('#Total proteins annoated by at least 1 GO term = ', total_prot_annot_GO_dir ,'\n'))
        cat(paste0('#Total proteins annoated directly by "Molecular_function" GO term = ', molecular_function_dataframe[2,2],'\t directly and implicitly = ',molecular_function_dataframe[1,2] ,'\n'))
        cat(paste0('#Total proteins annoated directly by "Biological_process" GO term = ', biological_process_dataframe[2,2],'\t directly and implicitly = ',biological_process_dataframe[1,2] ,'\n'))
        cat(paste0('#Total proteins annoated directly by "Cellular_component" GO term = ', cellular_component_dataframe[2,2],'\t directly and implicitly = ',cellular_component_dataframe[1,2] ,'\n'))

    sink()

