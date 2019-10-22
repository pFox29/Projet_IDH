#!/usr/bin/env Rscript
    library(tidyverse)
    library(igraph)
    library(optparse)
    library(psych)
    library(ggpubr)
    library(ggplot2)
    library(ggthemes)
    library(optparse)
#Param
    option_list = list(
        make_option(c("-i", "--input_path"), type="character", default="/home/guest/Desktop/Cheryn_ALI_Projet_IDH/Premiere_partie/DATA/", 
                help="input path", metavar="character"),
        
        make_option(c("-o", "--output_path"), type="character", default="/home/guest/Desktop/Cheryn_ALI_Projet_IDH/Premiere_partie/RESULTATS/", 
                help="output path", metavar="character")               
    ); 
    
    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);
#Paths
    input_path = opt$input_path 
    output_path = opt$output_path
    setwd(input_path)



#Functions
    ggboxplot_simple= function(vector, vector_name="", notch=TRUE,title="", sub="",ylab="",ylim=NULL, legend = TRUE, tick_angle=180, x_abline=NULL ){
        data = data.frame(group = vector_name, value = vector)
        plot.data <- rbind(data)

        plot = ggplot(plot.data) +
                geom_boxplot(aes(x=group, y=value, fill=group),notch = notch)+
                    ggtitle( title, subtitle = sub )+
                    scale_y_continuous(name = ylab, limits=ylim)+
                    
                    theme(  plot.title = element_text(hjust=0.5, vjust=0.5),
                            plot.subtitle = element_text(size = 8),
                            axis.title.x = element_blank(),
                            axis.text.x = element_blank()
                           
                            )  

        if(!(is.null(x_abline))){
            plot = plot + geom_vline(xintercept = x_abline, colour = "red")
        }
        
        if(legend == FALSE){
            plot= plot + theme(  plot.title = element_text(hjust=0.5, vjust=0.5),
                            plot.subtitle = element_text(size = 8),
                            axis.title.x = element_blank(),
                            legend.position = "none",
                            #axis.text.x = element_text(angle = tick_angle)
                            axis.text.x = element_blank()

                            ) 
        }
        return(plot)
        }


#Load dataFrames
    #Direct and implicit
        GOterm_prot_dir_imp_dataFrame = readRDS(file= paste0(input_path, "GOterm_prot_dir_imp_dataFrame.RDS" ))
        nbprot_per_goterm_dir_imp_dataFrame = readRDS(file= paste0(input_path, "nbprot_per_goterm_dir_imp_dataFrame.RDS" ))
        nbgoterm_per_prot_dir_imp_dataFrame = readRDS(file= paste0(input_path, "nbgoterm_per_prot_dir_imp_dataFrame.RDS" ))
        nbprot_nbgoterm_per_prot_dir_imp_dataFrame = readRDS(file= paste0(input_path, "nbprot_nbgoterm_per_prot_dir_imp_dataFrame.RDS" ))
    #Direct
        GOterm_prot_dir_dataFrame = readRDS(file= paste0(input_path, "GOterm_prot_dir_dataFrame.RDS" ))
        nbprot_per_goterm_dir_dataFrame = readRDS(file= paste0(input_path, "nbprot_per_goterm_dir_dataFrame.RDS" ))
        nbgoterm_per_prot_dir_dataFrame = readRDS(file= paste0(input_path, "nbgoterm_per_prot_dir_dataFrame.RDS" ))
        nbprot_nbgoterm_per_prot_dir_dataFrame = readRDS(file= paste0(input_path, "nbprot_nbgoterm_per_prot_dir_dataFrame.RDS" ))

        nbprot_nbgoterm_per_prot_dir_imp_dataFrame= dfOrder(nbprot_nbgoterm_per_prot_dir_imp_dataFrame,1)
        nbprot_nbgoterm_per_prot_dir_dataFrame= dfOrder(nbprot_nbgoterm_per_prot_dir_dataFrame,1)

#Major GO terms

    go_molecular_function = "GO:0003674"
    go_biological_process = "GO:0008150"
    go_cellular_component = "GO:0005575"

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

#nbr of goterms per proteines
    nbr_not_annotated_prot_dir = length(which(nbgoterm_per_prot_dir_dataFrame[,2]==0))
    nbgoterm_per_prot_dir_filt_dataFrame = nbgoterm_per_prot_dir_dataFrame[which(nbgoterm_per_prot_dir_dataFrame[,2]>0),]
    nbgoterm_per_prot_dir_filt_dataFrame = dfOrder(nbgoterm_per_prot_dir_filt_dataFrame,2)
   
    nbgoterm_per_prot_dir_boxplot = ggboxplot_simple(as.vector(nbgoterm_per_prot_dir_filt_dataFrame[,2]), "nbGo_per_prot", notch=FALSE,
        title="Nbr of direct annotation per protein", 
        sub=paste0("Nbr of protein = ", nrow(nbgoterm_per_prot_dir_filt_dataFrame), " (",nbr_not_annotated_prot_dir," non annotated proteines)  Mean = ",mean(nbgoterm_per_prot_dir_filt_dataFrame[,2]) ),
        ylim=c(0,180), legend = FALSE, tick_angle=0 )
   
    nbprot_per_goterm_dir_density = ggdensity(nbgoterm_per_prot_dir_filt_dataFrame, x = "nbGO",
        add = "mean", rug = TRUE,
        title = "Density Nbr of direct annotation per protein")

    nbr_not_annotated_prot_dir_imp = length(which(nbgoterm_per_prot_dir_imp_dataFrame[,2]==0))
    nbgoterm_per_prot_dir_imp_filt_dataFrame = nbgoterm_per_prot_dir_imp_dataFrame[which(nbgoterm_per_prot_dir_imp_dataFrame[,2]>0),]
    nbgoterm_per_prot_dir_imp_filt_dataFrame = dfOrder(nbgoterm_per_prot_dir_imp_filt_dataFrame,2)
    
    nbgoterm_per_prot_dir_imp_boxplot = ggboxplot_simple(as.vector(nbgoterm_per_prot_dir_imp_filt_dataFrame[,2]), "nbGo_per_prot", notch=FALSE,
        title="Nbr of direct and implicit annotation per protein", 
        sub=paste0("Nbr of protein = ", nrow(nbgoterm_per_prot_dir_imp_filt_dataFrame), " (",nbr_not_annotated_prot_dir_imp," non annotated proteines)  Mean = ", mean(nbgoterm_per_prot_dir_imp_filt_dataFrame[,2])),
        ylim=c(0,180), legend = FALSE, tick_angle=0 )
    nbprot_per_goterm_dir_imp_density = ggdensity(nbgoterm_per_prot_dir_imp_filt_dataFrame, x = "nbGO",
        add = "mean", rug = TRUE,
        title = "Density Nbr of direct and implicit annotation per protein")

    pdf(file=paste0(output_path,"nbgoterm_per_prot.pdf"))
        print(nbprot_per_goterm_dir_density)
        print(nbgoterm_per_prot_dir_boxplot)

        print(nbprot_per_goterm_dir_imp_density)
        print(nbgoterm_per_prot_dir_imp_boxplot)

    dev.off()
    
 
#nbr of per proteines per goterms

    nbrGO_whitout_data_dir = length(which(nbprot_per_goterm_dir_dataFrame[,3]==0))
    nbprot_per_goterm_dir_filt_dataFrame = nbprot_per_goterm_dir_dataFrame[which(nbprot_per_goterm_dir_dataFrame[,3]>0),] 
    nbprot_per_goterm_dir_filt_dataFrame = dfOrder(nbprot_per_goterm_dir_filt_dataFrame,3) 
    
    nbprot_per_goterm_dir_boxplot = ggboxplot_simple(as.vector(nbprot_per_goterm_dir_filt_dataFrame[,3]), "nbProt_per_GOterm", notch=FALSE,
        title="Nbr of protein directly annotated per term", 
        sub=paste0("Nbr of terms = ", nrow(nbprot_per_goterm_dir_filt_dataFrame), " (",nbrGO_whitout_data_dir," obselete term) Mean = ", mean(nbprot_per_goterm_dir_filt_dataFrame[,3]) ),
        ylim=c(0,25), legend = FALSE, tick_angle=0 )
    nbprot_per_goterm_dir_density = ggdensity(nbprot_per_goterm_dir_filt_dataFrame, x = "prot",
        add = "mean", rug = TRUE,
        title = "Density Nbr of protein directly annotated per term",
        color = "go_type", 
        palette = c("#00AFBB", "#E7B800","#F6A3A3"))
            
    nbrGO_whitout_data_dir_imp = length(which(nbprot_per_goterm_dir_imp_dataFrame[,3]==0))
    nbprot_per_goterm_dir_imp_filt_dataFrame = nbprot_per_goterm_dir_imp_dataFrame[which(nbprot_per_goterm_dir_imp_dataFrame[,3]>0),]
    nbprot_per_goterm_dir_imp_filt_dataFrame = dfOrder(nbprot_per_goterm_dir_imp_filt_dataFrame,3)
    
    nbprot_per_goterm_dir_imp_boxplot = ggboxplot_simple(as.vector(nbprot_per_goterm_dir_imp_filt_dataFrame[,3]), "nbGo_per_prot", notch=FALSE,
        title="Nbr of protein directly and implicitly annotated per term", 
        sub=paste0("Nbr of protein = ", nrow(nbprot_per_goterm_dir_imp_filt_dataFrame), " (",nbrGO_whitout_data_dir_imp," obselete term) Mean = ", mean(nbprot_per_goterm_dir_imp_filt_dataFrame[,3]) ),
        ylim=c(0,25), legend = FALSE, tick_angle=0 )
    nbprot_per_goterm_dir_imp_density = ggdensity(nbprot_per_goterm_dir_imp_filt_dataFrame, x = "prot",
        add = "mean", rug = TRUE,
        title = "Density Nbr of protein directly and implicitly annotated per term",
        color = "go_type", 
        palette = c("#00AFBB", "#E7B800","#F6A3A3"))

    pdf(file=paste0(output_path,"nbprot_per_goterm.pdf"))
        print(nbprot_per_goterm_dir_density)
        print(nbprot_per_goterm_dir_boxplot)

        print(nbprot_per_goterm_dir_imp_density)
        print(nbprot_per_goterm_dir_imp_boxplot)

    dev.off()
    
 
 
 
