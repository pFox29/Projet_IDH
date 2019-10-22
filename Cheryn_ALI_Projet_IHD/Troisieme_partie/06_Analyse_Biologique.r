#Cheryn ALI
#Octobre 2019
# Les genes co exprimés a travers differents type cellulaire sont-ils impliqué dans un meme preocessus biologique et/ou localisé preferentiellement sur le meme chromosome?

#Librairie
    library("tidyverse")
    library("gplots")
    library("ggplot2")
    library("ggdendro")
    library("limma")
    library("made4")
    library("cluster")

#Paths
    input_path = "/home/guest/Desktop/Cheryn_ALI_Projet_IHD/Troisieme_partie/DATA/"
    output_path = "/home/guest/Desktop/Cheryn_ALI_Projet_IHD/Troisieme_partie/RESULTATS/"
    setwd(input_path)

#Functions
    fish_pval = function(liste1, liste2, pop_length= NULL,test_side=NULL){
        commun = length(which(liste1 %in% liste2))
        liste1_only = length(liste1) - commun
        liste2_only = length(liste2) - commun

        if(is.null(pop_length)){
            rest_pop = 0
        }
        if(!(is.null(pop_length))){
            rest_pop = pop_length - (liste1_only +liste2_only +commun)
        }

        fish_mat = matrix(c(liste1_only,liste2_only,commun,rest_pop), nrow = 2)
       
        pval = fisher.test(fish_mat)$p.value
        
        if(!(is.null(test_side))){
            pval = fisher.test(fish_mat, alternative = test_side)$p.value
        }

        return(pval)
    }


#Data - Normalisation - Filtration
    raw_data_all =  read.csv(paste0(input_path,'DM_RH_RNA_Seq_FPKM_expression_matrix_for_DM_v4.03_13dec2013_desc.csv'),sep='\t', skipNul = TRUE,header=T) #http://solanaceae.plantbiology.msu.edu/pgsc_download.shtml # Données en FPKM
        
    raw_data = raw_data_all[,c(1:4,45:53)] # Recuperation des données sur les tissus

    for( i in 5:ncol(raw_data)){
        raw_data[,i] = as.numeric(sub(",", ".", sub(".", "", raw_data[,i], fixed=TRUE), fixed=TRUE))
        
    }
    dim(raw_data) #39000 genes

    #Filtrattion
        data_filt = raw_data[-which(rowSums(raw_data[,5:13])==0),] #filtration des genes exprimé dans aucune condition
        dim(data_filt) #30162 genes

        total_exp = rowSums(raw_data[,5:13])
        seuil_025 = quantile(total_exp,0.25)
        data_filt_qu025 = raw_data[-which(rowSums(raw_data[,5:13])<seuil_025),] #filtration des gens faiblement exprimé
        dim(data_filt_qu025) #29250 genes

        seuil_090 = quantile(total_exp,0.90)
        data_filt_qu025_qu090 = data_filt_qu025[-which(rowSums(data_filt_qu025[,5:13])>seuil_090),] #filtration des genes trop fortement exprimé
        dim(data_filt_qu025_qu090) #25350 genes

        SD = apply(data_filt_qu025_qu090[,5:13], 1, function(x) sd(x))
        seuil_sd = quantile(SD, 0.75)
        data_filt_qu025_qu090_sd75 = data_filt_qu025_qu090[which(SD>seuil_sd),] #filtration des genes les plus variables
        dim(data_filt_qu025_qu090_sd75) #6338 genes

       

    #Normalisation
        data_filt_qu025_qu090_sd75_norm = data_filt_qu025_qu090_sd75
        data_filt_qu025_qu090_sd75_norm[,5:13] = apply(data_filt_qu025_qu090_sd75[,5:13],2, function(x) x/quantile(x,0.75)) #normalisation par le 3eme quantile
        data_filt_qu025_qu090_sd75_norm_log = data_filt_qu025_qu090_sd75_norm
        data_filt_qu025_qu090_sd75_norm_log[,5:13] = log2(data_filt_qu025_qu090_sd75_norm[,5:13]+1) #normalisation log
        rownames(data_filt_qu025_qu090_sd75_norm_log) = data_filt_qu025_qu090_sd75_norm_log[,1]
            #plot(density(as.matrix(data_filt_qu025_qu090_sd75_norm_log[,5:13])))
        data_filt_norm = data_filt_qu025_qu090_sd75_norm_log[,5:13]
        rownames(data_filt_norm) = data_filt_qu025_qu090_sd75_norm_log[,1]
        tdata_filt_norm = as.data.frame(t(data_filt_norm))
    
#Visualisation des données
    ag='ward.D2' #Methode d'agglomeration
    corrMethods='pearson'#Methode de correlation
    
    titre=paste0('Heatmap Données (avec Dendrogramme sur les tissues), Ag=',ag,', Cor=',corrMethods )
    pdf(file=paste0(output_path, titre,'.pdf'))
        heatmap.2(as.matrix(data_filt_norm),
                    trace     = "none",
                    scale     = "row",
                    col       = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
                    hclustfun = function(x) hclust(x,method = ag), 
                    distfun = function(x) as.dist(1-cor( t(x), method=corrMethods)),
                    labRow = "",
                    dendrogram = "both",
                    Colv = T,
                    main=titre
                    )
    dev.off()

    distance_mat = as.dist(1-cor(tdata_filt_norm, method=corrMethods)) #matrice de distance
    dendro=hclust(distance_mat,method=ag)
    
        # titre=paste0('Dendrogramme Ag=',ag,', Cor=',corrMethods)
        # pdf(file=paste0(output_path, titre,'.pdf'))
        #     plot(dendro,main=titre,labels=F,xlab='Genes')
        #     rect.hclust(dendro, k=5,border="tomato")
        # dev.off()

#Clustering
    k =  5 #nombre de clusters
    clusters=cutree(dendro, k)
    
    pdf(file=paste0(output_path,'Clustering Donnees Ag=',ag,', Dist=1-cor(',corrMethods,'), K=',k,'.pdf'))
        nofclust.height <-  length(unique(as.vector(clusters)));
        hmcols <- rev(redgreen(2750))
        selcol <- colorRampPalette(brewer.pal(12,"Set3"))
        selcol2 <- colorRampPalette(brewer.pal(9,"Set1"))
        clustcol.height = selcol2(nofclust.height);
        data_filt_norm$class=clusters
        hm=heatmap.2(as.matrix(data_filt_norm[,-10]),
                    trace     = "none",
                    scale     = "row",
                    col       = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
                    hclustfun = function(x) hclust(x,method = ag),
                    distfun = function(x) as.dist(1-cor( t(x), method=corrMethods)),
                    labRow = "",
                    RowSideColors=clustcol.height[data_filt_norm$class],
                    Colv = T,
                    dendrogram='both')
        
        for (cl in 1:k){
            moncluster = subset(data_filt_norm,data_filt_norm$class==cl)
            plot(x=seq(1,length(data_filt_norm)-1),y=moncluster[1,1:(length(data_filt_norm)-1)],type="l",xlim=c(1,length(data_filt_norm)-1),col='lightgrey', ylim = c(min(moncluster),max(moncluster)), main=paste('Cluster: ',cl,'nb',sep=""))
            for (gene in 2:nrow(moncluster)){
            points(x=seq(1,length(data_filt_norm)-1),y=moncluster[gene,1:(length(data_filt_norm)-1)], type ='l', col ='lightgrey', pch = 16)
            }
            meanclust=apply(moncluster, 2, function(x) mean(x))
            points(x=seq(1,length(data_filt_norm)-1),y=meanclust[1:(length(data_filt_norm)-1)], type ='l', col ='tomato', pch = 16)
        
        }
    dev.off()

    data_filt_norm$class=clusters
    clusterings=list()
    for (i in seq(1,k)){
        moncluster = subset(rownames(data_filt_norm),data_filt_norm$class==i)
        clusterings = c(clusterings, list(moncluster))
    }
    # saveRDS(file= paste0(output_path, "clusterings_list.RDS" ),  clusterings)
    # saveRDS(file= paste0(output_path, "distance_mat.RDS" ),  distance_mat)
    # saveRDS(file= paste0(output_path, "data_filt_norm_DF.RDS" ), data_filt_norm)

#Listes de noms de genes
    list_cluster = list(
        GN_cluster_1 = rownames(data_filt_norm[which(data_filt_norm[,10] == 1),]), #1208
        GN_cluster_2 = rownames(data_filt_norm[which(data_filt_norm[,10] == 2),]), #1355
        GN_cluster_3 = rownames(data_filt_norm[which(data_filt_norm[,10] == 3),]), #1389
        GN_cluster_4 = rownames(data_filt_norm[which(data_filt_norm[,10] == 4),]), #901
        GN_cluster_5 = rownames(data_filt_norm[which(data_filt_norm[,10] == 5),]) #1485
        )

    list_chr = list(
        GN_chr01 = rownames(data_filt_qu025_qu090_sd75_norm_log[which(data_filt_qu025_qu090_sd75_norm_log[,2] == "chr01"),]), #766    
        GN_chr02 = rownames(data_filt_qu025_qu090_sd75_norm_log[which(data_filt_qu025_qu090_sd75_norm_log[,2] == "chr02"),]), #668
        GN_chr03 = rownames(data_filt_qu025_qu090_sd75_norm_log[which(data_filt_qu025_qu090_sd75_norm_log[,2] == "chr03"),]), #655
        GN_chr04 = rownames(data_filt_qu025_qu090_sd75_norm_log[which(data_filt_qu025_qu090_sd75_norm_log[,2] == "chr04"),]), #551
        GN_chr05 = rownames(data_filt_qu025_qu090_sd75_norm_log[which(data_filt_qu025_qu090_sd75_norm_log[,2] == "chr05"),]), #452
        GN_chr06 = rownames(data_filt_qu025_qu090_sd75_norm_log[which(data_filt_qu025_qu090_sd75_norm_log[,2] == "chr06"),]), #539
        GN_chr07 = rownames(data_filt_qu025_qu090_sd75_norm_log[which(data_filt_qu025_qu090_sd75_norm_log[,2] == "chr07"),]), #458
        GN_chr08 = rownames(data_filt_qu025_qu090_sd75_norm_log[which(data_filt_qu025_qu090_sd75_norm_log[,2] == "chr08"),]), #387
        GN_chr09 = rownames(data_filt_qu025_qu090_sd75_norm_log[which(data_filt_qu025_qu090_sd75_norm_log[,2] == "chr09"),]), #459
        GN_chr10 = rownames(data_filt_qu025_qu090_sd75_norm_log[which(data_filt_qu025_qu090_sd75_norm_log[,2] == "chr10"),]), #459                        
        GN_chr11 = rownames(data_filt_qu025_qu090_sd75_norm_log[which(data_filt_qu025_qu090_sd75_norm_log[,2] == "chr11"),]), #370                        
        GN_chr12 = rownames(data_filt_qu025_qu090_sd75_norm_log[which(data_filt_qu025_qu090_sd75_norm_log[,2] == "chr12"),]), #448                        
        GN_chr00 = rownames(data_filt_qu025_qu090_sd75_norm_log[which(data_filt_qu025_qu090_sd75_norm_log[,2] == "chr00"),]) #116
        )                            
    
    sink(paste0(output_path,"GN_clusters.txt"))
        for(cl in 1:length(list_cluster)){
            query = ""
            for(gn in 1:length(list_cluster[[cl]])){
                query = paste0(list_cluster[[cl]][[gn]], " ",query)
            } 
            cat(paste0(query, "\n"))

            }
    sink()

#Test d'enrichissement cluster - chromosome

    pop_len = nrow(data_filt_qu025_qu090_sd75_norm_log)

    fish_res = data.frame(
        cluster_1 = rep(5,13),
        cluster_2 = rep(5,13),
        cluster_3 = rep(5,13),
        cluster_4 = rep(5,13),
        cluster_5 = rep(5,13)
    )

    rownames(fish_res) = c("chr01","chr02","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12","chr00")

    commun_elements = data.frame(
        cluster_1 = rep(0,13),
        cluster_2 = rep(0,13),
        cluster_3 = rep(0,13),
        cluster_4 = rep(0,13),
        cluster_5 = rep(0,13)
    )
    rownames(commun_elements) = c("chr01","chr02","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12","chr00")

    for(cl in 1:length(list_cluster)){
        for(ch in 1: length(list_chr)){
            lcl = as.list(list_cluster[[cl]])
            lch = as.list(list_chr[[ch]])
            fish_res[ch,cl] = fish_pval(list_chr[[ch]],list_cluster[[cl]], pop_length = pop_len, test_side="greater")
            commun_elements[ch,cl] = length(which(list_chr[[ch]] %in% list_cluster[[cl]]))
        }
    }

    commun_elements_t <- as.table(as.matrix(commun_elements))

    pdf(file=paste0(output_path, 'fish_res.pdf'))
        textplot(fish_res[,1:3])
        textplot(fish_res[,4:5])
        balloonplot(t(commun_elements_t), main ="Commun elements", xlab ="", ylab="",
                    label = TRUE, show.margins = FALSE)
        textplot(paste0("p.value chi2: ",chisq.test(commun_elements)$p.value), cex=1)
    dev.off()

