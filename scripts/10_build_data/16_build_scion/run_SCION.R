#!/usr/bin/env Rscript
# VENDORED from MotrpacHumanPreSuspensionAnalysis/R/run_SCION.R.
#
# run_SCION() plus its nested SCION_infer() and the three internal helpers it reaches
# through the namespace - scion_data_processing(), RS.Get.Weight.Matrix() and RSGWM2() -
# copied so step 16 runs without the installed analysis package.
#
# CHANGES from upstream:
#   - make.names() on the matrix rownames is gone. It rewrote every rowname carrying a
#     hyphen, a space or a leading digit, while the differential-analysis table the
#     clustering keys on kept the original ids, so those features matched nothing, were
#     never clustered and never reached inference. The hyphen is in both isoform suffixes
#     and the "prot-ph" ome tag.
#   - the cluster lookup is a direct name match. Upstream rebuilds the sanitised names
#     with three regexes; scion_run_cmeans() now returns the matrix's own
#     feature..ome..tissue form.
#   - both data.frame() calls that build the edge tables pass check.names = FALSE. Their
#     columns are regulator feature ids that become the Regulator column, and the hub step
#     joins those back against myregdata's rownames - name repair severs that join, which
#     ends the run before write.table() and leaves a directory of per-cluster CSVs and no
#     network file.
#   - the DA table comes from the step-10 objects via .scion_load_da(), and the hub step's
#     feature-to-gene map from the step-07 object via .feature_to_gene().
#   - the five check_package_installation() calls are dropped. Stage 0 preflight checks
#     parallel, doParallel and randomForest when RUN_SCION=TRUE; pacman is not used.

run_SCION <- function (randomGroupCode = c("ADUResist", "ADUEndur"),
                       regulators,
                       targets,
                       permute = NULL,
                       dim = "col",
                       cluster = TRUE,
                       dir.name = "exp",
                       weightthreshold = 0,
                       normalize = FALSE,
                       num.cores = 1,
                       connect.hubs = TRUE,
                       verbose = TRUE) {
  message("Note: epigenetic targets are not currently supported")

  #create directory to save results
  if(!is.null(permute)){
    my.dir <- file.path(dir.name, as.character(permute))
    cat(paste("Starting permutation ",permute,"\n",sep=""))
  }else{
    my.dir <- dir.name
  }
  message(my.dir)
  if(!dir.exists(my.dir)) dir.create(my.dir, recursive = TRUE)

  cat("Processing data tables\n") #first need to pre-process the data
  tables = scion_data_processing(randomGroupCode = randomGroupCode,
                                 reg.mat = regulators,
                                 target.mat = targets,
                                 permute,
                                 dim)
  #Here instead of directly entering in the prefix at function start, I just
  #assign it according to the regulator syntax.
  #I just choose the first ome that shows up in the regulators in the first row
  prefix = unlist(strsplit(rownames(tables[["reg"]])[1], "..", fixed = TRUE))[2]
  if(grepl("metab", prefix)) prefix = "metab"

  #assign pre-processed data tables
  mytargetdata <- as.data.frame(tables$target)
  myregdata <- as.data.frame(tables$reg)
  # Rownames stay as the feature ids are. Sanitising them here is what severed the
  # join to the clustering and to the edge tables; see the header.

  #We need to manipulate the input structure to run properly.
  if(cluster){
    reg_tissue_omes = lapply(strsplit(rownames(myregdata), "\\.{3}|(?<!\\.)\\.\\.", perl=TRUE),
                             function(x) if (length(x) > 3) tail(x, 3) else x)

    target_tissue_omes = strsplit(rownames(mytargetdata), "\\.\\.")

    split_rownames_df = as.data.frame(do.call(rbind, c(reg_tissue_omes, target_tissue_omes)))
    colnames(split_rownames_df) = c("feature_id", "ome", "tissue")
    unique_pairs = split_rownames_df %>%
      dplyr::mutate(ome = dplyr::case_when(
        grepl("metab", ome) ~ "metab",
        TRUE ~ ome
      )) %>%
      dplyr::distinct(ome, tissue)
    #----building expected DA_list format for run_cmeans--------
    DA_Input_Cmeans = list()
    for(row in seq(nrow(unique_pairs))){
      selected_ome = unique_pairs[row,][["ome"]]
      selected_tissue = unique_pairs[row,][["tissue"]]
      selected_ome_mod = gsub("\\.", "-", selected_ome) #have to do this bc the rownames arent happy with "-"
      if(selected_tissue == "adipose" && (selected_ome_mod == "prot-pr" || selected_ome_mod == "prot-ph"))
        stop("Adipose Prot-ph, Prot-pr only have 2 timepoints. We don't recommend using these, it messes with the clustering")
      da_table = .scion_load_da(selected_tissue, selected_ome_mod)
      if(randomGroupCode == "ADUEndur") da_table = da_table %>%
        dplyr::filter(contrast_type == "exercise_with_controls") %>%
        dplyr::filter(contrast_category == "EE-CON")
      if(randomGroupCode == "ADUResist") da_table = da_table %>%
        dplyr::filter(contrast_type == "exercise_with_controls") %>%
        dplyr::filter(contrast_category == "RE-CON")

      filt_to_inputs = split_rownames_df %>%
        dplyr::mutate(ome = dplyr::case_when(
          grepl("metab", ome) ~ "metab",
          TRUE ~ ome
        )) %>%
        dplyr::filter(ome == selected_ome, tissue == selected_tissue) %>%
        dplyr::pull(feature_id)
      #i add "__" to make splitting the results and manipulating syntax a bit easier
      da_filt = da_table %>%
        dplyr::filter(feature_id %in% filt_to_inputs) %>%
        dplyr::mutate(feature_id = paste(feature_id, selected_tissue, sep = "__"))

      #expected format for c-means, cameraPR unnesting
      merged_name = paste(selected_tissue, selected_ome_mod, sep = ".")
      DA_Input_Cmeans[[merged_name]] = da_filt
    }
    clusters = scion_run_cmeans(DA_Input_Cmeans)

    cmeans_results = clusters[["cluster"]] #again we force everything to be one given tissue
    # A direct name match: scion_run_cmeans() returns names already in the
    # feature..ome..tissue form the regulator and target matrices carry.
    cluster_org = data.frame(rowname_format = names(cmeans_results),
                             cluster = as.integer(cmeans_results)) %>%
      tibble::column_to_rownames("rowname_format")
  }

  #infer a network using GENIE3 on each cluster
  SCION_infer <- function(mytargetdata,
                          myregdata,
                          clusterresults,
                          weightthreshold,
                          prefix,
                          permute=NULL,
                          normalize=FALSE,
                          verbose=F){

    finalnetwork <- data.frame(Regulator=character(),
                               Interaction=character(),
                               Target=character(),
                               Weight=double(),
                               Cluster=numeric(),
                               stringsAsFactors=FALSE)
    myhubs=NULL
    cat("Inferring cluster-specific networks.\n")
    for (i in 1:max(clusterresults$cluster)){
      #so here we need to reverse back into selecting the tissue
      mygenes = row.names(clusterresults)[clusterresults$cluster==i]
      clustertargetdata = mytargetdata[row.names(mytargetdata)%in%mygenes,]
      clusterregdata = myregdata[row.names(myregdata)%in%mygenes,]

      #we need at least one target and at least one regulator
      if (dim(clustertargetdata)[1]<1 | dim(clusterregdata)[1]<1){
        next
      }

      #infer the network
      if(verbose){
        cat(paste("Inferring network for cluster ",i,
                  " with ", dim(clusterregdata)[1], " regulators and ",
                  dim(clustertargetdata)[1], " targets.\n",sep=""))
      }
      network = RS.Get.Weight.Matrix(t(clustertargetdata),t(clusterregdata),normalize=normalize,num.cores=num.cores)
      #if network inference failed, move on
      if (is.null(network)){
        next
      }

      #now make a new network where we eliminate all the low confidence edges
      # check.names = FALSE: these columns are regulator feature ids and become the
      # Regulator column, which the hub step joins back against myregdata's rownames.
      trimmednet = data.frame(ifelse(network<weightthreshold,NaN,network),
                              check.names = FALSE)

      #translate the trimmed network into a table we can import into cytoscape
      networktable = data.frame(Regulator=character(), Interaction=character(), Target=character(), Weight=double(),Cluster=numeric(),
                                stringsAsFactors=FALSE)
      row=1
      for (j in 1:dim(trimmednet)[1]){
        for (k in 1:dim(trimmednet)[2]){
          #skip NaNs as these have no edge
          if (is.na(trimmednet[j,k])){
            next
          }else{
            networktable[row,] = cbind(colnames(trimmednet)[k],"regulates",rownames(trimmednet)[j],trimmednet[j,k],i)
            row = row+1
          }
        }
      }

      #save this network
      finalnetwork = rbind(finalnetwork,networktable)
      write.csv(finalnetwork,paste0(my.dir,"/network_cluster_",i,".csv"))

      #get the hub gene (highest outdegree) and save it to connect the clusters later
      #if there is a tie, we save both hubs
      myregs = unique(networktable$Regulator)
      for (j in 1:length(myregs)){
        numedges = sum(networktable$Regulator%in%myregs[j])
        if (j==1){
          hubedges = numedges
          hub = myregs[j]
        }else if (numedges>hubedges){
          hubedges = numedges
          hub = myregs[j]
        }else if (numedges==hubedges){
          hub = rbind(hub,myregs[j])
        }
      }
      if (!exists("myhubs")){
        myhubs = hub
      } else{
        myhubs = rbind(myhubs,hub)
      }
    }

    #connect the hubs for each cluster
    if (connect.hubs && exists("myhubs") && length(myhubs)>2){
      cat("Connecting hubs\n")
      write.csv(myhubs,paste0(my.dir,"/hubs.csv"))
      write.csv(mytargetdata,paste0(my.dir,"/mytargetdata.csv"))
      write.csv(myregdata,paste0(my.dir,"/myregdata.csv"))

      hubtargetdata = mytargetdata[row.names(mytargetdata)%in%myhubs,]
      hubregdata = myregdata[row.names(myregdata)%in%myhubs,]

      #Note: this pulls independent of tissue any features that match up to
      #gene symbols related to prot-ph
      if (dim(hubtargetdata)[1]==0 & prefix=='prot-ph'){
        #strip the PTM information so that we can get the targets
        genes <- unlist(strsplit(myhubs,'..', fixed = TRUE))
        genes <- genes[seq(1,length(genes),by=3)]

        any_ome_ids = .feature_to_gene() %>%
          dplyr::filter(feature_id %in% genes) %>%
          dplyr::left_join(.feature_to_gene(),
                           by = "gene_symbol") %>%
          dplyr::filter(!assay.y %in% c("epigen-methylcap-seq", "epigen-atac-seq")) %>%
          dplyr::pull(feature_id.y)
        #----
        hubtargetdata = mytargetdata %>%
          dplyr::filter(sub("\\.\\..*", "", rownames(.)) %in% any_ome_ids)

      }
      else if (dim(hubtargetdata)[1]==0 & prefix=="metab"){
        #use the regulator matrix AS the target matrix for these
        #good for metabolome
        hubtargetdata = hubregdata
      }
      network = RS.Get.Weight.Matrix(t(hubtargetdata),t(hubregdata),normalize=normalize)
      if (!is.null(network)){
        #now make a new network where we eliminate all the low confidence edges
        #There are some network values with <0 that do get trimmed
        trimmednet = data.frame(ifelse(network<weightthreshold,NaN,network),
                                check.names = FALSE)
        #translate the trimmed network into a table we can import into cytoscape
        networktable = data.frame(Regulator=character(), Interaction=character(), Target=character(), Weight=double(),Cluster=numeric(),
                                  stringsAsFactors=FALSE)
        row=1
        for (j in 1:dim(trimmednet)[1]){
          for (k in 1:dim(trimmednet)[2]){
            #skip NaNs as these have no edge
            if (is.na(trimmednet[j,k])){
              next
            }else{
              networktable[row,] = cbind(colnames(trimmednet)[k],"regulates",rownames(trimmednet)[j],trimmednet[j,k],"")
              row = row+1
            }
          }
        }
        #save this network
        finalnetwork = rbind(finalnetwork,networktable)
      }
    }

    #save the final network and the clustering information
    if(!is.null(permute)){
      cat(paste("Saving network for permutation ",permute,"\n",sep=""))
      write.table(finalnetwork,paste (my.dir,"/",prefix, '-SCION-network-permutation-',permute,'.tsv', sep=''),row.names=FALSE,quote=FALSE,sep='\t')
    }else{
      write.table(finalnetwork,paste (my.dir, "/", prefix, '-SCION-network-full.tsv', sep=''),row.names=FALSE,quote=FALSE,sep='\t')
    }
  }

  #infer the network
  #set seed so network inference is the same each time
  set.seed(0)
  SCION_infer(mytargetdata,
              myregdata,
              clusterresults = cluster_org,
              weightthreshold,
              prefix,
              permute,
              normalize,
              verbose)
  cat("Done\n")
}


scion_data_processing <- function(randomGroupCode,
                                  reg.mat,
                                  target.mat,
                                  permute=NULL,
                                  dim="row"){
  if(is.character(reg.mat) | is.character(target.mat))
    stop("It looks like your input is a character. Use .load_scion_matrixes to specify your input more specifically")

  find_shared_participants = .find_shared_scion_matrixes(randomGroupCode, reg.mat, target.mat)
  reg.mat = find_shared_participants[["scion_matrix_one"]]
  target.mat = find_shared_participants[["scion_matrix_two"]]

  #need to change the seed for each permutation. otherwise, the shuffling is exactly the same when running in parallel
  if(!is.null(permute)){
    set.seed(permute)
  }
  #permute data matrices if desired
  #target mat
  if(!is.null(permute)){
    if(dim=="col"){
      #permute the target matrix column-wise (by sample)
      target.names <- row.names(target.mat)
      target.mat <- apply(target.mat,2,sample)
      row.names(target.mat) <- target.names
    }else{
      #permute the target matrix row-wise (by gene)
      target.mat <- t(apply(target.mat,1,sample))
    }
  }
  #reg mat
  if(!is.null(permute)){
    if(dim=="col"){
      #permute the regulator matrix col-wise (by sample)
      reg.names <- row.names(reg.mat)
      reg.mat <- apply(reg.mat,2,sample)
      row.names(reg.mat) <- reg.names
    }else{
      #permute the regulator matrix row-wise (by site)
      reg.mat <- t(apply(reg.mat,1,sample))
    }
  }
  return(list("reg"=reg.mat,
              "target"=target.mat))
}

# target.matrix and input.matrix (TFs) have samples as rows and genes as columns
RS.Get.Weight.Matrix<- function(target.matrix,
                                input.matrix,
                                K="sqrt",
                                nb.trees=10000,
                                importance.measure="%IncMSE",
                                seed=NULL,
                                trace=TRUE,
                                normalize=TRUE,
                                num.cores=1, ...){
  #require(parallel)
  # set random number generator seed if seed is given
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # to be nice, report when parameter importance.measure is not correctly spelled
  if (importance.measure != "IncNodePurity" && importance.measure != "%IncMSE") {
    stop("Parameter importance.measure must be \"IncNodePurity\" or \"%IncMSE\"")
  }

  # normalize expression matrix
  target.matrix <- apply(target.matrix, 2, function(x) { (x - mean(x,na.rm = T)) / sd(x,na.rm = T) } )
  input.matrix <- apply(input.matrix, 2, function(x) { (x - mean(x,na.rm =T)) / sd(x,na.rm=T) } )
  input.matrix <- input.matrix[,!is.na(colSums(input.matrix))]
  # setup weight matrix
  num.samples <- dim(target.matrix)[1]
  num.targets <- dim(target.matrix)[2]
  num.inputs <- dim(input.matrix)[2]
  target.names <- colnames(target.matrix)
  input.names <- colnames(input.matrix)
  #print(input.names)
  #if no inputs or targets, return NULL
  if (is.null(num.inputs) | is.null(num.targets)){
    return(NULL)
  }

  weight.matrix <- matrix(0.0, nrow=num.targets, ncol=num.inputs)
  rownames(weight.matrix) <- target.names
  colnames(weight.matrix) <- input.names

  # set mtry
  if (is.numeric(K)) {
    mtry <- K
  } else if (K == "sqrt") {
    mtry <- round(sqrt(num.inputs))
  } else if (K == "all") {
    mtry <- num.inputs-1
  } else {
    stop("Parameter K must be \"sqrt\", or \"all\", or an integer")
  }
  # compute importances for every target gene
  names(target.names)<-target.names

  #parallelize if at least 3 cores, otherwise, don't
  if(num.cores>2){
    clst <- parallel::makeCluster(num.cores-1,type="FORK",outfile="log.txt")
    doParallel::registerDoParallel(clst)
    imList<-parallel::parLapply(cl=clst, X=target.names, function(x) RSGWM2(x,num.targets,target.names,input.matrix,target.matrix,trace,mtry,nb.trees,importance.measure,...))
    parallel::stopCluster(cl=clst)
  }else{
    imList<-lapply(target.names,function(x) RSGWM2(x,num.targets,target.names,input.matrix,target.matrix,trace,mtry,nb.trees,importance.measure,...))
  }

  for(nm in names(imList))
  {
    tcols<-names(imList[[nm]])
    weight.matrix[nm,tcols] <- imList[[nm]]
  }

  mynet <- weight.matrix/num.samples
  if(normalize==TRUE){
    mynet <- (mynet-min(mynet,na.rm=TRUE))/(max(mynet,na.rm=TRUE)-min(mynet,na.rm=TRUE))
  }
  return(mynet)
}

RSGWM2<-function(target.gene.name,
                 num.targets,
                 target.names,
                 input.matrix,
                 target.matrix,
                 trace,
                 mtry,
                 nb.trees,
                 importance.measure,
                 ...){
  target.gene.idx<-which(target.names==target.gene.name)
  if (trace)
  {
    flush.console()
  }
  #do NOT remove the target gene from the input matrix. Upstream did, to drop
  #autoregulation from the network, and it breaks a network with only 1 regulator.
  temp.input.matrix<-input.matrix

  x <- temp.input.matrix
  y <- target.matrix[,target.gene.name]

  rf <- randomForest::randomForest(x = x, y = y, mtry=mtry, ntree=nb.trees, keep.forest=F, importance=TRUE,...)

  im <- randomForest::importance(rf)[,importance.measure]

  return(im)
}
