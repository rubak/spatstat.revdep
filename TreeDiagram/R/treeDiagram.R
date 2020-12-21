# # required preinstall packages: ape, tree, ggplot2, cowplot
# # also see tree diagram paper (2020 Summer)
# packages <- c("ggplot2","ape","cowplot","tree","stringr","spatstat","stats")
# # check for missing packages
# installed_packages <- packages %in% rownames(installed.packages())
# if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages])
# }
# # Packages loading
# invisible(lapply(packages, library, character.only = TRUE))


# =====================get tree information function ==================================
## This function allows users to retrieve tree information to know how many density plots to draw and how to draw
#' @importFrom spatstat is.empty
#' @importFrom stringr str_remove
#' @keywords internal
tree_info <- function(dataset){

  # avoid empty dataset
  if(is.empty(dataset)){stop("argument is of an empty dataset")}

  # library(stringr)

  ## In the tree information, each row represents either the split node or the leaf(end point of split)
  ## However, we only want the information about where the splits are, so we firstly remove rows when var==<leaf>
  split_condition <- dataset[which(dataset[,1]!="<leaf>"),]

  ## Now we want to furthur clean this dataset. As the row index is not consecutive
  ## and is not in order, we firstly sort by row index
  ## This step is not necessary but it helps later when we plot multiple densityplot together
  split_condition <-split_condition[order(as.numeric(rownames(split_condition))),]

  ## also, in "splits" column, it is a nested dataframe of two columns, which are "cutleft" and "cutright"
  ## these two columns indicate what value is used of a cut; but now, we will just save the values instead of two columns
  split_condition$values <- str_remove(split_condition$splits[,1],"<")

  # then remove "splits" to prevent storing the ducplicated information
  split_condition$splits <- NULL
  # split_condition #debug

  ## add one more variable, "node_index"
  split_condition$node_index <- as.integer(rownames(split_condition))

  ## return tree information as a data frame
  return(split_condition)
}

# e.g.
# t1.pr <- tree::tree(predictor ~.,data = your.data)
# split_condition <- tree_info(t1.pr[[1]])


#======= function reads newick string and turns to information of each split of a tree =========
#' @importFrom ape read.tree
newickToTree <- function(string){
  # turn string to phylo object in order to read as tree
  cat(string, file = "ex.tre", sep = "\n")

  # library(ape)
  tree <- read.tree("ex.tre")

  # Debug: stop if user typed in string with wrong Newick's format
  if(length(which(tree$node.label=="")==TRUE)!=0){
    stop("Error: Incorrect Newick's format is detected. Please check if your parent node is at the end of the brackets of child nodes")
  }

  # number of tip nodes
  tipnodes_num <- length(tree$tip.label)
  # number of total nodes
  totalNnode <-  tipnodes_num+tree$Nnode

  # we contrust tree information with actual node number as we used in tree diagram
  tree.dat <- as.data.frame(tree$edge)

  # add variable and value used in each internal node
  tree.dat[which(tree$edge[,2] > tipnodes_num),3] <- tree$node.label[-1]

  # add variable and value used in each tip node (end node)
  tree.dat[which(tree$edge[,2] <= tipnodes_num),3] <- tree$tip.label

  # add information for root node
  tree.dat <- rbind(tree.dat, c(0,tipnodes_num+1,tree$node.label[1]))

  # add = or < or > to tell if the node is on both, left or right side of its parent node
  tree.dat[,4] <- unlist(stringr::str_extract_all(string=tree.dat$V3, pattern="(=|<|>)"))

  # abstract variable used in each node
  tree.dat$var <- unlist(sub("(=|<|>).*", "", tree.dat$V3))

  # abstract value used in each node
  tree.dat$values <- unlist(sub(".*(=|<|>)", "", tree.dat$V3))


  # ================================================================================================= #
  # algorithm for turning tree.dat to a list of node number that can be used for drawing tree diagram #
  # ================================================================================================= #

  # Column 1 returns node number for root node(length of tip.node plus 1) to the last internal
  # node (excluding tip nodes). If we only look at one parent node at a time from column 1, usually,
  # there are three situations:
  # a parent node may have 1) two children nodes, or 2)one left child node, or 3) one right child node.
  # In each situation, we will have, for:
  # 1): abstract rows in tree.dat when v1 (column 1)==parent node number, and then we will find two
  # children nodes for such parent node.

  # NOTE: in the first object return from read.tree function(tree$edge), the order of its second
  # column(each node) is recorded by reading each flow from left to right.
  # e.g. Here is a balance tree with seven nodes in tree diagram
  #              ----1-----
  #             |         |
  #          ---2---   ---3---
  #         |      |  |      |
  #         4      5  6      7
  # In tree.dat$V2, these nodes will be recorded as 1,2,4,5,3,6,7
  # Thus, if there are two children nodes, and let update_row be the two rows of these
  # two children nodes, the first element, tree.dat$V2[update_row][1] is always the left-side node

  # after knowing which side these two children nodes belong to, we assign generation = i, and
  # calculate left-side node number = tree diagram's parent node number*2 (see tree diagram paper).
  # For right-side node number, we also assign generation = i and calculate node number = tree diagram
  # parent node number*2+1; then move to next parent node (e.g. if root node = 4, then we move to 4+1)

  # 2)&3): abstract rows in tree.dat when v1== parent node number, and since there is only one node,
  # we then check which side does this node belongs to. This is done by checking when v2==parent node
  # number, what is the indicator symbol in column 4.
  # (e.g ">" means right-side node, and < is left-side node), then we can assign generation =i, and
  # calculate actual node number

  tree.dat$generation <- 0
  tree.dat$actualNodeNum <- 1
  # for ith parent node
  for (i in (tipnodes_num+1):totalNnode){
    # print(i) #debug

    # child node(s) is(are) at row(s):
    update_row <- which(tree.dat$V1==i)
    # print(c("update_row is(are):",update_row)) #debug

    # update generation for child node(s)
    tree.dat$generation[update_row] <- tree.dat$generation[which(tree.dat$V2==i)]+1
    # print(c("generation update:",tree.dat$generation)) #debug

    # If there are two children nodes, then do
    if (length(update_row)==2){

      # calculate actual node number as we used in tree diagram for left node
      tree.dat$actualNodeNum[update_row][1] <- tree.dat$actualNodeNum[which(tree.dat$V2==i)]*2
      # calculate actual node number as we used in tree diagram for right node
      tree.dat$actualNodeNum[update_row][2] <- tree.dat$actualNodeNum[which(tree.dat$V2==i)]*2+1

    }

    # if only one child node, then check which side of the child node is at
    # (check what does its parent node indicate: > means right; < means left)
    else{
      if ( tree.dat$V4[which(tree.dat$V2==i)] == "<"){
        # print("single node on the left") #debug
        # calculate actual node number as we used in tree diagram for left node
        tree.dat$actualNodeNum[update_row] <- tree.dat$actualNodeNum[which(tree.dat$V2==i)]*2
      }
      else{
        # print("single node on the right") #debug
        # calculate actual node number as we used in tree diagram for right node
        tree.dat$actualNodeNum[update_row] <- tree.dat$actualNodeNum[which(tree.dat$V2==i)]*2+1
      }
    }
  }

  # tree.dat
  tree.dat <- tree.dat[order(tree.dat$actualNodeNum),]


  # Debug: warning message if missing variable names or values of split nodes
  if(length(which(is.na(as.numeric(tree.dat$values))==TRUE))!=0){
    warning(paste0(c("missing value(s) at row ",paste0(which(is.na(as.numeric(tree.dat$values))),","))))
  }

  else if(length(which(tree.dat$var==""))!=0){
    warning(paste0(c("missing variable name(s) at row ",paste0(which((tree.dat$var)==""),","))))
  }

  # This step helps us make this data frame easier fitting in our functions that previously created
  tree.dat$node_index <- tree.dat$actualNodeNum

  return(tree.dat)
}

# ======== get a list of node numbers along each hierarchical path to a big list ===================
# Description: this process finds a list of decision node numbers along each hierarchy path
# from root node to the node at a split ( see tree diagram paper Algorithm 1), and all these
# lists are then collected as a nested list. e.g. If a tree has three splits nodes at node
# number 1,2,3(N_1,N_2,N_3), in the return object of get_nodes_list function, "ls_node_num" will be:
# [[1]]
# [1] 1
#
# [[1]]
# [2] 1 2
#
# [[1]]
# [3] 1 3

get_nodes_list <- function(all_nodes){

  ## make sure our input is read as interger
  node_num <- as.integer(all_nodes) #each split is given a node num
  ## also, we need to check node_num is not empty (empty means no split)
  if (length(node_num)==0){stop("there's no splitting node")}
  # print(node_num) #debug

  ls_node_num <- list()

  for(i in node_num){

    last_node_num <- i
    # print(i) #debug

    if (last_node_num == 1) {
      ls_node_num = c(list(i),ls_node_num)
      # print(paste("list of node when i = 1 is",i)) #debug
    }

    else{
      ls = last_node_num
      for(j in 1:floor(log2(last_node_num))){
        previous_node <- floor(last_node_num/2)
        if (previous_node >= 1 ){
          ls = c(ls,previous_node)
          last_node_num = previous_node
        }
        # print(ls) #debug
      }

      ls_node_num = c(ls_node_num,list(sort(unique(ls))))
      # print(ls_node_num) #debug
    }
  }

  return(ls_node_num)

}

# e.g.
# get_nodes_list <- function(a_list_of_node_numbers)

# ================ transfer node number to split conditions ================
## This function returns condition (e.g values and variables used )
# at each split to subtract subdata
get_subsets_cond <- function(list_node,tree_info_data){

  cond1 = cond2 <- NULL

  for (i in list_node){
    # print(i) #debug
    rownum <- which(tree_info_data$node_index==i)
    if((i %% 2) == 0){cond2 = cond1 }
    else{cond1=cond2}
    cond1 <- paste(cond1,as.character(tree_info_data$var[rownum]),"<",tree_info_data$values[rownum],"&")
    cond2 <- paste(cond2,as.character(tree_info_data$var[rownum]),">",tree_info_data$values[rownum],"&")

  }
  # remove the last character, "&", from the condition string
  cond1 <- substr(cond1,1,nchar(cond1)-1)
  cond2 <- substr(cond2,1,nchar(cond2)-1)

  # print(cond1) #debug
  # print(cond2) #debug
  return(list(cond1,cond2))
}

# e.g.1: find the condition used at the root node (ls_node_num[[1]])
# get_subsets_cond(ls_node_num[[1]],split_condition)

# e.g.2: find all conditions using a for loop
# t2_conditions <- list()
# for (i in 1:length(ls_node_num)){
#   t2_conditions <- c(t2_conditions,get_subsets_cond(ls_node_num[[i]]))
# }

# ===============  get estimate label function =========================
# if number of label 1/total number of samples >= threshold then est_label = 1; else est_label=2
# note: this function is only suitable for new informative plot function
# future task: make it more general(?)
get_est_label <- function(data,cat_var,lab1,est_lab1=1,est_lab2=2,threshold=0.5){

  # count number of label 1
  num_lab1 <- length(which(data[,cat_var]==lab1))
  # print(c("num_lab1 is:",num_lab1)) #debug

  # estimate label
  est_lab <- ifelse(num_lab1/dim(data)[1]>=threshold,est_lab1,est_lab2)

}

# ====================== ***NEW*** informative_plot function ==========================
# ===================(previously named as) new.double.density2 function ===============
# what changed? ---> 1) we replace the old function with this new one by
# adding a new feature for showing classification
# 2) removed x-axis and add y-axis

# # parameters:
# data:dataset
# cont_var: classification variable at a split (or normal vector in non axis-align cut)
# cat_var: true classification variable (which is the predition variable in tree method)
# filename: name for the output png
# classify.value: value used at each split (e.g. Glucose>91.5 is the split condition,
# so 91.5 is the value used for that split)

#' @import ggplot2
#' @importFrom stats density
# informative (double density) plot function for one split node
informative_plot <- function(data,cont_var,cat_var,filename=NA,classify.value=NULL){

  # # Libraries
  # library(ggplot2)

  #stop running function if categorical var is more than 2 levels
  if(length(levels(data[,cat_var])) > 2) { stop("categorical variable is not binary") }
  #stop running function if var is not factored
  if(is.null(levels(data[,cat_var]))) { stop("categorical variable is not factored") }
  # avoid emtpy dataset
  if((dim(data)[1])<1){stop("informative plot requires at least two data points to draw")}

  # relevel data set
  levels(data[,cat_var])<-c(1,2)

  # subset data according to each level
  subset1 <- data[data[,cat_var] == 1,]
  subset2 <- data[data[,cat_var] == 2,]

  # Plot:
  # if data point <2, return empty plot
  if ((nrow(subset1)<2)&(nrow(subset2)<2)){
    # print("both subset have only one data point") #debug
    max_density <- 0.5 # set maximum y scale value

    p<-ggplot()+theme_void() +
      coord_cartesian(ylim = c(-max_density, max_density))
    return(p)
  }

  # set up color scheme
  # true label = 1
  lab1 <- 1
  # estimate label
  est_lab1=1
  est_lab2=2
  # threshold for determine majority
  threshold = 0.5

  # determine majority and estimated label for left and right-side data
  left.est <- get_est_label(data[data[,cont_var]<=classify.value,],cat_var,lab1,est_lab1,est_lab2,threshold)
  right.est <- get_est_label(data[data[,cont_var]>classify.value,],cat_var,lab1,est_lab1,est_lab2,threshold)
  # print(c("left.est is:",left.est, " &right.est is:", right.est))#debug

  # set up color scheme in following order: 1)Dark blue (true positive), 2)light blue(false positive),
  # 3)light pink(false negative;bottom right), 4)dark pink (true negative;bottom left)
  color_scheme <- c("#699fb3","#9cbfcd","#cd9ca7","#b3697a")

  # determind color for left-side top density,left-side bottom density,
  # right-side top density, right-side bottom density plot
  if (left.est==1){
    ifelse(right.est==1,color_num<-c(1,3,1,3),color_num<-c(1,3,2,4))
  }
  else{
    ifelse(right.est==1,color_num<-c(2,4,1,3),color_num<-c(2,4,2,4))
  }

  color_for_plot <- color_scheme[color_num]

  #avoid R CMD check from getting NOTE
  x <- y <-list()
  # subset1 has only one data point --> draw bottom density
  if (nrow(subset1)<2) {
    # print("subset1 has only one data point") #debug

    # prepare for density plot:
    # density data
    bottom_density_dat <- with(density(subset2[,cont_var]), data.frame(x=x, y=-y))
    # find maximum y scale value
    max_density <- max(density(subset2[,cont_var])$y)
    # max y-axis value
    seg_y_end <- bottom_density_dat$y[length(which(bottom_density_dat$x <= classify.value))]

    # plot
    p<-ggplot()+
      #Bottom density plot
      geom_area(data=bottom_density_dat,aes(x=x,y=y), fill = color_for_plot[4]) +
      geom_area(data=subset(bottom_density_dat,x<=classify.value),aes(x=x,y=y), fill = color_for_plot[2])+
      #remove legend,label, axis grid, ticks value, original (x,y)-axis
      theme_void() +
      # new y-axis
      geom_segment(aes(x=classify.value,xend=classify.value,y=0,yend=seg_y_end)) +
      coord_cartesian(ylim = c(-max_density, max_density))
    return(p)
  }

  # subset2 has only one data point --> draw top density
  else if (nrow(subset2)<2){
    # print("subset2 has only one data point") #debug

    # prepare for density plot: density data
    top_density_dat <- with(density(subset1[,cont_var]), data.frame(x=x, y=y))
    # find maximum y scale value
    max_density <- max(density(subset1[,cont_var])$y)
    # max y-axis value
    seg_y_end <- top_density_dat$y[length(which(top_density_dat$x <= classify.value))]

    p<-ggplot()+
      #Top density plot
      geom_area(data=top_density_dat,aes(x=x,y=y), fill = color_for_plot[3])+#9ccdc1
      geom_area(data=subset(top_density_dat,x<=classify.value),aes(x=x,y=y), fill = color_for_plot[1])+#69b3a2
      #remove legend,label, axis grid, ticks value, original (x,y)-axis
      theme_void() +
      # new y-axis
      geom_segment(aes(x=classify.value,xend=classify.value,y=0,yend=seg_y_end)) +
      coord_cartesian(ylim = c(-max_density, max_density))
    return(p)
  }

  # otherwise draw both top and bottom densities
  else{
    # prepare for density plot: density data
    top_density_dat <- with(density(subset1[,cont_var]), data.frame(x=x, y=y))
    bottom_density_dat <- with(density(subset2[,cont_var]), data.frame(x=x, y=-y))
    # find maximum y scale value
    max_density <- max(c(max(density(subset1[,cont_var])$y),
                         max(density(subset2[,cont_var])$y)))
    # max y-axis value
    seg_y_end <- top_density_dat$y[length(which(top_density_dat$x <= classify.value))]
    seg_y_start <- bottom_density_dat$y[length(which(bottom_density_dat$x <= classify.value))]

    p<-ggplot()+
      #Top density plot
      geom_area(data=top_density_dat,aes(x=x,y=y), fill = color_for_plot[3])+#9ccdc1
      geom_area(data=subset(top_density_dat,x<=classify.value),aes(x=x,y=y), fill = color_for_plot[1])+#69b3a2
      #Bottom density plot
      geom_area(data=bottom_density_dat,aes(x=x,y=y), fill = color_for_plot[4]) +
      geom_area(data=subset(bottom_density_dat,x<=classify.value),aes(x=x,y=y), fill = color_for_plot[2])+
      #remove legend,label, axis grid, ticks value, original (x,y)-axis
      theme_void() +
      # new y-axis
      geom_segment(aes(x=classify.value,xend=classify.value,y=seg_y_start,yend=seg_y_end)) +
      coord_cartesian(ylim = c(-max_density, max_density))

    # print(p) #debug

    # save pic automatically
    # ggsave(paste0(filename,".png"),bg = "transparent", height = pic_height, width = pic_width)
    return(p)
  }
}


# =============== location of each sub-plot in the output plot's coordinate ==================
loca_plot <- function(node_list,mid_x,mid_y,width,length){

  list_len <- length(node_list)

  # if length is less than 1
  if(list_len<=0){stop("length is less than 1")}
  # if length is greater than or equal to 1
  x0=mid_x
  y0=mid_y

  if(list_len>1){
    x_list = node_list[seq(2, list_len, 2)]
    y_list = setdiff(node_list,x_list)[-1]

    # find center coordinate of the plot on x-aixs
    if(length(x_list) != 0){
      n=1
      for(i in x_list){
        if( (i%%2) == 0){ x0 = x0 - ((1/2^n)+(1/2^(n+1)))*(width/2);n=n+2
        # print(paste("even x0:",x0)) #debug
        }
        else{ x0 = x0 + ((1/2^n)+(1/2^(n+1)))*(width/2);n=n+2
        # print(paste("odd x0:",x0)) #debug
        }
      }
    }
    # find center coordinate for the plot on y-axis
    if(length(y_list) != 0){
      k=1
      for(j in y_list){
        if( (j%%2) == 0){ y0 = y0 - ((1/2^k)+(1/2^(k+1)))*(length/2);k=k+2
        # print(paste("even y0:",y0)) #debug
        }
        else{ y0 = y0 + ((1/2^k)+(1/2^(k+1)))*(length/2);k=k+2
        # print(paste("odd y0:",y0)) #debug
        }
      }
    }

    width = width/(2^list_len) ; #print(paste("width:",width)) #debug
    length = length/(2^(list_len-1)); #print(paste("length:",length)) #debug
    return(c(x0-width/2,y0-length/2,width,length,list_len,nodei=node_list[list_len]))
  }
  # if length is 1
  else{
    return(c(x0-width/4,y0-length/2,width/2,length,list_len,nodei=node_list[list_len]))
  }
}

# e.g.
# loca_plot(node_list,0.5,0.5,1,1)


# ==================== function to combine multiple informative plots ==============
# We gather most of separate steps(& functions) together into one to create our tree diagram.
#' @import ggplot2
#' @import cowplot
multi_densPlot <- function(data,conditions_string=NULL,ls_node_num=NULL,split_var,cat_var,filename,pic_height=10,pic_width=10){

  # avoid empty split node
  if(length(conditions_string)!=0){

    # collecting positions for all plots (find the bottom left coordinate to put our plot
    # as well as their width and length)
    # note1: I believe that in "cowplot", the fixed width and length for an output plot is 1 for both.
    # If we increase width and length, we will only see part of output plot, which is like zooming in a plot
    # for more detail, please check "cowplot" function in "Cowplot" package

    mid_x=0.5 #center of x-axis
    mid_y=0.5 #center of y-axis
    width=1   #length of x-axis (see note1)
    length=1  #length of y-axis (see note1)
    width0 <- width
    length0 <- length

    loca_list <- list()
    for(i in 1:length(ls_node_num)){
      node_list = ls_node_num[[i]]
      loca_list[[i]] <- loca_plot(node_list,mid_x,mid_y,width,length)
    }
    # print(loca_list) #debug

    # contain a list of picture for later combining into one plot
    pic_list <- list()

    # contain a list of information for partitioned data at each node
    lab_info <- list()

    # change default seting of coordinate in ggplot to be fixed.
    # we need this step in order to do the ratation later
    # library(cowplot)
    # library(ggplot2)

    # this step is to fixed coordinate system when we rotate picture
    cf <- coord_fixed()
    cf$default <- TRUE

    for(i in 1:length(conditions_string)){
      # print(i) #debug

      # collect partitioned data
      subset1 <- data[with(data,which(eval(parse(text=conditions_string[[i]][1])))),]
      subset2 <- data[with(data,which(eval(parse(text=conditions_string[[i]][2])))),]
      combined_subset <- data.frame(rbind(subset1,subset2))
      # print(dim(combined_subset)) #debug

      # choose the normal vector to draw density plot(in axis align cut, it is the split variable)
      cont_var <- as.character(split_var[i])
      # print(cont_var) #debug

      # find the split value to draw y-segment showing where the split occur
      split_value <- as.numeric(sub(".*<", "", conditions_string[[i]][1]))
      # print(split_value) #debug

      # size of tick labels
      nodei=loca_list[[i]][6]
      # print(nodei) #debug

      # prepare for labels (node num, classify.var, min, classify value, max) for the first three cuts
      lab_info[[i]] <- c(paste0("Node",nodei,":", cont_var,"(n=", dim(combined_subset)[1],")\n",
                                "min=",min(combined_subset[,cont_var]),"\n",
                                "classify.value=",split_value,"\n",
                                "max=",max(combined_subset[,cont_var])),split_value)
      # print(lab_info[[i]]) #debug

      # draw each informative plot
      p<-informative_plot(combined_subset,cont_var,cat_var,classify.value=split_value)
      # print(p) #debug

      # for the every 2nd layer, we need to rotate our plot to the left with 90 degree
      if(((length(stringr::str_extract_all(conditions_string[[i]][1],"&")[[1]])+1) %% 2) == 0){
        p<-p+scale_y_reverse()+cf+coord_flip()+
          theme(panel.background = element_rect(fill = "transparent",colour = NA))
      }

      # we put labels for the first three cuts
      if(loca_list[[i]][5]<3){
        p<- p+geom_label(aes(x=as.numeric(lab_info[[i]][2]),y=0,label = lab_info[[i]][1]),
                         vjust = 0.5,hjust=0.5,alpha=0.5)
      }

      # print(p) #debug

      pic_list[[i]] <- p

    }

    # initial an empty plot
    p0 <- ggplot() + theme(panel.background = element_rect(fill = "transparent", colour = NA),
                           plot.background = element_rect(fill = "transparent", colour = NA))
    p <- ggdraw() +draw_plot(p0,mid_x-mid_x,mid_y-mid_y,width0,length0)

    # recursively adding each double-density plot
    for (i in 1:length(ls_node_num)){

      p<- p+draw_plot(pic_list[[i]],
                      loca_list[[i]][1],loca_list[[i]][2],loca_list[[i]][3],loca_list[[i]][4])

    }
    # save plot
    ggsave(paste0(filename,".png"),p,height=pic_height,width=pic_width,bg = "transparent")
    # return(data.frame(cbind(pic_list,loca_list))) #debug
  }
}

# # e.g.
# split_var=split_condition$var
# cat_var="Classification"
# filename = "breast cancer tree diagram"
# multi_densPlot(cancer,conditions,ls_node_num,split_var,cat_var,pic_height=30,pic_width=30,filename=filename)



# ============= combine all steps into one function =================
# NOTE!!!! function name changed from "new_tree_diagram" to "treeDiagram"
# treedat can be either a tree string or tree information in decision tree
# maybe in the future we can have 1,2,3.. to indicate different types of tree information
# constructing by different tree functions (e.g 1=newick, 2=tree information in tree(),3=...)
# see treeDiagram.Rd for more description on parameters

#' @export
treeDiagram <- function(data,treedat,cat_var,filename,pic_height=10,pic_width=10){
  # get tree information
  if (is.data.frame(treedat)==FALSE){
    # print(is.data.frame(treedat)==FALSE) #debug
    split_condition <- newickToTree(treedat)
  }
  else{
    split_condition <- tree_info(treedat)
    # split_condition #debug
  }

  # get list of sublist which records nodes along each hierarchy path from root to each terminal nodes
  all_nodes <- split_condition$node_index
  ls_node_num <- get_nodes_list(all_nodes)
  # str(ls_node_num) #debug

  # transform these node numbers to actual split conditions in character string
  # using a for loop to get all splitting conditions
  conditions <- list()
  for (i in 1:length(ls_node_num)){
    conditions[[i]] <- get_subsets_cond(ls_node_num[[i]],split_condition)
  }
  # str(conditions) #debug

  split_var=split_condition$var
  multi_densPlot(data,conditions,ls_node_num,split_var,cat_var,filename,pic_height,pic_width)

}



#==========================================================================
# ================= diagram for tesselation random forest =================
#==========================================================================

get_node_list_trf <- function(CutInfo){
  numCut <- length(CutInfo)

  if(numCut<1){stop("number of cut is zero")}

  # one cut
  else if (numCut==1){return(list(bicode.list = 1,
                                  splitNodeNum.list =1))}
  # more than one cut
  else{
    bi_code_list <- c(1,rep(NA,numCut-1))
    splitNodeNum.list <- c(1,rep(NA,numCut-1))

    for(i in 2:numCut){

      # print(i) #debug

      left_prev <- sort(as.numeric(unlist(CutInfo[[i-1]][11]))) #data on the left-hand side of parent node
      right_prev <- sort(as.numeric(unlist(CutInfo[[i-1]][12])))#data on right side of parent node

      left_current <- as.numeric(unlist(CutInfo[[i]][11]))
      right_current <- as.numeric(unlist(CutInfo[[i]][12]))
      total <- sort(c(left_current,right_current)) #data on the child node

      # print(c(left_len,right_len,total_len)) #debug

      # if child node is on the left hand side of the parent node
      if(length(which((total==left_prev)==FALSE))==0){
        bi_code_list[i] <- paste0(bi_code_list[i-1],1)
        splitNodeNum.list[i] <- splitNodeNum.list[i-1]*2
        # print(c("A::left:",bi_code_list[i],splitNodeNum.list[i])) #debug
      }

      # if child node is on the right hand side of the parent node
      else if(length(which((total==right_prev)==FALSE))==0){
        bi_code_list[i] <- paste0(bi_code_list[i-1],0)
        splitNodeNum.list[i] <- splitNodeNum.list[i-1]*2+1
        # print(c("A::right:",bi_code_list[i],splitNodeNum.list[i])) #debug
      }

      # if child node is on neither left or right hand side of the parent node, we search its
      # previous ancient node
      else{
        j=i-2
        # print(c("A:",j)) #debug
        left_prev1 <- sort(as.numeric(unlist(CutInfo[[j]][11]))) #ata on the left-hand side of parent node
        right_prev1 <- sort(as.numeric(unlist(CutInfo[[j]][12])))#data on right side of parent node

        while((length(which((total==left_prev1)==FALSE))!=0)&
              (length(which((total==right_prev1)==FALSE))!=0)){

          j=j-1
          # print(c("B: in the while loop:",j)) #debug
          left_prev1 <- sort(as.numeric(unlist(CutInfo[[j]][11]))) #number of data on the left-hand side of parent node
          right_prev1 <- sort(as.numeric(unlist(CutInfo[[j]][12])))#num. of data on right side of parent node

          # print("B:left vs total") #debug
          # print(length(which((total==left_prev1)==FALSE))) #debug
          # print("B:right vs total") #debug
          # print(length(which((total==right_prev1)==FALSE))) #debug

          if(j==0){stop("B: Error: check while loop in get_node_num_trf function;something is wrong your provided dataset")} #prevent infinite loop
        }

        if(length(which((total==left_prev1)==FALSE))==0){
          bi_code_list[i] <- paste0(bi_code_list[j],1)
          splitNodeNum.list[i] <- splitNodeNum.list[j]*2
          # print(c("C:left:",bi_code_list[i],splitNodeNum.list[i])) #debug
        }
        else if(length(which((total==right_prev1)==FALSE))==0){
          bi_code_list[i] <- paste0(bi_code_list[j],0)
          splitNodeNum.list[i] <- splitNodeNum.list[j]*2+1
          # print(c("C:right:",bi_code_list[i],splitNodeNum.list[i])) #debug
        }

        else{
          stop("Err occured in the last else{} in get_node_num_trf function;
               we suppose not to see this error message")
        }
        }
      }

    return(list(bicode.list = bi_code_list,
                splitNodeNum.list =as.numeric(splitNodeNum.list)))
    }
}

# # e.g. binary partitioning for the breast cancer cancer using tessellation
# # random forest (one tree)
# length(Cut)
# test <- get_node_list_trf(Cut)


# Visualization for result from random tessellation random forest model

# @param Cutinfo A nested list. Return output from \code{tessellation random process} software
# @param response_var Binary categorical variable. Response variable name in the tessellation random forest method
# @param ls_node_num A nested list. A list of split node number along each hierarchical path from root-to-node.
# @param filename Character string. The name of output figure
# @param pic_height numeric output figure's height; default argument is set as 10.
# @param pic_width numeric output figure's width; default argument is set as 10.
#
# @return a tree diagram
#
# @references S. Ge, S. Wang, Y. W. Teh, L. Wang, and L. T. Elliott. \emph{Random tessellation forests}. In Proceedings of the 33rd Conference on Neural Information Processing Systems, 2019.
#
# @import ggplot2
# @import cowplot
# @import stringr
# @importFrom utils str
#
# @keywords internal
multi_densPlot_trf <- function(Cutinfo,response_var,ls_node_num=NULL,filename,pic_height=NULL,pic_width=NULL){

  if((length(Cutinfo)-1)!=0){

    # contain a list of picture for later combining into one plot
    pic_list <- list()

    # change default seting of coordinate in ggplot to be fixed.
    # we need this step in order to do the ratation later
    # library(cowplot)
    # library(ggplot2)
    cf <- coord_fixed()
    cf$default <- TRUE

    for(i in 1:(length(Cutinfo)-1)){
      subset1 <- data.frame(t.scale=unlist(Cutinfo[[i]][9]),label=response_var[as.numeric(unlist(Cutinfo[[i]][11]))])
      subset2 <- data.frame(t.scale=unlist(Cutinfo[[i]][10]),label=response_var[as.numeric(unlist(Cutinfo[[i]][12]))])
      colnames(subset1)[1] =colnames(subset2)[1]="t.scale"
      # print("before combining data") #debug
      combined_subset <- data.frame(rbind(subset1,subset2))
      combined_subset$label <- as.factor(combined_subset$label)

      # print(str(combined_subset)) #debug

      # plot picture separately
      cont_var <- "t.scale"
      cat_var <- "label"
      split_value <- as.numeric(unlist(Cutinfo[[i]][13]))
      p<-informative_plot(combined_subset,cont_var,cat_var,classify.value=split_value)
      # readline(prompt = "Pause. Press <Enter> to continue...") #debug
      # print(p) #debug

      # for the every 2nd layer, we need to rotate out plot
      if((length(ls_node_num[[i]]) %% 2) == 0){
        p<-p+scale_y_reverse()+cf +coord_flip()+theme(panel.background = element_rect(fill = "transparent",colour = NA))

      }

      pic_list[[i]] <- p
    }
    # return(pic_list)}} #debug
    # {{ #debug

    # after collecting all plots, now we need to figure out the bottom left coordinate to put our plot
    # as well as their width and length
    mid_x=0.5
    mid_y=0.5
    width=1
    length=1
    width0 <- width
    length0 <- length

    loca_list <- list()
    for(i in 1:length(ls_node_num)){
      node_list = ls_node_num[[i]]
      loca_list[[i]] <- loca_plot(node_list,mid_x,mid_y,width,length)
    }

    # print(loca_list) #debug

    # recursively adding each double-density plot
    # initial an empty plot
    p0 <- ggplot() + theme(panel.background = element_rect(fill = "transparent", colour = NA),
                           plot.background = element_rect(fill = "transparent", colour = NA))
    p <- ggdraw() +draw_plot(p0,mid_x-mid_x,mid_y-mid_y,width0,length0)
    for (i in 1:(length(ls_node_num)-1)){
      p<- p+draw_plot(pic_list[[i]],
                      loca_list[[i]][1],loca_list[[i]][2],loca_list[[i]][3],loca_list[[i]][4])
    }
    p
    ggsave(paste0(filename,".png"),p,height=pic_height,width=pic_width,bg = "transparent")

    # return(data.frame(cbind(pic_list,loca_list))) #debug
  }
}

# e.g.
# filename = "breast cancer tess19 visualization"
# multi_densPlot_trf(Cut,ls_node_num,filename,10,10)
