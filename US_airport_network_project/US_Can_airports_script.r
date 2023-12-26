library(CINNA)
library(igraph)
library(ggplot2)
library(fitdistrplus)
library(readxl)


p<- read.table("reachability_no_headers.txt", sep = "", header = F,
               stringsAsFactors = FALSE)
G<- graph_from_data_frame(p, directed=T, vertices=NULL )
G<-set.edge.attribute(G,"weight", value = abs(E(G)$V3))
print(p)

airport_data <- read_excel("reachability-meta.xlsx")
print(airport_data)

################################################################################
#                                                                              #
#                           NETWORK CHARACTERISTICS                            #
#                                                                              #
################################################################################

is.weighted(G)   # TRUE
is.directed(G)   # TRUE
is_simple(G)     # do loops or / and multiple edges TRUE
is_connected(G)  # TRUE
length(V(G))     # 456
vcount(G)        # order of the graph 456
ecount(G)        # size of the graph 71959
radius(G)        # 2
average.path.length(G) # 288.4665
diameter(G)      # 938
transitivity(G, type="global")
edge_density(G)  # 0.3468238

################################################################################
#                                                                              #
#                                 CENTRALITIES                                 #
#                                                                              #
################################################################################

pr_cent<-proper_centralities(G) # Shows all the centralities available for this graph. 
a<-pr_cent[c(8, 11, 14, 16, 27)]
calc_cent <- calculate_centralities(G, include  = a)

################# DATAFRAME FROM THE CENTRALITIES ######################
nodesg <- as_ids(V(G)) #length(nodesg)=456


cl_cent<-calculate_centralities(G, include = "Closeness Centrality (Freeman)")
V(G)$closeness<-cl_cent

eig_cent<-calculate_centralities(G, include = "eigenvector centralities")
V(G)$Eigenv<-eig_cent

Deg_Cent<-calculate_centralities(G, include = "Degree Centrality")
V(G)$Degree<-Deg_Cent

Betw_cent<-calculate_centralities(G, include = "Shortest-Paths Betweenness Centrality")
V(G)$Betweeness<-Betw_cent

Ecc_Centr<-calculate_centralities(G, include = "Eccentricity Centrality")
V(G)$Ecc<-Ecc_Centr

Ver_deg<-calculate_centralities(G, include = "Weighted Vertex Degree")
V(G)$Vertex_deg<-Ver_deg


componentg<-as.numeric(components(G)$membership)


inDegree<-centr_degree(
  G,
  mode = c("in"),
  loops = TRUE,
  normalized = TRUE)
V(G)$In_Degree<-inDegree

outDegree<-centr_degree(
  G,
  mode = c("out"),
  loops = TRUE,
  normalized = TRUE)
V(G)$Out_Degree<-outDegree

inDegree
outDegree
Deg_Cent
sorted_IN <- lapply(inDegree,sort,decreasing=TRUE)
sorted_OUT<-lapply(outDegree,sort,decreasing=TRUE)
sorted_IN
sorted_OUT
list3<-mapply("+", sorted_IN, sorted_OUT, SIMPLIFY = FALSE)
list3


centralitiesg<-data.frame(x1=nodesg, 
                          x2=componentg, 
                          x3=cl_cent, 
                          x4=eig_cent,
                          x5=Deg_Cent,
                          x6=Betw_cent,
                          x7=Ecc_Centr # έβγαλα το Katz_Cent και δούλεψε
                          )

dimnames(centralitiesg)[[2]]<-c("node", 
                                "component", 
                                "closeness", 
                                "Eigenvc", 
                                "Degree",
                                "Betweenness",
                                "Eccentricity_C")


head(centralitiesg, n=20)
str(centralitiesg)               ## variables in the data frame

gcdata.f<-centralitiesg



## Sort according to degree and get the first 20 rows
ord.deg.gcda.f <- gcdata.f[order(-gcdata.f$Degree), ][1:10,];  
ord.deg.gcda.f

ord.closeness.gcda.f <- gcdata.f[order(-gcdata.f$closeness), ][1:10,]; 
ord.closeness.gcda.f

ord.Eigenvc.gcda.f <- gcdata.f[order(-gcdata.f$Eigenvc), ][1:10,]; 
ord.Eigenvc.gcda.f

ord.Eccentricity_C.gcda.f <- gcdata.f[order(-gcdata.f$Eccentricity_C), ][1:10,]; 
ord.Eccentricity_C.gcda.f

ord.Betweenness.gcda.f <- gcdata.f[order(-gcdata.f$Betweenness), ][1:10,]; 
ord.Betweenness.gcda.f


deg <- degree(G)
hist(deg,   main="Histogram of node degree")


################################################################################
#                                                                              #
#                             HAVERSINE DISTANCE                               #
#                                                                              #
################################################################################

# Define haversine function calculating the distances in km. 
haversine <- function(lat1, lon1, lat2, lon2) {
  # Convert degrees to radians
  lat1 <- (lat1 * pi) / 180
  lon1 <- (lon1 * pi) / 180
  lat2 <- (lat2 * pi) / 180
  lon2 <- (lon2 * pi) / 180
  
  # Haversine formula
  dlon <- lon2 - lon1
  dlat <- lat2 - lat1
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  radius <- 6371  # Radius of the Earth in kilometers. Use 3956 for miles.
  distance <- radius * c
  
  return(distance)
}



# Create a data frame to store results
result_df <- data.frame(
  node1 = numeric(),
  node2 = numeric(),
  node1_name = character(),
  node2_name = character(),
  distance = numeric()
)


# Iterate over each row in 'p'
for (i in 1:nrow(p)) {
  node1_id <- p[i, 1]
  node2_id <- p[i, 2]
  
  # Get information for node1
  node1_info <- airport_data[airport_data$node_id == node1_id, c("latitude_2", "longitude_2", "node_id", "name")]
  node1_coords <- node1_info[c("latitude_2", "longitude_2")]
  
  # Get information for node2
  node2_info <- airport_data[airport_data$node_id == node2_id, c("latitude_2", "longitude_2", "node_id", "name")]
  node2_coords <- node2_info[c("latitude_2", "longitude_2")]
  
  # Calculate haversine distance
  distance <- haversine(node1_coords$latitude_2, node1_coords$longitude_2,
                        node2_coords$latitude_2, node2_coords$longitude_2)
  
  # Store the results
  result_df <- rbind(
    result_df,
    data.frame(
      node1 = node1_id,
      node2 = node2_id,
      node1_name = node1_info$name,
      node2_name = node2_info$name,
      distance = distance
    )
  )
}
# 'result_df' now contains the distances for each network connection
print(result_df)

# Sort the result_df dataframe by the distance
sorted_result_df_dec <- result_df[order(result_df$distance, decreasing = TRUE), ]
print(sorted_result_df_dec)

sorted_result_df_inc <- result_df[order(result_df$distance, decreasing = FALSE), ]
print(sorted_result_df_inc)
