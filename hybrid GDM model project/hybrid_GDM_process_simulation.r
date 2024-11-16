
################################################################################
#                                                                              #
#                          GDM process simulation                              #
#                       Script by Milan Tsompanoglou                           #
#                                                                              #
################################################################################

# Based on the model of the publication: 
# Dong, Q., Zhou, X., Martínez, L. A Hybrid Group Decision Making Framework for 
# Achieving Agreed Solutions Based on Stable Opinions. Inf. Sci. 2019, 490, 227–243. 



library(CINNA)
library(igraph)
library(ggplot2)
library(fitdistrplus)
library(readxl)
library(reshape2)

################################################################################
#                                                                              #
#                            CONSTANT PARAMETERS                               #
#                                                                              #
################################################################################

# n <- 10       # Number of agents. 
# ε <- 0.5      # Confidence threshold. 
# φ <- 0.05     # Consensus threshold. 
# p <- 2        # Persistence Degree. 



# Beginning of function. 
simulate_GDM <- function(n, ε, φ, p) { 

################################################################################
#                                                                              #
#                    INITIAL VARIABLE VALUES & DEFINITIONS                     #
#                                                                              #
################################################################################
  
  t <- 1        # Time iteration of the GDM process. 
  a_i <- list() # Parameter for opinion update of agent i. 
  a_j <- list() # Parameter for opinion update of agent j. 
  
  
  
  # Defining the Initial Opinion Profile (t=0). 
  O <- list() 
  #O[[t]] <- list(0.75, 0.9, 0.5, 0.3, 0.25) # Enter specific opinion values. 
  
  # The opinion profile can alternatively have random values. 
  O[[t]] <- list() 
  for (i in 1:n) { 
    O[[t]][[i]] <- runif(1) # Random generation of opinion values ranging from 0 to 1. 
  } 
  #print(O[[t]]) 
  
  
  
  # Defining Distance Matrix for the current time step t. 
  D <- list() 
  D[[t]] <- matrix(NA, nrow = n, ncol = n) # Initial Distance Matrix (t=0). 
  
  # Fill in the Distance Matrix with absolute differences. 
  for (i in 1:n) { 
    for (j in 1:n) { 
      D[[t]][i, j] <- abs(O[[t]][[i]] - O[[t]][[j]]) 
    } 
  } 
  #print(D[[t]]) 
  
  
  
  # Defining Adjacency Matrix. 
  A <- list() 
  graph <- erdos.renyi.game(n, p = 0.7, type = "gnp") # Random matrix formation with ER graph. 
  
  # p is the probability to have a connection. 
  A[[t]] <- as.matrix(get.adjacency(graph, sparse = FALSE))   # Initial Adjacency Matrix. 
  #print(A[[t]]) 
  
  
  
  # Defining Agents' weights as their node degree. 
  W <- list() 
  W[[t]] <- colSums(A[[t]]) 
  #print(W[[t]]) 
  
  
  
  # Check if the Stable State can be reached with these opinions (all elements of D are <= φ or >= ε). 
  Stable_State_Reached <- all(D[[t]] < φ | D[[t]] > ε) 
  #print(Stable_State_Reached) 
  
  # In the case it can't be reached, a new discussion round is initiated (t = t + 1). 
  if (Stable_State_Reached == FALSE) { 
    t <- t + 1 
  } 
  print(t) 
  
################################################################################
#                                                                              #
#                       OPINION UPDATE / NETWORK EVOLUTION                     #
#                                                                              #
################################################################################
  
  # Repeat the process until Stable_State_Reached is TRUE. 
  while (!Stable_State_Reached) { 
    
    O[[t]] <- O[[t-1]] 
    
    
    
    # Find all the agent pairs with opinion difference ranging (φ,ε). 
    Negotiation_agent_pairs <- list() 
    
    # Loop through the elements of the distance matrix. 
    for (i in 1:n) { 
      for (j in 1:n) { 
        # Check if the conditions are met. 
        if (D[[t-1]][i, j] < ε & D[[t-1]][i, j] > φ) { 
          # Store the agent pair. 
          Negotiation_agent_pairs <- append(Negotiation_agent_pairs, list(c(i, j))) 
        } 
      } 
    } 
    #print(Negotiation_agent_pairs) 
    
    
    
    # Choose a random agent pair. 
    rap <- sample(Negotiation_agent_pairs, 1) # rap for random agent pair. 
    #print(rap) 
    
    # Check if both i and j agent's weights are zero to avoid division by zero. 
    if (W[[t-1]][rap[[1]][2]] == 0 && W[[t-1]][rap[[1]][1]] == 0) { 
      a_i[[t-1]] <- 1 - 1 / (p * (1 + 1)) 
      a_j[[t-1]] <- 1 - 1 / (p * (1 + 1)) 
    } else { 
      # Check if the i agent's weight W[[t]][rap[[1]][2]] is zero to avoid division by zero. 
      if (W[[t-1]][rap[[1]][2]] == 0) { 
        a_i[[t-1]] <- 1 - 1 / (p * (1 + W[[t-1]][rap[[1]][1]])) 
      } else { 
        a_i[[t-1]] <- 1 - W[[t-1]][rap[[1]][2]] / (p * (W[[t-1]][rap[[1]][2]] + W[[t-1]][rap[[1]][1]])) 
      } 
      
      # Check if the j agent's weight W[[t]][rap[[1]][1]] is zero to avoid division by zero. 
      if (W[[t-1]][rap[[1]][1]] == 0) { 
        a_j[[t-1]] <- 1 - 1 / (p * (W[[t-1]][rap[[1]][2]] + 1)) 
      } else { 
        a_j[[t-1]] <- 1 - W[[t-1]][rap[[1]][1]] / (p * (W[[t-1]][rap[[1]][2]] + W[[t-1]][rap[[1]][1]])) 
      } 
    } 
    
    
    
    # Opinion update. 
    #O[[t]][[rap[[1]][1]]] <- list() 
    # rap[[1]][1]] the first number in here signifies the first and only element since its a list and the second number is the agent. 1 is for i and 2 for j. 
    O[[t]][[rap[[1]][1]]] <- O[[t-1]][[rap[[1]][1]]] * a_i[[t-1]] + (1 - a_i[[t-1]]) * O[[t-1]][[rap[[1]][2]]] 
    #print(O[[t]][[rap[[1]][1]]]) 
    #O[[t]][[rap[[1]][2]]] <- list() 
    O[[t]][[rap[[1]][2]]] <- O[[t-1]][[rap[[1]][2]]] * a_j[[t-1]] + (1 - a_j[[t-1]]) * O[[t-1]][[rap[[1]][1]]] 
    #print(O[[t]][[rap[[1]][2]]]) 
    #print(O[[t]]) 
    #print(O) 
    
    
    
    # Defining Distance Matrix for the current t. 
    D[[t]] <- matrix(NA, nrow = n, ncol = n) # Initial Distance Matrix (t=0). 
    
    # Fill in the Distance Matrix with absolute differences. 
    for (i in 1:n) { 
      for (j in 1:n) { 
        D[[t]][i, j] <- abs(O[[t]][[i]] - O[[t]][[j]]) 
      } 
    } 
    #print(D[[t]]) 
    
    
    
    A[[t]] <- as.matrix(get.adjacency(graph, sparse = FALSE))   # Adjacency Matrix for this t. 
    # Adjacency Matrix update. 
    for (i in 1:n) { 
      for (j in 1:n) { 
        if (D[[t]][i, j] > ε) { 
          # If the distance is greater than ε, set the connection to 0. 
          A[[t]][i, j] <- 0 
        } else { 
          # If the distance is less than or equal to ε, set the connection to 1. 
          A[[t]][i, j] <- 1 
        } 
      } 
    } 
    #print(A[[t]]) 
    # Set the main diagonal elements to 0. 
    diag(A[[t]]) <- 0 
    #print(A[[t]]) 
    
    
    
    # Agents' weights update. 
    W[[t]] <- colSums(A[[t]]) 
    #print(W[[t]]) 
    
    
    
    # Check if the Stable State can be reached with these opinions (all elements of D are <= φ or >- ε). 
    Stable_State_Reached <- all(D[[t]] <= φ | D[[t]] >= ε) 
    #print(Stable_State_Reached) 
    
    # In the case it can't be reached, a new discussion round is initiated (t=t+1). 
    if (!Stable_State_Reached) { 
      t <- t + 1 
    } 
    print(t) 
    
    
    
    # In order to prevent infinite loop. 
    if (t == 200) { 
      Stable_State_Reached <- TRUE 
    } 
  } 
  
  
  
  #print(W[[t]]) 
  #print(A[[t]]) 
  #print(D[[t]]) 
  #print(O) 
  #print(t) # This time is the programming time iteration. 
  T_final <- t - 1 
  
  print(paste("The true GDM time iteration for agreement is T =", T_final)) 
  
################################################################################
#                                                                              #
#                  OPINION AGGREGATION FOR SELECTION PROCESS                   #
#                                                                              #
################################################################################
  
  # Calculate the sum of weights. 
  total_weight <- sum(W[[t]]) 
  #print(W[t]) 
  #print(total_weight) 
  
  
  
  # Convert the lists to matrices. 
  O_matrix <- matrix(unlist(O[[t]]), nrow = 1) 
  #print(O_matrix) 
  W_matrix <- matrix(unlist(W[[t]]), nrow = 1) 
  #print(W_matrix) 
  
  
  
  # Calculate the aggregated opinion of the GDM. 
  Agg_O <- O_matrix %*% t(W_matrix) / total_weight 
  
  
  
  print(paste("The final solution of the GDM problem is the aggregated opinion value Agg_O =", Agg_O)) 
  print(paste("It was calculated after", T_final, "time iterations."))
  
################################################################################
#                                                                              #
#                                VISUALIZATIONS                                #
#                                                                              #
################################################################################
  
  # Convert the O list to a data frame for ggplot. 
  O_df <- data.frame(t = 1:T_final, do.call(rbind, lapply(O[1:T_final], unlist))) 
  
  # Reshape the data frame for easier plotting. 
  O_melted <- melt(O_df, id.vars = "t", variable.name = "Agent", value.name = "Opinion") 
  
  # Plot the step plot. 
  plot <- ggplot(O_melted, aes(x = t, y = Opinion, group = Agent, color = Agent)) +
    geom_step(size = 1) +
    labs(title = "Opinions of Agents Over Time", x = "Time (t)", y = "Opinion") +
    theme_minimal() 
  
  print(plot) 
  
  
  # Values returned by the function. 
  return(list(O = O, D = D, A = A, W = W, T_final = T_final, Agg_O = Agg_O, plot = plot)) 
} 

################################################################################
#                                                                              #
#                                RUN THE SCRIPT                                #
#                                                                              #
################################################################################

# Run the script once with the following parameters. 
simulation_results <- simulate_GDM(n = 10, ε = 0.5, φ = 0.05, p = 2) 



# Access the results from the returned list. 
O <- simulation_results$O 
D <- simulation_results$D 
A <- simulation_results$A 
W <- simulation_results$W 
T_final <- simulation_results$T_final 
Agg_O <- simulation_results$Agg_O 
plot <- simulation_results$plot 



print(paste("The final solution of the GDM problem is the aggregated opinion value Agg_O =", Agg_O)) 
print(paste("It was calculated after", T_final, "time iterations.")) 

################################################################################
#                                                                              #
#                                  STATISTICS                                  #
#                                                                              #
################################################################################

# Run the function with different φ values. 
n <- 10 
results_list <- list() 
time_iterations <- numeric(n)  # Initialize as a numeric vector. 

for (i in 1:n) { 
  φ_value <- 0.005 * i 
  results <- simulate_GDM(n, ε = 0.5, φ = φ_value, p = 2) 
  time_iterations[i] <- results$T_final 
  results_list[[paste0("φ_", φ_value)]] <- results 
} 

print(time_iterations) 



# Create a scatter plot of the required time iterations as a function of Consensus threshold φ. 
# First create a data frame for the scatter plot. 
scatter_data <- data.frame(φ = seq(0.005, 0.05, by = 0.005), T_final = time_iterations) 

# Create the scatter plot. 
ggplot(scatter_data, aes(x = φ, y = T_final)) +
  geom_point(color = "red", size = 3) +
  labs(title = "Effect of φ on Time Iterations", 
       x = "φ", 
       y = "Time Iterations") +
  theme_minimal() 


################################################################################


# Run the function with different number of agents n. 
results_list_n <- list() 
time_iterations_n <- numeric(5)  # Initialize as a numeric vector. 

for (i in 1:10) { 
  n_value <- 2 * i 
  results <- simulate_GDM(n = n_value, ε = 0.5, φ = 0.05, p = 2) 
  time_iterations_n[i] <- results$T_final 
  results_list[[paste0("n_", n_value)]] <- results 
} 

print(time_iterations_n) 



# Create a scatter plot of the required time iterations as a function of Consensus threshold φ. 
# First create a data frame for the scatter plot. 
scatter_data <- data.frame(n = seq(2, 20, by = 2), T_final = time_iterations_n) 

# Create the scatter plot. 
ggplot(scatter_data, aes(x = n, y = T_final)) +
  geom_point(color = "red", size = 3) +
  labs(title = "Effect of n on Time Iterations",
       x = "Number of agents",
       y = "Time Iterations") +
  theme_minimal() 







