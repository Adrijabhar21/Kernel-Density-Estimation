rm(list=ls())

#required libraries
library(ggplot2)
library(patchwork)
library(tidyverse)
library(dplyr)
library(ggplot2)

#-------------------------------------------------------------------------------

#function for mixture distribution
rmix<-function(n, a_1, b_1, a_2, b_2)
{
  u=runif(n)
  z=rnorm(n,a_1, b_1)*(u<=0.5)+rnorm(n, a_2, b_2)*(u>0.5)
  return(z)
}
dmix<-function(x, a_1, b_1, a_2, b_2)
{
  y=(0.5)*(((1/(sqrt(2*pi)*b_1)) * exp(-(0.5/(b_1)^2)*((x-a_1)^2)))+
             (1/(sqrt(2*pi)*b_2)) * exp(-(0.5/(b_2)^2)*((x-a_2)^2)))
  return(y)
}
dmix_st <-function(x)
{
  y=(0.5)*(((1/(sqrt(2*pi)*1)) * exp(-(0.5/(1)^2)*((x+2)^2)))+
             (1/(sqrt(2*pi)*1)) * exp(-(0.5/(1)^2)*((x-2)^2)))
  return(y)
}

#-----------------------------------------------------------------------------
#function for kernels
naive<-function(x)
{
  y=0.5*(abs(x)<=1)
  return(y)
}
epanechnikov<-function(x)
{
  y=0.75*(1-x^2)*(abs(x)<=1)
  return(y)
}
tricube<-function(x)
{
  y=(35/32)*((1-abs(x)^3)^3)*(abs(x)<=1)
  return(y)
}
gaussian<-function(x)
{
  y=(1/(sqrt(2*pi)))*exp(-0.5*(x^2))
  return(y)
}
logistic<-function(x)
{
  y=1/(exp(x)+2+exp(-x))
  return(y)
}
cosine<-function(x)
{
  y=(pi/4)*cos((pi/2)*x)*(abs(x)<=1)
  return(y)
}
exponential<-function(x)
{
  y=exp(-x)*(x>0)
  return(y)
}
cauchy<-function(x)
{
  y=1/(pi*((1+x)^2))
  return(y)
}

#--------------------------------------------------------------------------------------------------
#plots for different kernel on same graph
#we have to manually change  f, g and also the number of parameters in the 
#argument of the function as per our need


onepara_est_cont <- function(f, g, h, n, a) {
  z <- f(n, a)
  x <- seq(min(min(z), -3), max(max(z), 3), length.out = 5000)
  r <- matrix(0, nrow = length(x), ncol = 4)
  
  # Define the kernel functions
  kernels <- list(
    Epanechnikov = function(u) ifelse(abs(u) <= 1, 3/4 * (1 - u^2), 0),
    Naive = function(u) 1/2 * (abs(u) <= 1),
    Gaussian = function(u) dnorm(u),
    Exponential = function(u) ifelse(u >= 0, exp(-u), 0)
  )
  
  for (j in 1:4) {
    for (i in 1:length(x)) {
      exp_values <- (x[i] - z) / h
      r[i, j] <- mean(kernels[[j]](exp_values)) / h
    }
  }
  
  plot(x, g(x, a), type = "l", col = "black", lwd = 2,
       xlab = "X", ylab = "Density",
       main = paste("Kernel Density Estimation (n =", n, ", h =", h, ")"),
       ylim = c(0, max(g(x, a)) + 0.1))
  
  lines(x, dnorm(x, mean = a, sd = 1), col = "blue", lwd = 2)
  
  colors <- c("red", "green", "purple", "orange")
  legend_titles <- c("Epanechnikov", "Naive", "Gaussian", "Exponential")
  
  for (j in 1:4) {
    lines(x, r[, j], col = colors[j], lwd = 2)
  }
  
  abline(h = 0, lty = 3, col = "black")
  legend("topright", legend = c(legend_titles, "Actual Density"), col
         = c(colors, "blue"), lty = 1, bty = "n", lwd = 2, y.intersp = 0.7,
         xjust = 1, yjust = 1)  # Adjusted xjust and yjust values
  
}

# Call the function with appropriate parameters
onepara_est_cont(rexp, dexp, h = 0.5, n = 100, a = 1)




#----------------------------------------------------------------------------------------------
#plots of estimated density and actual density for different sample sizes for one
#parameter continuous
generate_and_arrange_plots <- function(kernel_func, kernel_name, f, g, h, n_values, a_1) {
  plots_list <- lapply(n_values, function(n) {
    z <- f(n, a_1)
    x <- seq(min(min(z),-3), max(max(z),3), by = 0.1)
    r <- sapply(x, function(xi) mean(kernel_func((xi - z) / h)) / h)
    df <- data.frame(x, estimated_value = r, true_value = g(x, a_1))
    
    ggplot(df, aes(x = x)) +
      geom_line(aes(y = estimated_value, color = "Estimated Density")) +
      geom_line(aes(y = true_value, color = "Actual Density")) +
      scale_color_manual(values = c("Actual Density" = "blue", "Estimated Density" = "red")) +
      labs(title = paste("Sample Size n =", n, "| h =", h), subtitle = kernel_name) +
      theme_minimal() +
      theme(legend.title = element_blank())
  })
  
  plot_layout <- reduce(plots_list, `+`) +
    plot_layout(guides = 'collect') +
    plot_annotation(title = paste("Kernel Density Estimation Using", kernel_name, "Kernel"),
                    theme = theme(plot.title = element_text(hjust = 0.5, color = "dark blue", face = "bold", size = 14)))
  
  return(plot_layout)
}

n_values <- c(20, 35, 50, 73, 100,150)
plot_layout <- generate_and_arrange_plots(
  kernel_func = naive,
  kernel_name = "Box",
  f = rexp,
  g = dexp,
  h = 0.45,
  n_values = n_values,
  a_1 = 1
)

plot_layout

#---------------------------------------------------------------------------------------------------
#plots of estimated density and actual density for different sample sizes for two
#parameter continuous
#we have to manually change the title
generate_and_arrange_plots <- function(kernel_func, kernel_name, f, g, h, n_values, a_1, b_1) {
  plots_list <- lapply(n_values, function(n) {
    z <- f(n, a_1, b_1)
    x <- seq(min(min(z),-3), max(max(z),3), by = 0.1)
    r <- sapply(x, function(xi) mean(kernel_func((xi - z) / h)) / h)
    df <- data.frame(x, estimated_value = r, true_value = g(x, a_1, b_1))
    
    ggplot(df, aes(x = x)) +
      geom_line(aes(y = estimated_value, color = "Estimated Density")) +
      geom_line(aes(y = true_value, color = "Actual Density")) +
      scale_color_manual(values = c("Actual Density" = "blue", "Estimated Density" = "red")) +
      labs(title = paste("Sample Size n =", n, "| h =", h), subtitle = kernel_name) +
      theme_minimal() +
      theme(legend.title = element_blank())
  })
  
  plot_layout <- reduce(plots_list, `+`) +
    plot_layout(guides = 'collect') +
    plot_annotation(title = paste("Kernel Density Estimation Using", kernel_name, "Kernel"),
                    theme = theme(plot.title = element_text(hjust = 0.5, color = "dark blue", face = "bold", size = 14)))
  
  return(plot_layout)
}

n_values <- c(20, 35, 50, 73, 100,150)
plot_layout <- generate_and_arrange_plots(
  kernel_func = naive,
  kernel_name = "Box",
  f = rnorm,
  g = dnorm,
  h = 0.7,
  n_values = n_values,
  a_1 = 0,
  b_1 = 1
)

plot_layout
#-------------------------------------------------------------------------------------------------------
#plots of estimated density and actual density for different sample sizes for one
#parameter discrete
#we have to manually change the title
generate_and_arrange_plots <- function(kernel_func, kernel_name, f, g, h, n_values, a_1) {
  plots_list <- lapply(n_values, function(n) {
    z <- f(n, a_1)
    x <- seq(min(min(z),-3), max(max(z),3), by = 0.1)
    r <- sapply(x, function(xi) mean(kernel_func((xi - z) / h)) / h)
    df <- data.frame(x, estimated_value = r, true_value = g(x, a_1))
    
    ggplot(df, aes(x = x)) +
      geom_line(aes(y = estimated_value, color = "Estimated Density")) +
      geom_col(aes(y = true_value, fill = "Actual PMF"), width = 0.03) +
      scale_fill_manual(values = c("Actual PMF" = "blue"), guide = guide_legend(override.aes = list(color = "blue"))) +
      scale_color_manual(values = c("Estimated Density" = "red")) +
      labs(title = paste("Sample Size n =", n, "| h =", h), subtitle = kernel_name) +
      theme_minimal() +
      theme(legend.title = element_blank())
  })
  
  plot_layout <- reduce(plots_list, `+`) +
    plot_layout(guides = 'collect') +
    plot_annotation(title = paste("Kernel Density Estimation Using", kernel_name, "Kernel"),
                    theme = theme(plot.title = element_text(hjust = 0.5, color = "dark blue", face = "bold", size = 14)))
  
  return(plot_layout)
}

n_values <- c(20, 35, 50, 73, 100,150)
plot_layout <- generate_and_arrange_plots(
  kernel_func = naive,
  kernel_name = "Box",
  f = rpois,
  g = dpois,
  h = 0.45,
  n_values = n_values,
  a_1 = 1
)

plot_layout
#--------------------------------------------------------------------------------------------------
#plots of estimated density and actual density for different sample sizes for two
#parameter discrete
#we have to manually change the title
generate_and_arrange_plots <- function(kernel_func, kernel_name, f, g, h, n_values, a_1, b_1) {
  plots_list <- lapply(n_values, function(n) {
    z <- f(n, a_1, b_1)
    x <- seq(min(min(z),-3), max(max(z),3), by = 0.1)
    r <- sapply(x, function(xi) mean(kernel_func((xi - z) / h)) / h)
    df <- data.frame(x, estimated_value = r, true_value = g(x, a_1, b_1))
    
    ggplot(df, aes(x = x)) +
      geom_line(aes(y = estimated_value, color = "Estimated Density")) +
      geom_col(aes(y = true_value, fill = "Actual PMF"), width = 0.03) +
      scale_fill_manual(values = c("Actual PMF" = "blue"), guide = guide_legend(override.aes = list(color = "blue"))) +
      scale_color_manual(values = c("Estimated Density" = "red")) +
      labs(title = paste("Sample Size n =", n, "| h =", h), subtitle = kernel_name) +
      theme_minimal() +
      theme(legend.title = element_blank())
  })
  
  plot_layout <- reduce(plots_list, `+`) +
    plot_layout(guides = 'collect') +
    plot_annotation(title = paste("Kernel Density Estimation Using", kernel_name, "Kernel"),
                    theme = theme(plot.title = element_text(hjust = 0.5, color = "dark blue", face = "bold", size = 14)))
  
  return(plot_layout)
}

n_values <- c(20, 35, 50, 73, 100,150)
plot_layout <- generate_and_arrange_plots(
  kernel_func = naive,
  kernel_name = "Box",
  f = rbinom,
  g = dbinom,
  h = 0.7,
  n_values = n_values,
  a_1 = 10,
  b_1 = 0.2
)

plot_layout
#-------------------------------------------------------------------------------------------------------
#plots of estimated density and actual density for different sample sizes for mixture
#of two normals
#we have to manually change the title
generate_and_arrange_plots <- function(kernel_func, kernel_name, f, g, h, n_values, a_1, b_1, a_2, b_2) {
  plots_list <- lapply(n_values, function(n) {
    z <- f(n, a_1, b_1, a_2, b_2)
    x <- seq(min(min(z),-3), max(max(z),3), by = 0.1)
    r <- sapply(x, function(xi) mean(kernel_func((xi - z) / h)) / h)
    df <- data.frame(x, estimated_value = r, true_value = g(x, a_1, b_1, a_2, b_2))
    
    ggplot(df, aes(x = x)) +
      geom_line(aes(y = estimated_value, color = "Estimated Density")) +
      geom_line(aes(y = true_value, color = "Actual Density")) +
      scale_color_manual(values = c("Actual Density" = "blue", "Estimated Density" = "red")) +
      labs(title = paste("Sample Size n =", n, "| h =", h), subtitle = kernel_name) +
      theme_minimal() +
      theme(legend.title = element_blank())
  })
  
  plot_layout <- reduce(plots_list, `+`) +
    plot_layout(guides = 'collect') +
    plot_annotation(title = paste("Kernel Density Estimation Using", kernel_name, "Kernel"),
                    theme = theme(plot.title = element_text(hjust = 0.5, color = "dark blue", face = "bold", size = 14)))
  
  return(plot_layout)
}

n_values <- c(20, 35, 50, 73, 100,150)
plot_layout <- generate_and_arrange_plots(
  kernel_func = naive,
  kernel_name = "Box",
  f = rmix,
  g = dmix,
  h = 0.7,
  n_values = n_values,
  a_1 = -2,
  b_1 = 1,
  a_2 = 2,
  b_2 = 1
)

plot_layout
#---------------------------------------------------------------------------------------------------
#Glivenko-cantelli lemma 
#we have to manually change the kernel, f, g and also the number of parameters in the 
#argument of the function as per our need
library(ggplot2)
library(gridExtra)

n_values <- c(10, 40, 73, 100, 200, 300)

dif <- function(n, f, g, kernel, h, a, b) {
  z <- f(n, a, b)
  x <- seq(from = min(z), to = max(z), length.out = 10000)
  est <- sapply(x, function(y) mean(kernel((z - y) / h)) / h)
  act <- g(x, a, b)
  maxdif <- round(max(abs(est - act)), 4)
  
  df <- data.frame(x = x, est = est, act = act)
  
  p <- ggplot(df, aes(x = x)) +
    geom_line(aes(y = act, color = "Actual Density"), linetype = "solid") +
    geom_line(aes(y = est, color = "Estimated Density"), linetype = "dashed") +
    scale_color_manual(values = c("Actual Density" = "red", "Estimated Density" = "blue"),
                       labels = c("Actual Density", "Estimated Density")) +
    ylim(0, 0.5) +
    labs(x = "X", y = "Density", title = paste("Max Diff:", maxdif, "n =", n)) +
    theme_minimal() +
    theme(legend.position = "top")
  
  return(p)
}

# Plotting with ggplot2
plots <- lapply(n_values, function(n) dif(n, rcauchy, dcauchy, epanechnikov, 0.5, 0, 1))

# Arrange plots in a grid
grid.arrange(grobs = plots, ncol = 3)

#------------------------------------------------------------------------------------------
#Asymptotic normality
#we have to manually change the title and the distribution from which the sample is taken
#and have to set the directory where we want the plots to be saved
# we also have to manually change the quantiles
library(ggplot2)
library(gridExtra)
setwd("C:/Users/srijani/Desktop/density estimation")

# Define the kernel functions
# Define the kernel functions
gaussian_kernel <- function(u) dnorm(u)
logistic_kernel <- function(u) dlogis(u)
naive_kernel <- function(u) ifelse(abs(u) <= 1, 0.5, 0)
tricube_kernel <- function(u) ifelse(abs(u) <= 1, 35/32 * (1 - abs(u)^3)^3, 0)
cosine_kernel <- function(u) ifelse(abs(u) <= 1, pi/4 * cos(pi/2 * u), 0)
epanechnikov_kernel <- function(u) ifelse(abs(u) <= 1, 3/4 * (1 - u^2), 0)

# Set the values for n and calculate h_n for each
n_values <- c(5, 10,20, 30, 40, 50,73,100,120, 200, 350, 500)
h_values <- n_values^(-1/5)

# Function to simulate and plot for a given kernel
simulate_and_plot <- function(kernel_func, kernel_name) {
  for (i in seq_along(n_values)) {
    n <- n_values[i]
    h <- h_values[i]
    
    # Generate random samples
    x_i <- rlogis(n,0,1)
    
    # Compute the 0.1 quantile of the sample
    x_q <- qlogis(0.1, 0, 1)
    
    # For each sample, estimate the density at x_q using the kernel function
    f_n_x <- replicate(1000, {
      sample_x_i <- rlogis(n,0,1)          # New sample for each replication
      mean(kernel_func((x_q - sample_x_i) / h)) / h
    })
    
    
    df<-D(D(expression(exp(-x)/(1+exp(-x))),"x"),"x")
    mode(df)
    x<-x_q
    #E_f_n_x <- dweibull(x_q,0.5,1)+ (0.5 * h^2 * D(D(expression(dweibull), "x_q"), "x_q"))
    E_f_n_x <- dlogis(x_q,0,1)+ (0.5 * h^2 * eval(df))
    kernel_int <- integrate(function(u) kernel_func(u)^2, lower = -Inf, upper = Inf)$value
    var_f_n_x <- (1 / (n * h)) * dlogis(x_q,0,1) * kernel_int
    
    # Standardize f_n(x)
    standardized_f_n_x <- (f_n_x-E_f_n_x)/sqrt( var_f_n_x)
    
    # Create a dataframe for plotting
    data_to_plot <- data.frame(standardized_f_n_x = standardized_f_n_x)
    
    # Plot the histogram with the density of N(0,1) overlaid
    p <- ggplot(data_to_plot, aes(x = standardized_f_n_x)) +
      geom_histogram(aes(y = ..density..), bins = 30, color = "black", fill = "skyblue") +
      stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "red", size = 1) +
      labs(title = paste("Kernel:", kernel_name, "| n =", n, "| h_n =", round(h, 3)),
           x = "Standardized f_n(x_q)", y = "Density")
    
    
    filename <- paste0("plot_", kernel_name, "sample_logistic",  "_n", n, "_hn", round(h, 3), ".png")
    ggsave(filename = filename, plot = p, width = 8, height = 6)
  }
}

# Run simulation and save plots for each kernel
simulate_and_plot(gaussian_kernel, "Gaussian")
simulate_and_plot(naive_kernel, "Box")
simulate_and_plot(tricube_kernel, "Tricube")
simulate_and_plot(epanechnikov_kernel, "Epanechnikov")