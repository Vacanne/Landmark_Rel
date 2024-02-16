datasss<-structure(c(643.6778, 643.5982, NA, 643.2662, 643.3656, NA,
                  644.3629, 642.4032, 644.2238, 641.9081, 643.918, 644.818,
                  643.0176, 645.1885, 644.8896, 643.6229, 643.3475, 643.7109,
                  645.4201, 643.849, 643.4063, 643.9002, 644.869, 644.1398,
                  643.1982, 645.2147, 644.0666, 644.0027, 644.3825, 644.1253,
                  534.8515, 535.4345, 535.4598, 537.0303, 534.2815, 536.4774,
                  535.1373, 535.3429, 532.8413, 535.7005, 534.6512, 534.9217,
                  536.3685, 532.6456, 534.5621, 536.1483, 535.1335, 535.5629,
                  534.6126, 534.8484, 535.144, 535.7275, 534.2598, 535.1761,
                  533.6064, 536.3735, 534.6211, 535.5147, 535.5664, 532.6782,
                  747.7311, 746.9847, 747.9334, 747.3719, 747.0898, 747.6721,
                  748.8754, 747.2877, 747.3629, 746.4112, 747.4311, 747.7401,
                  745.7479, 748.7068, 748.0264, 746.917, 746.6313, 746.352,
                  747.9826, 748.8281, 747.2261, 747.3404, 749.1582, 747.4016,
                  747.3727, 747.7921, 748.0741, 747.8167, 748.4421, 748.1684,
                  455.4268, 454.1801, 454.5107, 455.3041, 454.7709, 454.4288,
                  455.5102, 455.8622, 452.8014, 456.3649, 457.122, 454.4475,
                  454.5527, 454.49, 455.3574, 453.3824, 453.366, 454.1397,
                  454.7946, 454.2844, 454.955, 454.8379, 454.9831, 452.4035,
                  454.0239, 454.5158, 454.679, 453.9102, 454.7414, 454.1674,
                  642.017, 641.9899, 641.5111, 642.162, 640.7423, 641.21,
                  641.7984, 641.5366, 642.0789, 643.2813, 642.0298, 642.3489,
                  640.8125, 640.5547, 643.6272, 641.3697, 641.5671, 641.1236,
                  641.2855, 641.5975, 640.3547, 641.4822, 641.7491, 642.0568,
                  642.319, 642.1502, 641.6577, 641.6475, 641.8089, 643.1683,
                  824.7686, 824.5463, 824.9609, 825.2234, 824.6971, 824.2282,
                  825.8532, 825.7691, 826.3475, 823.6111, 823.737, 823.3856,
                  824.7205, 822.5188, 824.3462, 826.8957, 825.8131, 824.0807,
                  823.0208, 826.0652, 823.9313, 823.5996, 824.9055, 825.11,
                  824.5842, 824.1129, 824.509, 823.7596, 824.7395, 826.6186), .Dim=as.integer(c(3,2,30)))

k=3
n=30
# Create a dataframe from the "datasss" list
datasss_df <- as.data.frame(datasss)
# Assign column names to the dataframe
colnames(datasss_df) <- c("x", "y")
# Add a column for the landmark and unit
n <- nrow(datasss_df)
datasss_df$landmark <- rep(1:n, each = 1)
datasss_df$unit <- rep(1:3, each = 1, times = 1)




F_approach_missing_landmark <- function(data, k, n) {
  # Step I: Automatically detect the missing landmark by checking for NA values in the x and y columns
  missing_landmark <- data[is.na(data$x) | is.na(data$y),]
  m <- missing_landmark$landmark
  # Step II: Automatically detect the i-th and j-th landmarks by checking which landmarks are not missing in the same sample of m
  reference_landmarks <- data[data$landmark %in% c(i, j) & data$unit == missing_landmark$unit & !is.na(data$x) & !is.na(data$y),]
  i <- reference_landmarks[1,]$landmark
  j <- reference_landmarks[2,]$landmark
  # Step III: Transform the data set into Bookstein coordinates by taking the i-th and j-th landmarks as reference
  data_bookstein <- transform_to_bookstein(data, reference_landmarks)
   # Step IV: In units where the m.th landmark is not missing, calculate the distances between the i.th and j.th landmarks and the m.th landmark using Euclidean distances
  data_bookstein <- data_bookstein[!is.na(data_bookstein$x_m) & !is.na(data_bookstein$y_m),]
  data_bookstein$distance_i_m <- sqrt((data_bookstein$x_i - data_bookstein$x_m)^2 + (data_bookstein$y_i - data_bookstein$y_m)^2)
  data_bookstein$distance_j_m <- sqrt((data_bookstein$x_j - data_bookstein$x_m)^2 + (data_bookstein$y_j - data_bookstein$y_m)^2)
  # Step V: The mean and standard error of the i-m (d2) and j-m (d3) distances are calculated according to the distances between the landmark calculated for each unit.
  mean_i_m <- mean(data_bookstein$distance_i_m)
  sd_i_m <- sd(data_bookstein$distance_i_m)
  mean_j_m <- mean(data_bookstein$distance_j_m)
  sd_j_m <- sd(data_bookstein$distance_j_m)
  # Step VI: Calculate 95% confidence intervals for d2 and d3 distances
  conf_int_i_m <- confint(lm(distance_i_m ~ 1, data = data_bookstein), level = 0.95)
  conf_int_j_m <- confint(lm(distance_j_m ~ 1, data = data_bookstein), level = 0.95)
   # Step VII: Accept lower limit values of confidence intervals as beginning of iteration and upper limit as end of iteration
  begin_iteration_i_m <- conf_int_i_m[1,1]
  end_iteration_i_m <- conf_int_i_m[2,1]
  begin_iteration_j_m <- conf_int_j_m[1,1]
  end_iteration_j_m <- conf_int_j_m[2,1]
  # Step VIII: Determine iteration coefficient
  iteration_coefficient <- estimate_iteration_coefficient(data_bookstein)
  # Step IX: Estimate the coordinates of the missing m-th landmark using the circle equation
  x_i <- data[data$landmark == i,]$x
  y_i <- data[data$landmark == i,]$y
  x_j <- data[data$landmark == j,]$x
  y_j <- data[data$landmark == j,]$y
  # Use the estimated iteration coefficient k, x and y coordinates of the i-th and j-th landmarks to solve for the x and y coordinates of the missing m-th landmark
  x_m <- (iteration_coefficient*x_i + x_j)/(1+iteration_coefficient)
  y_m <- (iteration_coefficient*y_i + y_j)/(1+iteration_coefficient)
  
  # Step X: Calculate F statistics for the predicted m-th landmark coordinates
  F_x_m <- (x_m - mean_i_m)/sd_i_m
  F_y_m <- (y_m - mean_j_m)/sd_j_m
  
  # Step XI: Repeat the iteration until the d2 and d3 distances reach the upper limit value
  # Step XII: Calculate Min(F) and Max(F) statistics considering all iterations
  F_x_m_all <- c()
  F_y_m_all <- c()
  while (F_x_m < end_iteration_i_m && F_y_m < end_iteration_j_m) {
    iteration_coefficient <- estimate_iteration_coefficient(data_bookstein)
    x_m <- (iteration_coefficient*x_i + x_j)/(1+iteration_coefficient)
    y_m <- (iteration_coefficient*y_i + y_j)/(1+iteration_coefficient)
    F_x_m <- (x_m - mean_i_m)/sd_i_m
    F_y_m <- (y_m - mean_j_m)/sd_j_m
    F_x_m_all <- c(F_x_m_all, F_x_m)
    F_y_m_all <- c(F_y_m_all, F_y_m)
  }
  min_F <- min(c(F_x_m_all, F_y_m_all))
  max_F <- max(c(F_x_m_all, F_y_m_all))
  
  # Step XIII: According to Min(F) and Max(F) statistics, the corresponding x and y coordinates are considered as missing landmark coordinates
  estimated_missing_landmark <- c(x_m, y_m)
  if (F_x_m == min_F) {
    estimated_missing_landmark[1] <- x_m
  }
  if (F_y_m == min_F) {
    estimated_missing_landmark[2] <- y_m
  }
  
  # Step XIV: Transform the coordinates of all landmarks from Bookstein coordinates to their original coordinates
  original_coordinates <- transform_from_bookstein(data_bookstein, reference_landmarks)
  print(min_F)
  print(max_F)
  print(estimated_missing_landmark)
  return(list(estimated_missing_landmark = estimated_missing_landmark, min_F = min_F, max_F = max_F, original_coordinates = original_coordinates))

}

F_approach_missing_landmark(datasss_df,3,30)