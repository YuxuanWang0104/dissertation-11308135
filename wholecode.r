# 定义参数
N <-  68350000 # 总体规模
confidence_level <- 0.95  # 置信水平
E <- 0.01  # 允许误差
p <- 0.5  # 假设的总体比例，用0.5，这样会给出最大的样本量，提供最保守的估计。

# 计算z值
Z <- qnorm(1 - (1 - confidence_level) / 2)

# 计算样本量
n <- (N * Z^2 * p * (1 - p)) / ((N - 1) * E^2 + Z^2 * p * (1 - p))

# 输出结果
n

# 生成初始数据集函数，添加阈值调整
generate_data_with_beta0 <- function(n, beta_0, age_mean, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age, beta_sex, beta_cpd, beta_pcs, beta_cps) {
  #特征变量的分布
  #age <- rnorm(n, mean = age_mean, sd = age_sd)
  #age <- pmin(pmax(age, 18), 95)
    
  age_log_mean <- log(age_mean^2 / sqrt(age_sd^2 + age_mean^2))
  age_log_sd <- sqrt(log(1 + (age_sd^2 / age_mean^2)))
  
  # 生成年龄数据
  age <- rlnorm(n, meanlog = age_log_mean, sdlog = age_log_sd)
  age <- pmin(pmax(age, 18), 95)
  
  # 调整年龄分布，使75岁及以上患者占10%
  num_75_or_over <- round(n * 0.10)
  indices_75_or_over <- sample(1:n, num_75_or_over)
  age[indices_75_or_over] <- runif(num_75_or_over, min = 75, max = 94)
  
  sex <- rbinom(n, size = 1, prob = sex_prob)
  
  cpd <- rbinom(n, size = 1, prob = cpd_prob)
  
  pcs <- rbinom(n, size = 1, prob = pcs_prob)
    
  cps <- rbinom(n, size = 1, prob = cps_prob)
  
  # 计算线性预测变量
  lp <- beta_0 + beta_age * age + beta_sex * sex + beta_cpd * cpd + beta_pcs * pcs + beta_cps * cps
  # 计算概率，使用logistic函数
  p <- 1 / (1 + exp(-lp))
  
  mortality <- rbinom(n, size = 1, prob = p)
  
  data <- data.frame(
    age = age,
    sex = sex,
    chronic_pulmonary_disease = cpd,
    previous_cardiac_surgery = pcs,
    critical_preoperative_state = cps,
    mortality = mortality
  )
  
  return(data)
}

# 估算 beta_0 的函数
estimate_beta0 <- function(target_mortality_rate, n, age_mean, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age, beta_sex, beta_cpd, beta_pcs, beta_cps, tolerance = 0.001) {
  beta_0 <- 0  # 初始假设的 beta_0
  step <- 1  # 步长
  difference <- tolerance + 1  # 初始化差异值
  
  while (abs(difference) > tolerance) {
    data <- generate_data_with_beta0(n, beta_0, age_mean, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age, beta_sex, beta_cpd, beta_pcs, beta_cps)
    observed_rate <- mean(data$mortality)
    difference <- target_mortality_rate - observed_rate
    
    # 调整 beta_0
    beta_0 <- beta_0 + step * difference
  }
  
  return(beta_0)
}

# 使用估算的 beta_0 生成最终数据集
generate_final_data <- function(n, target_mortality_rate, age_mean, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age, beta_sex, beta_cpd, beta_pcs, beta_cps) {
  beta_0 <- estimate_beta0(target_mortality_rate, n, age_mean, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age, beta_sex, beta_cpd, beta_pcs, beta_cps)
  data <- generate_data_with_beta0(n, beta_0, age_mean, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age, beta_sex, beta_cpd, beta_pcs, beta_cps)
  return(data)
}

set.seed(13)
# 验证数据集
data_val <- generate_final_data(n = 3000, target_mortality_rate = 0.2, age_mean = 64.6, age_sd = 12.5, sex_prob = 0.309, cpd_prob = 0.107, pcs_prob = 0.282, cps_prob = 0.241, 
                            beta_age = 0.0285181, beta_sex = 0.2196434, beta_cpd = 0.1886564, beta_pcs = 1.118599, beta_cps = 1.086517)

# 查看生成的数据
head(data_val)

mean(data_val$mortality)
mean(data_val$age)
mean(data_val$sex)
mean(data_val$chronic_pulmonary_disease)
mean(data_val$critical_preoperative_state)
mean(data_val$previous_cardiac_surgery)


train_model <- function(data) {
  model <- glm(mortality ~ age + sex + chronic_pulmonary_disease + previous_cardiac_surgery + critical_preoperative_state, data = data, family = binomial)
  return(model)
}

library(pROC)
library(rms)
library(DescTools)

# 定义计算性能指标的函数
calculate_performance_metrics <- function(model) {

  # 计算AUC
  pred <- predict(model,newdata = data_val,type = "response")
  auc <- roc(data_val$mortality, pred)$auc
  
  # 计算O/E
  o_e <- sum(data_val$mortality) / sum(pred)
  
    
  # 存储预测结果和实际标签
  temp_df <- data.frame(pred = pred, mortality = data_val$mortality)
  
  # 绘制校准曲线并提取斜率
  calibration <- val.prob(temp_df$pred, temp_df$mortality, pl = FALSE, smooth = FALSE)
  c_slope <- calibration["Slope"]
  
  # 计算Brier Score
  brier_score <- mean((pred - data_val$mortality)^2)
  
 
  return(list(
    auc = auc,
    o_e = o_e,
    c_slope = c_slope,
    brier_score = brier_score
    
    #pred = pred  # 返回预测值以便绘制校准曲线
  ))
}


# 定义一个函数来生成beta Age变化的数据集
generate_shifted_data_age <- function(time_periods, n, target_mortality_rate, age_mean, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age_start, beta_age_end, beta_sex, beta_cpd, beta_pcs, beta_cps) {
  data_list <- list()
  beta_age_values <- seq(beta_age_start, beta_age_end, length.out = time_periods)
  
  for (i in 1:time_periods) {
    beta_age <- beta_age_values[i]
    data <- generate_final_data(n/time_periods, target_mortality_rate, age_mean, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age, beta_sex, beta_cpd, beta_pcs, beta_cps)
    #data$time_period <- i
    #data$beta <- beta_age
    data_list[[i]] <- data
  }
  
  final_data <- do.call(rbind, data_list)
  return(final_data)
}

# 定义一个函数来生成beta cpd变化的数据集
generate_shifted_data_cpd <- function(time_periods, n, target_mortality_rate, age_mean, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age, beta_cps, beta_sex, beta_cpd_start, beta_cpd_end, beta_pcs) {
  data_list <- list()
  beta_cpd_values <- seq(beta_cpd_start, beta_cpd_end, length.out = time_periods)
  
  for (i in 1:time_periods) {
    beta_cpd <- beta_cpd_values[i]
    data <- generate_final_data(n / time_periods, target_mortality_rate, age_mean, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age, beta_sex, beta_cpd, beta_pcs, beta_cps)
    #data$time_period <- i
    data_list[[i]] <- data
  }
  
  final_data <- do.call(rbind, data_list)
  return(final_data)
}

#定义一个函数来生成cpd和age共同变化的数据集
# 定义一个函数来生成逐渐变化的数据集
generate_shifted_data_both <- function(time_periods, n, target_mortality_rate, age_mean, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age_start, beta_age_end, beta_cps, beta_sex, beta_cpd_start, beta_cpd_end, beta_pcs) {
  data_list <- list()
  beta_age_values <- seq(beta_age_start, beta_age_end, length.out = time_periods)
  beta_cpd_values <- seq(beta_cpd_start, beta_cpd_end, length.out = time_periods)
  
  for (i in 1:time_periods) {
    beta_age <- beta_age_values[i]
    beta_cpd <- beta_cpd_values[i]
    data <- generate_final_data(n / time_periods, target_mortality_rate, age_mean, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age, beta_sex, beta_cpd, beta_pcs, beta_cps)
    #data$time_period <- i
    data_list[[i]] <- data
  }
  
  final_data <- do.call(rbind, data_list)
  return(final_data)
}


generate_event_rate_shifted_data <- function(time_periods, n, target_mortality_start, target_mortality_end, age_mean, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age, beta_cps, beta_sex, beta_cpd, beta_pcs) {
  data_list <- list()
  target_mortality_values <- seq(target_mortality_start, target_mortality_end, length.out = time_periods)
  
  for (i in 1:time_periods) {
    target_mortality_rate <- target_mortality_values[i]
    data <- generate_final_data(n / time_periods, target_mortality_rate, age_mean, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age, beta_sex, beta_cpd, beta_pcs, beta_cps)
    #data$time_period <- i
    data_list[[i]] <- data
  }
  
  final_data <- do.call(rbind, data_list)
  return(final_data)
}

#协变量偏移
generate_covariate_shifted_data_age <- function(time_periods, n, target_mortality_rate, age_mean_start, age_mean_end, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age, beta_cps, beta_sex, beta_cpd, beta_pcs) {
  data_list <- list()
  age_mean_values <- seq(age_mean_start, age_mean_end, length.out = time_periods)
  
  for (i in 1:time_periods) {
    age_mean <- age_mean_values[i]
    data <- generate_final_data(n / time_periods, target_mortality_rate, age_mean, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age, beta_sex, beta_cpd, beta_pcs, beta_cps)
    #data$time_period <- i
    data_list[[i]] <- data
  }
  
  final_data <- do.call(rbind, data_list)
  return(final_data)
}

generate_covariate_shifted_data_cpd <- function(time_periods, n, target_mortality_rate, age_mean, age_sd, sex_prob, cpd_prob_start, cpd_prob_end, pcs_prob, cps_prob, beta_age, beta_cps, beta_sex, beta_cpd, beta_pcs) {
  data_list <- list()
  cpd_prob_values <- seq(cpd_prob_start, cpd_prob_end, length.out = time_periods)
  
  for (i in 1:time_periods) {
    cpd_prob <- cpd_prob_values[i]
    data <- generate_final_data(n / time_periods, target_mortality_rate, age_mean, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age, beta_sex, beta_cpd, beta_pcs, beta_cps)
    #data$time_period <- i
    data_list[[i]] <- data
  }
  
  final_data <- do.call(rbind, data_list)
  return(final_data)
}

generate_covariate_shifted_data_both <- function(time_periods, n, target_mortality_rate, age_mean_start, age_mean_end, age_sd, sex_prob, cpd_prob_start, cpd_prob_end, pcs_prob, cps_prob, beta_age, beta_cps, beta_sex, beta_cpd, beta_pcs) {
  data_list <- list()
  age_mean_values <- seq(age_mean_start, age_mean_end, length.out = time_periods)
  cpd_prob_values <- seq(cpd_prob_start, cpd_prob_end, length.out = time_periods)
  
  for (i in 1:time_periods) {
    age_mean <- age_mean_values[i]
    cpd_prob <- cpd_prob_values[i]
    data <- generate_final_data(n / time_periods, target_mortality_rate, age_mean, age_sd, sex_prob, cpd_prob, pcs_prob, cps_prob, beta_age, beta_sex, beta_cpd, beta_pcs, beta_cps)
    data$time_period <- i
    data_list[[i]] <- data
  }
  
  final_data <- do.call(rbind, data_list)
  return(final_data)
}



#可视化


create_error_bar_plots <- function(data_frame, data_frame_name) {
  # 动态获取 experiment 列的唯一值，并按正确顺序排列
  unique_experiments <- unique(data_frame$experiment)
  
  # 确保 experiment 列作为因子并按唯一值排序
  data_frame$experiment <- factor(data_frame$experiment, levels = unique_experiments)
  
  # 创建误差条图和折线图的函数
  create_error_bar_plot <- function(metric_median, metric_25, metric_75, metric_name, metric_color) {
    # 检查数据框中是否存在指定的列
    if (!all(c(metric_median, metric_25, metric_75) %in% colnames(data_frame))) {
      stop(paste("数据框中找不到指定的列:", metric_median, metric_25, metric_75))
    }

    ggplot(data_frame, aes(x = experiment, y = !!sym(metric_median))) +
      geom_line(aes(group = 1, color = metric_name), linewidth = 1) +  # 使用 linewidth 代替 size
      geom_point(aes(color = metric_name), size = 2) +
      geom_errorbar(aes(ymin = !!sym(metric_25), ymax = !!sym(metric_75), color = metric_name), width = 0.2) +
      scale_color_manual(values = metric_color) +
      theme_minimal() +
      labs(title = paste("Error Bar Plot for", metric_name, "in", data_frame_name),
           x = "Experiment", y = metric_name) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
      #expand_limits(y = 0.2)  # 确保 y 轴从 0 开始
  }

  # 创建各个性能指标的图表
  plot_auc <- create_error_bar_plot("auc_median", "auc_25", "auc_75", "AUC", "blue")
  plot_oe <- create_error_bar_plot("o_e_median", "o_e_25", "o_e_75", "O/E", "green")
  plot_cslope <- create_error_bar_plot("c_slope_median", "c_slope_25", "c_slope_75", "C-Slope", "red")
  plot_brier_score <- create_error_bar_plot("brier_score_median", "brier_score_25", "brier_score_75", "Brier Score", "purple")

  # 打印图表
  print(plot_auc)
  print(plot_oe)
  print(plot_cslope)
  print(plot_brier_score)
}
