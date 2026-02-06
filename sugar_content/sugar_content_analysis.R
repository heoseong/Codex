# ============================================================================
# 기후변화 시나리오에 따른 사과 당도(Sugar Content) 변화 예측 분석
# ============================================================================
# 과피색(CIRG) 분석 코드를 기반으로 작성
# scenario_predictions.csv의 기후 데이터를 재활용하여 분석 시간 단축
# ============================================================================

# 필요한 패키지 설치 및 로드
required_packages <- c("tidyverse", "lubridate", "caret", "randomForest",
                       "gbm", "e1071", "glmnet", "ggplot2", "corrplot",
                       "gridExtra", "scales", "Metrics", "car")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

sapply(required_packages, install_if_missing)
lapply(required_packages, library, character.only = TRUE)

# ============================================================================
# 1. 기본 설정
# ============================================================================

# 작업 디렉토리 설정
BASE_PATH <- "G:/"
OUTPUT_PATH <- "G:/sugar_content/"
setwd(BASE_PATH)

# 지역 좌표 정보
region_coords <- data.frame(
  Region = c("Gunwi", "Cheongsong", "Yeongju", "Jangsu",
             "Geochang", "Chungju", "Pocheon", "Hwaseong"),
  lat = c(36.24, 36.43, 36.81, 35.65, 35.67, 36.97, 38.11, 37.20),
  lon = c(128.57, 129.16, 128.62, 127.52, 127.91, 127.99, 127.20, 126.82)
)

# 기후변화 시나리오 목록
scenarios <- c("SSP1-26", "SSP2-45", "SSP3-70", "SSP5-85")

# 분석 대상 변수
TARGET_VAR <- "Sugar.content"
TARGET_VAR_KR <- "당도"
TARGET_UNIT <- "°Brix"

# ============================================================================
# 2. 데이터 로드 및 전처리
# ============================================================================

cat("========================================\n")
cat("데이터 로드 및 전처리 시작\n")
cat("========================================\n")

# 사과 품질 데이터 로드
apple_data <- read.csv(file.path(BASE_PATH, "apple_fruit_quality.csv"), fileEncoding = "UTF-8-BOM")

# 날짜 형식 변환
apple_data$Full.bloom <- as.Date(apple_data$Full.bloom)
apple_data$Ripening.date <- as.Date(apple_data$Ripening.date)

# 연도 추출
apple_data$Year <- year(apple_data$Ripening.date)

# DAFB (Days After Full Bloom) 계산
apple_data$DAFB <- as.numeric(apple_data$Ripening.date - apple_data$Full.bloom)

cat("\n사과 품질 데이터 요약:\n")
cat("- 총 관측 수:", nrow(apple_data), "\n")
cat("- 연도 범위:", min(apple_data$Year), "-", max(apple_data$Year), "\n")
cat("- 지역 수:", length(unique(apple_data$Region)), "\n")

# 당도 값 확인
cat("\n", TARGET_VAR_KR, " 통계:\n", sep = "")
print(summary(apple_data[[TARGET_VAR]]))

# 기상 데이터 로드
weather_data <- read.csv(file.path(BASE_PATH, "weather_data_all_updated.csv"), fileEncoding = "UTF-8-BOM")
weather_data$date <- as.Date(weather_data$date)

# 일교차 계산
weather_data$temp_range <- weather_data$max_temp - weather_data$min_temp

cat("\n기상 데이터 요약:\n")
cat("- 총 관측 수:", nrow(weather_data), "\n")
cat("- 날짜 범위:", as.character(min(weather_data$date)), "-", as.character(max(weather_data$date)), "\n")

# ============================================================================
# 3. 기상 데이터와 품질 데이터 결합
# ============================================================================

# 각 지역/연도별로 만개일부터 수확일까지의 기상 데이터 평균 계산
calculate_season_weather <- function(apple_row, weather_df) {
  region <- apple_row$Region
  full_bloom <- apple_row$Full.bloom
  ripening <- apple_row$Ripening.date

  # 해당 지역, 기간의 기상 데이터 필터링
  season_weather <- weather_df %>%
    filter(Region == region,
           date >= full_bloom,
           date <= ripening)

  if(nrow(season_weather) == 0) {
    return(data.frame(
      avg_temp_mean = NA,
      max_temp_mean = NA,
      min_temp_mean = NA,
      solar_radiation_mean = NA,
      temp_range_mean = NA
    ))
  }

  return(data.frame(
    avg_temp_mean = mean(season_weather$avg_temp, na.rm = TRUE),
    max_temp_mean = mean(season_weather$max_temp, na.rm = TRUE),
    min_temp_mean = mean(season_weather$min_temp, na.rm = TRUE),
    solar_radiation_mean = mean(season_weather$solar_radiation, na.rm = TRUE),
    temp_range_mean = mean(season_weather$temp_range, na.rm = TRUE)
  ))
}

# 기상 데이터 계산 및 결합
cat("\n기상 데이터 결합 중...\n")

weather_summary <- do.call(rbind, lapply(1:nrow(apple_data), function(i) {
  calculate_season_weather(apple_data[i,], weather_data)
}))

# 데이터 결합
analysis_data <- cbind(apple_data, weather_summary)

# 결측치 확인 및 제거
cat("\n결측치 확인:\n")
print(colSums(is.na(analysis_data[, c(TARGET_VAR, "avg_temp_mean", "max_temp_mean",
                                       "min_temp_mean", "solar_radiation_mean", "temp_range_mean")])))

# 결측치 제거
analysis_data_clean <- analysis_data %>%
  filter(!is.na(!!sym(TARGET_VAR)) & !is.na(avg_temp_mean) & !is.na(max_temp_mean) &
         !is.na(min_temp_mean) & !is.na(solar_radiation_mean) & !is.na(temp_range_mean))

cat("\n분석에 사용할 데이터 수:", nrow(analysis_data_clean), "\n")

# ============================================================================
# 4. 상관관계 분석
# ============================================================================

cat("\n========================================\n")
cat("상관관계 분석\n")
cat("========================================\n")

# 분석에 사용할 변수 선택
corr_vars <- c(TARGET_VAR, "avg_temp_mean", "max_temp_mean", "min_temp_mean",
               "solar_radiation_mean", "temp_range_mean", "DAFB")

corr_data <- analysis_data_clean[, corr_vars]
corr_matrix <- cor(corr_data, use = "complete.obs")

cat("\n상관계수 행렬:\n")
print(round(corr_matrix, 3))

# 상관관계 그래프 저장
png(file.path(OUTPUT_PATH, "correlation_matrix.png"), width = 800, height = 600)
corrplot(corr_matrix, method = "color", type = "upper",
         addCoef.col = "black", number.cex = 0.8,
         tl.col = "black", tl.srt = 45,
         title = paste0("기후요소와 ", TARGET_VAR_KR, " 간의 상관관계"),
         mar = c(0, 0, 2, 0))
dev.off()

cat("\n상관관계 그래프 저장: correlation_matrix.png\n")

# 산점도 그래프 생성
scatter_plots <- list()

scatter_plots[[1]] <- ggplot(analysis_data_clean, aes(x = avg_temp_mean, y = !!sym(TARGET_VAR))) +
  geom_point(aes(color = Region), alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(title = paste0("평균기온 vs ", TARGET_VAR_KR), x = "평균기온 (°C)", y = paste0(TARGET_VAR_KR, " (", TARGET_UNIT, ")")) +
  theme_minimal()

scatter_plots[[2]] <- ggplot(analysis_data_clean, aes(x = max_temp_mean, y = !!sym(TARGET_VAR))) +
  geom_point(aes(color = Region), alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(title = paste0("최고기온 vs ", TARGET_VAR_KR), x = "최고기온 (°C)", y = paste0(TARGET_VAR_KR, " (", TARGET_UNIT, ")")) +
  theme_minimal()

scatter_plots[[3]] <- ggplot(analysis_data_clean, aes(x = min_temp_mean, y = !!sym(TARGET_VAR))) +
  geom_point(aes(color = Region), alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(title = paste0("최저기온 vs ", TARGET_VAR_KR), x = "최저기온 (°C)", y = paste0(TARGET_VAR_KR, " (", TARGET_UNIT, ")")) +
  theme_minimal()

scatter_plots[[4]] <- ggplot(analysis_data_clean, aes(x = solar_radiation_mean, y = !!sym(TARGET_VAR))) +
  geom_point(aes(color = Region), alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(title = paste0("일사량 vs ", TARGET_VAR_KR), x = "일사량 (MJ/m²)", y = paste0(TARGET_VAR_KR, " (", TARGET_UNIT, ")")) +
  theme_minimal()

scatter_plots[[5]] <- ggplot(analysis_data_clean, aes(x = temp_range_mean, y = !!sym(TARGET_VAR))) +
  geom_point(aes(color = Region), alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(title = paste0("일교차 vs ", TARGET_VAR_KR), x = "일교차 (°C)", y = paste0(TARGET_VAR_KR, " (", TARGET_UNIT, ")")) +
  theme_minimal()

scatter_plots[[6]] <- ggplot(analysis_data_clean, aes(x = DAFB, y = !!sym(TARGET_VAR))) +
  geom_point(aes(color = Region), alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  labs(title = paste0("DAFB vs ", TARGET_VAR_KR), x = "DAFB (일)", y = paste0(TARGET_VAR_KR, " (", TARGET_UNIT, ")")) +
  theme_minimal()

# 산점도 저장
png(file.path(OUTPUT_PATH, "scatter_plots.png"), width = 1200, height = 800)
grid.arrange(grobs = scatter_plots, ncol = 3)
dev.off()

cat("산점도 그래프 저장: scatter_plots.png\n")

# ============================================================================
# 5. 이상치 제거 (Cook's Distance 기반)
# ============================================================================

cat("\n========================================\n")
cat("이상치 제거 (Cook's Distance)\n")
cat("========================================\n")

# 단순회귀모델로 Cook's Distance 계산
simple_formula <- as.formula(paste(TARGET_VAR, "~ avg_temp_mean"))
simple_model <- lm(simple_formula, data = analysis_data_clean)
cooks_d <- cooks.distance(simple_model)

# 이상치 기준: Cook's D > 4/n
threshold <- 4 / nrow(analysis_data_clean)
outliers <- which(cooks_d > threshold)

cat("\n이상치 기준 (4/n):", round(threshold, 4), "\n")
cat("검출된 이상치 수:", length(outliers), "\n")

if(length(outliers) > 0) {
  cat("이상치 데이터:\n")
  print(analysis_data_clean[outliers, c("Region", "Year", TARGET_VAR, "avg_temp_mean")])
}

# 이상치 제거
if(length(outliers) > 0) {
  analysis_data_final <- analysis_data_clean[-outliers, ]
} else {
  analysis_data_final <- analysis_data_clean
}
cat("\n이상치 제거 후 데이터 수:", nrow(analysis_data_final), "\n")

# ============================================================================
# 6. 다중회귀분석
# ============================================================================

cat("\n========================================\n")
cat("다중회귀분석\n")
cat("========================================\n")

# 전체 변수 다중회귀모델
mlr_formula <- as.formula(paste(TARGET_VAR, "~ avg_temp_mean + max_temp_mean + min_temp_mean +
                                 solar_radiation_mean + temp_range_mean + DAFB"))
mlr_model <- lm(mlr_formula, data = analysis_data_final)

cat("\n다중회귀분석 결과:\n")
print(summary(mlr_model))

# VIF 확인 (다중공선성)
cat("\nVIF (Variance Inflation Factor):\n")
tryCatch({
  vif_values <- car::vif(mlr_model)
  print(vif_values)
}, error = function(e) {
  cat("VIF 계산 불가 (변수 수 부족 또는 공선성 문제)\n")
})

# 단계적 변수 선택
mlr_model_step <- step(mlr_model, direction = "both", trace = 0)
cat("\n단계적 변수선택 후 모델:\n")
print(summary(mlr_model_step))

# 최종 회귀모델 선택
final_mlr_model <- mlr_model_step

# 회귀계수 저장
cat("\n최종 다중회귀모델 계수:\n")
print(coef(final_mlr_model))

# ============================================================================
# 7. 머신러닝 모델 비교
# ============================================================================

cat("\n========================================\n")
cat("머신러닝 모델 비교\n")
cat("========================================\n")

# 학습/테스트 데이터 분할
set.seed(123)
train_idx <- createDataPartition(analysis_data_final[[TARGET_VAR]], p = 0.8, list = FALSE)
train_data <- analysis_data_final[train_idx, ]
test_data <- analysis_data_final[-train_idx, ]

# 교차검증 설정
train_control <- trainControl(method = "cv", number = 5)

# ML 공식
ml_formula <- as.formula(paste(TARGET_VAR, "~ avg_temp_mean + max_temp_mean + min_temp_mean +
                                solar_radiation_mean + temp_range_mean + DAFB"))

# 1. Random Forest
cat("\n1. Random Forest 모델 학습 중...\n")
rf_model <- train(
  ml_formula,
  data = train_data,
  method = "rf",
  trControl = train_control,
  tuneLength = 5
)

# 2. Gradient Boosting Machine
cat("2. GBM 모델 학습 중...\n")
gbm_grid <- expand.grid(
  interaction.depth = c(1, 3),
  n.trees = c(50, 100),
  shrinkage = 0.1,
  n.minobsinnode = min(5, floor(nrow(train_data) / 10))
)
gbm_model <- tryCatch({
  train(
    ml_formula,
    data = train_data,
    method = "gbm",
    trControl = train_control,
    tuneGrid = gbm_grid,
    verbose = FALSE
  )
}, error = function(e) {
  cat("  GBM 에러 발생, 기본 설정으로 재시도...\n")
  train(
    ml_formula,
    data = train_data,
    method = "gbm",
    trControl = trainControl(method = "cv", number = 3),
    tuneGrid = expand.grid(
      interaction.depth = 1,
      n.trees = 50,
      shrinkage = 0.1,
      n.minobsinnode = 3
    ),
    verbose = FALSE
  )
})

# 3. Support Vector Machine
cat("3. SVM 모델 학습 중...\n")
svm_model <- train(
  ml_formula,
  data = train_data,
  method = "svmRadial",
  trControl = train_control,
  tuneLength = 5,
  preProcess = c("center", "scale")
)

# 4. Elastic Net
cat("4. Elastic Net 모델 학습 중...\n")
enet_model <- tryCatch({
  train(
    ml_formula,
    data = train_data,
    method = "glmnet",
    trControl = train_control,
    tuneGrid = expand.grid(alpha = c(0, 0.5, 1), lambda = 10^seq(-3, 1, length = 10)),
    preProcess = c("center", "scale")
  )
}, error = function(e) {
  cat("  Elastic Net 에러, 기본 선형회귀 사용\n")
  train(
    ml_formula,
    data = train_data,
    method = "lm",
    trControl = train_control
  )
})

# 테스트 데이터 예측
predictions <- data.frame(
  Actual = test_data[[TARGET_VAR]]
)

predictions$MLR <- predict(final_mlr_model, newdata = test_data)
predictions$RF <- predict(rf_model, newdata = test_data)
predictions$GBM <- predict(gbm_model, newdata = test_data)
predictions$SVM <- predict(svm_model, newdata = test_data)
predictions$ENet <- predict(enet_model, newdata = test_data)

# 평가지표 계산
evaluate_model <- function(actual, predicted) {
  rmse_val <- sqrt(mean((actual - predicted)^2))
  mae_val <- mean(abs(actual - predicted))
  r2_val <- 1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
  return(c(RMSE = rmse_val, MAE = mae_val, R2 = r2_val))
}

model_results <- data.frame(
  Model = c("MLR", "RF", "GBM", "SVM", "ENet")
)

for(model_name in model_results$Model) {
  metrics <- evaluate_model(predictions$Actual, predictions[[model_name]])
  model_results[model_results$Model == model_name, c("RMSE", "MAE", "R2")] <- metrics
}

cat("\n모델 성능 비교 (테스트 데이터):\n")
print(model_results[order(-model_results$R2), ])

# 최고 성능 모델 선택
best_model_name <- model_results$Model[which.max(model_results$R2)]
cat("\n최고 성능 모델:", best_model_name, "\n")

if(best_model_name == "RF") {
  best_ml_model <- rf_model
} else if(best_model_name == "GBM") {
  best_ml_model <- gbm_model
} else if(best_model_name == "SVM") {
  best_ml_model <- svm_model
} else if(best_model_name == "ENet") {
  best_ml_model <- enet_model
} else {
  best_ml_model <- rf_model
}

# 변수 중요도 (Random Forest)
cat("\nRandom Forest 변수 중요도:\n")
print(varImp(rf_model))

# ============================================================================
# 8. 실측값 vs 예측값 비교 그래프
# ============================================================================

cat("\n========================================\n")
cat("실측값 vs 예측값 비교\n")
cat("========================================\n")

# 전체 데이터에 대한 예측
analysis_data_final$MLR_pred <- predict(final_mlr_model, newdata = analysis_data_final)
analysis_data_final$ML_pred <- predict(best_ml_model, newdata = analysis_data_final)

# 전체 데이터 평가지표
mlr_metrics_all <- evaluate_model(analysis_data_final[[TARGET_VAR]], analysis_data_final$MLR_pred)
ml_metrics_all <- evaluate_model(analysis_data_final[[TARGET_VAR]], analysis_data_final$ML_pred)

cat("\n전체 데이터 평가지표:\n")
cat("\n다중회귀모델:\n")
cat("  R² =", round(mlr_metrics_all["R2"], 4), "\n")
cat("  RMSE =", round(mlr_metrics_all["RMSE"], 4), "\n")
cat("  MAE =", round(mlr_metrics_all["MAE"], 4), "\n")

cat("\n머신러닝모델 (", best_model_name, "):\n")
cat("  R² =", round(ml_metrics_all["R2"], 4), "\n")
cat("  RMSE =", round(ml_metrics_all["RMSE"], 4), "\n")
cat("  MAE =", round(ml_metrics_all["MAE"], 4), "\n")

# 실측값 vs 예측값 그래프
p1 <- ggplot(analysis_data_final, aes(x = !!sym(TARGET_VAR), y = MLR_pred)) +
  geom_point(aes(color = Region), alpha = 0.7, size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = paste0("다중회귀모델: 실측값 vs 예측값\nR² = ",
                      round(mlr_metrics_all["R2"], 3),
                      ", RMSE = ", round(mlr_metrics_all["RMSE"], 3)),
       x = paste0("실측값 (", TARGET_VAR_KR, ")"), y = paste0("예측값 (", TARGET_VAR_KR, ")")) +
  theme_minimal() +
  coord_equal()

p2 <- ggplot(analysis_data_final, aes(x = !!sym(TARGET_VAR), y = ML_pred)) +
  geom_point(aes(color = Region), alpha = 0.7, size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = paste0("머신러닝모델 (", best_model_name, "): 실측값 vs 예측값\nR² = ",
                      round(ml_metrics_all["R2"], 3),
                      ", RMSE = ", round(ml_metrics_all["RMSE"], 3)),
       x = paste0("실측값 (", TARGET_VAR_KR, ")"), y = paste0("예측값 (", TARGET_VAR_KR, ")")) +
  theme_minimal() +
  coord_equal()

png(file.path(OUTPUT_PATH, "observed_vs_predicted.png"), width = 1200, height = 500)
grid.arrange(p1, p2, ncol = 2)
dev.off()

cat("\n실측값 vs 예측값 그래프 저장: observed_vs_predicted.png\n")

# ============================================================================
# 9. 기후변화 시나리오 데이터 활용 (scenario_predictions.csv 재활용)
# ============================================================================

cat("\n========================================\n")
cat("기후변화 시나리오 데이터 처리 (기존 데이터 재활용)\n")
cat("========================================\n")

# scenario_predictions.csv 로드 (CIRG 분석에서 생성된 기후 데이터)
scenario_climate_data <- read.csv(file.path(BASE_PATH, "scenario_predictions.csv"))

cat("\n기존 시나리오 기후 데이터 로드 완료\n")
cat("- 레코드 수:", nrow(scenario_climate_data), "\n")
cat("- 시나리오:", unique(scenario_climate_data$Scenario), "\n")
cat("- 연도 범위:", min(scenario_climate_data$Year), "-", max(scenario_climate_data$Year), "\n")

# 기후 데이터만 추출 (CIRG 예측값 제외)
all_scenario_data <- scenario_climate_data %>%
  select(Scenario, Year, Region, avg_temp_mean, max_temp_mean, min_temp_mean,
         solar_radiation_mean, temp_range_mean, DAFB)

# 당도 예측
cat("\n당도 예측 중...\n")
all_scenario_data$Sugar_MLR <- predict(final_mlr_model, newdata = all_scenario_data)
all_scenario_data$Sugar_ML <- predict(best_ml_model, newdata = all_scenario_data)

cat("예측 완료:", nrow(all_scenario_data), "레코드\n")

# 결과 저장
write.csv(all_scenario_data, file.path(OUTPUT_PATH, "scenario_predictions_sugar.csv"), row.names = FALSE)
cat("예측 결과 저장: scenario_predictions_sugar.csv\n")

# ============================================================================
# 10. 시각화: 시나리오 비교 그래프 (전체 평균 추이)
# ============================================================================

cat("\n========================================\n")
cat("시각화: 시나리오 비교 그래프\n")
cat("========================================\n")

# 연도별 시나리오별 평균 계산
yearly_summary <- all_scenario_data %>%
  group_by(Scenario, Year) %>%
  summarize(
    Sugar_MLR_mean = mean(Sugar_MLR, na.rm = TRUE),
    Sugar_MLR_sd = sd(Sugar_MLR, na.rm = TRUE),
    Sugar_ML_mean = mean(Sugar_ML, na.rm = TRUE),
    Sugar_ML_sd = sd(Sugar_ML, na.rm = TRUE),
    .groups = "drop"
  )

# 다중회귀모델 그래프
p_mlr_scenario <- ggplot(yearly_summary, aes(x = Year)) +
  geom_ribbon(aes(ymin = Sugar_MLR_mean - Sugar_MLR_sd,
                  ymax = Sugar_MLR_mean + Sugar_MLR_sd,
                  fill = Scenario), alpha = 0.2) +
  geom_line(aes(y = Sugar_MLR_mean, color = Scenario), size = 1) +
  scale_color_manual(values = c("SSP1-26" = "blue", "SSP2-45" = "green",
                                "SSP3-70" = "orange", "SSP5-85" = "red")) +
  scale_fill_manual(values = c("SSP1-26" = "blue", "SSP2-45" = "green",
                               "SSP3-70" = "orange", "SSP5-85" = "red")) +
  labs(title = paste0("기후변화 시나리오별 ", TARGET_VAR_KR, " 변화 예측 (다중회귀모델)"),
       subtitle = "음영: ±1 표준편차",
       x = "연도", y = paste0(TARGET_VAR_KR, " (", TARGET_UNIT, ")")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# 머신러닝모델 그래프
p_ml_scenario <- ggplot(yearly_summary, aes(x = Year)) +
  geom_ribbon(aes(ymin = Sugar_ML_mean - Sugar_ML_sd,
                  ymax = Sugar_ML_mean + Sugar_ML_sd,
                  fill = Scenario), alpha = 0.2) +
  geom_line(aes(y = Sugar_ML_mean, color = Scenario), size = 1) +
  scale_color_manual(values = c("SSP1-26" = "blue", "SSP2-45" = "green",
                                "SSP3-70" = "orange", "SSP5-85" = "red")) +
  scale_fill_manual(values = c("SSP1-26" = "blue", "SSP2-45" = "green",
                               "SSP3-70" = "orange", "SSP5-85" = "red")) +
  labs(title = paste0("기후변화 시나리오별 ", TARGET_VAR_KR, " 변화 예측 (", best_model_name, ")"),
       subtitle = "음영: ±1 표준편차",
       x = "연도", y = paste0(TARGET_VAR_KR, " (", TARGET_UNIT, ")")) +
  theme_minimal() +
  theme(legend.position = "bottom")

png(file.path(OUTPUT_PATH, "scenario_comparison_mlr.png"), width = 1000, height = 600)
print(p_mlr_scenario)
dev.off()

png(file.path(OUTPUT_PATH, "scenario_comparison_ml.png"), width = 1000, height = 600)
print(p_ml_scenario)
dev.off()

cat("시나리오 비교 그래프 저장: scenario_comparison_mlr.png, scenario_comparison_ml.png\n")

# ============================================================================
# 11. 시각화: 지역별 그래프
# ============================================================================

cat("\n시각화: 지역별 그래프\n")

# 지역별 연도별 평균
regional_summary <- all_scenario_data %>%
  group_by(Scenario, Year, Region) %>%
  summarize(
    Sugar_MLR = mean(Sugar_MLR, na.rm = TRUE),
    Sugar_ML = mean(Sugar_ML, na.rm = TRUE),
    .groups = "drop"
  )

# 다중회귀모델 지역별 그래프
regional_plots_mlr <- list()
for(region in unique(regional_summary$Region)) {
  region_data <- regional_summary %>% filter(Region == region)

  p <- ggplot(region_data, aes(x = Year, y = Sugar_MLR, color = Scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("SSP1-26" = "blue", "SSP2-45" = "green",
                                  "SSP3-70" = "orange", "SSP5-85" = "red")) +
    labs(title = paste(region, "- 다중회귀모델"),
         x = "연도", y = paste0(TARGET_VAR_KR, " (", TARGET_UNIT, ")")) +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 10))

  regional_plots_mlr[[region]] <- p
}

png(file.path(OUTPUT_PATH, "regional_sugar_mlr.png"), width = 1600, height = 1200)
grid.arrange(grobs = regional_plots_mlr, ncol = 3)
dev.off()

# 머신러닝모델 지역별 그래프
regional_plots_ml <- list()
for(region in unique(regional_summary$Region)) {
  region_data <- regional_summary %>% filter(Region == region)

  p <- ggplot(region_data, aes(x = Year, y = Sugar_ML, color = Scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("SSP1-26" = "blue", "SSP2-45" = "green",
                                  "SSP3-70" = "orange", "SSP5-85" = "red")) +
    labs(title = paste(region, "-", best_model_name),
         x = "연도", y = paste0(TARGET_VAR_KR, " (", TARGET_UNIT, ")")) +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 10))

  regional_plots_ml[[region]] <- p
}

png(file.path(OUTPUT_PATH, "regional_sugar_ml.png"), width = 1600, height = 1200)
grid.arrange(grobs = regional_plots_ml, ncol = 3)
dev.off()

cat("지역별 그래프 저장: regional_sugar_mlr.png, regional_sugar_ml.png\n")

# ============================================================================
# 12. 시각화: 구간별 막대그래프
# ============================================================================

cat("\n시각화: 구간별 막대그래프\n")

# 구간 정의
all_scenario_data$Period <- cut(all_scenario_data$Year,
                                breaks = c(2020, 2040, 2060, 2080, 2100),
                                labels = c("2021-2040", "2041-2060",
                                           "2061-2080", "2081-2100"))

# 구간별 평균 계산
period_summary <- all_scenario_data %>%
  group_by(Scenario, Region, Period) %>%
  summarize(
    Sugar_MLR = mean(Sugar_MLR, na.rm = TRUE),
    Sugar_ML = mean(Sugar_ML, na.rm = TRUE),
    .groups = "drop"
  )

# 다중회귀모델 막대그래프
period_plots_mlr <- list()
for(region in unique(period_summary$Region)) {
  region_data <- period_summary %>% filter(Region == region)

  p <- ggplot(region_data, aes(x = Period, y = Sugar_MLR, fill = Scenario)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("SSP1-26" = "blue", "SSP2-45" = "green",
                                 "SSP3-70" = "orange", "SSP5-85" = "red")) +
    labs(title = paste(region, "- 다중회귀모델"),
         x = "기간", y = paste0("평균 ", TARGET_VAR_KR, " (", TARGET_UNIT, ")")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          plot.title = element_text(size = 10))

  period_plots_mlr[[region]] <- p
}

png(file.path(OUTPUT_PATH, "period_barplot_mlr.png"), width = 1600, height = 1200)
grid.arrange(grobs = period_plots_mlr, ncol = 3)
dev.off()

# 머신러닝모델 막대그래프
period_plots_ml <- list()
for(region in unique(period_summary$Region)) {
  region_data <- period_summary %>% filter(Region == region)

  p <- ggplot(region_data, aes(x = Period, y = Sugar_ML, fill = Scenario)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("SSP1-26" = "blue", "SSP2-45" = "green",
                                 "SSP3-70" = "orange", "SSP5-85" = "red")) +
    labs(title = paste(region, "-", best_model_name),
         x = "기간", y = paste0("평균 ", TARGET_VAR_KR, " (", TARGET_UNIT, ")")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          plot.title = element_text(size = 10))

  period_plots_ml[[region]] <- p
}

png(file.path(OUTPUT_PATH, "period_barplot_ml.png"), width = 1600, height = 1200)
grid.arrange(grobs = period_plots_ml, ncol = 3)
dev.off()

cat("구간별 막대그래프 저장: period_barplot_mlr.png, period_barplot_ml.png\n")

# ============================================================================
# 13. 최종 평가지표 요약
# ============================================================================

cat("\n========================================\n")
cat("최종 분석 결과 요약\n")
cat("========================================\n")

# 시나리오별 변화 요약
scenario_change_summary <- all_scenario_data %>%
  group_by(Scenario) %>%
  summarize(
    Sugar_MLR_2021_2040 = mean(Sugar_MLR[Year >= 2021 & Year <= 2040], na.rm = TRUE),
    Sugar_MLR_2081_2100 = mean(Sugar_MLR[Year >= 2081 & Year <= 2100], na.rm = TRUE),
    Sugar_ML_2021_2040 = mean(Sugar_ML[Year >= 2021 & Year <= 2040], na.rm = TRUE),
    Sugar_ML_2081_2100 = mean(Sugar_ML[Year >= 2081 & Year <= 2100], na.rm = TRUE),
    .groups = "drop"
  )

scenario_change_summary$MLR_Change <- scenario_change_summary$Sugar_MLR_2081_2100 -
  scenario_change_summary$Sugar_MLR_2021_2040
scenario_change_summary$ML_Change <- scenario_change_summary$Sugar_ML_2081_2100 -
  scenario_change_summary$Sugar_ML_2021_2040

cat("\n시나리오별", TARGET_VAR_KR, "변화 (2021-2040 vs 2081-2100):\n")
print(scenario_change_summary)

# 지역별 변화 요약
regional_change_summary <- all_scenario_data %>%
  group_by(Region, Scenario) %>%
  summarize(
    Sugar_MLR_2021_2040 = mean(Sugar_MLR[Year >= 2021 & Year <= 2040], na.rm = TRUE),
    Sugar_MLR_2081_2100 = mean(Sugar_MLR[Year >= 2081 & Year <= 2100], na.rm = TRUE),
    Sugar_ML_2021_2040 = mean(Sugar_ML[Year >= 2021 & Year <= 2040], na.rm = TRUE),
    Sugar_ML_2081_2100 = mean(Sugar_ML[Year >= 2081 & Year <= 2100], na.rm = TRUE),
    .groups = "drop"
  )

regional_change_summary$MLR_Change <- regional_change_summary$Sugar_MLR_2081_2100 -
  regional_change_summary$Sugar_MLR_2021_2040
regional_change_summary$ML_Change <- regional_change_summary$Sugar_ML_2081_2100 -
  regional_change_summary$Sugar_ML_2021_2040

cat("\n지역별-시나리오별", TARGET_VAR_KR, "변화:\n")
print(regional_change_summary)

# 결과 저장
write.csv(scenario_change_summary, file.path(OUTPUT_PATH, "scenario_change_summary.csv"), row.names = FALSE)
write.csv(regional_change_summary, file.path(OUTPUT_PATH, "regional_change_summary.csv"), row.names = FALSE)

cat("\n변화량 요약 저장: scenario_change_summary.csv, regional_change_summary.csv\n")

# ============================================================================
# 14. 모델 저장
# ============================================================================

cat("\n========================================\n")
cat("모델 저장\n")
cat("========================================\n")

# 모델 저장
saveRDS(final_mlr_model, file.path(OUTPUT_PATH, "model_mlr.rds"))
saveRDS(best_ml_model, file.path(OUTPUT_PATH, "model_ml_best.rds"))
saveRDS(rf_model, file.path(OUTPUT_PATH, "model_rf.rds"))

cat("모델 저장 완료:\n")
cat("- model_mlr.rds (다중회귀모델)\n")
cat("- model_ml_best.rds (최적 머신러닝모델)\n")
cat("- model_rf.rds (Random Forest 모델)\n")

cat("\n========================================\n")
cat(TARGET_VAR_KR, "분석 완료!\n")
cat("========================================\n")

# 생성된 파일 목록
cat("\n생성된 파일 (", OUTPUT_PATH, "):\n", sep = "")
cat("1. correlation_matrix.png - 상관관계 행렬\n")
cat("2. scatter_plots.png - 산점도 그래프\n")
cat("3. observed_vs_predicted.png - 실측값 vs 예측값 비교\n")
cat("4. scenario_comparison_mlr.png - 시나리오 비교 (다중회귀)\n")
cat("5. scenario_comparison_ml.png - 시나리오 비교 (머신러닝)\n")
cat("6. regional_sugar_mlr.png - 지역별 (다중회귀)\n")
cat("7. regional_sugar_ml.png - 지역별 (머신러닝)\n")
cat("8. period_barplot_mlr.png - 구간별 막대그래프 (다중회귀)\n")
cat("9. period_barplot_ml.png - 구간별 막대그래프 (머신러닝)\n")
cat("10. scenario_predictions_sugar.csv - 시나리오 예측 결과\n")
cat("11. scenario_change_summary.csv - 시나리오별 변화 요약\n")
cat("12. regional_change_summary.csv - 지역별 변화 요약\n")
cat("13. model_mlr.rds - 다중회귀 모델\n")
cat("14. model_ml_best.rds - 최적 머신러닝 모델\n")
cat("15. model_rf.rds - Random Forest 모델\n")
