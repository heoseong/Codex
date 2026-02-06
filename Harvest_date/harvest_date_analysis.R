# ============================================================================
# 기후변화 시나리오에 따른 사과 수확시기 예측 (RSM 기반)
# ============================================================================
# 당도와 과피색(CIRG)을 기반으로 Response Surface Methodology 적용
# 당도는 높을수록, CIRG는 5에 가까울수록 수확 적기
# ============================================================================

# 필요한 패키지 설치 및 로드
required_packages <- c("tidyverse", "lubridate", "rsm", "ggplot2",
                       "gridExtra", "scales", "sf",
                       "rnaturalearth", "rnaturalearthdata")

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
OUTPUT_PATH <- "G:/Harvest_date/"
setwd(BASE_PATH)

# 기후변화 시나리오 목록
scenarios <- c("SSP1-26", "SSP2-45", "SSP3-70", "SSP5-85")

# ============================================================================
# 2. CIRG 계산 함수
# ============================================================================

calculate_CIRG <- function(L_star, a_star, b_star) {
  h_rad <- atan2(b_star, a_star)
  h_deg <- h_rad * 180 / pi
  h_deg <- ifelse(h_deg < 0, h_deg + 360, h_deg)
  C_star <- sqrt(a_star^2 + b_star^2)
  CIRG <- (180 - h_deg) / (L_star + C_star)
  return(CIRG)
}

# ============================================================================
# 3. 데이터 로드 및 전처리
# ============================================================================

cat("========================================\n")
cat("데이터 로드 및 전처리 시작\n")
cat("========================================\n")

apple_data <- read.csv(file.path(BASE_PATH, "apple_fruit_quality.csv"),
                       fileEncoding = "UTF-8-BOM")

apple_data$Full.bloom <- as.Date(apple_data$Full.bloom)
apple_data$Ripening.date <- as.Date(apple_data$Ripening.date)

apple_data$Year <- year(apple_data$Ripening.date)
apple_data$DAFB <- as.numeric(apple_data$Ripening.date - apple_data$Full.bloom)
apple_data$CIRG <- calculate_CIRG(apple_data$L., apple_data$a., apple_data$b.)

analysis_data <- apple_data %>%
  select(Region, Year, Full.bloom, Ripening.date, DAFB, Sugar.content, CIRG) %>%
  filter(!is.na(Sugar.content) & !is.na(CIRG) & !is.na(DAFB))

cat("\n분석에 사용할 데이터 수:", nrow(analysis_data), "\n")

# ============================================================================
# 4. RSM 모델 적합 (DAFB ~ Sugar + CIRG)
# ============================================================================

cat("\n========================================\n")
cat("RSM 모델 적합\n")
cat("========================================\n")

rsm_model <- rsm(
  DAFB ~ FO(Sugar.content, CIRG) + TWI(Sugar.content, CIRG) + PQ(Sugar.content, CIRG),
  data = analysis_data
)

cat("\nRSM 모델 요약:\n")
print(summary(rsm_model))

# ============================================================================
# 5. 수확 적기 기준(당도 높음, CIRG=5 근접) 정의
# ============================================================================

sugar_range <- range(analysis_data$Sugar.content, na.rm = TRUE)
cirg_range <- range(analysis_data$CIRG, na.rm = TRUE)
max_cirg_dev <- max(abs(cirg_range - 5))

sugar_desirability <- function(x) {
  scaled <- (x - sugar_range[1]) / diff(sugar_range)
  pmax(pmin(scaled, 1), 0)
}

cirg_desirability <- function(x) {
  scaled <- 1 - (abs(x - 5) / max_cirg_dev)
  pmax(pmin(scaled, 1), 0)
}

overall_desirability <- function(sugar, cirg) {
  sqrt(sugar_desirability(sugar) * cirg_desirability(cirg))
}

# 최적 조합 탐색
sugar_grid <- seq(sugar_range[1], sugar_range[2], length.out = 100)
cirg_grid <- seq(cirg_range[1], cirg_range[2], length.out = 100)
opt_grid <- expand.grid(Sugar.content = sugar_grid, CIRG = cirg_grid)
opt_grid$desirability <- overall_desirability(opt_grid$Sugar.content, opt_grid$CIRG)

best_row <- opt_grid[which.max(opt_grid$desirability), ]
optimal_dafb <- predict(rsm_model, newdata = best_row)

cat("\n최적 수확 조건(그리드 탐색):\n")
cat("- 당도:", round(best_row$Sugar.content, 2), "\n")
cat("- CIRG:", round(best_row$CIRG, 2), "\n")
cat("- 예측 DAFB:", round(optimal_dafb, 1), "\n")

# ============================================================================
# 6. 기후변화 시나리오 예측 데이터 결합
# ============================================================================

cat("\n========================================\n")
cat("시나리오 예측 데이터 결합\n")
cat("========================================\n")

sugar_scenario <- read.csv(file.path(BASE_PATH, "sugar_content", "scenario_predictions_sugar.csv"))
cirg_scenario <- read.csv(file.path(BASE_PATH, "Fruit_color", "scenario_predictions.csv"))

scenario_data <- cirg_scenario %>%
  select(Scenario, Year, Region, DAFB, CIRG_MLR, CIRG_ML) %>%
  inner_join(
    sugar_scenario %>%
      select(Scenario, Year, Region, DAFB, Sugar_MLR, Sugar_ML),
    by = c("Scenario", "Year", "Region", "DAFB")
  )

cat("\n결합된 시나리오 데이터:", nrow(scenario_data), "레코드\n")

# ============================================================================
# 7. 수확시기 예측 및 날짜 계산
# ============================================================================

cat("\n========================================\n")
cat("수확시기 예측\n")
cat("========================================\n")

bloom_doy_by_region <- analysis_data %>%
  group_by(Region) %>%
  summarize(mean_bloom_doy = round(mean(yday(Full.bloom), na.rm = TRUE)), .groups = "drop")

scenario_data <- scenario_data %>%
  left_join(bloom_doy_by_region, by = "Region") %>%
  mutate(
    Harvest_DAFB_MLR = predict(rsm_model, newdata = data.frame(
      Sugar.content = Sugar_MLR,
      CIRG = CIRG_MLR
    )),
    Harvest_DAFB_ML = predict(rsm_model, newdata = data.frame(
      Sugar.content = Sugar_ML,
      CIRG = CIRG_ML
    )),
    Desirability_MLR = overall_desirability(Sugar_MLR, CIRG_MLR),
    Desirability_ML = overall_desirability(Sugar_ML, CIRG_ML),
    Full_bloom_est = as.Date(paste0(Year, "-01-01")) + days(mean_bloom_doy - 1),
    Harvest_date_MLR = Full_bloom_est + round(pmax(Harvest_DAFB_MLR, 0)),
    Harvest_date_ML = Full_bloom_est + round(pmax(Harvest_DAFB_ML, 0)),
    Harvest_doy_MLR = yday(Harvest_date_MLR),
    Harvest_doy_ML = yday(Harvest_date_ML)
  )

write.csv(scenario_data, file.path(OUTPUT_PATH, "harvest_date_predictions.csv"), row.names = FALSE)
cat("\n예측 결과 저장: harvest_date_predictions.csv\n")

# ============================================================================
# 8. 시각화: 시나리오별 수확시기 변화
# ============================================================================

cat("\n========================================\n")
cat("시각화: 시나리오별 수확시기 변화\n")
cat("========================================\n")

harvest_summary <- scenario_data %>%
  group_by(Scenario, Year) %>%
  summarize(
    Harvest_doy_MLR_mean = mean(Harvest_doy_MLR, na.rm = TRUE),
    Harvest_doy_MLR_sd = sd(Harvest_doy_MLR, na.rm = TRUE),
    Harvest_doy_ML_mean = mean(Harvest_doy_ML, na.rm = TRUE),
    Harvest_doy_ML_sd = sd(Harvest_doy_ML, na.rm = TRUE),
    .groups = "drop"
  )

p_mlr <- ggplot(harvest_summary, aes(x = Year)) +
  geom_ribbon(aes(ymin = Harvest_doy_MLR_mean - Harvest_doy_MLR_sd,
                  ymax = Harvest_doy_MLR_mean + Harvest_doy_MLR_sd,
                  fill = Scenario), alpha = 0.2) +
  geom_line(aes(y = Harvest_doy_MLR_mean, color = Scenario), size = 1) +
  scale_color_manual(values = c("SSP1-26" = "blue", "SSP2-45" = "green",
                                "SSP3-70" = "orange", "SSP5-85" = "red")) +
  scale_fill_manual(values = c("SSP1-26" = "blue", "SSP2-45" = "green",
                               "SSP3-70" = "orange", "SSP5-85" = "red")) +
  labs(title = "기후변화 시나리오별 수확시기 예측 (다중회귀 기반 RSM)",
       subtitle = "음영: ±1 표준편차",
       x = "연도", y = "수확시기 (연중일, DOY)") +
  theme_minimal() +
  theme(legend.position = "bottom")

p_ml <- ggplot(harvest_summary, aes(x = Year)) +
  geom_ribbon(aes(ymin = Harvest_doy_ML_mean - Harvest_doy_ML_sd,
                  ymax = Harvest_doy_ML_mean + Harvest_doy_ML_sd,
                  fill = Scenario), alpha = 0.2) +
  geom_line(aes(y = Harvest_doy_ML_mean, color = Scenario), size = 1) +
  scale_color_manual(values = c("SSP1-26" = "blue", "SSP2-45" = "green",
                                "SSP3-70" = "orange", "SSP5-85" = "red")) +
  scale_fill_manual(values = c("SSP1-26" = "blue", "SSP2-45" = "green",
                               "SSP3-70" = "orange", "SSP5-85" = "red")) +
  labs(title = "기후변화 시나리오별 수확시기 예측 (머신러닝 기반 RSM)",
       subtitle = "음영: ±1 표준편차",
       x = "연도", y = "수확시기 (연중일, DOY)") +
  theme_minimal() +
  theme(legend.position = "bottom")

png(file.path(OUTPUT_PATH, "scenario_comparison_harvest_mlr.png"), width = 1000, height = 600)
print(p_mlr)
dev.off()

png(file.path(OUTPUT_PATH, "scenario_comparison_harvest_ml.png"), width = 1000, height = 600)
print(p_ml)
dev.off()

cat("시나리오 비교 그래프 저장: scenario_comparison_harvest_mlr.png, scenario_comparison_harvest_ml.png\n")

# ============================================================================
# 9. 시각화: 지역별 수확시기 변화
# ============================================================================

cat("\n시각화: 지역별 수확시기 변화\n")

regional_summary <- scenario_data %>%
  group_by(Scenario, Year, Region) %>%
  summarize(
    Harvest_doy_MLR = mean(Harvest_doy_MLR, na.rm = TRUE),
    Harvest_doy_ML = mean(Harvest_doy_ML, na.rm = TRUE),
    .groups = "drop"
  )

regional_plots_mlr <- list()
regional_plots_ml <- list()

for(region in unique(regional_summary$Region)) {
  region_data <- regional_summary %>% filter(Region == region)

  p_mlr_region <- ggplot(region_data, aes(x = Year, y = Harvest_doy_MLR, color = Scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("SSP1-26" = "blue", "SSP2-45" = "green",
                                  "SSP3-70" = "orange", "SSP5-85" = "red")) +
    labs(title = paste(region, "- 다중회귀 기반 RSM"),
         x = "연도", y = "수확시기 (연중일, DOY)") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 10))

  p_ml_region <- ggplot(region_data, aes(x = Year, y = Harvest_doy_ML, color = Scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("SSP1-26" = "blue", "SSP2-45" = "green",
                                  "SSP3-70" = "orange", "SSP5-85" = "red")) +
    labs(title = paste(region, "- 머신러닝 기반 RSM"),
         x = "연도", y = "수확시기 (연중일, DOY)") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 10))

  regional_plots_mlr[[region]] <- p_mlr_region
  regional_plots_ml[[region]] <- p_ml_region
}

png(file.path(OUTPUT_PATH, "regional_harvest_mlr.png"), width = 1600, height = 1200)
grid.arrange(grobs = regional_plots_mlr, ncol = 3)
dev.off()

png(file.path(OUTPUT_PATH, "regional_harvest_ml.png"), width = 1600, height = 1200)
grid.arrange(grobs = regional_plots_ml, ncol = 3)
dev.off()

cat("지역별 그래프 저장: regional_harvest_mlr.png, regional_harvest_ml.png\n")

# ============================================================================
# 10. 시각화: 남한 상세지도 기반 수확시기 변화
# ============================================================================

cat("\n========================================\n")
cat("시각화: 남한 상세지도 기반 수확시기 변화\n")
cat("========================================\n")

region_coords <- data.frame(
  Region = c("Gunwi", "Cheongsong", "Yeongju", "Jangsu",
             "Geochang", "Chungju", "Pocheon", "Hwaseong"),
  lat = c(36.24, 36.43, 36.81, 35.65, 35.67, 36.97, 38.11, 37.20),
  lon = c(128.57, 129.16, 128.62, 127.52, 127.91, 127.99, 127.20, 126.82)
)

korea_map <- rnaturalearth::ne_states(country = "South Korea", returnclass = "sf")

map_base <- scenario_data %>%
  filter(Year >= 2071, Year <= 2100) %>%
  group_by(Scenario, Region) %>%
  summarize(
    Harvest_doy_MLR = mean(Harvest_doy_MLR, na.rm = TRUE),
    Harvest_doy_ML = mean(Harvest_doy_ML, na.rm = TRUE),
    .groups = "drop"
  )

map_summary <- map_base %>%
  left_join(region_coords, by = "Region") %>%
  mutate(model = "MLR", Harvest_doy = Harvest_doy_MLR) %>%
  select(Scenario, Region, lat, lon, model, Harvest_doy) %>%
  bind_rows(
    map_base %>%
      left_join(region_coords, by = "Region") %>%
      mutate(model = "ML", Harvest_doy = Harvest_doy_ML) %>%
      select(Scenario, Region, lat, lon, model, Harvest_doy)
  )

p_map <- ggplot() +
  geom_sf(data = korea_map, fill = "gray95", color = "gray70") +
  geom_point(
    data = map_summary,
    aes(x = lon, y = lat, color = Harvest_doy),
    size = 3
  ) +
  scale_color_viridis_c(option = "plasma", name = "수확시기 (DOY)") +
  coord_sf(xlim = c(125.5, 131.0), ylim = c(33.0, 39.5)) +
  facet_grid(model ~ Scenario) +
  labs(
    title = "기후변화 시나리오별 수확시기(2071-2100 평균)",
    subtitle = "남한 상세지도 기반 지역별 예측"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

png(file.path(OUTPUT_PATH, "map_harvest_doy_2071_2100.png"), width = 1600, height = 900)
print(p_map)
dev.off()

cat("지도 그래프 저장: map_harvest_doy_2071_2100.png\n")

# ============================================================================
# 11. 결과 요약
# ============================================================================

cat("\n========================================\n")
cat("결과 요약\n")
cat("========================================\n")
cat("1. harvest_date_predictions.csv - 시나리오별 수확시기 예측 결과\n")
cat("2. scenario_comparison_harvest_mlr.png - 시나리오 비교 (MLR 기반)\n")
cat("3. scenario_comparison_harvest_ml.png - 시나리오 비교 (ML 기반)\n")
cat("4. regional_harvest_mlr.png - 지역별 수확시기 (MLR 기반)\n")
cat("5. regional_harvest_ml.png - 지역별 수확시기 (ML 기반)\n")
cat("6. map_harvest_doy_2071_2100.png - 남한 상세지도 기반 수확시기\n")
