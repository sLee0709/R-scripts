####修改源代码####

add_column_if_missing<-function (df, ...) 
{
  column_values <- list(...)
  for (column_name in names(column_values)) {
    default_val <- rep(column_values[[column_name]], nrow(df))
    if (column_name %in% colnames(df)) {
      df[[column_name]] <- ifelse(is.na(df[[column_name]]), 
                                  default_val, df[[column_name]])
    }
    else {
      df[[column_name]] <- default_val
    }
  }
  df
}

compute_bounds<-function (row_pos, column_pos, segment_data, rect_data, circle_data, 
          funkyrect_data, pie_data, text_data) 
{
  suppressWarnings({
    minimum_x <- min(column_pos$xmin, segment_data$x, segment_data$xend, 
                     rect_data$xmin, circle_data$x - circle_data$r, funkyrect_data$x - 
                       funkyrect_data$r, pie_data$xmin, text_data$xmin, 
                     na.rm = TRUE)
    maximum_x <- max(column_pos$xmax, segment_data$x, segment_data$xend, 
                     rect_data$xmax, circle_data$x + circle_data$r, funkyrect_data$x + 
                       funkyrect_data$r, pie_data$xmax, text_data$xmax, 
                     na.rm = TRUE)
    minimum_y <- min(row_pos$ymin, segment_data$y, segment_data$yend, 
                     rect_data$ymin, circle_data$y - circle_data$r, funkyrect_data$y - 
                       funkyrect_data$r, pie_data$ymin, text_data$ymin, 
                     na.rm = TRUE)
    maximum_y <- max(row_pos$ymax, segment_data$y, segment_data$yend, 
                     rect_data$ymax, circle_data$y + circle_data$r, funkyrect_data$y + 
                       funkyrect_data$r, pie_data$ymax, text_data$ymax, 
                     na.rm = TRUE)
  })
  list(minimum_x = minimum_x, maximum_x = maximum_x, minimum_y = minimum_y, 
       maximum_y = maximum_y)
}

compose_ggplot <- function (geom_positions, position_args) 
{
  g <- ggplot() + coord_equal(expand = FALSE) + scale_alpha_identity() + 
    scale_colour_identity() + scale_fill_identity() + scale_size_identity() + 
    scale_linewidth_identity() + scale_linetype_identity() + 
    cowplot::theme_nothing()
  row_pos <- (geom_positions$row_pos %||% tibble(colour_background = logical(0))) %>% 
    filter(.data$colour_background)
  if (nrow(row_pos) > 0) {
    g <- g + geom_rect(aes(xmin = min(geom_positions$column_pos$xmin) - 
                             0.25, xmax = max(geom_positions$column_pos$xmax) + 
                             0.25, ymin = .data$ymin - (geom_positions$viz_params$row_space/2), 
                           ymax = .data$ymax + (geom_positions$viz_params$row_space/2)), 
                       row_pos, fill = "#DDDDDD")
  }
  if (nrow(geom_positions$segment_data %||% tibble()) > 0) {
    # 这里修改了源代码，把源代码中ticks隐藏了，也就是lintype为solid的给过滤掉了。
    # 过滤 linetype 为 'dashed' 的行
    filtered_segment_data <- geom_positions$segment_data %>%
      filter(linetype == "dashed")
    
    # 检查过滤后的数据行数
    if (nrow(filtered_segment_data) > 0) {
      # 如果存在符合条件的数据，则继续添加列并绘图
      filtered_segment_data <- filtered_segment_data %>% 
        add_column_if_missing(size = 0.5, colour = "black", linetype = "solid")
      
      g <- g + geom_segment(aes(x = .data$x, xend = .data$xend, 
                                y = .data$y, yend = .data$yend, linewidth = .data$size, 
                                colour = .data$colour, linetype = .data$linetype), 
                            filtered_segment_data)
      
    }
  }
  if (nrow(geom_positions$rect_data %||% tibble()) > 0) {
    # 处理 column_id 不为空的部分
    rect_data_non_na <- geom_positions$rect_data %>%
      filter(!is.na(.data$column_id)) %>%  # 筛选 column_id 不为空的行
      add_column_if_missing(alpha = 1, border = TRUE, border_colour = "black") %>%
      mutate(
        border_colour = ifelse(.data$border, .data$border_colour, NA_character_)
      )
    
    # 处理 column_id 为空的部分
    rect_data_na <- geom_positions$rect_data %>%
      filter(is.na(.data$column_id)) %>%  # 筛选 column_id 为空的行
      add_column_if_missing(alpha = 1, border = FALSE, border_colour = NA) %>%
      mutate(
        fill = "gray",  # 示例逻辑：设置默认填充颜色
        xmin = .data$xmin,  # 示例逻辑：调整 xmin 位置
        xmax = .data$xmax,  # 示例逻辑：调整 xmax 位置
        ymin = .data$ymin-0.5,  # 示例逻辑：调整 ymin 位置
        ymax = .data$ymax+0.5   # 示例逻辑：调整 ymax 位置
      )
    
    # 合并两部分数据
    geom_positions$rect_data <- bind_rows(rect_data_non_na, rect_data_na)
    
    # 绘图
    g <- g + geom_rect(aes(xmin = .data$xmin, xmax = .data$xmax, 
                           ymin = .data$ymin, ymax = .data$ymax, 
                           fill = .data$colour, 
                           colour = .data$border_colour, 
                           alpha = .data$alpha), 
                       geom_positions$rect_data, linewidth = 0.25)
    
  }
  if (nrow(geom_positions$circle_data %||% tibble()) > 0) {
    g <- g + ggforce::geom_circle(aes(x0 = .data$x0, y0 = .data$y0, 
                                      fill = .data$colour, r = .data$r), geom_positions$circle_data, 
                                  linewidth = 0.25)
  }
  if (nrow(geom_positions$funkyrect_data %||% tibble()) > 
      0) {
    g <- g + geom_rounded_rect(aes(xmin = .data$xmin, xmax = .data$xmax, 
                                   ymin = .data$ymin, ymax = .data$ymax, radius = .data$corner_size, 
                                   fill = .data$colour), geom_positions$funkyrect_data, 
                               size = 0.25, colour = "black")
  }
  if (nrow(geom_positions$pie_data %||% tibble()) > 0) {
    g <- g + ggforce::geom_arc_bar(aes(x0 = .data$x0, y0 = .data$y0, 
                                       r0 = .data$r0, r = .data$r, start = .data$rad_start, 
                                       end = .data$rad_end, fill = .data$colour), data = geom_positions$pie_data, 
                                   linewidth = 0.25)
  }
  if (nrow(geom_positions$img_data %||% tibble()) > 0) {
    if (!requireNamespace("magick", quietly = TRUE)) {
      cli_alert_warning("Package `magick` is required to draw images. Skipping columns with geom == \"image\".")
    }
    else {
      for (r in seq_len(nrow(geom_positions$img_data))) {
        image <- geom_positions$img_data[[r, "path"]]
        if (!inherits(image, "magick-image")) {
          if (is.character(image)) {
            assert_that(file.exists(image), msg = paste0("Image '", 
                                                         image, "' does not exist."))
          }
          image <- magick::image_read(image)
        }
        g <- g + cowplot::draw_image(image = image, 
                                     x = geom_positions$img_data[[r, "xmin"]], 
                                     y = geom_positions$img_data[[r, "ymin"]], 
                                     width = geom_positions$img_data[[r, "width"]], 
                                     height = geom_positions$img_data[[r, "height"]])
      }
    }
  }
  if (nrow(geom_positions$text_data %||% tibble()) > 0) {
    # 1. 处理 geom 不为空的行
    geom_text_non_na <- geom_positions$text_data %>%
      filter(!is.na(.data$geom)) %>%  # 筛选 geom 列不为空的行
      add_column_if_missing(hjust = 0.5, vjust = 0.5, 
                            size = 4, fontface = "plain", colour = "black", 
                            lineheight = 1, angle = 0) %>% 
      mutate(
        angle2 = .data$angle / 360 * 2 * pi, 
        cosa = cos(.data$angle2) %>% round(2), 
        sina = sin(.data$angle2) %>% round(2), 
        alphax = ifelse(.data$cosa < 0, 1 - .data$hjust, .data$hjust) * abs(.data$cosa) + 
          ifelse(.data$sina > 0, 1 - .data$vjust, .data$vjust) * abs(.data$sina), 
        alphay = ifelse(.data$sina < 0, 1 - .data$hjust, .data$hjust) * abs(.data$sina) + 
          ifelse(.data$cosa < 0, 1 - .data$vjust, .data$vjust) * abs(.data$cosa), 
        x = (1 - .data$alphax) * .data$xmin + .data$alphax * .data$xmax, 
        y = (1 - .data$alphay) * .data$ymin + .data$alphay * .data$ymax
      ) %>% 
      filter(.data$label_value != "")
    
    # 2. 处理 geom 为空且 hjust == 0 的行
    geom_text_na_hjust_0 <- geom_positions$text_data %>%
      filter(is.na(.data$geom) & .data$hjust == 0) %>%  # 筛选 geom 为空且 hjust == 0 的行
      mutate(
        hjust = 1,
        vjust = 0.5,
        size = 4,
        fontface = "plain",
        colour = "black",
        angle = 90,
        angle2 = angle / 360 * 2 * pi,
        label_value = ifelse(is.na(.data$label_value), paste0("Default_", row_number()), .data$label_value),
        x = (.data$xmin + .data$xmax) / 2,  # 示例逻辑，可以根据需求调整
        y = (.data$ymin + .data$ymax) / 2 - 20   # 向下移动
      )
    
    # 3. 处理 geom 为空且 hjust != 0 的行
    geom_text_na_hjust_not_0 <- geom_positions$text_data %>%
      filter(is.na(.data$geom) & .data$hjust != 0) %>%  # 筛选 geom 为空且 hjust != 0 的行
      mutate(
        hjust = .data$hjust,  # 保留原始 hjust
        vjust = 0.3,
        size = 4,
        fontface = "plain",
        colour = "black",
        angle = 0,
        angle2 = angle / 360 * 2 * pi,
        label_value = ifelse(is.na(.data$label_value), paste0("Default_", row_number()), .data$label_value),
        x = (.data$xmin + .data$xmax) / 2,  # 示例逻辑，可以根据需求调整
        y = (.data$ymin + .data$ymax) / 2   # 向上移动
      )
    
    # 合并三部分数据
    geom_positions$text_data <- bind_rows(geom_text_non_na, geom_text_na_hjust_0, geom_text_na_hjust_not_0)
    
    # 绘图
    g <- g + geom_text(aes(x = .data$x, y = .data$y, label = .data$label_value, 
                           colour = .data$colour, hjust = .data$hjust, vjust = .data$vjust, 
                           size = .data$size, fontface = .data$fontface, angle = .data$angle), 
                       data = geom_positions$text_data)
  }
  
  
  if (is.null(geom_positions$bounds)) {
    geom_positions$bounds <- compute_bounds(row_pos = geom_positions$row_pos, 
                                            column_pos = geom_positions$column_pos, segment_data = geom_positions$segment_data, 
                                            rect_data = geom_positions$rect_data, circle_data = geom_positions$circle_data, 
                                            funkyrect_data = geom_positions$funkyrect_data, 
                                            pie_data = geom_positions$pie_data, text_data = geom_positions$text_data)
  }
  minimum_x <- geom_positions$bounds$minimum_x - (position_args$expand_xmin %||% 
                                                    0)
  maximum_x <- geom_positions$bounds$maximum_x + (position_args$expand_xmax %||% 
                                                    0)
  minimum_y <- geom_positions$bounds$minimum_y - (position_args$expand_ymin %||% 
                                                    0)
  maximum_y <- geom_positions$bounds$maximum_y + (position_args$expand_ymax %||% 
                                                    0)
  g <- g + expand_limits(x = c(minimum_x, maximum_x), y = c(minimum_y, 
                                                            maximum_y))
  g$minimum_x <- minimum_x
  g$maximum_x <- maximum_x
  g$minimum_y <- minimum_y
  g$maximum_y <- maximum_y
  g$width <- (maximum_x - minimum_x)/4
  g$height <- (maximum_y - minimum_y)/4
  g
}

assignInNamespace("compose_ggplot", compose_ggplot, ns = "funkyheatmap")
