#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
#test
#install.packages(c("shinylive", "httpuv"))
library(pacman)
p_load(shiny, tidyverse, terra, tidyterra, ggnewscale, shinylive, httpuv, tmap, sf, reactlog, stars, DT) #lwgeom, 
reactlog::reactlog_enable()

#read in data
looking_good <- rast("pansharpened_terra.tif")

the_map <- tm_shape(looking_good) + 
  tm_rgb(col = tm_vars(c(1,2,3), multivariate = TRUE),
         col.scale = tm_scale_rgb(stretch = FALSE, max_color_value = 1)
  )+
  tm_basemap(server = "Esri.WorldGrayCanvas", group = NA, alpha = NA)

points1 <- st_read("./w3_stic_locs_snap.shp") %>% 
  select(pk) %>% 
  rename("value" = pk) %>% 
  mutate("dataset" = "pk")
points2 <- st_read("./depths.shp") %>% 
  select(depth) %>% 
  rename("value" = depth) %>% 
  mutate("dataset" = "depth")

#read in elevation raster
elevation <- (rast("./dem.tif"))

tmap_mode("view")

# Define your working CRS: NAD83 / UTM Zone 19N
utm_crs <- 26919


ui <- fluidPage(
  h3("Interactive Depth-to-Bedrock Profile Analysis"),
  p("Click twice on map to define a rectangle.  Rotate with the rotation angle slider. Click Run Analysis to plot cross-section of surface elevation, depth measurements, and proportion of time flowing of sensor locations below map."),
  actionButton("reset", "Reset Selection"),
  sliderInput("angle", "Rotation Angle (degrees):", min = 0, max = 360, value = 0),
  actionButton("run_analysis", "Run Analysis"),
  tmapOutput("map", height = "600px"),
  plotOutput("profile_plot", height = "300px"),
  actionButton("remove_selected", "Remove Selected Points"),
  actionButton("reset_points", "Reset Removed Points"),
  DT::dataTableOutput("data_table")
  
)

server <- function(input, output, session) {
  click_coords <- reactiveVal(data.frame(lon = numeric(0), lat = numeric(0)))
  
  observeEvent(input$map_click, {
    coord <- input$map_click
    pt <- st_sfc(st_point(c(coord$lng, coord$lat)), crs = 4326) |>
      st_transform(utm_crs)
    coords <- rbind(click_coords(), st_coordinates(pt))
    if (nrow(coords) > 2) coords <- coords[2:3, ]
    click_coords(coords)
  })
  
  observeEvent(input$reset, {
    click_coords(data.frame(lon = numeric(0), lat = numeric(0)))
    analysis_data(NULL)
  })
  
  observeEvent(input$reset_points, {
    excluded_rows(integer(0))
  })
  
  
  get_rotated_rectangle <- reactive({
    coords <- click_coords()
    if (nrow(coords) < 2) return(NULL)
    
    x_range <- range(coords[, 1])
    y_range <- range(coords[, 2])
    
    rect_coords <- rbind(
      c(x_range[1], y_range[1]),
      c(x_range[1], y_range[2]),
      c(x_range[2], y_range[2]),
      c(x_range[2], y_range[1]),
      c(x_range[1], y_range[1])
    )
    
    rect_poly <- st_polygon(list(rect_coords)) |> st_sfc(crs = utm_crs)
    center <- st_centroid(rect_poly)
    
    angle_rad <- input$angle * pi / 180
    rot_mat <- matrix(c(cos(angle_rad), sin(angle_rad), -sin(angle_rad), cos(angle_rad)), nrow = 2)
    
    coords_mat <- st_coordinates(rect_poly)[, 1:2]
    coords_centered <- sweep(coords_mat, 2, st_coordinates(center))
    coords_rotated <- (coords_centered %*% rot_mat) + matrix(st_coordinates(center), nrow(coords_mat), 2, byrow = TRUE)
    
    rotated_poly <- st_polygon(list(coords_rotated)) |> st_sfc(crs = utm_crs)
    st_sf(geometry = rotated_poly)
  })
  
  output$map <- renderTmap({
    rect <- get_rotated_rectangle()
    center <- centerline_geom()
    
    tm <- the_map
    tm <- tm+tm_shape(st_transform(points1, 4326)) + tm_dots(fill = "value", size = 0.6, shape = 24,
                                                             tm_scale_continuous(values = "brewer.blues"),
                                                             fill.legend = tm_legend(title = "Proportion of Time Flowing"),
                                                             group = "Sensor locations")
    tm <- tm + tm_shape(st_transform(points2, 4326)) + tm_dots(fill = "value", size = 0.5, shape = 21,
                                                               fill.scale = tm_scale_continuous(values = "brewer.greys"),
                                                               fill.legend = tm_legend(title = "Depth (m)"),
                                                               group = "Seismic Depths")
    if (!is.null(rect)) {
      tm <- tm + tm_shape(st_transform(rect, 4326)) + tm_borders(col = "red", lwd = 2)
    } 
    if (!is.null(center)) {
      tm <- tm + tm_shape(st_transform(center, 4326)) + tm_lines(col = "black", lwd = 2)
    }
    tm
  })
  
  # Helper: measure projected distance along a centerline
  measure_along_line <- function(line, point) {
    point <- st_transform(point, st_crs(line))
    start <- st_startpoint(line)
    dist <- st_distance(st_startpoint(line), point)
    as.numeric(dist)
    # line_pts <- st_line_sample(line, sample = seq(0, 1, length.out = 1000))
    # dists <- st_distance(line_pts, point)
    #as.numeric(dists) #commented these out
    # min_index <- which.min(dists)
    # total_length <- st_length(line)
    # rel_pos <- (min_index - 1) / (length(line_pts) - 1)
    # as.numeric(rel_pos * total_length)
  }
  
  #reactive values
  analysis_data <- reactiveVal(NULL)
  centerline_geom <- reactiveVal(NULL)
  excluded_rows <- reactiveVal(integer(0))
  
  observeEvent(input$remove_selected, {
    selected <- input$data_table_rows_selected
    current <- excluded_rows()
    
    if (length(selected) > 0) {
      new_excluded <- unique(c(current, selected))
      excluded_rows(new_excluded)
      
      # Clear selection
      DT::selectRows(DT::dataTableProxy("data_table"), NULL)
    }
  })
  
  
  
  
  observeEvent(input$run_analysis, {
    rect <- get_rotated_rectangle()
    if (is.null(rect)) return()
    
    coords <- st_coordinates(rect)[1:4, 1:2]
    edge_lengths <- c(
      sqrt(sum((coords[1, ] - coords[2, ])^2)),
      sqrt(sum((coords[2, ] - coords[3, ])^2))
    )
    
    if (edge_lengths[1] < edge_lengths[2]) { #swapped >
      p1 <- (coords[1, ] + coords[2, ]) / 2
      p2 <- (coords[3, ] + coords[4, ]) / 2
    } else {
      p1 <- (coords[2, ] + coords[3, ]) / 2
      p2 <- (coords[4, ] + coords[1, ]) / 2
    }
    
    centerline <- st_linestring(rbind(p1, p2)) |> st_sfc(crs = utm_crs)
    centerline_geom(centerline)
    
    #centerline <- st_cast(c(p1, p2),"LINESTRING")|> st_sfc(crs = utm_crs)
    sampled_pts <- st_line_sample(centerline, n = 50) |> st_cast("POINT")
    sampled_elev <- terra::extract(elevation, vect(sampled_pts))[, 2]
    
    centerline_profile <- tibble(
      dataset = "Centerline",
      distance_m = seq(0, as.numeric(st_length(centerline)), length.out = length(sampled_pts)),
      elevation = sampled_elev,
      value = NA,
      distance_from_line = 0
    )
    
    
    clip_and_project <- function(points, dataset_name) {
      clipped <- suppressWarnings(st_intersection(points, rect))
      if (nrow(clipped) == 0) return(NULL)
      
      clipped <- st_transform(clipped, st_crs(centerline))  # ensure match
      projected <- st_nearest_points(clipped, centerline)
      snapped <- st_cast(projected, "POINT")[seq(2, length(projected) * 2, by = 2)]
      snapped <- st_set_crs(snapped, st_crs(centerline))
      
      distance_from <- st_distance(clipped, centerline)
      

      dists <- sapply(1:length(snapped), function(i) {
        measure_along_line(centerline, snapped[i])
      })
      
      elevations_ex <- extract(elevation, vect(snapped))
      
      tibble(
        dataset = dataset_name,
        distance_m = dists,
        elevation = elevations_ex$dem,
        value = clipped$value,
        distance_from_line = as.numeric(distance_from)
      )
    }
    
    data1 <- clip_and_project(points1, "pk")  
      
    data2 <- clip_and_project(points2, "depth")
    
    if (!is.null(data1) && !is.null(data2)) {
      analysis_data(bind_rows(data1, data2, centerline_profile))
    } else {
      analysis_data(NULL)
    }
  })
  
  output$profile_plot <- renderPlot({
    df <- analysis_data()
    if (is.null(df)) return()
    excluded <- excluded_rows()
    if (length(excluded) > 0) {
      df <- df[-excluded, ]
    }
    
    ggplot(df, aes(x = distance_m)) +
      geom_line(filter(df, dataset == "Centerline"), mapping = aes(x = distance_m, y = elevation))+
      geom_line(filter(df, dataset == "depth"), mapping = aes(x = distance_m, y = elevation - value), lty = 2)+
      geom_point(filter(df, dataset == "pk"), mapping = aes(x = distance_m, y = elevation, fill = value), pch = 24, size = 3)+
      geom_point(filter(df, dataset == "depth"), mapping = aes(x = distance_m, y = elevation - value), pch = 21, size = 3) +
      scale_fill_distiller("brewer.blues", name = "Proportion of Time Flowing", values = c(1, 0))+
      #geom_smooth(se = FALSE) +
      labs(x = "Distance along centerline (m)", y = "Elevation (m)") +
      theme_classic()
  })
  
  
  output$data_table <- DT::renderDataTable({
    df <- analysis_data()
    if (is.null(df)) return()
    
    df <- df |> filter(dataset != "Centerline")
      #arrange(dataset, distance_m)
    
    excluded <- excluded_rows()
    if (length(excluded) > 0) {
      df <- df[-excluded, ]
    }
    
    DT::datatable(
      df,
      selection = list(mode = "multiple", selected = NULL, target = 'row'),
      options = list(pageLength = 10,
                     searching = FALSE,
                     lengthChange = FALSE)
    ) %>% formatRound(columns = c("distance_m","elevation", "value", "distance_from_line"), digits = 1)
  })
  
  
}

shinyApp(ui, server)

#shinylive::export(appdir = "seismicMap", destdir = "docs")
#runStaticServer("docs")
