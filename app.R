options(shiny.maxRequestSize=30*1024^2) #increase size limit of uploaded files to 30MB

list_of_packages = c('tidyverse','DT', 
                     'shiny', 'shinyjs', 'shinyWidgets',
                     'sf', 'spdep', 'leaflet', 'leaflet.providers',
                     'scales', 'RColorBrewer', 'viridis')

lapply(list_of_packages, 
       function(x) if(!require(x,character.only = TRUE)) install.packages(x))

#general packages
library(tidyverse)
library(DT)
#shiny packages
library(shiny)
library(shinyjs)
library(shinyWidgets)
#mapping packages
library(sf)
library(spdep)
library(leaflet)
library(leaflet.providers)
#data viz packages
library(scales)
library(RColorBrewer)
library(viridis)

#load in shapefiles for plotting
uk.areas.sf = read_sf("shapefiles/Areas.shp")
uk.districts.sf = read_sf("shapefiles/Districts_small.shp")

#create centroids for hotspot analysis
centroids <- uk.districts.sf %>%
  st_make_valid() %>%
  st_centroid() %>%
  st_coordinates() %>%
  as.data.frame() %>%
  rename(centroid_long = X, centroid_lat = Y) %>%
  mutate(name = uk.districts.sf$name) %>%
  dplyr::select(name, centroid_long, centroid_lat)

#function for calculating % use of variant for choropleth maps
calcVals <- function(data, location, variable, level) {
  #calculate number of responses by the location type of interest
  num_responses <- data %>%
    rename('var' = !! variable) %>%
    filter(!is.na(var) & var != '') %>%
    group_by(!! location) %>%
    summarise(total = n())
  
  #tabulate responses in each region, filter down to variant of interest,
  #then join with num_responses dataset and calculate proportion:
  data %>%
    rename('var' = !! variable) %>%
    #mutate(var = case_when(                  #only needed if counting 'more than ones'
    #  str_detect(var, level) ~ level,
    #  TRUE ~ var)) %>%
    count(!! location, var) %>%
    filter(var == level) %>%             #if only want to count people with that sole response
    #filter(str_detect(var, level)) %>%    #if counting people who answered 'more than one'
    dplyr::select(!! location, n) %>%
    right_join(num_responses) %>%
    mutate(value = case_when(
      n > 0 ~ round((n / total)*100, 2),
      is.na(n) ~ 0)) %>%
    dplyr::select(-n) %>%
    ungroup() %>%
    rename('name' = !! location)
}


###
##### SPATIAL WEIGHTS MATRIX
###

#best method: uses reciprocal distance, AND takes itself into account (sets distance as half of next nearest)
calculateSWM <- function(data, k) {
  data %<>%
    as.data.frame() %>%
    left_join(centroids) %>%
    filter(!is.na(value) & !is.na(centroid_long)) %>%
    dplyr::select(name, centroid_long, centroid_lat, value)
  
  #create the nearest neighbour weights matrix
  xy = as.matrix(data.frame(data$centroid_long, data$centroid_lat))
  
  #find nearest neighbours
  neighbours = knn2nb(knearneigh(xy, k=k-1, longlat=TRUE))
  neighbours = include.self(neighbours)
  
  #calculate distances between neighbours
  dsts = nbdists(neighbours, xy, longlat=TRUE) 
  
  #replace 0s in the distance matrix with the next lowest distance from each set of neighbours
  #(needed when including self as a neighbour because itself has a distance of 0 and you can't do reciprocal)
  for (x in 1:length(dsts)) {
    dsts[[x]] = replace(dsts[[x]], dsts[[x]] == 0, min(dsts[[x]][dsts[[x]] > 0]))
  }
  
  #weighting of the neighbors by reciprocal distance (1/x) (closer neighbour = higher value = more weight)
  recipdsts = lapply(dsts, function(x) 1/x) 
  
  #return this in a format for spdep
  nb2listw(neighbours, glist=recipdsts) 
}

#create a matrix of all postcode district centroids
xy <- centroids %>% dplyr::select(-name) %>% as.matrix()

#use this to calculate the 250 nearest neighbours for each postcode district 
neighbours <- xy %>% knearneigh(k=250, longlat=TRUE) %>% knn2nb
dsts <- nbdists(neighbours, xy, longlat=TRUE)

interpolateMissing <- function(getis_data, no_neighbours) {
  #postcode districts for which we have no responses
  missing_vals <- uk.districts.sf %>%
    filter(!name %in% getis_data$name) %>%
    left_join(centroids) %>%
    as.data.frame() %>%
    dplyr::select(name, centroid_long, centroid_lat) %>%
    mutate(value = NA) #%>%
    #rename(pc.district = name)
  
  #for each postcode district with a missing value, ...
  for (x in 1:nrow(missing_vals)) {
    index <- which(centroids$name %in% missing_vals[x, ]$name) #gets index of missing pc in list
    nn_names <- centroids$name[neighbours[[index]]] #finds its nearest neighbours
    nn_dsts <- dsts[[index]] #extracts distances of these nearest neighbours
    
    #get a final list of the k nearest neighbours, their values, and their distances
    nn_df <- data.frame(nn_names, nn_dsts) %>%
      left_join(getis_data, by=c('nn_names'='name')) %>%
      arrange(nn_dsts) %>%
      filter(!is.na(value)) %>%
      head(n = no_neighbours)
    
    #interpolate by calculating a weighted mean (weighted by reciprocal distance)
    est.val <- weighted.mean(nn_df$value, 1/nn_df$nn_dsts)
    
    #add to the missing_vals dataframe
    missing_vals[x, ]$value = est.val
  }
  
  missing_vals %<>%
    mutate(total = 0) %>%
    mutate(value.raw = NA) %>%
    dplyr::select(name, total, value, value.raw)
  
  #combine the original z-score dataframe with the newly-interpolated values
  rbind(getis_data, missing_vals)
}



#setup categorical colour palettes for point maps
pointPal <- data.frame(
  val = c('rainbow', 'Dark2', 'Accent', 'Set3', 'Pastel1')
)
pointPal$img = c(
  sprintf("<img src='rainbow_pal.png' width=100px height=20px><div class='jhr'>%s</div></img>", pointPal$val[1]),
  sprintf("<img src='dark2_pal.png' width=100px height=20px><div class='jhr'>%s</div></img>", pointPal$val[2]),
  sprintf("<img src='accent_pal.png' width=100px height=20px><div class='jhr'>%s</div></img>", pointPal$val[3]),
  sprintf("<img src='set3_pal.png' width=100px height=20px><div class='jhr'>%s</div></img>", pointPal$val[4]),
  sprintf("<img src='pastel1_pal.png' width=100px height=20px><div class='jhr'>%s</div></img>", pointPal$val[5])
)

#setup continuous colour palettes for choropleth maps
polyPal <- data.frame(
  val = c('YlOrRd', 'Blues', 'GnBu', 'viridis', 'magma')
)
polyPal$img = c(
  sprintf("<img src='ylorrd_pal.png' width=100px height=20px><div class='jhr'>%s</div></img>", polyPal$val[1]),
  sprintf("<img src='blues_pal.png' width=100px height=20px><div class='jhr'>%s</div></img>", polyPal$val[2]),
  sprintf("<img src='gnbu_pal.png' width=100px height=20px><div class='jhr'>%s</div></img>", polyPal$val[3]),
  sprintf("<img src='viridis_pal.png' width=100px height=20px><div class='jhr'>%s</div></img>", polyPal$val[4]),
  sprintf("<img src='magma_pal.png' width=100px height=20px><div class='jhr'>%s</div></img>", polyPal$val[5])
)


############
######## UI
############

ui <- fluidPage(
  useShinyjs(),
  tags$head(tags$style(
    HTML("
         #titlePanel {
           color: white;
           background: #0B517F;
           border-style: solid;
           border-size: 2px;
           border-radius: 20px;
           width: 99%;
           height: auto;
           padding: 10px 20px 20px 20px;
           margin: 10px 10px 10px 10px;
           display: inline-block;
         }
         
         #subtext {
           color: grey;
           font-style: italic;
           font-size: 12px;
         }
         
         hr {
            border-top: 1px solid #b4b4b4;
         }
         
         #titlePanel a {
            color: #a0e9f9;
            font-weight: bold;
         }
         
         #bannerLogo {
            display: inline-block;
            float: right;
         }
         
         .fa {
            vertical-align: middle;
            padding: 10px 10px 10px 0px;
         }
         
         .jhr{
            display: inline;
            vertical-align: middle;
            padding-left: 10px;
         }
          
          .colname {
          color: #FDDAFF;
          font-weight: bold;
          }
         "))),
  
  tags$div(id="titlePanel",
           #img(id='bannerLogo', src = 'map_plotter_logo.png', height = 165),
           titlePanel(windowTitle = 'Dialect Map Plotter', 
                      title = HTML("<b>Dialect Map Plotter</b>")),
           HTML('<p>Upload a spreadsheet of dialect survey responses in <b>.csv</b> format and this tool will generate interactive dialect maps plots made using <a href="https://rstudio.github.io/leaflet/" target="_blank">leaflet</a> and <a href="https://shiny.rstudio.com/" target="_blank">Shiny</a> in <a href="https://www.r-project.org/about.html" target="_blank">R</a>.</p> <p>For full functionality, the spreadsheet you upload should contain columns labelled <span class="colname">age</span>, <span class="colname">dob</span> (date of birth), <span class="colname">sex</span>, <span class="colname">longitude</span>, <span class="colname">latitude</span>, <span class="colname">pc.district</span> (postcode district), <span class="colname">pc.area</span> (postcode area), and then whatever columns correspond to your dialect features. Click <a href="dialect_survey_demo.csv" target="_blank">here</a> to download a demo dataset to use with the app, which demonstrates how your own spreadsheet should be formatted.</p> <p style="float:left;">Developed by <b><a href="https://www.gbailey.uk/" target="_blank">George Bailey</a></b> (2022) at the University of York.</p>'),
  ),
  
  sidebarLayout(
    sidebarPanel(width = 4,
      tags$div(id = "dataSelectionPanel",
        fileInput(inputId = "filedata", label = "1. Upload dialect survey responses (csv format)", 
                  accept = c(".csv")),
        uiOutput('pickFeature'), uiOutput('pickMapType'),
        uiOutput('pickLevel'), 
        uiOutput('pickPointLevels'), 
        uiOutput('pickGender'), uiOutput('pickAge'), uiOutput('pickDoB'), 
        uiOutput('pickPointPal'),
        uiOutput('pickPolyPal')
      ),
      
      conditionalPanel(
        condition = "input.varSel != '' & ((input.mapType == 'points' & input.pointLvlSel != '') | (input.mapType != 'points' & input.lvlSel != ''))",
        uiOutput('plotButtonUI'))
    ),
    
    mainPanel(width = 8,
              tabsetPanel(
                tabPanel("Map", leafletOutput(outputId = "map", width = "100%", height = "700")),
                tabPanel("Data", DTOutput(outputId = "table"))))
  )
)

############
#### SERVER
############

server <- function(input, output) {
  
  mapSelection <- reactiveValues(type = 'NA',
                                 polyPal = 'NA',
                                 pointPal = 'NA"')
  
  observeEvent(input$plotButton, {
    mapSelection$type = input$mapType
    mapSelection$polyPal = input$polyPalSel
    mapSelection$pointPal = input$pointPalSel
  })
  
  data <- reactive({
    req(input$filedata)
    read.csv(input$filedata$datapath)
  })
  
  # RENDER UI ELEMENTS FOR DATA SELECTION
  
    output$pickFeature <- renderUI({
    feature_list <- data() %>%
      select(bread:last_col()) %>%
      colnames()
    
    selectInput("varSel", "2. Select variable to plot:",
                choices = c('Choose one' = '', feature_list),
                selectize = T)
  })
    
  output$pickMapType <- renderUI({
      req(data())
      radioButtons(
        inputId = "mapType", "3. Choose map type", 
        choices = list('points' = 'points', 
                       'choropleth (by postcode area)' = 'pc.area',
                       'choropleth (by postcode district)' = 'pc.district', 
                       'smoothed choropleth (by postcode district)' = 'pc.district.getis'))
    })
  
  output$pickLevel <- renderUI({
    req(input$varSel)
    if (input$mapType != 'points') {
      levels <- data()[[input$varSel]] %>% table() %>% sort(decreasing = TRUE) %>% names()
      
      selectInput("lvlSel", HTML('4. Select response to plot<br><span id="subtext">(ordered here by overall frequency)</span>'),
                  choices = c('Choose one' = '', levels), selectize = T)
    }
  })
  
  output$pickPointLevels <- renderUI({
    req(input$varSel)
    if (input$mapType == 'points') {
      levels <- data()[[input$varSel]] %>% table() %>% sort(decreasing = TRUE) %>% names()
      
      selectInput("pointLvlSel", HTML('4. Select responses to plot<br><span id="subtext">(select as many as you want; ordered here by overall frequency)</span>'),
                  choices = c('Choose responses' = '', levels),
                  multiple = T, selectize = T)
    }
  })
  
  output$pickAge <- renderUI({
    req(data())
    sliderInput(
      inputId = "ageSel", label = "Filter by age?",
      min = min(data()$age), max = max(data()$age), value = c(min(data()$age), max(data()$age)))
  })
  
  output$pickDoB <- renderUI({
    req(data())
    sliderInput(
      inputId = "dobSel", label = "Filter by date of birth?",
      min = min(data()$dob), max = max(data()$dob), value = c(min(data()$dob), max(data()$dob)),
      sep = "")
  })
  
  output$pickGender <- renderUI({
    req(data())
    selectInput(
      inputId = "sexSel", label = HTML(paste(h4('Filtering the data'), "Filter by gender?")),
      choices = c("no", unique(data()$sex)), selected = 'no')
  })
  
  output$pickPointPal <- renderUI({
    req(data())
    if (input$mapType == 'points') {
      pickerInput(inputId = "pointPalSel",
                  label = HTML(paste(hr(), "Choose a colour palette")),
                  choices = pointPal$val,
                  choicesOpt = list(content = pointPal$img))
    }
  })
  
  output$pickPolyPal <- renderUI({
    req(data())
    if (input$mapType != 'points') {
      pickerInput(inputId = "polyPalSel",
                  label = HTML(paste(hr(), "Choose a colour palette")),
                  choices = polyPal$val,
                  choicesOpt = list(content = polyPal$img))
    }
  })
  
  output$plotButtonUI <- renderUI({
    req(data())
    actionButton(inputId = 'plotButton', label = HTML('Plot/update map<br><span style="font-size: 12px;">(this might take a while depending<br>on the data and style of map)</span>'),
                 icon = icon("map", style = "font-size: 30px;"),
                 style="color: #fff; background-color: #009961; font-size: 17px; width: 100%;")
  })
  
  # UPDATE MAPPING DATA BASED ON SELECTIONS
  
  mapData <- eventReactive(input$plotButton, {
    
    if (input$sexSel == 'no') {
      tempData <- data() %>%
        filter(age >= input$ageSel[1] & age <= input$ageSel[2]) %>%
        filter(dob >= input$dobSel[1] & dob <= input$dobSel[2])
    } else {
      tempData <- data() %>%
        filter(sex == input$sexSel) %>%
        filter(age >= input$ageSel[1] & age <= input$ageSel[2]) %>%
        filter(dob >= input$dobSel[1] & dob <= input$dobSel[2])
    }
    
    if (input$mapType == 'points') {
      tempData %>%
        mutate(vartoplot = tempData[[input$varSel]]) %>%
        filter(vartoplot %in% input$pointLvlSel)
    } else if (input$mapType == 'pc.area') {
      vals <- calcVals(tempData, quo(pc.area), input$varSel, input$lvlSel)
      uk.areas.sf <- left_join(uk.areas.sf, vals)
      uk.areas.sf
    } else if (input$mapType == 'pc.district') {
      vals <- calcVals(tempData, quo(pc.district), input$varSel, input$lvlSel)
      uk.districts.sf <- left_join(uk.districts.sf, vals)
      uk.districts.sf
    } else if (input$mapType == 'pc.district.getis') {
      vals <- calcVals(tempData, quo(pc.district), input$varSel, input$lvlSel) %>%
        filter(name %in% centroids$name)
      swm <- calculateSWM(vals, 25)
      vals <- vals %>%
        mutate(value.raw = value) %>%
        mutate(value = localG(value, swm)) %>%
        interpolateMissing(30)
      uk.districts.sf <- left_join(uk.districts.sf, vals)
      uk.districts.sf
    }

  })
  
  output$map <- renderLeaflet({
    if (is.null(data())) { return(NULL) }

    if (mapSelection$type != 'points') {

      #CHOROPLETH MAP

      pal <- colorBin(mapSelection$polyPal, domain = mapData()$value, bins = 7)

      if (mapSelection$type == 'pc.district.getis') {
        labels <- sprintf("<strong>%s</strong><br/>z = %g (N = %s)",
                          mapData()$name, mapData()$value, mapData()$total) %>%
          lapply(htmltools::HTML)
      } else {
        labels <- sprintf("<strong>%s</strong><br/>%g%% (N = %s)",
                          mapData()$name, mapData()$value, mapData()$total) %>%
          lapply(htmltools::HTML)
      }
      
      l <- leaflet(mapData()) %>%
        addProviderTiles("CartoDB.Positron") %>%
        setView(-4.174805, 54.757916, zoom = 6) %>%
        addPolygons(
          fillColor = ~pal(value),
          weight = 1,
          color = "white",
          fillOpacity = 0.7,
          highlightOptions = highlightOptions(
            weight = 3,
            color = "#666",
            fillOpacity = 0.7,
            bringToFront = TRUE),
          label = labels,
          labelOptions = labelOptions(
            style = list("font-weight" = "normal", padding = "3px 8px"),
            direction = "auto",
            textsize = "15px")
        ) %>%
        leaflet::addLegend(
          pal = pal, values = ~value,
          opacity = 0.7, title = NULL
        )

    } else {

      #POINT MAP

      if (mapSelection$pointPal == 'rainbow') {
        pal <- colorFactor(rainbow(length(unique(mapData()$vartoplot))), 
                           mapData()$vartoplot)
      } else {
        pal <- colorFactor(brewer.pal(length(unique(mapData()$vartoplot)), mapSelection$pointPal), 
                           mapData()$vartoplot)
      }

      l <- leaflet(mapData()) %>%
        addProviderTiles("CartoDB.Positron", group="Light mode") %>%
        addProviderTiles("CartoDB.DarkMatter", group="Dark mode") %>%
        setView(-4.174805, 54.757916, zoom = 6) %>%
        addCircleMarkers(
          lng = ~longitude,
          lat = ~latitude,
          stroke = FALSE, fillOpacity = .7, radius = 3,
          color = ~pal(vartoplot),
        ) %>%
        leaflet::addLegend(
          pal = pal, values = ~vartoplot,
          opacity = 0.7, title = NULL
        ) %>%
        addLayersControl(
          baseGroups = c("Light mode", "Dark mode"),
          options = layersControlOptions(collapsed = FALSE)
        )
    }


  })
  
  output$table <- renderDT(data())
}

shinyApp(ui, server)