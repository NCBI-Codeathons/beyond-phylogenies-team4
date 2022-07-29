#SHINY APP

library(shiny)
## Use R Shiny package for interactive plotting of distribution data
ui <- fluidPage(
  titlePanel(title=h4("Cluster trait distributions", align="center")),
  sidebarPanel( 
    sliderInput("date_spec", "Date:", min = min_date, max = max_date, step=7, value=c(min_date, max_date))),
  selectInput(inputId = "cluster", label = "cluster", 
              choice = cluster),
  
  mainPanel(plotOutput("plot2")))

server <- function(input,output){
  
  dat <- reactive({
    test <- traits[traits$DATE %in% seq(from=min(input$date_spec),to=max(input$date_spec),by=1),]
    #     print(test)
    test
  })
  
  output$plot2<-renderPlot({
    ggplot(dat()) +
      geom_density(aes(x=DATE, fill=trait), alpha=0.8) +
      facet_wrap(cluster_id~field, scales="free_y", nrow=length(unique(traits$cluster_id)), ncol=length(unique(traits$field))) +
      theme(text = element_text(size=20), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20) ) +
      theme_minimal()
  }, height = 2000,width = 2000)
}


shinyApp(ui, server)
