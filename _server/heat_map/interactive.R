
#----------------------------------------------------------------------
# interactive.R handles user interactions with data plots
#----------------------------------------------------------------------

# enable interactive brush selection on the heat map to set the viewport region
observeEvent(input$heatMapBrush, {
    reportProgress('heatMapBrush')
    minStackI <- round(input$heatMapBrush$ymin, 0)
    maxStackI <- round(input$heatMapBrush$ymax, 0)
    brushSelection_$cells <- viewport_$stack[minStackI:maxStackI]
    xScalar <- if(input$chrom == "all") 1 else 1e6
    updateTextInput(session, 'start', value = round(input$heatMapBrush$xmin * xScalar, 0))
    updateTextInput(session, 'end',   value = round(input$heatMapBrush$xmax * xScalar, 0) - 1)  
})

