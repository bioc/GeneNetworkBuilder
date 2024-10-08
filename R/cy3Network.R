#' visualize the network by Cy3
#' @description
#' Using RCy3 to visualize the network.
#' 
#' @param gR an object of \link[graph:graphNEL-class]{graphNEL}
#' @param ... parameters will be passed to \link[RCy3]{createNetworkFromGraph}
#' @param stringify Run STRINGify or not
#' @param species if stringify is TRUE, teh species will be used to retreive
#' the nodes and edges properties.
#' @param style The default style when create the network
#' @param widths The link width range.
#' @return The network SUID.
#' @importFrom RCy3 createNetworkFromGraph getVisualStyleNames setVisualStyle
#' setNodeSizeMapping setNodeBorderColorBypass setNodeColorBypass 
#' setEdgeLineWidthMapping commandsRun getNetworkSuid getTableColumns
#' lockNodeDimensions
#' @export
#' @examples
#' data("ce.miRNA.map")
#' data("example.data")
#' data("ce.interactionmap")
#' data("ce.IDsMap")
#' sifNetwork<-buildNetwork(example.data$ce.bind, ce.interactionmap, level=2)
#' cifNetwork<-filterNetwork(rootgene=ce.IDsMap["DAF-16"], sifNetwork=sifNetwork, 
#'   exprsData=uniqueExprsData(example.data$ce.exprData, "Max", condenseName='logFC'),
#'   mergeBy="symbols",
#'   miRNAlist=as.character(ce.miRNA.map[ , 1]), tolerance=1)
#' gR<-polishNetwork(cifNetwork)
#' if(interactive()){
#'   cy3Network(gR)
#' }
cy3Network <- function(gR = graphNEL(), ...,
                       stringify=FALSE, species = 'Homo sapiens',
                       style='Marquee', widths = c(0.25, 5)){
  stopifnot(is(gR, 'graphNEL'))
  availableStyles <- getVisualStyleNames()
  style <- match.arg(style, choices = availableStyles)
  network <- createNetworkFromGraph(gR, ...)
  if(stringify){
    ## doc https://github.com/RBVI/stringApp/blob/master/doc_automation.md#compound-query
    string.cmd = paste0('string stringify column="name" networkNoGui="current" species="', species, '"')
    commandsRun(string.cmd)
    network <- getNetworkSuid()
    stringstyle <- "STRING - From graph"
    if(!stringstyle %in% getVisualStyleNames()){
      stop('Default string style is not in list.',
           'Please select correct one from the output of getVisualStyleNames()')
    }else{
      style <- stringstyle
    }
  }else{
    setVisualStyle(style, network = network)
    lockNodeDimensions(TRUE, style.name = style)
    setNodeSizeMapping('size', style.name = style, network = network)
  }
  nodeData <- getTableColumns(table = 'node', network = network)
  borderColor <- graph::nodeRenderInfo(gR, 'col')
  if(length(borderColor)>0){
    if(stringify){
      N <- nodeData$name[match(names(borderColor), nodeData$`query term`)]
    }else{
      N <- names(borderColor)
    }
    borderColor[is.na(borderColor)] <- '#000000'
    setNodeBorderColorBypass(N[!is.na(N)], 
                             borderColor[!is.na(N)],
                             network = network)
  }
  color <- graph::nodeRenderInfo(gR, 'fill')
  if(length(color)>0){
    if(stringify){
      N <- nodeData$name[match(names(color), nodeData$`query term`)]
    }else{
      N <- names(color)
    }
    color[is.na(color)] <- '#FFFFFF'
    setNodeColorBypass(N[!is.na(N)],
                       new.colors=color[!is.na(N)],
                       network = network)
  }
  setEdgeLineWidthMapping('weight',
                          style.name = style,
                          widths = widths,
                          network = network)
  return(network)
}
