#Author: Mostafa Rahnama 
#Date: 3/28/20
#### Input: delta file generated with nucmmer 

library(dplyr)
library(magrittr)
library(GenomicRanges)
library(knitr)
library(ggplot2)
library(tidyr)
library(shiny)

# change to specify infile on request?

file = ("~/FHC121vCBS133.delta.delta")
name <- strsplit(file, "/")
name <- noquote(name)
name <- tail(name[[1]], n =1)
name <- paste(name[1])
name <- gsub(".delta", "", name)

#a function to read a delta file
readDelta <- function(deltafile){
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  lines = lines[-1]
  lines.l = strsplit(lines, ' ')
  lines.len = lapply(lines.l, length) %>% as.numeric
  lines.l = lines.l[lines.len != 1]
  lines.len = lines.len[lines.len != 1]
  head.pos = which(lines.len == 4)
  head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
  mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
  res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}


# reading dalta file
mumgp = readDelta(file)


mumgp %>% head %>% kable

#Filter contigs with poor alignments
filterMum <- function(df, minl=10000){   #, flanks=1e4
  df = df %>% filter(abs(re-rs)>minl) %>% group_by(qid, rid) # %>%
    #summarize(qsL=min(qs)-flanks, qeL=max(qe)+flanks, rs=median(rs)) %>%
   # ungroup %>% arrange(desc(rs)) %>%
   # mutate(qid=factor(qid, levels=unique(qid))) %>% select(-rs)
  #merge(df, coord) %>% filter(qs>qsL, qe<qeL) %>%
  #  mutate(qid=factor(qid, levels=levels(coord$qid))) %>% select(-qsL, -qeL)
}

mumgp.filt = filterMum(mumgp) #, minl=1e4
mumgp.filt %>% head %>% kable

#############Diagonalize
#For each contig, I compute the major strand (strand with most bases aligned) and flip if necessary. 
#The contigs are also ordered based on the reference region with most bases and 
# the weighted means of the start position in this matched reference region.

diagMum <- function(df){
  ## Find best qid order
  rid.o = df %>% group_by(qid, rid) %>% summarize(base=sum(abs(qe-qs)),
                                                  rs=weighted.mean(rs, abs(qe-qs))) %>%
    ungroup %>% arrange(desc(base)) %>% group_by(qid) %>% do(head(., 1)) %>%
    ungroup %>% arrange(desc(rid), desc(rs)) %>%
    mutate(qid=factor(qid, levels=unique(qid)))
  ## Find best qid strand
  major.strand = df %>% group_by(qid) %>%
    summarize(major.strand=ifelse(sum(sign(qe-qs)*abs(qe-qs))>0, '+', '-'),
              maxQ=max(c(qe, qs)))
  merge(df, major.strand) %>% mutate(qs=ifelse(major.strand=='-', maxQ-qs, qs),
                                     qe=ifelse(major.strand=='-', maxQ-qe, qe),
                                     qid=factor(qid, levels=levels(rid.o$qid)))
}

mumgp.filt.diag = diagMum(mumgp.filt)
head(mumgp.filt.diag)
#mumgp.filt.diag$rid <- gsub("chr", "Chr", mumgp.filt.diag$rid)

#mumgp.filt.diag <- subset(mumgp.filt.diag, qid %in% longer300$V1)   #filter

#mumgp.filt.diag$qid = factor(mumgp.filt.diag$qid, levels = c("MiniChr2", "MiniChr1","Chr7","Chr6","Chr5","Chr4","Chr3","Chr2", "Chr1"))


mumgp.filt.diag$qid = factor(mumgp.filt.diag$qid, levels = unique(mumgp.filt.diag$qid))
#mumgp.filt.diag <- filter(mumgp.filt.diag, rid != "scaf3" )
#mumgp.filt.diag <- filter(mumgp.filt.diag, rid != "scaf4" )
#mumgp.filt.diag <- filter(mumgp.filt.diag, rid != "scaf5")

rchromoList <- unique(mumgp.filt.diag$rid)
qchromoList <- unique(mumgp.filt.diag$qid)
chromoLengths <- c()
for(i in 1:length(rchromoList)) {
  chr =  filter (mumgp.filt.diag, rchromoList[i] == rid)
  len = max(chr$re)
  chromoLengths[i] = len
  }


ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
#        fileInput('file1', 'Select the XXX.delta file',
#                 accept=c('text/csv','text/tab-separated-values,text/plain','.delta')),
        
       checkboxGroupInput(inputId = "inputqChromo", "Chromosome of Query to plot:",
                          inline = FALSE, choices = qchromoList, selected = qchromoList),       
       checkboxGroupInput(inputId = "inputrChromo", "Chromosome of Reference to plot:",
                          inline = FALSE, choices = rchromoList, selected = rchromoList),

    sliderInput("xaxis", "Select chromosome region range:",
                min = 0, max = max(chromoLengths), value = c(1,max(chromoLengths)))),
  
    mainPanel(plotOutput("Windows"), width = 8)
  )
)


server <- function(input, output, session){
  df <- reactive(mumgp.filt.diag %>%
                   filter (rid %in% input$inputrChromo) %>%
                  filter(qid %in% input$inputqChromo) )
  region <- reactive({input$region}   )
  
  output$Windows <- renderPlot({
       
        ggplot(data = df(), 
               aes(x=rs/1000000, xend=re/1000000, y=qs/1000000, yend=qe/1000000, colour=strand)) +
      geom_segment() + 
      geom_point(alpha=.5) + theme_bw() + 
      labs(subtitle = "") + #, caption = "CD156 assembly"
      facet_grid(qid~rid, scales='free', space='free', switch='y') +
      scale_y_continuous(breaks=c(0,0.5,1, 2, 4,6,8), position="right") +
      scale_x_continuous(breaks=c(0, 2, 4,6,8)) + 
      #scale_fill_discrete(guide = guide_legend(keywidth = 10, keyheight = 10)) +
      labs(title = paste(name,"\n"), x = "B71 Chromosome position (Mb)\n", y= "LpKY97 Chromosome position (Mb)\n", color = "Strand")+
            theme( plot.margin=unit(c(0, 0, 0, 3), units="line"), # t, r, b, l Dimensions of each margin 
                   plot.title = element_text(hjust = 1.5),
                           strip.text.y.left=element_text(angle=0, size=20,face="bold", margin = margin(0,0.5,0,1, "cm")),
                           strip.text.x=element_text(angle=90, size=20,face="bold"), 
                           strip.background=element_blank(),
                           legend.position="bottom", #c(.99,.01)
                           #legend.key = element_rect(fill = "white"),
                           #legend.justification=c(1,0),
                           legend.title=element_text(face="bold", size=25), 
                           legend.text=element_text(face="bold", size=25),
                           legend.key.size = unit(5,"cm"),
                           legend.key.width = unit(5, "cm"),
                           #legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="white"),
                           panel.spacing.y = unit(2, "mm"),
                           #panel.margin = unit(0.3, "cm"),
                           #plot.margin = margin(2, 2, 2, 2, "cm"),
                           axis.title.x = element_text(size = 40,face="bold", margin = margin(1, 20, 20, 20), vjust = -3),
                           axis.title.y = element_text(size = 40,face="bold",  vjust = -1),
                           axis.text.y=element_text(size=rel(1.5), colour = "black",face="bold"),
                           axis.text.x=element_text(size=rel(1.5), colour = "black", face="bold"), #angle = 90,
                           axis.ticks.length = unit(0.3,"cm"),
                           axis.ticks = element_line(linewidth=1, colour = "black"),
                           plot.subtitle = element_text(size = 22, face = "bold", hjust = 0.5)
                           #plot.caption = element_text(size = 22, face = "bold", hjust = 0.5, vjust = 1, angle = 90),
             ) + xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set1')
    
    p <- ggplot(data = df(), 
           aes(x=rs/1000000, xend=re/1000000, y=qs/1000000, yend=qe/1000000, colour=strand)) +
      geom_segment() + 
      geom_point(alpha=.5) + theme_bw() + 
      labs(subtitle = "") + #, caption = "CD156 assembly"
      facet_grid(qid~rid, scales='free', space='free', switch='y') +
      scale_y_continuous(breaks=c(0,0.5,1, 2, 4,6,8), position="right") +
      scale_x_continuous(breaks=c(0, 2, 4,6,8)) + 
      #scale_fill_discrete(guide = guide_legend(keywidth = 10, keyheight = 10)) +
      labs(title = paste(name,"\n"), x = "B71 Chromosome position (Mb)\n", y= "LpKY97 Chromosome position (Mb)\n", color = "Strand")+
      theme( plot.margin=unit(c(0, 0, 0, 3), units="line"), # t, r, b, l Dimensions of each margin 
             plot.title = element_text(hjust = 1.5),
             strip.text.y.left=element_text(angle=0, size=20,face="bold", margin = margin(0,0.5,0,1, "cm")),
             strip.text.x=element_text(angle=90, size=20,face="bold"), 
             strip.background=element_blank(),
             legend.position="bottom", #c(.99,.01)
             #legend.key = element_rect(fill = "white"),
             #legend.justification=c(1,0),
             legend.title=element_text(face="bold", size=25), 
             legend.text=element_text(face="bold", size=25),
             legend.key.size = unit(5,"cm"),
             legend.key.width = unit(5, "cm"),
             #legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="white"),
             panel.spacing.y = unit(2, "mm"),
             #panel.margin = unit(0.3, "cm"),
             #plot.margin = margin(2, 2, 2, 2, "cm"),
             axis.title.x = element_text(size = 40,face="bold", margin = margin(1, 20, 20, 20), vjust = -3),
             axis.title.y = element_text(size = 40,face="bold",  vjust = -1),
             axis.text.y=element_text(size=rel(1.5), colour = "black",face="bold"),
             axis.text.x=element_text(size=rel(1.5), colour = "black", face="bold"), #angle = 90,
             axis.ticks.length = unit(0.3,"cm"),
             axis.ticks = element_line(linewidth=1, colour = "black"),
             plot.subtitle = element_text(size = 22, face = "bold", hjust = 0.5)
             #plot.caption = element_text(size = 22, face = "bold", hjust = 0.5, vjust = 1, angle = 90),
      ) + xlab('FHC121') + ylab('CBS133791') + scale_colour_brewer(palette='Set1')
    
    pdf("~/mummerPlotX.pdf", 30, 30)
    print(p)
    dev.off()
    
  }, height = 900 )
  
}

shinyApp(ui = ui, server = server)








