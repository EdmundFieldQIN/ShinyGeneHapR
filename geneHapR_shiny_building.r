library(shiny)
library(shinydashboard)
library(colourpicker)
library(DT)
library(geneHapR)
library(GenomicRanges)
options(shiny.maxRequestSize = 100000 * 1024^2)
# windows下加载
# source("D:/R/R_Workspace/workspace/ShinyAPP/haplotype_shinyapp/www/vis_LDheatmap.R")
# source("D:/R/R_Workspace/workspace/ShinyAPP/haplotype_shinyapp/www/HaploWorldMap.r")
# source("D:/R/R_Workspace/workspace/ShinyAPP/haplotype_shinyapp/www/zzz.R")
# source("D:/R/R_Workspace/workspace/ShinyAPP/haplotype_shinyapp/www/plotHapTable.r")
# linux下加载
# source("www/vis_LDheatmap.R")
# source("www/HaploWorldMap.r")
# source("www/zzz.R")
# source("www/hap2vcf.R")
# source("www/plotHapTable.r")
# tbtools插件加载
source(file.path(find.package("geneHapR"), "extdata", "www","vis_LDheatmap.R"))
source(file.path(find.package("geneHapR"), "extdata", "www","HaploWorldMap.r"))
# source(file.path(find.package("geneHapR"), "extdata", "www","zzz.R"))
# source(file.path(find.package("geneHapR"), "extdata", "www", "plotHapTable.r"))
# ========================
# 工具函数(获取bcftools)
# ========================
get_tools <- function() {
  sys <- Sys.info()[["sysname"]]
  
  if (grepl("Windows", sys, ignore.case = TRUE)) {
    base <- R.home()
    bcftools <- file.path(base, "extract_software", "cygwin64", "usr", "local", "bin", "bcftools.exe")
    tabix    <- file.path(base, "extract_software", "cygwin64", "usr", "local", "bin", "tabix.exe")
    ## 把可执行文件打包到一个包下,tbtools插件查找
    # base <- find.package("geneHapR")
    # bcftools <- file.path(base, "extdata","extract_software", "cygwin64", "usr", "local", "bin", "bcftools.exe")
    # tabix    <- file.path(base, "extdata","extract_software", "cygwin64", "usr", "local", "bin", "tabix.exe")
    if (!file.exists(bcftools)) stop("bcftools 未找到")
    if (!file.exists(tabix))    stop("tabix 未找到")
    return(list(bcftools = bcftools, tabix = tabix, is_windows = TRUE))
  } else {
    return(list(bcftools = "bcftools", tabix = "tabix", is_windows = FALSE))
  }
}
to_cygwin_path <- function(path) {
  path <- normalizePath(path, winslash = "/")
  path <- sub("^([A-Za-z]):", "/cygdrive/\\L\\1", path, perl = TRUE)
  return(path)
}
# ========================
# 自动index（large_vcf用）
# ========================
ensure_index <- function(tools, vcf, threads = 1) {
  index_file <- paste0(vcf, ".tbi")
  if (file.exists(index_file)) return(list(skipped = TRUE, elapsed = 0))
  
  vcf_use <- vcf
  if (tools$is_windows) vcf_use <- to_cygwin_path(vcf)
  
  t0 <- proc.time()[["elapsed"]]
  system2(
    tools$bcftools,
    args = c("index", "--tbi", "--threads", threads, vcf_use),
    stdout = TRUE, stderr = TRUE
  )
  elapsed <- round(proc.time()[["elapsed"]] - t0, 1)
  
  if (!file.exists(index_file)) stop("index 创建失败")
  return(list(skipped = FALSE, elapsed = elapsed))
}
# ========================
# 读取contig（largevcf用）
# ========================
get_contigs <- function(tools, vcf) {
  vcf_use <- vcf
  if (tools$is_windows) vcf_use <- to_cygwin_path(vcf)
  
  res <- tryCatch({
    system2(tools$bcftools, c("view", "-h", vcf_use), stdout = TRUE)
  }, error = function(e) NULL)
  
  if (is.null(res)) return(NULL)
  contigs <- grep("^##contig=", res, value = TRUE)
  if (length(contigs) == 0) return(NULL)
  sub('.*ID=([^,>]+).*', '\\1', contigs)
}
# ==================================
# 判断字符串x是否以某个后缀surfix结尾。
# ==================================
is.surfix <- function(x, surfix) {
  x <- tolower(x)
  surfix <- tolower(surfix)
  return(any(endsWith(x, surfix)))
}
ui <- dashboardPage(
  # 定义home部分的文字，也是标签页名字
  dashboardHeader(title = "Shiny GeneHapR"),
  # 定义左侧导航栏的内容，基本代码结构:
  dashboardSidebar(
    # sidebarSearchForm(textId="searchText",buttonId="searchButton",label="Search",icon = shiny::icon("search")),
    sidebarMenu(
      menuItem("Main Page",tabName = "Main_Page",icon=icon("file-import",style ="color: #F2406D;" ),badgeLabel="Welcome",badgeColor="green"),
      menuItem("Import Data",tabName="Import",icon=icon("file-import",style ="color: #F2406D;" ),badgeLabel="Step One",badgeColor="maroon"),
      menuItem("Filter Genotypes",tabName="Filter",icon=icon("filter",style ="color: #F2406D;" ),badgeLabel="Step Two",badgeColor="maroon"),
      menuItem("Haplotype",tabName="Haplotype",icon=icon("table-cells",style ="color: #F2406D;" ),badgeLabel="Step Three",badgeColor="maroon"),
      menuItem("Visualize",tabName="Visualize",icon=icon("images",style ="color: #F2406D;"),badgeLabel="Step Final",badgeColor="maroon")
    )
  ),
  # 皮肤颜色
  skin = "purple",
  
  # 控件
  dashboardBody(
    tabItems(
      #### 首页Page ####
      tabItem(tabName = "Main_Page",
              fluidRow(
                # 介绍本人单位信息，以及这个geneHapR包的Github地址
                column(width = 4,
                box(title = "Welcome to use geneHapR Shiny App",status = "primary", solidHeader = TRUE,width = 12,
                    p("让你某些单倍型分析及可视化更加顺畅!:D",style = "font-family: 'times'; font-size:12pt;color:grey"),
                    br(),
                    strong("欢迎各位使用，请不要将本工具用于商业用途，谢谢！", style = "font-family: 'times'; font-size:16pt"),
                    br(),
                    strong("🙇🙇🙇", style = "font-family: 'times'; font-size:16pt"),
                    hr(),
                    p("本工具使用Shiny开发，单倍型数据分析与可视化功能使用geneHapR包实现。geneHapR包提供了一系列函数用于单倍型数据的处理和可视化。具体请访问：", a("Github Code.", href="https://gitee.com/zhangrenl/genehapr",style = "font-family: 'times'; font-size:14pt",target = "_blank"),style = "font-family: 'times'; font-size:14pt"),
                    br(),
                    p("TBtools插件版本基于R4.2.1开发，部分包功能不及新版本完善。如果有能力单独运行R脚本，可以下载本脚本新版并在本地运行：", a("Github Code.", href="https://github.com/EdmundFieldQIN/ShinyGeneHapR",style = "font-family: 'times'; font-size:14pt",target = "_blank"),style = "font-family: 'times'; font-size:14pt"),
                    br(),
                    p("有任何问题或建议请反馈给：", a("Github Issues.", href="https://github.com/EdmundFieldQIN/ShinyGeneHapR/issues",style = "font-family: 'times'; font-size:14pt",target = "_blank"), style = "font-family: 'times'; font-size:14pt"),
                    hr(),
                    hr(),
                    p(em("Email: 871729982@qq.com",style = "font-family: 'times'; font-size:12pt;color:grey")),
                    br(),
                    p(em("地址:华南农业大学15号楼岭南现代农业科学与技术广东省实验室",style = "font-family: 'times'; font-size:12pt;color:grey")),
                    br(),
                    p(em("水稻基因挖掘与利用研究团队",style = "font-family: 'times'; font-size:12pt;color:grey"))
                )
                ),
                # 介绍这个包的功能和使用方法
                column(width = 8,
                    box(title = "Introduction of geneHapR Shiny App",status = "primary", solidHeader = TRUE,width = 12,
                        fluidRow(
                            box(title = strong("Welcome to geneHapR Shiny App ！",style = "font-family: 'times'; font-size:16pt"),status = "danger", solidHeader = F,width = 12,
                                p(span("功能介绍：",style="font-weight:bold;color:#F2406D;font-family: 'SimSun'"),
                                    "geneHapR Shiny App是一个基于Shiny框架开发的单倍型数据分析和可视化工具。\n\t\t它提供了从数据导入、过滤、单倍型识别到结果可视化的一站式解决方案，适用于不同格式的基因型数据和注释文件。用户可以通过交互式界面轻松完成单倍型分析，并生成高质量的图表展示结果。",style = "font-family: 'times'; font-size:14pt"),
                                # 分点介绍这个工具的使用方法，风格统一style = "font-family: 'times'; font-size:14pt"
                                p(span("数据导入：",style="font-weight:bold;color:#F2406D;font-family: 'SimSun'"),"支持VCF、Plink、Fasta和Geno格式的基因型数据，以及GFF和BED格式的注释文件。用户可以根据自己的数据格式选择相应的导入方式，并在界面上查看导入结果。",style = "font-family: 'times'; font-size:14pt"),
                                p(span("数据过滤：",style="font-weight:bold;color:#F2406D;font-family: 'SimSun'"),"提供基于位置、类型或两者结合的过滤选项，用户可以根据需要筛选特定区域或类型的变异进行单倍型分析。",style = "font-family: 'times'; font-size:14pt"),
                                p(span("单倍型分析：",style="font-weight:bold;color:#F2406D;font-family: 'SimSun'"),"用户可以设置单倍型命名前缀，选择是否移除杂合基因型和未知基因型的个体，并运行单倍型识别和结果总结功能。分析结果可以在界面上查看，并支持下载。",style = "font-family: 'times'; font-size:14pt"),
                            )
                        ),
                        fluidRow(
                            box(title = strong("Visualization of geneHapR Shiny App ！",style = "font-family: 'times'; font-size:16pt"),status = "danger", solidHeader = F,width = 12,
                                # 分点介绍所有的可视化功能，风格统一style = "font-family: 'times'; font-size:14pt"
                                p(span("可视化一：",style="font-weight:bold;color:#F2406D;font-family: 'SimSun'"),"单倍型表格展示",style = "font-family: 'times'; font-size:14pt"),
                                p(span("可视化二：",style="font-weight:bold;color:#F2406D;font-family: 'SimSun'"),"基因模型上的变异展示",style = "font-family: 'times'; font-size:14pt"),
                                p(span("可视化三：",style="font-weight:bold;color:#F2406D;font-family: 'SimSun'"),"单倍型Network",style = "font-family: 'times'; font-size:14pt"),
                                p(span("可视化四：",style="font-weight:bold;color:#F2406D;font-family: 'SimSun'"),"单倍型与表型结合分析",style = "font-family: 'times'; font-size:14pt"),
                                p(span("可视化五：",style="font-weight:bold;color:#F2406D;font-family: 'SimSun'"),"单倍型位点效应分析",style = "font-family: 'times'; font-size:14pt"),
                                p(span("可视化六：",style="font-weight:bold;color:#F2406D;font-family: 'SimSun'"),"单倍型LD Heatmap",style = "font-family: 'times'; font-size:14pt"),
                                p(span("可视化七：",style="font-weight:bold;color:#F2406D;font-family: 'SimSun'"),"单倍型地理分布图",style = "font-family: 'times'; font-size:14pt"),
                            )
                        ),
                        fluidRow(
                            box(title = strong("Usage of geneHapR Shiny App ！",style = "font-family: 'times'; font-size:16pt"),status = "danger", solidHeader = F,width = 12,
                            # 使用方法：按照左侧导航栏的步骤依次进行数据导入、过滤、单倍型分析和可视化。每个步骤都有详细的设置选项和结果展示区域。风格统一style = "font-family: 'times'; font-size:14pt"
                                p(span("使用方法：",style="font-weight:bold;color:#F2406D;font-family: 'SimSun'"),
                                    "按照左侧导航栏的步骤依次进行数据导入、过滤、单倍型分析和可视化。每个步骤都有详细的设置选项和结果展示区域。",style = "font-family: 'times'; font-size:14pt"),
                                )
                        ),
                        
                        hr(),
                        div(div(p("Copyright Wgroup 2026. All rights reserved. Designed by EdmundFieldQIN.")),style = "font-family: 'times'; font-size:14pt"),
                        hr()
                    )

                )
                
              )
      ),
      
      #### 导入数据page ####
      tabItem(tabName="Import",
              # Boxes need to be put in a row (or column)
              fluidRow(
                tabBox(title = p("All File Import",style = "font-weight:bold;color : #5e58ca"),id = "tabset1",selected = p("Geno Data",style = "font-weight:bold;color : #5e58ca"),width = 12,
                       ####---- 基因型数据读入 ----####
                       tabPanel(p("Geno Data",style = "font-weight:bold;color : #5e58ca"),
                                # 数据读取页面 #
                                fluidRow(box(title = "Genotype File Input",status="primary", solidHeader=TRUE,
                                             selectInput("genotype_format", "Select your genotypes file format:",
                                                         choices = c("VCF" = "vcf_format",
                                                                     "Plink" = "plink_format",
                                                                     "Fasta" = "fasta_format",
                                                                     "Geno" = "geno_format"
                                                         ),
                                                         selected = "VCF_Format"),
                                             conditionalPanel(condition = "input.genotype_format == 'vcf_format'",
                                                              fileInput("vcf_file", "Choose VCF File",
                                                                        accept = c(".vcf",".vcf.gz")),
                                                              actionButton("vcf_import", "Import VCF File",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;")
                                             ),
                                             conditionalPanel(condition = "input.genotype_format == 'plink_format'",
                                                              fileInput("plink_file1", "Choose Ped File",accept = ".ped"),
                                                              fileInput("plink_file2", "Choose Map File",accept = ".map"),
                                                              radioButtons('plink_file_sep', '数据分隔符号', 
                                                                           c(Tab = '\t', Comma = ',', Space = ' '), 
                                                                           selected = '\t', inline = TRUE),
                                                              actionButton("plink_import", "Import Plink Files",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;")
                                             ),
                                             conditionalPanel(condition = "input.genotype_format == 'fasta_format'",
                                                              fileInput("fasta_file", "Choose Fasta File",accept = c(".fa",".fasta")),
                                                              actionButton("fasta_import", "Import Fasta File",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;")
                                             ),
                                             conditionalPanel(condition = "input.genotype_format == 'geno_format'",
                                                              fileInput("geno_file", "Choose Geno File",accept = ".geno"),
                                                              checkboxInput('geno_file_header', '第一行作为表头', TRUE),
                                                              radioButtons('geno_file_sep', '数据分隔符号', 
                                                                           c(Tab = '\t', Comma = ',', Semicolon = ';'), 
                                                                           selected = ',', inline = TRUE),
                                                              actionButton("geno_import", "Import Geno File",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;")
                                             )
                                ),
                                
                                # 数据展示页面 #
                                box(title = "File Import Show",status="danger", solidHeader = TRUE,
                                    # 当genotype_format等于vcf_format或plink_format时显示下面的内容
                                    conditionalPanel(condition = "input.genotype_format == 'vcf_format' || input.genotype_format == 'plink_format'", 
                                                     div(verbatimTextOutput("vcf_plink_import_info")
                                                     )
                                    ),
                                    # 当genotype_format等于fasta_format或Geno_format时显示下面的内容
                                    conditionalPanel(condition = "input.genotype_format == 'fasta_format' || input.genotype_format == 'geno_format'", 
                                                     div(style = "overflow-y: auto;overflow-x: auto;",
                                                         DT::DTOutput("geno_fasta_import_info")
                                                     )
                                    )
                                )
                                )
                       ), # tabPanel "Geno Data"
                       
                       ####---- 注释导入 ----####
                       tabPanel("Annotation",
                                fluidRow(box(title = "Annotation File Input",status="primary", solidHeader=TRUE,
                                             selectInput("annotation_format", "Select your annotation file format:",
                                                         choices = c("GFF" = "gff_format","BED" = "bed_format"),
                                                         selected = "gff_format"),
                                             conditionalPanel(condition = "input.annotation_format == 'gff_format'",
                                                              selectInput("gff_version", "Select GFF version:",
                                                                          choices = c("GFF3" = "gff3", "GFF2" = "gff2","GFF1" = "gff1","GFF" = "gff","GTF"= "gtf","GVF"= "gvf"),
                                                                          selected = "gff")
                                             ),
                                             
                                             
                                             fileInput("gff_file", "Choose Annotation File",
                                                       accept = c(".gff","gff1","gff2",".gff3",".gtf",".gvf",".bed",".bed6",".bed4")),
                                             
                                             checkboxInput("add_pro",label = "Add_Promoter",value = TRUE),
                                             numericInput("add_pro_length",label = "Length from ATG:(bp)",value = 2000),
                                             
                                             actionButton("gff_import", "Import Annotation File",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;")
                                ),
                                box(title = "Annotation Import Show",status="danger", solidHeader=TRUE,
                                    column(width = 12,verbatimTextOutput("annotation_import_info"))
                                )
                                )
                       ),
                       
                       
                       ####---- 表型信息导入 ----####
                       tabPanel("Phenotype",
                                fluidRow(
                                  box(title = "Phenotype File Input",status="primary", solidHeader=TRUE,
                                      fileInput("pheno_file", "Choose Phenotype File",
                                                accept = c(".csv",".tsv",".txt")),
                                      checkboxInput('pheno_file_header', '第一行作为表头', TRUE),
                                      radioButtons('pheno_file_sep', '数据分隔符号', 
                                                   c(Tab = '\t', Comma = ','), 
                                                   selected = ',', inline = TRUE),
                                      actionButton("pheno_import", "Import Phenotype File",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;")
                                  ),
                                  box(title = "Phenotype Import Show",status="danger", solidHeader=TRUE,
                                      div(style = "overflow-y: auto;overflow-x: auto;",
                                          DT::DTOutput("phenotype_import_info"))
                                  )
                                )
                       ),
                       
                       ####---- 基因组信息导入 ----#### 
                       tabPanel("Accession Info",
                                fluidRow(box(title = "Accession Info File Input",status="primary", solidHeader=TRUE,
                                             fileInput("accinfo_file", "Choose Accession Info File",
                                                       accept = c(".csv",".tsv",".txt")),
                                             checkboxInput('accinfo_file_header', '第一行作为表头', TRUE),
                                             radioButtons('accinfo_file_sep', '数据分隔符号', 
                                                          c(Tab = '\t', Comma = ','), 
                                                          selected = ',', inline = TRUE),
                                             actionButton("accinfo_import", "Import Accession Info File",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;")
                                ),
                                box(title = "Accession Info Import Show",status="danger", solidHeader=TRUE,
                                    div(style = "overflow-y: auto;overflow-x: auto;",
                                        DT::DTOutput("accession_import_info"))
                                )
                                )
                       )
                )    # tabBox
              ),    # fluidRow
              
              fluidRow(
                box(title = "All Import File Status Boxes",status="primary", solidHeader=TRUE,width = 12,
                    # 动态生成四个大类的状态卡片
                    uiOutput("statusBoxes")
                )
              )  # fluidRow
      ),
      
      
      #### 数据过滤page ####
      tabItem(tabName="Filter",
              fluidRow(
                tabBox(title = p("Filter Data",style = "font-weight:bold;color : #5e58ca"),id = "tabset2",selected = p("Filter Genotypes Data",style = "font-weight:bold;color : #5e58ca"),width = 12,
                       tabPanel(p("Filter Genotypes Data",style = "font-weight:bold;color : #5e58ca"),
                                fluidRow(
                                  ## 过滤参数 ##
                                  box(title = "Filter Setting",status="primary", solidHeader=TRUE,
                                      fluidRow(column(width = 6,
                                                      radioButtons("filter_mode", label = "Filter variants by",
                                                                   selected = "none",
                                                                   choices = c("Position" = "POS",
                                                                               "Type" = "type",
                                                                               "Both of above" = "both",
                                                                               "Do not filter" = "none"
                                                                   )
                                                      ),
                                                      selectInput("filter_type",label = "Type",
                                                                  choices = c("CDS"="CDS","UTR"="UTR"),
                                                                  selected = "CDS")
                                      ),
                                      column(width = 6,
                                             selectInput("filter_chr", label = "Chromosome",
                                                         choices = list()),
                                             numericInput("filter_start", label = "Start", value = 0),
                                             numericInput("filter_end", label = "End",value = 0)
                                      )
                                      ),
                                      actionButton("filter_geno","Start Filter",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;")
                                  ),
                                  
                                  ## 过滤后数据展示 ##
                                  box(title = "After Filter Show",status="danger", solidHeader=TRUE,
                                      div(style = "overflow-y: auto;overflow-x: auto;",
                                          verbatimTextOutput("after_filter_vcf_plink_show"),
                                          DTOutput("after_filter_geno_show")
                                      )
                                  )
                                )
                       ),
                       tabPanel("Extract Large VCF",
                                fluidRow(box(title = "Extract Large VCF Setting",status="primary", solidHeader=TRUE,
                                             fileInput("large_vcf", "input VCF.gz"), 
                                             numericInput("large_vcf_threads", "threads", value = 4, max = 20,min = 1, step = 1),
                                             selectInput("large_vcf_chr", "Chromosome", choices = NULL),
                                             numericInput("large_vcf_start", "Start", 1),
                                             numericInput("large_vcf_end", "End", 10000),
                                             checkboxInput("large_vcf_header", "keep vcf header", TRUE),
                                             fluidRow(
                                               column(width = 6,
                                                      actionButton("large_vcf_run", "Run Filter",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;")),
                                               column(width = 6,actionButton("large_vcf_clean", "清理临时文件"))
                                             )
                                ),# box
                                box(title = "After Filter Large Genotypes Show",status="danger", solidHeader=TRUE,
                                    div(style = "overflow-y: auto;overflow-x: auto;",
                                        verbatimTextOutput("large_vcf_log"),
                                        downloadButton("large_vcf_download", "下载结果",style = "color: white; background-color: #1abc9c; border-color: #16a085;")
                                    )
                                )
                                )
                       ) # tabPanel
                )
              )
      ),
      
      #### 单倍型分析page ####
      tabItem(tabName="Haplotype",
              fluidRow(
                tabBox(title = p("Haplotype Identification",style = "font-weight:bold;color : #5e58ca"),id = "tabset3",selected = p("Identification",style = "font-weight:bold;color : #5e58ca"),width = 12,
                       tabPanel(p("Identification",style = "font-weight:bold;color : #5e58ca"),
                                fluidRow(
                                  column(width = 12,
                                         box(title = "Identification setting",status="primary", solidHeader=TRUE,collapsible = T,width = 12,
                                             textInput("hap_prefix","Prefix of haplotype names",value = "Hap"),
                                             checkboxInput("remove_hetero_ornot","Remove heterozygote individuals",value = TRUE),
                                             checkboxInput("remove_unknown_ornot","Remove genotype unknown individuals",value = TRUE),
                                             fluidRow(
                                               column(width = 6,
                                                      actionButton("haplo_result_run", "Haplotype Find",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;")),
                                               column(width = 6,actionButton("haplo_summary_run", "Summary Result",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;"))
                                             )
                                         )),
                                ),
                                fluidRow(
                                  column(width = 6,
                                         box(title = "Haplotype result",status="danger", solidHeader=TRUE,width = 12,
                                             div(style = "overflow-y: auto;overflow-x: auto;",
                                                 DT::DTOutput("haplo_result_show"),
                                                 downloadButton("haplo_result_download","Download Result",style = "color: white; background-color: #1abc9c; border-color: #16a085;")
                                             )
                                         )),
                                  column(width = 6,
                                         box(title = "Summary of haplotype result",status="danger", solidHeader=TRUE,width = 12,
                                             div(style = "overflow-y: auto;overflow-x: auto;",
                                                 DT::DTOutput("haplo_summary_show"),
                                                 downloadButton("haplo_summary_download","Download Summary",style = "color: white; background-color: #1abc9c; border-color: #16a085;")
                                             )
                                         ))
                                )
                                
                       )
                )
              )
      ),
      
      #### 可视化page ####
      tabItem(tabName="Visualize",
              fluidRow(
                tabBox(title = p("Visualize Input",style = "font-weight:bold;color : #5e58ca"),id = "tabset4",selected = p("Active Data",style = "font-weight:bold;color : #5e58ca"),width = 12,
                       tabPanel(p("Active Data",style = "font-weight:bold;color : #5e58ca"),
                                fluidRow(
                                  column(width = 3,
                                         box(title = "Data Setting",status = "primary",solidHeader = T,width = 12,
                                             checkboxInput("import_own_data",em("Want to import HapResult you have ?",style = "font-weight:bold;font-size:12pt"),value = F),
                                             conditionalPanel(condition = "input.import_own_data == true",
                                                              fileInput("own_hapres",label = "Input your HapResult"),
                                                              radioButtons("own_hapres_sep","HapResult File Separated by :",c(Tab = "\t",Comma= ","),selected = ',', inline = TRUE),
                                                              actionButton("own_hapres_import",label = "Import your HapResult",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;")
                                             ),
                                             hr(),
                                             checkboxInput("filter_hap_ornot",em("Want to filter HapResult?",style = "font-weight:bold;font-size:12pt"),value = F),
                                             conditionalPanel(condition = "input.filter_hap_ornot == true",
                                                              selectInput("rm_mode",label = "Choose remove mode:",
                                                                          choices = c("Position" = "position",
                                                                                      "Accession" = "accession",
                                                                                      "Haplotype" = "haplotype",
                                                                                      "Frequent" = "freq"),
                                                                          selected = "freq"),
                                                              uiOutput("rm_condition"),
                                                              actionButton("filter_hap","Filter HapResult",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;")
                                             ),
                                             hr(),
                                             checkboxInput("setatgas0_ornot",em("Set HapResult Sites ATG as 0 ?",style = "font-weight:bold;font-size:12pt"),value = F),
                                             conditionalPanel(condition = "input.setatgas0_ornot == true",
                                                              p("Only use in HapTable and GeneModel",style = "font-weight:bold;font-size:12pt ;color:red"),
                                                              textInput("atgas0_geneid","Input Gene Id:",placeholder = "Must be included in Annotation file"),
                                                              htmlOutput("atgas0_geneid_info"),
                                                              selectInput("atgas0_chr", label = "Chromosome",choices = list()),
                                                              numericInput("atgas0_start", label = "Start", value = 0),
                                                              numericInput("atgas0_end", label = "End",value = 0),
                                                              actionButton("set_atgas0","Set ATG as 0",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;")
                                             ),
                                         )
                                  ),
                                  column(width = 9,
                                         box(title = "Current Active HapResult",status="danger", solidHeader=TRUE,width = 12,
                                             verbatimTextOutput("active_hapresult_info"),
                                             DT::DTOutput("current_HapResult"),
                                             downloadButton("current_hapresult_download",style = "color: white; background-color: #1abc9c; border-color: #16a085;")
                                         ),
                                         box(title = "Current Active HapSummary",status="danger", solidHeader=TRUE,width = 12,
                                             verbatimTextOutput("active_hapsummary_info"),
                                             DT::DTOutput("current_hapsummary"),
                                             downloadButton("current_hapsummary_download",style = "color: white; background-color: #1abc9c; border-color: #16a085;")
                                         ),
                                  )
                                )        
                       ),
                       tabPanel("Haplo Table",
                                fluidRow(
                                  column(width = 8,
                                         box(title = "Haplotype Table Plot",status = "primary",solidHeader = T,width = 12,
                                             div(style = "overflow-y: auto;overflow-x: auto;",
                                                 selectInput("hap_table_input","Select import data type:",choices = c("Origin"="ori"),selected = "ori"),
                                                 plotOutput("hap_table_plotout",height = "80%" ) %>% shinycssloaders::withSpinner()#添加加载动画
                                             )
                                         )),
                                  column(width = 4,box(title = "HaploTable Setting",status = "danger",solidHeader = T,width = 12,
                                                       fluidRow(
                                                         column(width = 6,
                                                                textInput("hap_table_prefix","Prefix of Haplotype Names",value = "Hap"),
                                                                textInput("hap_table_genename", "Gene name", value = ""),
                                                                colourInput("hap_table_bg_color",label = "Cell background color",value = "grey90",allowTransparent = T),
                                                                sliderInput("hap_table_angle", label = "Angle of x-axis label",
                                                                            min = 0, max = 90, step = 45, value = 45)),
                                                         
                                                         
                                                         column(width = 6,
                                                                selectInput("hap_table_infotag", "Tag name in INFO", choices = list()),
                                                                textInput("hap_table_tagname", "Tag name in IMG", value = ""),
                                                                textInput("hap_table_tagsplit", "Tag split", value = "|"),
                                                                numericInput("hap_table_tagfield", "Tag field", value = "1"))
                                                       ),
                                                       fluidRow(
                                                         column(width = 12,
                                                                actionButton("hap_table_plot","Start plot HapTable",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;"),
                                                                htmlOutput("hap_table_plot_info"))
                                                       ),
                                                       fluidRow(
                                                         column(width = 6,sliderInput("hap_table_plotout_width",label = "Width",min = 240,max = 1920,value = 960,step = 10),numericInput("hap_table_down_width","Out Width",min = 1,max = 20,value = 9,step = 1)),
                                                         column(width = 6,sliderInput("hap_table_plotout_height",label = "Height",min = 240,max = 1920,value = 720,step = 10),numericInput("hap_table_down_height","Out Height",min = 1,max = 20,value = 7,step = 1))
                                                       ),
                                                       radioButtons("hap_table_out_type","Plot Output Type",choices = c( 'PDF' = 'pdf', "PNG" = 'png','JPEG' = 'jpeg'),inline = T),
                                                       downloadButton("hap_table_download","Download Plot",style = "color: white; background-color: #1abc9c; border-color: #16a085;")
                                  ))
                                )
                       ),
                       tabPanel("Variants on Gene Model",
                                fluidRow(
                                  column(width = 8,
                                         box(title = "Gene Model Show",status = "primary",solidHeader = T,width = 12,
                                             div(style = "overflow-y: auto;overflow-x: auto;",
                                                 selectInput("genemodel_input","Select import data type:",choices = c("Origin"="ori"),selected = "ori"),
                                                 plotOutput("genemodel_plotout",height = "80%" )
                                             )
                                         )),
                                  column(width = 4,
                                         box(title = "Gene Model Setting",status = "danger",solidHeader = T,width = 12,collapsible = T,
                                             fluidRow(
                                               column(width = 6,selectInput("genemodel_chr","Chromosome",choices = list()),
                                                      numericInput("genemodel_start",label = "Length from model start",value = -200,max = 0)),
                                               column(width = 6,selectInput("genemodel_type", label = "Variants type",
                                                                            choices = list("circle" = "circle", "pie" = "pie", "pin" = "pin", "flag" = "flag"),
                                                                            selected = "pin"),
                                                      numericInput("genemodel_width",label = "Model length:",value = 2000,min = 100)
                                               )
                                             ),
                                             sliderInput("genemodel_cex",label = "Size of label", min = 0, max = 3, value = 0.8, step = 0.05),
                                             sliderInput("genemodel_cdsh",label = "Size of CDS Height",min = 0.02,max = 0.2,value = 0.06,step = 0.01),
                                             sliderInput("genemodel_utrvscds",label = "The ratio of UTR/CDS",min = 0.4,max = 1.6,value = 0.66),
                                             checkboxGroupInput("geneElement","Gene Element",choiceNames = NULL,choiceValues = NULL,selected = NULL),
                                             fluidRow(
                                               column(width = 12,
                                                      actionButton("genemodel_plot","Start plot GeneMo",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;"),
                                                      htmlOutput("genemodel_plot_info"))
                                             ),
                                             fluidRow(
                                               column(width = 6,sliderInput("genemodel_plotout_width",label = "Width",min = 240,max = 1920,value = 960,step = 10),numericInput("genemodel_down_width","Out Width",min = 1,max = 20,value = 9,step = 1)),
                                               column(width = 6,sliderInput("genemodel_plotout_height",label = "Height",min = 240,max = 1920,value = 720,step = 10),numericInput("genemodel_down_height","Out Height",min = 1,max = 20,value = 7,step = 1))
                                             ),
                                             radioButtons("genemodel_out_type","Plot Output Type",choices = c( 'PDF' = 'pdf', "PNG" = 'png','JPEG' = 'jpeg'),inline = T),
                                             downloadButton("genemodel_download","Download Plot",style = "color: white; background-color: #1abc9c; border-color: #16a085;")
                                         ))
                                )
                       ),
                       tabPanel("Haplo Network",
                                fluidRow(
                                  column(width = 8,
                                         box(title = "Haplotype Network Show",status = "primary",solidHeader = T,width = 12,
                                             div(style = "overflow-y: auto;overflow-x: auto;",
                                                 plotOutput("hapnet_plotout",height = "80%" )
                                             )
                                         )
                                  ),
                                  column(width = 4,
                                         box(title = "Haplotype Network Setting",status = "danger",solidHeader = T,width = 12,collapsible = T,
                                             fluidRow(
                                               column(width = 6,
                                                      selectInput("plothapnet_group", "Group by", choices = list("none"="none")),
                                                      numericInput("plothapnet_collink", "Link color", value = 1),
                                                      selectInput("plothapnet_labels","Display haplotype names",choices = c("TRUE"="TRUE","FALSE"="FALSE"),selected = "TRUE"),
                                                      sliderInput("plothapnet_xlim", "X limit range",min = -200, max = 200, value =c(-1,1)),
                                                      sliderInput("plothapnet_pie.lim", "pie size range",min = 0, max = 1, value =c(0.1,1), step = 0.01),
                                                      actionButton("hapnet_plot","Plot Haplo Network",width = "95%",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;")
                                               ),
                                               column(width = 6,
                                                      selectInput("plothapnet_showmutation", label = "Mutation symbol",selected = "line",choices = list("none" = 0, "line" = 1,"dot" = 2, "number of mutantions" = 3)),
                                                      numericInput("plothapnet_linkwidth", "Line width", value = 1),
                                                      selectInput("plothapnet_scale", label = "Scale",selected = "log2",choices = list("none" = 1, "log10" = "log10", "log2" = "log2")),
                                                      sliderInput("plothapnet_ylim", "Y limit range",min = -200, max = 200, value =c(-1,1)),
                                                      sliderInput("plothapnet_pie.lim.folder", "2^x",min = 0, max = 4, value = 0, step = 1),
                                                      sliderInput("plothapnet_cex",label = "Size of label", min = 0, max = 3, value = 0.8, step = 0.05)
                                               )
                                             )
                                         ),
                                         box(title = "Legend and output Setting",status = "danger",solidHeader = T,width = 12,collapsible = T,
                                             fluidRow(
                                               column(width = 6,
                                                      checkboxInput("plothapnet_legend_size", "Show size legend", value = TRUE),
                                                      sliderInput("plothapnet_x", "x", min = -200, max = 200, value = -10),
                                                      # hapnet_out_type
                                                      radioButtons("hapnet_out_type","Plot Output Type",choices = c( 'PDF' = 'pdf', "PNG" = 'png','JPEG' = 'jpeg'),inline = T),
                                               ),
                                               column(width = 6,
                                                      checkboxInput("plothapnet_legend_color", "Show color legend", value = TRUE),
                                                      sliderInput("plothapnet_y", "y", min = -200, max = 200, value = -10),
                                                      sliderInput("plothapnet_cexlegend",label = "Size of legend", min = 0, max = 3, value = 0.8, step = 0.05)
                                               )
                                             ),downloadButton("hapnet_download","Download Plot",style = "color: white; background-color: #1abc9c; border-color: #16a085;")
                                         )
                                  )
                                )
                       ),
                       tabPanel("Haplo Phenotype",
                                fluidRow(
                                  column(width = 8,
                                         box(title = "Haplotype Vs Phenotype Show",status = "primary",solidHeader = T,width = 12,
                                             div(style = "overflow-y: auto;overflow-x: auto;",
                                                 plotOutput("hapvspheno_plotout",height = "80%"),
                                                 verbatimTextOutput("hapvspheno_testout")
                                             )
                                         )
                                  ),
                                  column(width = 4,
                                         box(title = "HapVsPheno Setting",status = "danger",solidHeader = T,width = 12,
                                             fluidRow(
                                               column(width = 6,
                                                      selectInput("hapvspheno_pheno_name","Pheno Name",choices = list("first pheno" = "1", "second pheno" = 2),selected = "2"),
                                                      textInput("hapvspheno_title","Title",value = " ")
                                               ),
                                               column(width = 6,
                                                      numericInput("hapvspheno_minacc","Minimum number of accession",min = 0,value = 5),
                                                      numericInput("hapvspheno_angle","Angle of x labels",min = 0,max = 89,value = 45)
                                               )
                                             ),
                                             fluidRow(
                                               column(width = 12,actionButton("hapvspheno_plot","Start plot GeneMo",width = "100%",icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;"))
                                             ),
                                             fluidRow(
                                               column(width = 6,sliderInput("hapvspheno_plotout_width",label = "Width",min = 240,max = 1920,value = 960,step = 10),numericInput("hapvspheno_down_width","Out Width",min = 1,max = 20,value = 10,step = 1)),
                                               column(width = 6,sliderInput("hapvspheno_plotout_height",label = "Height",min = 240,max = 1920,value = 480,step = 10),numericInput("hapvspheno_down_height","Out Height",min = 1,max = 20,value = 5,step = 1))
                                             ),
                                             radioButtons("hapvspheno_out_type","Plot Output Type",choices = c( 'PDF' = 'pdf', "PNG" = 'png','JPEG' = 'jpeg'),inline = T),
                                             downloadButton("hapvspheno_download","Download Plot",style = "color: white; background-color: #1abc9c; border-color: #16a085;")
                                         )
                                  )
                                )
                       ),
                       tabPanel("Haplo SiteEff",
                                fluidRow(
                                  column(width = 8,
                                         fluidRow(
                                           column(width = 12,
                                                  box(title = "Haplotype Effect",status = "primary",solidHeader = T,width = 12,
                                                      plotOutput("HapEFF_plotout",height = "80%")
                                                  )
                                           )
                                         ),
                                         fluidRow(
                                           column(width = 6,box(title = "Haplo Site Effect",status = "primary",solidHeader = T,width = 12,collapsible = T,
                                                                div(DT::DTOutput("hapeff_eff_output"),
                                                                downloadButton("hapeff_eff_output_download",style = "color: white; background-color: #1abc9c; border-color: #16a085;")))),
                                           column(width = 6,box(title = "Haplo Site Pvalue",status = "primary",solidHeader = T,width = 12,collapsible = T,
                                                                div(DT::DTOutput("hapeff_p_output"),
                                                                downloadButton("hapeff_p_output_download",style = "color: white; background-color: #1abc9c; border-color: #16a085;"))))
                                         )
                                  ),
                                  column(width = 4,
                                         box(title = "HapEff Setting",status = "danger",solidHeader = T,width = 12,
                                             fluidRow(column(width = 6,
                                                             selectInput("hapeff_chr","Chromesome",choices = list(),selected = NULL),
                                                             numericInput("hapeff_start","Start",value = 0),
                                                             sliderInput("hapeff_cex","Point Size",min = 0.1,max = 1.5,step = 0.05,value = 0.5),
                                                             sliderInput("hapeff_cdsh","CDS Height",min = 0.1,max = 1.5,value = 0.5,step = 0.1),
                                                             sliderInput("hapeff_line_color","Line Color",min = 1,max = 9,value = 1,step = 1)
                                                             
                                             ),
                                             column(width = 6,
                                                    selectInput("hapeff_y","Y axis is:",choices = c("Pvalue"="pvalue","Effect" = "effect"),selected = "effect"),
                                                    numericInput("hapeff_end","End",value = 0),
                                                    sliderInput("hapeff_legend_cex","Legend size",min = 0.1,max = 1.5,step = 0.1,value = 0.8),
                                                    sliderInput("hapeff_point_type","Point Type",min = 1,max = 25,value = 20,step = 1),
                                                    sliderInput("hapeff_line_type","Line Type",min = 1,max = 6,value = 1,step = 1)
                                             )                                           
                                             ),
                                             fluidRow(column(width = 12,
                                                             textInput("hapeff_title","Title",value = " "),
                                                             checkboxGroupInput("hapeff_show_type","showType",
                                                                                choiceNames = c("five_prime_UTR", "CDS", "three_prime_UTR"),
                                                                                choiceValues = c("five_prime_UTR", "CDS", "three_prime_UTR"),
                                                                                selected = c("five_prime_UTR", "CDS", "three_prime_UTR"),inline = T),
                                                             actionButton("hapeff_plot","Start draw Haplo Eff Plot",width = "100%",
                                                                          icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;")
                                             )),
                                             fluidRow(
                                               column(width = 6,sliderInput("hapeff_plotout_width",label = "Width",min = 240,max = 1920,value = 960,step = 10),numericInput("hapeff_down_width","Out Width",min = 1,max = 20,value = 10,step = 1)),
                                               column(width = 6,sliderInput("hapeff_plotout_height",label = "Height",min = 240,max = 1920,value = 480,step = 10),numericInput("hapeff_down_height","Out Height",min = 1,max = 20,value = 5,step = 1))
                                             ),
                                             radioButtons("hapeff_out_type","Plot Output Type",choices = c( 'PDF' = 'pdf', "PNG" = 'png','JPEG' = 'jpeg'),inline = T),
                                             downloadButton("hapeff_download","Download Plot",style = "color: white; background-color: #1abc9c; border-color: #16a085;")
                                         )
                                  )
                                )
                       ),
                       tabPanel("LD Heatmap",
                                fluidRow(
                                  column(width = 8,
                                         box(title = "LD Heatmap Show",status = "primary",solidHeader = T,width = 12,
                                             div(style = "overflow-y: auto;overflow-x: auto;",
                                                 plotOutput("ldheatmap_plotout",height = "80%" )
                                             )
                                         )
                                  ),
                                  column(width = 4,
                                         box(title = "LD Heatmap Setting",status = "danger",solidHeader = T,width = 12,collapsible = T,
                                             fluidRow(
                                               column(width = 6,
                                                      selectInput("ldheatmap_chr","Chromosome",choices = list(),selected = NULL),
                                                      numericInput("ldheatmap_start","Start",value = 0),
                                                      colourInput("ldheatmap_snpcolor",label = "SNP Color",value = "black",allowTransparent = T)
                                               ),
                                               column(width = 6,
                                                      textInput("ldheatmap_title","Title",value = "Pairwise LD Heatmap"),
                                                      numericInput("ldheatmap_end","End",value = 0),
                                                      sliderInput("ldheatmap_cex",label = "Size of label", min = 0, max = 1, value = 0.8, step = 0.05)
                                               )
                                             ),
                                             fluidRow(
                                               column(width = 4,
                                                      colourInput("ldheatmap_lowcolor",label = "Low Color",value = "grey90",allowTransparent = T),
                                                      selectInput("ldheatmap_method", label = "LD method",choices = list("r2" = "r", "D'" = "D"),selected = "r")
                                               ),
                                               column(width = 4,
                                                      colourInput("ldheatmap_midcolor",label = "Mid Color",value = "#0084FF",allowTransparent = T),
                                                      selectInput("ldheatmap_distances",label = "Distances",choices = list("physical distance" = "physical", "genetic distance" = "genetic"),selected = "physical")
                                               ),
                                               column(width = 4,
                                                      colourInput("ldheatmap_highcolor",label = "High Color",value = "#F2406D",allowTransparent = T),
                                                      selectInput("ldheatmap_addmap", label = "Show GeneModel",choices = list("Yes" = "TRUE", "No" = "FALSE"),selected = "FALSE")
                                               )
                                             ),
                                             fluidRow(
                                               column(width = 6,
                                                      # textInput("ldheatmap_geneid", "Gene ID", value = ""),
                                                      selectInput("ldheatmap_text", "Add LD Value", choices = list("yes" = "TRUE", "no" = "FALSE"), selected = "FALSE"),
                                                      sliderInput("ldheatmap_maploc",label = "Gene Model Location",min = 0.05, max = 0.5, value = 0.15, step = 0.05)
                                               ),
                                               column(width = 6,
                                                      colourInput("ldheatmap_genemodelcolor",label = "Gene Model Color",value = "#F2406D",allowTransparent = T),
                                                      sliderInput("ldheatmap_mapheight",label = "Gene Model Height",min = 0.01, max = 0.2, value = 0.01, step = 0.01)
                                               )
                                             ),
                                             fluidRow(
                                               column(width = 12,
                                                      actionButton("ldheatmap_plot","Start plot LD Heatmap",width = "100%",
                                                                   icon = icon("play"),style = "color: white; background-color: #3498db; border-color: #2980b9;")
                                               ))
                                         ),
                                         box(title = "OutPut setting",status = "danger",solidHeader = T,width = 12,collapsible = T,
                                             fluidRow(
                                               column(width = 6,sliderInput("ldheatmap_plotout_width",label = "Width",min = 240,max = 1920,value = 960,step = 10),numericInput("ldheatmap_down_width","Out Width",min = 1,max = 20,value = 10,step = 1)),
                                               column(width = 6,sliderInput("ldheatmap_plotout_height",label = "Height",min = 240,max = 1920,value = 960,step = 10),numericInput("ldheatmap_down_height","Out Height",min = 1,max = 20,value = 10,step = 1))
                                             ),
                                             radioButtons("ldheatmap_out_type","Plot Output Type",choices = c( 'PDF' = 'pdf', "PNG" = 'png','JPEG' = 'jpeg'),inline = T),
                                             downloadButton("ldheatmap_download","Download Plot",style = "color: white; background-color: #1abc9c; border-color: #16a085;")
                                         )
                                  )
                                )
                       ),
                       tabPanel("GEO Distribution",
                                fluidRow(column(width = 8,
                                                box(title = "Geographical Distribution Show",status = "primary",solidHeader = T,width = 12,
                                                    div(style = "overflow-y: auto;overflow-x: auto;",
                                                        plotOutput("geodis_plotout",height = "80%" )
                                                    )
                                                ),
                                                box(title = "Geographical Map Setting",status = "danger",solidHeader = T,width = 12,collapsible = T,
                                                    fluidRow(
                                                      column(width = 3,colourInput("geodis_mapfill",label = "Regions fill color",value = "grey60",allowTransparent = T)),
                                                      column(width = 3,colourInput("geodis_mapborder",label = "Regions border color",value = "black",allowTransparent = T)),
                                                      column(width = 3,colourInput("geodis_mapbackground",label = "Map background color",value = "white",allowTransparent = T)),
                                                      column(width = 3,sliderInput("geodis_mapborderwidth",label = "Regions border width",min = 0, max = 1, value = 0.2, step = 0.05))
                                                    ),
                                                    fluidRow(
                                                        column(width = 3,
                                                            selectInput("geodis_exclude","Exclude Regions ?",choices = list("Yes" = "TRUE", "No" = "FALSE"),selected = "TRUE")
                                                        ),
                                                        column(width = 9,
                                                            selectInput("geodis_excluderegions", "Exclude Regions", choices = list("Antarctica" = "Antarctica"), selected = "Antarctica", multiple = TRUE)
                                                        )
                                                    ),
                                                    # 根据用户选择的haplotype数量动态生成颜色选择器
                                                    fluidRow(
                                                        column(width = 3,
                                                            selectInput("geodis_hapcolor","Customize Haplotype Colors ?",choices = list("Yes" = "TRUE", "No" = "FALSE"),selected = "FALSE")
                                                        ),
                                                        column(width = 9,
                                                            uiOutput("geodis_hapcolor_ui")
                                                        )
                                                    )
                                                )
                                        ),
                                        column(width = 4,
                                               box(title = "Legend and output Setting",status = "danger",solidHeader = T,width = 12,collapsible = T,
                                                   fluidRow(
                                                       column(width = 12,
                                                              selectInput("geodis_hapnames", "Select Haplotype Names (less than 3)", choices = list(), selected = NULL, multiple = TRUE)
                                                       )
                                                   ),
                                                   fluidRow(
                                                       column(width = 6,
                                                              selectInput("geodis_loncol", "Longitude column", choices = list(), selected = NULL),
                                                              numericInput("geodis_pielogbase", label = "Pie log base", value = 5),
                                                              sliderInput("geodis_pielinewidth", label = "Pie line width", min = 0, max = 1, value = 0, step = 0.05),
                                                              numericInput("geodis_pielegend_posix", label = "Pie legend x", value = -150,min =-180,max = 180,step = 10),
                                                              numericInput("geodis_pielegendtitle_posix", label = "Pie legend title x", value = -140,min =-180,max = 180,step = 10),
                                                              sliderInput("geodis_pielegendtitle_size", label = "Pie legend title size", min = 0.5, max = 10, value = 5, step = 0.5),
                                                              hr(),
                                                              actionButton("geodis_plot", "Start plot", width = "95%", icon = icon("play"), style = "color: white; background-color: #3498db; border-color: #2980b9;")
                                                       ),
                                                       column(width = 6,
                                                              selectInput("geodis_latcol", "Latitude column", choices = list(), selected = NULL),
                                                              colourInput("geodis_pielinecolor", label = "Pie line color", value = "white", allowTransparent = T),
                                                              sliderInput("geodis_pielegendsize", label = "Pie legend size", min = 0.5, max = 5, value = 3, step = 0.5),
                                                              numericInput("geodis_pielegend_posiy", label = "Pie legend y", value = -30,min =-90,max = 90,step = 10),
                                                              numericInput("geodis_pielegendtitle_posiy", label = "Pie legend title y", value = -20,min =-90,max = 90,step = 10),
                                                              sliderInput("geodis_haplegend_size", label = "Hap legend size", min = 0.5, max = 10, value = 10, step = 0.5),
                                                              selectInput("geodis_haplegend_posi", "Hap legend position", choices = list("right" = "right", "left" = "left", "top" = "top", "bottom" = "bottom"), selected = "bottom"),
                                                       )
                                                   ),
                                                   # 图片下载面板
                                                   fluidRow(
                                                         column(width = 6,
                                                                  sliderInput("geodis_plotout_width",label = "Width",min = 240,max = 1920,value = 960,step = 10),
                                                                  numericInput("geodis_down_width","Out Width",min = 1,max = 20,value = 10,step = 1),
                                                                  radioButtons("geodis_out_type","Plot Output Type",choices = c( 'PDF' = 'pdf', "PNG" = 'png','JPEG' = 'jpeg'),inline = T)
                                                         ),
                                                         column(width = 6,
                                                                  sliderInput("geodis_plotout_height",label = "Height",min = 240,max = 1920,value = 480,step = 10),
                                                                  numericInput("geodis_down_height","Out Height",min = 1,max = 20,value = 5,step = 1),
                                                                  hr(),
                                                                  downloadButton("geodis_download","Download Plot",style = "color: white; background-color: #1abc9c; border-color: #16a085;")
                                                         )
                                                   )
                                               )
                                        )
                                ),
                       )
                )
              )
      )  # tabItem 可视化page
    )    # tabItems
  )    # dashboardBody
)    # dashboardPage
server <- function(input, output, session) {
  
  ##########---------- 数据输入----------##########
  #### 重置话初始数据
  vcf <- seq <- p_link <- geno.gt  <- pheno <- gff <-
    accinfo <- hapResult <- hapSummary <- NULL
  # 追踪当前激活的文件类型（NULL 表示未导入）
  active_genotype_type <- reactiveVal(NULL)
  active_annotation_type <- reactiveVal(NULL)
  active_phenotype <- reactiveVal(FALSE)
  active_accession <- reactiveVal(FALSE)
  
  #### VCF文件输入
  vcf <- eventReactive(input$vcf_import, {
    req(input$vcf_file)
    ## 判断文件后缀是否为vcf或vcf.gz，若是则调用geneHapR::import_vcf函数进行读入，否则返回NULL
    vcf <- if(is.surfix(input$vcf_file$name, c("vcf","vcf.gz"))) {
      geneHapR::import_vcf(input$vcf_file$datapath)
    } else {return(NULL)}
    
    if(!is.null(vcf)){
      chrs <- as.list(unique(vcf@fix[,1]))
      names(chrs) <- as.vector(chrs)
      updateSelectInput(session, "filter_chr",
                        choices = chrs,
                        selected = chrs[[1]])
      start <- min(as.numeric(vcf@fix[,2]))
      end <- max(as.numeric(vcf@fix[,2]))
      updateNumericInput(session,"filter_start",
                         value = start,
                         min = start,
                         max = end)
      updateNumericInput(session,"filter_end",
                         value = end,
                         min = start,
                         max = end)
    }
    return(vcf)
  })
  
  observeEvent(input$vcf_import,{
    Sys.sleep(0.1)
    req(vcf())
    output$vcf_plink_import_info <- renderPrint({ vcf() })
    # 更新状态
    active_genotype_type("VCF")
  })
  
  
  #### plink文件输入
  plink <- eventReactive(input$plink_import, {
    req(input$plink_file1)
    req(input$plink_file2)
    p.link <- if(is.surfix(input$plink_file2$name, "map") & is.surfix(input$plink_file1$name, "ped")){
      geneHapR::import_plink.pedmap(pedfile = input$plink_file1$datapath,
                                    mapfile = input$plink_file2$datapath,
                                    sep_ped = input$plink_file_sep,
                                    sep_map = input$plink_file_sep) 
    } else {return(NULL)}
    
    if(!is.null(p.link)){
      chrs <- as.list(unique(p.link$map[,1]))
      names(chrs) <- as.vector(chrs)
      updateSelectInput(session, "filter_chr",
                        choices = chrs,
                        selected = chrs[[1]])
      start <- min(as.numeric(p.link$map[,4]))
      end <- max(as.numeric(p.link$map[,4]))
      updateNumericInput(session,"filter_start",
                         value = start,
                         min = start,
                         max = end)
      updateNumericInput(session,"filter_end",
                         value = end,
                         min = start,
                         max = end)
    }
    return(p.link)
  })
  
  observeEvent(input$plink_import, {
    Sys.sleep(0.1)
    req(plink())
    plink_object <- plink()
    output$vcf_plink_import_info <- renderPrint({
      cat("Indeividuals: ", nrow(plink_object$ped),"\n")
      cat("Markers: ", nrow(plink_object$map))
    })
    active_genotype_type("Plink")
  })
  
  
  #### fasta文件输入
  fasta <- eventReactive(input$fasta_import, {
    req(input$fasta_file)
    seqs <- if(is.surfix(input$fasta_file$name, c("fasta", "fa","FASTA","FA"))){
      geneHapR::import_seqs(input$fasta_file$datapath)
    } else {
      return(NULL)
    }
    
    if(!is.null(seqs)){
      withProgress(message = "序列比对中，请稍候...", value = 0, {
        setProgress(0.1, detail = "正在调用 muscle 进行多序列比对")
        seqs <- muscle::muscle(seqs)
        setProgress(0.9, detail = "格式转换中 (DNAStringSet)")
        seqs <- as(seqs, "DNAStringSet")
        setProgress(1,   detail = "完成")
      })
    }
    
    return(seqs)
  })
  
  
  observeEvent(input$fasta_import,{
    Sys.sleep(0.1)
    req(fasta())
    geno_fasta_import_info <- data.frame("seq_name" = names(fasta()),
                                         "seq_length" = Biostrings::nchar(fasta()))
    output$geno_fasta_import_info <- renderDT({geno_fasta_import_info},options = list(pageLength = 5))
    active_genotype_type("Fasta")
  })
  
  
  #### geno文件输入
  geno <- eventReactive(input$geno_import, {
    req(input$geno_file)
    if(is.surfix(input$geno_file$name, "geno")){
      # 读取该行及其后内容为新的数据框
      geno_df <- data.table::fread(input$geno_file$datapath,
                                   header = input$geno_file_header,
                                   sep = input$geno_file_sep,
                                   stringsAsFactors = FALSE,
                                   check.names = FALSE  # ?? 禁用 make.names，防止乱码列名
      )
      
    } else {
      geno_df <- NULL
    }
    
    
    if(!is.null(geno_df)){
      geno_df <- as.data.frame(geno_df)
      chrs <- as.list(unique(geno_df[,1]))
      names(chrs) <- as.vector(chrs)
      updateSelectInput(session, "filter_chr",
                        choices = chrs,
                        selected = chrs[[1]])
      start <- min(as.numeric(geno_df[,2]))
      end <- max(as.numeric(geno_df[,2]))
      updateNumericInput(session,"filter_start",
                         value = start,
                         min = start,
                         max = end)
      updateNumericInput(session,"filter_end",
                         value = end,
                         min = start,
                         max = end)
    }
    
    return(geno_df)
  })
  
  observeEvent(input$geno_import,{
    Sys.sleep(0.1)
    req(geno())
    geno_fasta_import_info <- DT::datatable(geno(), options = list(pageLength = 5))
    output$geno_fasta_import_info <- DT::renderDT({geno_fasta_import_info})
    active_genotype_type("Geno")
  })
  
  
  #### 注释文件输入
  gff <- eventReactive(input$gff_import,{
    
    req(input$gff_file)
    
    gff <- if(is.surfix(input$gff_file$name, c("gff","gff1","gff2","gff3","gtf","gvf","bed","bed6","bed4"))){
      switch(input$annotation_format,
             "gff_format" = geneHapR::import_gff(input$gff_file$datapath, format = input$gff_version),
             "bed_format" = geneHapR::import_bed(input$gff_file$datapath))
    } else {return(NULL)}
    
    if (input$add_pro) {gff <- addPromoter(gff,PromoterLength = input$add_pro_length)}
    
    if(!is.null(gff)){
      
      type_choices <- as.list(as.vector(unique(gff$type)))
      names(type_choices) <- as.vector(type_choices)
      updateSelectInput(session,
                        "filter_type",
                        choices = type_choices,
                        selected = type_choices[[1]]
      )
      updateCheckboxGroupInput(session,
                               "geneElement", inline = TRUE,
                               choiceNames = as.vector(unique(gff$type)),
                               choiceValues = as.vector(unique(gff$type)),
                               selected = as.vector(unique(gff$type)[2])
      )
      # choices <- as.vector(unique(gff$type))
    }
    
    return(gff)
  })
  
  # import annotation file
  observeEvent(input$gff_import,{
    Sys.sleep(0.1)
    output$annotation_import_info <- renderPrint(gff())
    
    # 更新激活类型（GFF 细化到子格式，BED 直接显示 BED）
    if (!is.null(gff())) {
      label <- if (input$annotation_format == "gff_format") {
        toupper(input$gff_version)   # "GFF3" / "GTF"
      } else {
        "BED"
      }
      active_annotation_type(label)
    }
  })
  
  observeEvent(input$atgas0_geneid,{
    req(gff())
    req(input$atgas0_geneid)
    # 筛选出 type == "gene" 的所有特征
    genes_gr <- gff()[gff()$type == "gene"]
    # 根据 ID 匹配目标基因
    target_gr <- genes_gr[genes_gr$ID == input$atgas0_geneid]
    # 检查是否找到
    if (length(target_gr) == 0) {
      output$atgas0_geneid_info <- renderPrint(p("This gene id is not in your gff file",style = "font-weight:bold;font-size:12pt ;color:red"))
    }
    else {
      output$atgas0_geneid_info <- renderPrint(p("Detected gene id in your gff file",style = "font-weight:bold;font-size:12pt ;color:green"))
      # 提取信息
      atgas0_chr <- as.list(as.vector(as.character(target_gr@seqnames)))
      names(atgas0_chr) <- as.vector(atgas0_chr)
      updateSelectInput(session,"atgas0_chr",choices = atgas0_chr,selected = atgas0_chr[[1]])
      updateNumericInput(session,"atgas0_start",value = as.numeric(target_gr@ranges@start))
      updateNumericInput(session,"atgas0_end",value = as.numeric(target_gr@ranges@start) + as.numeric(target_gr@ranges@width)-1)
      # 刷新genemodel的长宽
      updateNumericInput(session,"genemodel_width",value = as.numeric(target_gr@ranges@width)+200)
      updateNumericInput(session,"genemodel_start",value = -200)
    }
  })
  #### 表型文件输入
  phenotypes <- eventReactive(input$pheno_import,{
    req(input$pheno_file)
    pheno_df <- if(is.surfix(input$pheno_file$name, c("txt","csv","tsv"))){
      import_AccINFO(input$pheno_file$datapath, 
                     sep = input$pheno_file_sep,                      # 分隔符号，默认为制表符"\t"
                     na.strings = "NA")
    } else {return(NULL)}
    return(pheno_df)
    
  })
  
  observeEvent(phenotypes(),{
    Sys.sleep(0.1)
    req(phenotypes())
    phenotype_import_info <- datatable(phenotypes(), options = list(pageLength = 5))
    output$phenotype_import_info <- renderDT({phenotype_import_info})
    
    choices <- as.list(names(phenotypes()))
    names(choices) <- names(phenotypes())
    updateSelectInput(session, "hapvspheno_pheno_name", choices = choices)
    active_phenotype(TRUE)
  })
  
  
  
  #### 基因组信息输入
  accession <- eventReactive(input$accinfo_import,{
    req(input$accinfo_file)
    accinfo_df <- if(any(is.surfix(input$accinfo_file$name, c("txt","csv","tsv")))){
      import_AccINFO(input$accinfo_file$datapath, 
                     sep = input$accinfo_file_sep,                      # 分隔符号，默认为制表符"\t"
                     na.strings = "NA")
      
    } else {return(NULL)}
    return(accinfo_df)
  })
  
  
  
  
  observeEvent(input$accinfo_import,{
    Sys.sleep(0.1)
    req(accession())
    
    accession_import_info <- datatable(as.data.frame(accession()), options = list(pageLength = 5))
    output$accession_import_info <- DT::renderDT({accession_import_info})
    
    choices <- as.list(names(accession()))
    choices <- c(choices, "none" = "none")
    names(choices) <- as.vector(choices)
    updateSelectInput(session, "plothapnet_group",
                      choices = choices, selected = "none")
    updateSelectInput(session, "geodis_loncol",
                      choices = choices, selected = "none")
    updateSelectInput(session, "geodis_latcol",
                      choices = choices, selected = "none")
    active_accession(TRUE)
    
  })
  
  ##########---------- 大型基因组数据过滤 ----------##########
  tools <- get_tools()
  
  log_val <- reactiveVal("")
  result_file <- reactiveVal(NULL)
  
  append_log <- function(txt) {
    log_val(paste0(log_val(), txt, "\n"))
  }
  
  ## 读取contig确定染色体
  observeEvent(input$large_vcf, {
    req(input$large_vcf)
    
    append_log("读取 contig...")
    
    contigs <- get_contigs(tools, input$large_vcf$datapath)
    
    if (!is.null(contigs)) {
      updateSelectInput(session, "large_vcf_chr", choices = contigs, selected = contigs[1])
      append_log(paste("检测到染色体:", paste(head(contigs), collapse = ", ")))
    } else {
      append_log("未检测到 contig")
    }
  })
  
  ## 开始过滤
  observeEvent(input$large_vcf_run, {
    req(input$large_vcf, input$large_vcf_chr)
    # 锁定按钮
    shinyjs::disable("large_vcf_run")
    on.exit(shinyjs::enable("large_vcf_run"))   # 无论成功失败都最终解锁
    
    infile <- input$large_vcf$datapath
    outfile <- tempfile(fileext = ".vcf")
    
    infile_use <- infile
    outfile_use <- outfile
    
    if (tools$is_windows) {
      infile_use <- to_cygwin_path(infile)
      outfile_use <- to_cygwin_path(outfile)
    }
    
    region <- paste0(input$large_vcf_chr, ":", input$large_vcf_start, "-", input$large_vcf_end)
    
    append_log(paste("Region:", region))
    append_log("检查 index...")
    
    tryCatch({
      withProgress(message = "Building VCF index (.tbi)...", value = 0, {
        setProgress(0.2, detail = "Use bcftools index, first index will being long Time")
        ensure_index(tools, infile, threads = input$large_vcf_threads)
        setProgress(1,   detail = "Index build success")
      })
      append_log("index OK")
    }, error = function(e) {
      append_log(paste("index Fail:", e$message))
      return(NULL)
    })
    
    
    # 构建参数
    if (input$large_vcf_header){
      args <- c("view","-r",region,"-o",outfile_use, infile_use)
    }else{
      args <- c("view","-H","-r",region,"-o",outfile_use, infile_use)
    }
    
    append_log("运行 bcftools...")
    
    res <- tryCatch({
      system2(tools$bcftools, args = args, stdout = TRUE, stderr = TRUE)
    }, error = function(e) {
      paste("ERROR:", e$message)
    })
    
    append_log(paste(res, collapse = "\n"))
    
    if (file.exists(outfile)) {
      result_file(outfile)
      append_log("完成")
    } else {
      append_log("失败：未生成输出文件")
    }
  })
  
  output$large_vcf_log <- renderText({
    log_val()
  })
  
  output$large_vcf_download <- downloadHandler(
    filename = function() "result.vcf",
    content = function(file) {
      req(result_file())
      file.copy(result_file(), file)
    }
  )
  
  observeEvent(input$large_vcf_clean, {
    f <- input$large_vcf$datapath
    
    if (!is.null(f) && file.exists(f)) {
      file.remove(f)
      append_log("输入 VCF 文件已清理")
      
      tbi <- paste0(f, ".tbi")
      if (file.exists(tbi)) {
        file.remove(tbi)
        append_log("索引文件已清理")
      }
    } else {
      append_log("无可清理的输入文件")
    }
    
    gc()
    append_log("GC 完成")
  })
  
  
  ##########---------- 基因型数据过滤 ----------##########
  ## 过滤vcf
  filtered_vcf <- eventReactive(input$filter_geno,{
    req(vcf())   # 防止未导入就运行
    
    if (input$filter_mode == "none"){
      vcf_data <- vcf()
    } else if(input$filter_mode == "POS" || is.null(active_annotation_type())){  
      vcf_data <- filter_vcf(vcf(),mode = "POS",
                             Chr = input$filter_chr,
                             start = input$filter_start,
                             end = input$filter_end
      )
    }else if (input$filter_mode %in% c("type","both")){
      vcf_data <- filter_vcf(vcf(),mode = input$filter_mode,
                             gff = gff(),
                             type = input$filter_type,
                             Chr = input$filter_chr,
                             start = input$filter_start,
                             end = input$filter_end
      )
    }
    
    return(vcf_data)
  })
  
  ## 过滤plink
  filtered_plink <- eventReactive(input$filter_geno,{
    req(plink())
    if (input$filter_mode == "none"){
      plink_data <- plink()
      
    } else if(input$filter_mode == "POS" || is.null(active_annotation_type())){  
      plink_data <- filter_plink.pedmap(plink(),mode = "POS",
                                        Chr = input$filter_chr,
                                        start = input$filter_start,
                                        end = input$filter_end
      )
    }else if (input$filter_mode %in% c("type","both")){
      plink_data <- filter_plink.pedmap(plink(),mode = input$filter_mode,
                                        gff = gff(),
                                        type = input$filter_type,
                                        Chr = input$filter_chr,
                                        start = input$filter_start,
                                        end = input$filter_end
      )
    }
    
    return(plink_data)
  })
  
  ## 过滤geno
  filtered_geno <- eventReactive(input$filter_geno,{
    req(geno())
    if (input$filter_mode == "none"){
      geno_data <- geno()
      
    } else if(input$filter_mode == "POS" || is.null(active_annotation_type())){  
      geno_data <- filter_table(geno(),mode = "POS",
                                Chr = input$filter_chr,
                                start = input$filter_start,
                                end = input$filter_end
      )
    }else if (input$filter_mode %in% c("type","both")){
      geno_data <- filter_table(geno(),mode = input$filter_mode,
                                gff = gff(),
                                type = input$filter_type,
                                Chr = input$filter_chr,
                                start = input$filter_start,
                                end = input$filter_end
      )
    }
    return(geno_data)
  })
  
  
  # VCF 摘要
  output$after_filter_vcf_plink_show <- renderPrint({
    req(active_genotype_type())
    switch (active_genotype_type(),
            "VCF" = {
              req(filtered_vcf())  # vcf() 是 eventReactive 定义的对象
              if (input$filter_mode == "none") {
                cat("filter mode: ",input$filter_mode,"\n")
                print(filtered_vcf())
              }
              
              if (input$filter_mode == "POS") {
                cat("filter mode: ",input$filter_mode,"\n")
                cat("filter position: ",
                    input$filter_chr, ":",
                    input$filter_start, "-",
                    input$filter_end, "\n")
                print(filtered_vcf())
                
              }
              
              if (input$filter_mode %in% c("type","both")) {
                cat("filter mode: ",input$filter_mode,"\n")
                
                if(is.null(active_annotation_type()) || is.null(input$filter_type)){
                  # if(is.null(gff()) || is.null(input$filter_type)){
                  return("Please input annotation or cancel filter by type")
                }
                
                print(filtered_vcf())
              }
            },
            "Plink" = {
              req(filtered_plink())  # vcf() 是 eventReactive 定义的对象
              # filter_vcf()       # 打印 VCF 对象信息
              if (input$filter_mode == "none") {
                cat("filter mode: ",input$filter_mode,"\n")
                cat("Individuals: ", nrow(filtered_plink()$ped), "\n")
                cat("Markers: ", nrow(filtered_plink()$map))
                # print(filtered_plink())
              }
              
              if (input$filter_mode == "POS") {
                cat("filter mode: ",input$filter_mode,"\n")
                cat("filter position: ",
                    input$filter_chr, ":",
                    input$filter_start, "-",
                    input$filter_end, "\n")
                cat("Individuals: ", nrow(filtered_plink()$ped), "\n")
                cat("Markers: ", nrow(filtered_plink()$map))
                # print(filtered_plink())
                
              }
              
              if (input$filter_mode %in% c("type","both")) {
                cat("filter mode: ",input$filter_mode,"\n")
                
                if(is.null(active_annotation_type()) || is.null(input$filter_type)){
                  # if(is.null(gff()) || is.null(input$filter_type)){
                  return("Please input annotation or cancel filter by type")
                }
                cat("Individuals: ", nrow(filtered_plink()$ped), "\n")
                cat("Markers: ", nrow(filtered_plink()$map))
                # print(filtered_plink())
              }
            },
            "Fasta" = {
              req(fasta())
              cat("Fsata data can not filter")
              cat("Fasta data: use table view below","\n")
              # print(fasta())
              # 将对象转换为字符并去除颜色码
              out <- capture.output(print(fasta()))
              cat(crayon::strip_style(out), sep = "\n")
            },
            "Geno" = {
              cat("Genotype data: use table view below")
            }
    )# swith end
  })
  
  observeEvent(input$filter_geno,{
    req(filtered_geno())
    if(active_genotype_type() != "Geno"){
      output$after_filter_geno_show <- renderDT({
        NULL
      })
    }else if (active_genotype_type() == "Geno" && input$filter_mode %in% c("type","both") && {is.null(active_annotation_type()) || is.null(input$filter_type)}){
      output$after_filter_vcf_plink_show <- renderPrint("Please input annotation or cancel filter by type")
      output$after_filter_geno_show <- renderDT({
        NULL
      })
    }else if (active_genotype_type() == "Geno"){
      output$after_filter_vcf_plink_show <- renderPrint({NULL})
      output$after_filter_geno_show <- renderDT({
        DT::datatable(filtered_geno(), options = list(pageLength = 5))
      })
    }
    
  })
  
  # # Fasta 表
  # output$fasta_table <- DT::renderDT({
  #   req(fasta())
  #   df <- data.frame(seq_name = names(fasta()),
  #                    seq_length = Biostrings::nchar(fasta()))
  #   DT::datatable(df, options = list(pageLength = 5))
  # })
  # 
  
  ##########---------- haplotype鉴定 ----------##########
  HapResult <- eventReactive(input$haplo_result_run,{
    req(active_genotype_type())
    HapRes <- switch (active_genotype_type(),
                      "VCF" = {
                        req(filtered_vcf())  # vcf() 是 eventReactive 定义的对象
                        vcf2hap(filtered_vcf(),hapPrefix = input$hap_prefix,
                                hetero_remove = input$remove_hetero_ornot,
                                na_drop = input$remove_unknown_ornot
                        )
                      },
                      "Plink" = {
                        req(filtered_plink())  # vcf() 是 eventReactive 定义的对象
                        plink.pedmap2hap(filtered_plink(),hapPrefix = input$hap_prefix,
                                         hetero_remove = input$remove_hetero_ornot,
                                         na_drop = input$remove_unknown_ornot
                                         
                        )
                      },
                      "Fasta" = {
                        req(fasta())
                        seqs2hap(fasta(),
                                 Ref = names(fasta())[1],
                                 # Ref = fasta()@unmasked@ranges@NAMES[1], # Mutiplealignment format use it
                                 hapPrefix = input$hap_prefix,
                                 hetero_remove = input$remove_hetero_ornot,
                                 na_drop = input$remove_unknown_ornot
                        )
                      },
                      "Geno" = {
                        req(filtered_geno())
                        table2hap(filtered_geno(),hapPrefix = input$hap_prefix,
                                  hetero_remove = input$remove_hetero_ornot,
                                  na_drop = input$remove_unknown_ornot  
                        )
                      }
    )# swith end
    
    
    return(HapRes)
  })
  
  HapSummary <- eventReactive(input$haplo_summary_run,{
    req(HapResult())
    HapSum <- hap_summary(HapResult(),hapPrefix = input$hap_prefix)
    
    # if(!is.null(HapSum)){
    #   # 刷新genemodel的长宽
    #   chr <- unique(HapSum[1, c(2:(ncol(HapSum) - 2))])
    #   chr_list <- as.list(chr)
    #   names(chr_list) <- chr
    #   updateSelectInput(session, "genemodel_chr",
    #                     choices = chr_list)
    #   updateSelectInput(session, "hapeff_chr",
    #                     choices = chr_list)
    #   updateSelectInput(session, "ldheatmap_chr",
    #                     choices = chr_list)
    #   pos <- as.numeric(HapSum[2,c(2:(ncol(HapSum) - 2))])
    #   updateNumericInput(session, "ldheatmap_start", value = min(pos))
    #   updateNumericInput(session, "ldheatmap_end", value = max(pos))
    #   updateNumericInput(session,"genemodel_width",value = max(pos))
    #   updateNumericInput(session,"genemodel_start",value = min(pos))
    #   updateNumericInput(session,"hapeff_end",value = max(pos))
    #   updateNumericInput(session,"hapeff_start",value = min(pos))
      
    # }
    
    return(HapSum)
  })
  
  observeEvent(input$haplo_result_run,{
    req(HapResult())
    HapResultTab <- DT::datatable(HapResult(), 
                                  rownames = F,
                                  options = list(pageLength = 10,scrollX = TRUE),
                                  escape = FALSE) %>%
      DT::formatStyle(
        columns = names(HapResult()),
        `white-space` = 'nowrap',
        `overflow` = 'hidden',
        `text-overflow` = 'ellipsis',
        `max-width` = '200px'
      )
    output$haplo_result_show <- DT::renderDT({HapResultTab},
                                             server = TRUE
    )
  })
  
  observeEvent(input$haplo_summary_run,{
    req(HapSummary())
    HapSummaryTab <- DT::datatable(HapSummary(),
                                   rownames = F,
                                   options = list(pageLength = 10,scrollX = TRUE),
                                   escape = FALSE) %>% 
      DT::formatStyle(
        columns = names(HapSummary()),
        `white-space` = 'nowrap',
        `overflow` = 'hidden',
        `text-overflow` = 'ellipsis',
        `max-width` = '200px'
      )
    output$haplo_summary_show <- DT::renderDT({HapSummaryTab},
                                              server = TRUE
    )
  })
  
  output$haplo_result_download <- downloadHandler(
    filename = function(){
      paste0(Sys.Date(),"_hapResult.csv")
    },
    content = function(file){
      write.hap(HapResult(),file,sep = ",")
    }
  )
  output$haplo_summary_download <- downloadHandler(
    filename = function(){
      paste0(Sys.Date(),"_hapSummary.csv")
    },
    content = function(file){
      write.hap(HapSummary(),file,sep = ",")
    }
  )
  
  ##########---------- 数据处理页 ----------##########
  ## 更新filtermode
  output$rm_condition <- renderUI({
    
    if (input$filter_hap_ornot){
      switch(input$rm_mode,
             "freq" = {
               numericInput("freq_number","Filter \"freq\" less than :",value = 5,min = 0,step = 1)
             },"position" = {
               textAreaInput("rm_position","Input remove variants site num perline",placeholder = "1234\n1266\n1277\n1388",rows = 5)
             },"accession" = {
               textAreaInput("rm_accession","Input remove accession name perline",placeholder = "B001\nB003\nCX011\nSample_23",rows = 5)
             }, "haplotype" = {
               textAreaInput("rm_hap","Input remove Haplotype name perline",placeholder = "Hap001\nHap002\nHap12\nHap20",rows = 5)
             }
      )
    }else{
      return(NULL)
    }
  })
  
  ## 可视化状态管理 
  
  active_hapresult <- reactiveVal(NULL)
  atgas0_hapresult <- reactiveVal(NULL)
  display_hapresult <- reactiveVal(NULL)
  gff0 <- reactiveVal(NULL)
  # ── 关键修改：监听Step3结果，自动初始化（用户未勾选import时生效）──
  observeEvent(HapResult(), {
    if (!isTRUE(input$import_own_data)) {
      active_hapresult(HapResult())
      display_hapresult(HapResult())
      # output$active_hapresult_info <- renderPrint("HapResult from Step3 is now active.")
      showNotification("HapResult from Step3 is now active.", 
                       type = "message", duration = 3)
    }
  })
  
  # ── 取消勾选import → 回退到Step3结果 ──
  observeEvent(input$import_own_data, ignoreInit = TRUE, {
    if (!input$import_own_data) {
      # ↓ 直接用Step3的reactive，不再需要hapresult_from_step3中间变量
      if (!is.null(HapResult())) {
        active_hapresult(HapResult())
        display_hapresult(HapResult())
        showNotification("Reverted to Step3 HapResult.", type = "message")
      } else {
        showNotification("Step3 HapResult not available yet. Please run Step3 first.", 
                         type = "warning")
      }
    }
  })
  
  # ── 导入用户自己的 HapResult ──
  observeEvent(input$own_hapres_import, {
    req(input$own_hapres)
    tryCatch({
      hapres <- geneHapR::import_hap(
        input$own_hapres$datapath,
        # type = "hapresult",
        sep = input$own_hapres_sep
      )
      if (!inherits(hapres, "hapResult")) {output$active_hapresult_info <- renderPrint("Not a valid HapResult (must be hapResult).")} 
      else {
        active_hapresult(hapres)
        display_hapresult(hapres)
        showNotification("Your HapResult imported successfully!", type = "message")
      }
    }, error = function(e) {
      showNotification(paste("Import failed:", e$message), type = "error", duration = 8)
    })
  })
  
  # ── Filter ──
  observeEvent(input$filter_hap, {
    req(active_hapresult())
    result <- tryCatch({
      switch(input$rm_mode,
             "freq" = {
               req(input$freq_number)
               filter_hap(active_hapresult(),rm.mode = input$rm_mode,freq.min = input$freq_number)
             },
             "position" = {
               req(input$rm_position)
               pos <- as.numeric(trimws(unlist(strsplit(input$rm_position, "\n"))))
               filter_hap(active_hapresult(),rm.mode = input$rm_mode,position.rm = pos)
             },
             "accession" = {
               req(input$rm_accession)
               acc <- trimws(unlist(strsplit(input$rm_accession, "\n")))
               filter_hap(active_hapresult(),rm.mode = input$rm_mode,accession.rm = acc)
             },
             "haplotype" = {
               req(input$rm_hap)
               hap <- trimws(unlist(strsplit(input$rm_hap, "\n")))
               filter_hap(active_hapresult(),rm.mode = input$rm_mode,haplotype.rm = hap)
             }
      )
    }, error = function(e) {
      showNotification(paste("Filter failed:", e$message), type = "error", duration = 8)
      NULL
    })
    if (!is.null(result)) {
      active_hapresult(result)   # ← filter结果成为新工作基础
      display_hapresult(result)
      showNotification("Filter applied.", type = "message")
    }
  })
  # ================================================================
  # Set ATG as 0（独立分支，不覆盖 active_hapresult）
  # ================================================================
  
  observeEvent(input$set_atgas0, {
    req(active_hapresult(), gff(),input$atgas0_geneid)
    
    result <- tryCatch({
      geneHapR::hapSetATGas0(
        active_hapresult(),
        gff = gff(),
        geneID = input$atgas0_geneid,
        Chr = input$atgas0_chr,
        POS = c(input$atgas0_start,input$atgas0_end)
      )
    }, error = function(e) {
      showNotification(paste("Haplo Set ATG as 0 failed:", e$message), 
                       type = "error", duration = 8)
      NULL
    })
    gff0 <- tryCatch({
      geneHapR::gffSetATGas0(
        active_hapresult(),
        gff = gff(),
        geneID = input$atgas0_geneid,
        Chr = input$atgas0_chr,
        POS = c(input$atgas0_start,input$atgas0_end)
      )
    }, error = function(e) {
      showNotification(paste("GFF Set ATG as 0 failed:", e$message), 
                       type = "error", duration = 8)
      NULL
    })
    
    if (!is.null(result)) {
      atgas0_hapresult(result)   # 存入独立object，不动active
      display_hapresult(result)  # 但显示最新产物
      showNotification("Haplotype ATG set as 0 successfully! Stored as separate object.", 
                       type = "message", duration = 3)
    }
    if (!is.null(gff0)) {
      gff0(gff0)
      showNotification("GFF ATG set as 0 successfully! Stored as separate object.", 
                       type = "message", duration = 3)
    }
    
  })
  
  
  
  
  # ================================================================
  # Output：信息栏 + 表格
  # ================================================================
  
  # 显示当前Active HapResult的来源说明
  output$active_hapresult_info <- renderText({
    d <- display_hapresult()
    if (is.null(d)) return("No HapResult loaded yet.")
    
    # 判断当前显示的是哪个object（辅助用户理解）
    src_tag <- if (identical(d, atgas0_hapresult())) {
      "[Source: ATG-as-0 Result]"
    } else if (identical(d, active_hapresult())) {
      "[Source: Active HapResult (filtered/imported)]"
    } else {
      "[Source: Unknown]"
    }
    
    paste0(
      src_tag, "\n",
      "Haplotypes: ", length(unique(d$Hap))-4, " | ",
      "Sites: ",      sum(grepl("^[0-9]", colnames(d))), " | ",
      "Accessions: ", length(unique(d$Accession))-1
    )
  })
  
  output$current_HapResult <- DT::renderDT({
    req(display_hapresult())
    DT::datatable(
      display_hapresult(),
      options = list(scrollX = TRUE, pageLength = 7),
      rownames = FALSE
    ) %>% 
      DT::formatStyle(
        columns = names(display_hapresult()),
        `white-space` = 'nowrap',
        `overflow` = 'hidden',
        `text-overflow` = 'ellipsis',
        `max-width` = '200px'
      )
  },server = T)
  
  output$current_hapresult_download <- downloadHandler(
    filename = function(){
      paste0(Sys.Date(),{if (identical(display_hapresult(), atgas0_hapresult())) {
        "_ATG_As0_"
      } else if (identical(display_hapresult(), active_hapresult())) {
        "_Filtered_"
      } else {
        "_"
      }},"HapResult.csv")
    },
    content = function(file){
      write.hap(display_hapresult(),file,sep = ",")
    }
  )
  
  
  # HapSummary
  display_hap_summary <- reactive({
    req(display_hapresult())
    tryCatch(
      geneHapR::hap_summary(display_hapresult(),
                            # hapPrefix = input$hap_prefix
                            hapPrefix = "Hap"
      ),
      error = function(e) NULL
    )
  })
  
  output$active_hapsummary_info <- renderText({
    s <- display_hap_summary()
    if (is.null(s)) return("HapSummary not available.")
    paste0("Total haplotypes: ", nrow(s)-4)
  })
  
  output$current_hapsummary <- DT::renderDT({
    req(display_hap_summary())
    DT::datatable(
      display_hap_summary(),
      options = list(scrollX = TRUE, pageLength = 7),
      rownames = FALSE
    ) %>% 
      DT::formatStyle(
        columns = names(display_hap_summary()),
        `white-space` = 'nowrap',
        `overflow` = 'hidden',
        `text-overflow` = 'ellipsis',
        `max-width` = '200px'
      )
  },server = T)
  
  output$current_hapsummary_download <- downloadHandler(
    filename = function(){
      paste0(Sys.Date(),{if (identical(display_hapresult(), atgas0_hapresult())) {
        "_ATG_As0_"
      } else if (identical(display_hapresult(), active_hapresult())) {
        "_Filtered_"
      } else {
        "_"
      }},"HapSummary.csv")
    },
    content = function(file){
      write.hap(display_hap_summary(),file,sep = ",")
    }
  )
  
  ##########----------画图页----------##########
  ##### Visulization
  ##### plotHapTable
  observeEvent(display_hap_summary(),{
    req(atgas0_hapresult())
    if (!is.null(atgas0_hapresult())){
      updateSelectInput(session,"hap_table_input",choices = c("Origin"="ori","ATGas0" = "atgas0"),selected = "atgas0")
      updateSelectInput(session,"genemodel_input",choices = c("Origin"="ori","ATGas0" = "atgas0"),selected = "atgas0")
    }
  })
  
  # 读取到 HapSummary 后，自动更新相关输入选项（染色体、位置、Hap名称等），提升用户体验
  observeEvent(display_hap_summary(),{
    req(active_hapresult(),display_hap_summary())
    if(inherits(active_hapresult(), "hapResult")){
      
      # haplotable的
      info <- active_hapresult()[3, -c(1, ncol(active_hapresult()))] %>%
        data.frame() %>% t()
      info <- info[,1]
      info <- sapply(info, function(x) strsplit(x, "=")[[1]][1]) %>%
        unlist() %>% unique()
      if(length(info) >= 1){
        infos <- as.list(c(info,"none"))
        names(infos) <- c(info, "none")
        updateSelectInput(session, "hap_table_infotag",
                          choices = infos,
                          selected = "none")
      }
    }
    HapSum <- display_hap_summary()
    if (!is.null(HapSum)) {
        # 刷新genemodel的长宽
        chr <- unique(HapSum[1, c(2:(ncol(HapSum) - 2))])
        chr_list <- as.list(chr)
        names(chr_list) <- chr
        updateSelectInput(session, "genemodel_chr", choices = chr_list)
        updateSelectInput(session, "hapeff_chr", choices = chr_list)
        updateSelectInput(session, "ldheatmap_chr", choices = chr_list)
        
        pos <- as.numeric(HapSum[2, c(2:(ncol(HapSum) - 2))])
        updateNumericInput(session, "ldheatmap_start", value = min(pos))
        updateNumericInput(session, "ldheatmap_end", value = max(pos))
        updateNumericInput(session, "genemodel_width", value = max(pos))
        updateNumericInput(session, "genemodel_start", value = min(pos))
        updateNumericInput(session, "hapeff_end", value = max(pos))
        updateNumericInput(session, "hapeff_start", value = min(pos))
        
        hap_names_list <- as.list(HapSum[,"Hap"][5:nrow(HapSum)])
        names(hap_names_list) <- HapSum[,"Hap"][5:nrow(HapSum)]
        updateSelectInput(session, "geodis_hapnames",
                          choices = hap_names_list,
                          selected = hap_names_list[1:3])
    }
  })
  
  
  
  HapTablePlot <- eventReactive(input$hap_table_plot,{
    switch (input$hap_table_input,
            "ori" = {
              req(active_hapresult())
              hapsum <- hap_summary(active_hapresult(),hapPrefix = input$hap_table_prefix)
            },
            "atgas0" = {
              req(atgas0_hapresult())
              hapsum <- hap_summary(atgas0_hapresult(),hapPrefix = input$hap_table_prefix)
            }
    )
    
    if(input$hap_table_infotag == "none") {
      infotag <- tagsplit <- tagfield <- tagname <- NULL
    } else {
      infotag <- input$hap_table_infotag
      tagsplit <- input$hap_table_tagsplit
      tagfield <- input$hap_table_tagfield
      tagname <- input$hap_table_tagname
    }
    
    haptableplot <- geneHapR::plotHapTable(hapSummary = hapsum,
                                           hapPrefix = input$hap_table_prefix,
                                           geneName = input$hap_table_genename,
                                           INFO_tag = infotag,
                                           tag_split = tagsplit,
                                           tag_field = tagfield,
                                           tag_name = tagname,ALLELE.color = input$hap_table_bg_color,
                                           angle = input$hap_table_angle)
    output$hap_table_plot_info <- renderPrint(p("Press Button again to refresh plot", style = "font-weight:bold;font-size:12pt ;color:red"))
    
    return(haptableplot)
  })
  
  # 图片输出框
  output$hap_table_plotout <- renderPlot({
    req(HapTablePlot())
    HapTablePlot()
  },
  width = function() input$hap_table_plotout_width,
  height = function() input$hap_table_plotout_height,
  res = 96)
  # 图片下载框
  output$hap_table_download <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(),"_",input$hap_table_input, '_HapTable.', input$hap_table_out_type)
    },
    content = function(file) {
      if (input$hap_table_out_type == "pdf") {
        ggsave(filename = file, plot = HapTablePlot(),
               width = as.numeric(input$hap_table_down_width),
               height = as.numeric(input$hap_table_down_height), dpi = 300,
               device = cairo_pdf)
      } else {
        ggsave(filename = file, plot = HapTablePlot(),
               width = as.numeric(input$hap_table_down_width),
               height = as.numeric(input$hap_table_down_height), dpi = 300)
      }
      
      
    }
  )
  
  ##### Visulization
  ##### genemodel
  GeneModelPlot <- eventReactive(input$genemodel_plot,{
    switch (input$genemodel_input,
            "ori" = {
              req(active_hapresult(),gff())
              hapsum <- hap_summary(active_hapresult())
              gff1 <- gff()
            },
            "atgas0" = {
              req(atgas0_hapresult(),gff0())
              hapsum <- hap_summary(atgas0_hapresult())
              gff1 <- gff0()
            }
    )
    
    genemodelplot <- displayVarOnGeneModel(gff = gff1,          # Annotations
                                           hap = hapsum,        # Haplotype result
                                           Chr = input$genemodel_chr,        # Chromosome name
                                           startPOS = input$genemodel_start,        # Start position of gene model
                                           endPOS = input$genemodel_width,
                                           cex = input$genemodel_cex,type = input$genemodel_type,
                                           CDS_h = input$genemodel_cdsh,
                                           fiveUTR_h = input$genemodel_cdsh*input$genemodel_utrvscds, 
                                           threeUTR_h = input$genemodel_cdsh*input$genemodel_utrvscds,
                                           geneElement = input$geneElement
    )
    output$genemodel_plot_info <- renderPrint(p("Press Button again to refresh plot", style = "font-weight:bold;font-size:12pt ;color:red"))
    
    return(genemodelplot)
  })
  
  # 图片输出框
  output$genemodel_plotout <- renderPlot({
    req(GeneModelPlot())
    GeneModelPlot()
  },
  width = function() input$genemodel_plotout_width,
  height = function() input$genemodel_plotout_height,
  res = 96)
  # 图片下载框
  output$genemodel_download <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(),"_",input$genemodel_input, '_GeneModel.', input$genemodel_out_type)
    },
    content = function(file) {
      # ggsave(filename = file, plot = GeneModelPlot(),
      #        width = as.numeric(input$genemodel_down_width),
      #        height = as.numeric(input$genemodel_down_height), dpi = 300)
      # 根据当前输入值重新生成图形数据（不依赖 eventReactive）
      hapsum <- switch(input$genemodel_input,
                       "ori" = {
                         req(active_hapresult(), gff())
                         hap_summary(active_hapresult())
                       },
                       "atgas0" = {
                         req(atgas0_hapresult(), gff0())
                         hap_summary(atgas0_hapresult())
                       })
      gff1 <- switch(input$genemodel_input,
                     "ori" = gff(),
                     "atgas0" = gff0())
      
      
      switch(input$genemodel_out_type,
             "pdf" = {
               pdf(file, width = as.numeric(input$genemodel_down_width), height = as.numeric(input$genemodel_down_height), useDingbats = FALSE)
             },
             "png" = {
               png(file, width = as.numeric(input$genemodel_down_width), height = as.numeric(input$genemodel_down_height), units = "in", res = 300)
             },
             "jpeg" = {
               jpeg(file, width = as.numeric(input$genemodel_down_width), height = as.numeric(input$genemodel_down_height), units = "in", res = 300, quality = 90)
             }
             
      )
      
      # ? 关键三步
      grid::grid.newpage()
      grid::popViewport(0)
      
      g <- grid::grid.grabExpr({
        displayVarOnGeneModel(
          gff = gff1,
          hap = hapsum,
          Chr = input$genemodel_chr,
          startPOS = input$genemodel_start,
          endPOS = input$genemodel_width,
          cex = input$genemodel_cex,
          type = input$genemodel_type,
          CDS_h = input$genemodel_cdsh,
          fiveUTR_h = input$genemodel_cdsh * input$genemodel_utrvscds,
          threeUTR_h = input$genemodel_cdsh * input$genemodel_utrvscds,
          geneElement = input$geneElement
        )
      })
      
      grid::grid.draw(g)
      dev.off()
      
    })
  
  
  
  ##### Visulization
  ##### Haplo Network
  # plothapnet
  legend_hapnet <- reactive(
    {if(any(input$plothapnet_legend_size,
            input$plothapnet_legend_color))
      c(input$plothapnet_x, input$plothapnet_y) else
        FALSE
    })
  
  HapNet <- reactive({
    req(active_hapresult())
    if (input$plothapnet_group != "none") {
      geneHapR::get_hapNet(
        hap_summary(active_hapresult()),
        AccINFO = accession(),
        groupName = input$plothapnet_group
      )
    } else {
      geneHapR::get_hapNet(hap_summary(active_hapresult()))
    }
  })
  
  
  
  
  observeEvent(input$hapnet_plot,{
    req(HapNet())
    
    # if(input$plothapnet_group != "none"){
    #   accinfo <- accession()
    #   message("accinfo")
    #   hapNet <- geneHapR::get_hapNet(hap_summary(active_hapresult()), AccINFO = accinfo,
    #                                  groupName = input$plothapnet_group)
    # } else {
    #   hapNet <- geneHapR::get_hapNet(hap_summary(active_hapresult()))
    # }
    # message("ploting")
    hapNet <- HapNet()
    # message("lenged: ", legend)
    xlim <- input$plothapnet_xlim
    ylim <- input$plothapnet_ylim
    f <- tempfile()
    png(f)
    geneHapR::plotHapNet(hapNet)
    m <- round(par("usr"))
    dev.off()
    unlink(f)
    updateSliderInput(session, "plothapnet_xlim",
                      min = m[1]*3,
                      max = m[2]*3,
                      value = c(m[1],m[2]))
    updateSliderInput(session, "plothapnet_ylim",
                      min = m[3]*3,
                      max = m[4]*3,
                      value = c(m[3],m[4]))
    updateSliderInput(session, "plothapnet_x",
                      min = m[1]*3,
                      max = m[2]*3,
                      value = m[1] + m[2])
    updateSliderInput(session, "plothapnet_y",
                      min = m[3]*3,
                      max = m[4]*3,
                      value = m[3] + m[4])
    
    
    output$hapnet_plotout <- renderPlot({
      req(input$plothapnet_showmutation)
      pie_lim <- input$plothapnet_pie.lim * (2 ^ input$plothapnet_pie.lim.folder)
      geneHapR::plotHapNet(hapNet, scale = input$plothapnet_scale,
                           show.mutation = as.numeric(input$plothapnet_showmutation),
                           labels.cex = input$plothapnet_cex,
                           labels = as.logical(input$plothapnet_labels),
                           pie.lim = pie_lim,
                           cex.legend = input$plothapnet_cexlegend,
                           col.link = input$plothapnet_collink,
                           link.width = input$plothapnat_linkwidth,
                           legend = legend_hapnet(),
                           show_color_legend = input$plothapnet_legend_color,
                           show_size_legend = input$plothapnet_legend_size,
                           xlim = input$plothapnet_xlim,
                           ylim = input$plothapnet_ylim)
    }, width = 960, height = 720, res = 96)
    
    
  })
  
  
  
  
  
  # 图片下载框
  output$hapnet_download <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), '_HaploNetwork.', input$hapnet_out_type)
    },
    content = function(file) {
      # 根据当前输入值重新生成图形数据（不依赖 eventReactive）
      switch(input$hapnet_out_type,
             "pdf" = {
               pdf(file, width = 9.6, height = 7.2, useDingbats = FALSE)
             },
             "png" = {
               png(file, width = 9.6, height =7.2, units = "in", res = 300)
             },
             "jpeg" = {
               jpeg(file, width = 9.6, height =7.2, units = "in", res = 300, quality = 90)
             }
             
      )
      
      # ? 关键：重置 base 图形环境
      # graphics::plot.new()
      # （可选但推荐）控制边距，防止“出界”
      par(mar = c(5, 5, 2, 2))
      # 计算参数（与你 renderPlot 一致）
      pie_lim <- input$plothapnet_pie.lim * (2 ^ input$plothapnet_pie.lim.folder)
      # 绘图（直接画）
      geneHapR::plotHapNet(
        HapNet(),
        scale = input$plothapnet_scale,
        show.mutation = input$plothapnet_showmutation,
        labels.cex = input$plothapnet_cex,
        labels = as.logical(input$plothapnet_labels),
        pie.lim = pie_lim,
        cex.legend = input$plothapnet_cexlegend,
        col.link = input$plothapnet_collink,
        link.width = input$plothapnat_linkwidth,
        legend = legend_hapnet(),
        show_color_legend = input$plothapnet_legend_color,
        show_size_legend = input$plothapnet_legend_size,
        xlim = input$plothapnet_xlim,
        ylim = input$plothapnet_ylim
      )
      
      # 关闭设备
      dev.off()
      
    })
  
  ##### Visulization
  ##### HapVSPhone
  
  HapVsPheno <- eventReactive(input$hapvspheno_plot,{
    req(phenotypes(),active_hapresult())
    result <- try(geneHapR::hapVsPheno(hap = active_hapresult(), phenotypes(),
                                       mergeFigs = T,
                                       minAcc = input$hapvspheno_minacc,
                                       phenoName = input$hapvspheno_pheno_name,
                                       title = input$hapvspheno_title,
                                       angle = input$hapvspheno_angle
    ))
    return(result)
  })
  
  observeEvent(HapVsPheno(),{
    req(HapVsPheno())
    if(inherits(HapVsPheno(), "try-error")){
      output$hapvspheno_plotout <- renderPlot({
        plot.new()
        title(as.character(HapVsPheno()))
      })
    } else{
      output$hapvspheno_plotout <- renderPlot({HapVsPheno()$figs},
                                              width = function() input$hapvspheno_plotout_width,
                                              height = function() input$hapvspheno_plotout_height,
                                              res = 96)
      output$hapvspheno_testout <- renderPrint({HapVsPheno()$T.Result})
    }
  })
  # 图片输出框
  
  # 图片下载框
  output$hapvspheno_download <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), '_Haplotypes_Vs_',input$hapvspheno_pheno_name,'_Phenotypes.', input$hapvspheno_out_type)
    },
    content = function(file) {
      if (input$hapvspheno_out_type == "pdf") {
        ggsave(filename = file, plot = HapVsPheno()$figs,
               width = as.numeric(input$hapvspheno_down_width),
               height = as.numeric(input$hapvspheno_down_height), dpi = 300,
               device = cairo_pdf)
      } else {
        ggsave(filename = file, plot = HapVsPheno()$figs,
               width = as.numeric(input$hapvspheno_down_width),
               height = as.numeric(input$hapvspheno_down_height), dpi = 300)
      }
      
      
    }
  )
  
  
  ##### Visulization
  ##### Haplo Effect
  HapEFF <- eventReactive(input$hapeff_plot,{
    req(active_hapresult(),phenotypes())
    EFF <- siteEFF(active_hapresult(), phenotypes())
    return(EFF)
    
  })
  
  observeEvent(input$hapeff_plot,{
    req(HapEFF(),gff())
    output$HapEFF_plotout <- renderPlot({
      plotEFF(HapEFF(), gff = gff(),
              Chr = input$hapeff_chr, 
              start = input$hapeff_start, 
              end = input$hapeff_end,
              showType = input$hapeff_show_type,# c("five_prime_UTR", "CDS", "three_prime_UTR"), # see help(plotEFF)
              y = input$hapeff_y,                      # the means of y axis, one of effect or pvalue
              ylab = input$hapeff_y,                  # label of y axis
              cex = input$hapeff_cex,                         # Cex
              legend.cex = input$hapeff_legend_cex,                  # legend size
              main = input$hapeff_title,                     # main title
              CDS.height = input$hapeff_cdsh,                    # controls the height of CDS, heights of others will be half of that
              markMutants = TRUE,                # mark mutants by short lines
              mutants.col = input$hapeff_line_color, # 1-9
              mutants.type = input$hapeff_line_type, # parameters for appearance of mutants 1-6
              pch = input$hapeff_point_type)         
    },width = function() input$hapeff_plotout_width,
    height = function() input$hapeff_plotout_height,
    res = 96)
    
    output$hapeff_eff_output <- DT::renderDT({
      DT::datatable(
        as.data.frame(HapEFF()$EFF),
        options = list(scrollX = TRUE, pageLength = 7),
        rownames = FALSE
      )
    },server = T)
    
    output$hapeff_p_output <- DT::renderDT({
      DT::datatable(
        as.data.frame(HapEFF()$p),
        options = list(scrollX = TRUE, pageLength = 7),
        rownames = FALSE
      )
    },server = T)
  })
  
  output$hapeff_eff_output_download <- downloadHandler(
    filename = function(){
      paste0(Sys.Date(), '_HaploEFF_Effect.csv')
    },content = function(file){
      write.csv(as.data.frame(HapEFF()$EFF),file,quote = F,)
    })
  output$hapeff_p_output_download <- downloadHandler(
    filename = function(){
      paste0(Sys.Date(), '_HaploEFF_Pvalue.csv')
    },content = function(file){
      write.csv(as.data.frame(HapEFF()$p),file,quote = F,)
    })
  
  # 图片下载框
  output$hapeff_download <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(),'_HaploEffPlot.', input$hapeff_out_type)
    },
    content = function(file) {
      switch(input$hapeff_out_type,
             "pdf" = {
               pdf(file, width = as.numeric(input$hapeff_down_width), height = as.numeric(input$hapeff_down_height), useDingbats = FALSE)
             },
             "png" = {
               png(file, width = as.numeric(input$hapeff_down_width), height = as.numeric(input$hapeff_down_height), units = "in", res = 300)
             },
             "jpeg" = {
               jpeg(file, width = as.numeric(input$hapeff_down_width), height = as.numeric(input$hapeff_down_height), units = "in", res = 300, quality = 90)
             }
             
      )
      
      plotEFF(HapEFF(), gff = gff(),
              Chr = input$hapeff_chr, 
              start = input$hapeff_start, 
              end = input$hapeff_end,
              showType = input$hapeff_show_type,#c("five_prime_UTR", "CDS", "three_prime_UTR"), # see help(plotEFF)
              y = input$hapeff_y,                      # the means of y axis, one of effect or pvalue
              ylab = input$hapeff_y,                  # label of y axis
              cex = input$hapeff_cex,                         # Cex
              legend.cex = input$hapeff_legend_cex,                  # legend size
              main = input$hapeff_title,                     # main title
              CDS.height = input$hapeff_cdsh,                    # controls the height of CDS, heights of others will be half of that
              markMutants = TRUE,                # mark mutants by short lines
              mutants.col = input$hapeff_line_color, # 1-9
              mutants.type = input$hapeff_line_type, # parameters for appearance of mutants 1-6
              pch = input$hapeff_point_type)
      dev.off()
      
    })
  
  ##### Visulization
  ##### LD_heatmap
  observeEvent(input$ldheatmap_plot, {
    req(active_hapresult())
    output$ldheatmap_plotout <- renderPlot({
      plot_LDheatmap2(
        hap = active_hapresult(),
        add.map = as.logical(input$ldheatmap_addmap),
        gff = {
          if (as.logical(input$ldheatmap_addmap)) gff() else NULL
        },
        Chr = input$ldheatmap_chr, # 染色体名称，用于定位基因区域
        start = input$ldheatmap_start, # 基因区域的起始位置（bp）
        end = input$ldheatmap_end, # 基因区域的终止位置（bp）
        # geneID = input$ldheatmap_geneid, # 目标基因ID，用于从GFF中筛选相关注释
        map.height = input$ldheatmap_mapheight, # 基因模型图的高度比例（占图形总高度的比例）
        colorLegend = TRUE, # 是否显示LD值的颜色图例
        geneMapLocation = input$ldheatmap_maploc, # 基因模型在图形中的垂直位置（从顶部算起的比例）
        geneMapLabelX = NULL, # 基因模型标签的X坐标，NULL表示自动定位
        geneMapLabelY = NULL, # 基因模型标签的Y坐标，NULL表示自动定位
        distances = input$ldheatmap_distances, # 距离度量："physical"为物理距离(bp)，"genetic"为遗传距离(cM)
        LDmeasure = input$ldheatmap_method, # LD度量方式："r"表示r²，"D'"表示|D'|
        title = input$ldheatmap_title, # 图表的标题文字
        SNP.name = T, # 是否在热图边缘显示SNP名称
        color = c(input$ldheatmap_highcolor, input$ldheatmap_midcolor, input$ldheatmap_lowcolor), # 自定义颜色梯度，NULL使用默认红-蓝-灰渐变
        color_gmodel = input$ldheatmap_genemodelcolor, # 基因模型的颜色
        color_snp = input$ldheatmap_snpcolor, # SNP标记点的颜色
        color_snpname = input$ldheatmap_snpcolor, # SNP名称标签的颜色
        cex_snpname = input$ldheatmap_cex, # SNP名称标签的字体大小倍数
        snpmarks_height = NULL, # SNP标记线的高度，NULL表示自动计算
        newpage = T, # 是否在新页面中绘制图形（FALSE表示在当前页面添加）
        name = "ldheatmap", # 图形对象的名称，用于后续引用
        vp.name = NULL, # 视图端口的名称，用于指定图形绘制位置
        pop = FALSE, # 是否在绘制后弹出视图端口（FALSE保持图形对象）
        text = as.logical(input$ldheatmap_text) # 是否在热图单元格中显示LD数值（FALSE仅显示颜色）
      )
    },width = function() input$ldheatmap_plotout_width,
    height = function() input$ldheatmap_plotout_height,
    res = 96)
  })
  # 图片下载框
  output$ldheatmap_download <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), '_LDHeatmap.', input$ldheatmap_out_type)
    },
    content = function(file) {
      switch(
        input$ldheatmap_out_type,
        "pdf" = {
          pdf(file, width = as.numeric(input$ldheatmap_down_width), height = as.numeric(input$ldheatmap_down_height), useDingbats = TRUE)
        },
        "png" = {
          png(file, width = as.numeric(input$ldheatmap_down_width), height = as.numeric(input$ldheatmap_down_height), units = "in", res = 300)
        },
        "jpeg" = {
          jpeg(
            file,
            width = as.numeric(input$ldheatmap_down_width),
            height = as.numeric(input$ldheatmap_down_height),
            units = "in",
            res = 300,
            quality = 90
          )
        }
      )
      
      plot_LDheatmap2(
        hap = active_hapresult(),
        add.map = as.logical(input$ldheatmap_addmap),
        gff = {
          if (as.logical(input$ldheatmap_addmap)) gff() else NULL
        },
        Chr = input$ldheatmap_chr, # 染色体名称，用于定位基因区域
        start = input$ldheatmap_start, # 基因区域的起始位置（bp）
        end = input$ldheatmap_end, # 基因区域的终止位置（bp）
        # geneID = input$ldheatmap_geneid, # 目标基因ID，用于从GFF中筛选相关注释
        map.height = input$ldheatmap_mapheight, # 基因模型图的高度比例（占图形总高度的比例）
        colorLegend = TRUE, # 是否显示LD值的颜色图例
        geneMapLocation = input$ldheatmap_maploc, # 基因模型在图形中的垂直位置（从顶部算起的比例）
        geneMapLabelX = NULL, # 基因模型标签的X坐标，NULL表示自动定位
        geneMapLabelY = NULL, # 基因模型标签的Y坐标，NULL表示自动定位
        distances = input$ldheatmap_distances, # 距离度量："physical"为物理距离(bp)，"genetic"为遗传距离(cM)
        LDmeasure = input$ldheatmap_method, # LD度量方式："r"表示r²，"D'"表示|D'|
        title = input$ldheatmap_title, # 图表的标题文字
        SNP.name = T, # 是否在热图边缘显示SNP名称
        color = c(input$ldheatmap_highcolor, input$ldheatmap_midcolor, input$ldheatmap_lowcolor), # 自定义颜色梯度，NULL使用默认红-蓝-灰渐变
        color_gmodel = input$ldheatmap_genemodelcolor, # 基因模型的颜色
        color_snp = input$ldheatmap_snpcolor, # SNP标记点的颜色
        color_snpname = input$ldheatmap_snpcolor, # SNP名称标签的颜色
        cex_snpname = input$ldheatmap_cex, # SNP名称标签的字体大小倍数
        snpmarks_height = NULL, # SNP标记线的高度，NULL表示自动计算
        newpage = T, # 是否在新页面中绘制图形（FALSE表示在当前页面添加）
        name = "ldheatmap", # 图形对象的名称，用于后续引用
        vp.name = NULL, # 视图端口的名称，用于指定图形绘制位置
        pop = FALSE, # 是否在绘制后弹出视图端口（FALSE保持图形对象）
        text = as.logical(input$ldheatmap_text) # 是否在热图单元格中显示LD数值（FALSE仅显示颜色）
      )
      dev.off()
    }
  )
  

  ##### Visulization
  ##### Geographical Distribution
  # 当用户选择geodis_exclude选项为TRUE时，读取世界地图的所有地区信息，并返回给geodis_excluderegions输入
    observeEvent(input$geodis_exclude, {
        if (input$geodis_exclude) {
        world_map <- ggplot2::map_data("world")
        regions <- unique(world_map$region)
        regions_list <- as.list(regions)
        names(regions_list) <- regions
        updateSelectizeInput(session, "geodis_excluderegions", choices = regions_list, select ="Antarctica",server = TRUE)
        } else {
        updateSelectizeInput(session, "geodis_excluderegions", choices = NULL, server = TRUE)
        }
    })
  # 当geodis_hapnames输入发生变化且geodis_hapcolor为TRUE时，
  # 根据当前的geodis_hapnames选择更新与之对应数量的colourInput，返回到geodis_hapcolor_ui输出
  # 每个colourInput的id格式为"geodis_hapcolor_i"，其中i是haplotype名称在geodis_hapnames中的索引
  # 返回的colorInput横向排列，每行不超过3个
    # observeEvent(input$geodis_hapcolor, {
    #     hap_names <- input$geodis_hapnames
    #     if (length(hap_names) > 0 && input$geodis_hapcolor == "TRUE") {
    #     color_inputs <- lapply(seq_along(hap_names), function(i) {
    #         colourInput(
    #         inputId = paste0("geodis_hapcolor_", i),
    #         label = paste("Color for", hap_names[i]),
    #         value = RColorBrewer::brewer.pal(8, "Set2")[i %% 8 + 1]
    #         )
    #     })
    #     output$geodis_hapcolor_ui <- renderUI({
    #         do.call(tagList, color_inputs)
    #     })
    #     } else {
    #     output$geodis_hapcolor_ui <- renderUI({})
    #     }
    # })
  observeEvent(
    list(input$geodis_hapcolor, input$geodis_hapnames),
    {
        hap_names <- input$geodis_hapnames

        if (length(hap_names) > 0 && input$geodis_hapcolor == "TRUE") {
            color_inputs <- lapply(seq_along(hap_names), function(i) {
                colourInput(
                    inputId = paste0("geodis_hapcolor_", i),
                    label = paste("Color for", hap_names[i]),
                    value = RColorBrewer::brewer.pal(8, "Set2")[(i - 1) %% 8 + 1]
                )
            })

            # 每3个一行
            rows <- split(color_inputs, ceiling(seq_along(color_inputs) / 3))

            output$geodis_hapcolor_ui <- renderUI({
                tagList(
                    lapply(rows, function(one_row) {
                        fluidRow(
                            lapply(one_row, function(x) {
                                column(width = 4, x) # 12/3 = 4
                            })
                        )
                    })
                )
            })
        } else {
            output$geodis_hapcolor_ui <- renderUI(NULL)
        }
    },
    ignoreInit = FALSE
)
  
  # 绘图
  Geo_dis_plot <- eventReactive(input$geodis_plot, {
    req(active_hapresult(), accession())

    # 如果用户选择了为不同的haplotype指定颜色，则从输入中提取这些颜色值；否则，使用默认颜色方案
    HapColor <- if (input$geodis_hapcolor == "TRUE") {
        sapply(seq_along(input$geodis_hapnames), function(i) {
        input[[paste0("geodis_hapcolor_", i)]]
        })
    } else {
        NULL
    }

    # 绘图
    plot <- HaploWorldMap(
        hap = active_hapresult(),
        AccINFO = accession(),
        LON.col = input$geodis_loncol,
        LAT.col = input$geodis_latcol,
        hapNames = input$geodis_hapnames,
        hapColor = HapColor,
        pie.log.base = input$geodis_pielogbase,
        pie.line.color = input$geodis_pielinecolor,
        pie.line.width = input$geodis_pielinewidth,
        pie.legend.size = input$geodis_pielegendsize,
        pie.legend.posi = c(input$geodis_pielegend_posix, input$geodis_pielegend_posiy),
        pie.legend.title.size = input$geodis_pielegendtitle_size,
        pie.legend.title.posi = c(input$geodis_pielegendtitle_posix, input$geodis_pielegendtitle_posiy),
        hap.legend.size = input$geodis_haplegend_size,
        hap.legend.posi = input$geodis_haplegend_posi,
        map.fill.color = input$geodis_mapfill,
        map.border.color = input$geodis_mapborder,
        map.border.width = input$geodis_mapborderwidth,
        map.background.color = input$geodis_mapbackground,
        exclude.regions = as.logical(input$geodis_exclude),
        exclude.regions.names = input$geodis_excluderegions
    )
    return(plot)
  })

    output$geodis_plotout <- renderPlot({
        req(Geo_dis_plot())
        Geo_dis_plot()
    }, 
    width = function() input$geodis_plotout_width, 
    height = function() input$geodis_plotout_height, 
    res = 96)
  
    # 图片下载框
    output$geodis_download <- downloadHandler(
        filename = function() {
        paste0(Sys.Date(), '_Geographical_Distribution.', input$geodis_out_type)
        },
        content = function(file) {
        if (input$geodis_out_type == "pdf") {
            ggsave(filename = file, plot = Geo_dis_plot(),
            width = as.numeric(input$geodis_down_width),
            height = as.numeric(input$geodis_down_height), dpi = 300,
            device = cairo_pdf)
        } else {
            ggsave(filename = file, plot = Geo_dis_plot(),
            width = as.numeric(input$geodis_down_width),
            height = as.numeric(input$geodis_down_height), dpi = 300)
        }
        }
    )
  
  

  
  
  

  
  ##########---------- 补充函数 ----------##########
  
  # ========================
  # 更新状态卡片（数据输入状态栏）
  # ========================
  output$statusBoxes <- renderUI({
    
    # 辅助函数：生成单个 infoBox
    # subtitle 仅在有激活信息时显示基因型、注释
    make_box <- function(title, is_ready, icon_name, active_label = NULL) {
      infoBox(
        title = title,
        value = if (is_ready) "Imported" else "Waiting Import",
        subtitle = if (!is.null(active_label) && is_ready)
          paste0("Current Activated：", active_label)
        else NULL,
        icon  = icon(icon_name),
        color = if (is_ready) "green" else "navy",
        fill  = TRUE,
        width = 3
      )
    }
    
    tagList(
      make_box("Genotype(need)",!is.null(active_genotype_type()),"dna",active_genotype_type()),
      make_box("Annotation",!is.null(active_annotation_type()), "file-code",active_annotation_type()),
      make_box("Phenotype",active_phenotype(),"file-medical", NULL),
      make_box("Accession",active_accession(),"address-book", NULL)
    )
  })
  
}
shinyApp(ui = ui, server = server)
