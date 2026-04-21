library(shiny)
options(shiny.maxRequestSize = 100000 * 1024^2)

# ========================
# 工具函数
# ========================

get_tools <- function() {
  sys <- Sys.info()[["sysname"]]
  if (grepl("Windows", sys, ignore.case = TRUE)) {
    R.base <- R.home()
    bcftools <- file.path(
      R.base,
      "extract_software",
      "cygwin64",
      "usr",
      "local",
      "bin",
      "bcftools.exe"
    )
    tabix <- file.path(
      R.base,
      "extract_software",
      "cygwin64",
      "usr",
      "local",
      "bin",
      "tabix.exe"
    )
    if (!file.exists(bcftools)) {
      stop("bcftools 未找到")
    }
    if (!file.exists(tabix)) {
      stop("tabix 未找到")
    }

    return(list(
      bcftools = bcftools,
      tabix = tabix,
      is_windows = TRUE
    ))
  } else {
    return(list(
      bcftools = "bcftools",
      tabix = "tabix",
      is_windows = FALSE
    ))
  }
}

to_cygwin_path <- function(path) {
  path <- normalizePath(path, winslash = "/")
  path <- sub("^([A-Za-z]):", "/cygdrive/\\L\\1", path, perl = TRUE)
  return(path)
}

# ========================
# 自动 index
# ========================

# ensure_index <- function(tools, vcf) {
#   index_file <- paste0(vcf, ".tbi")
#
#   if (file.exists(index_file)) return(TRUE)
#   vcf_use <- vcf
#   if (tools$is_windows) vcf_use <- to_cygwin_path(vcf)
#
#   system2(
#     tools$tabix,
#     args = c("-p", "vcf", vcf_use),
#     stdout = TRUE,
#     stderr = TRUE
#   )
#
#   if (!file.exists(index_file)) {
#     stop("index 创建失败")
#   }
# }

ensure_index <- function(tools, vcf, threads = 4) {
  index_file <- paste0(vcf, ".tbi")

  if (file.exists(index_file)) {
    return(TRUE)
  }

  vcf_use <- vcf
  if (tools$is_windows) {
    vcf_use <- to_cygwin_path(vcf)
  }

  system2(
    tools$bcftools,
    args = c(
      "index",
      "--tbi", # 生成 .tbi 格式索引，与原 tabix 一致
      "--threads",
      threads,
      vcf_use
    ),
    stdout = TRUE,
    stderr = TRUE
  )

  if (!file.exists(index_file)) {
    stop("index 创建失败")
  }
}

# ========================
# 读取 contig
# ========================

get_contigs <- function(tools, vcf) {
  vcf_use <- vcf
  if (tools$is_windows) {
    vcf_use <- to_cygwin_path(vcf)
  }

  res <- tryCatch(
    {
      system2(tools$bcftools, c("view", "-h", vcf_use), stdout = TRUE)
    },
    error = function(e) NULL
  )

  if (is.null(res)) {
    return(NULL)
  }

  contigs <- grep("^##contig=", res, value = TRUE)
  if (length(contigs) == 0) {
    return(NULL)
  }

  sub('.*ID=([^,>]+).*', '\\1', contigs)
}

# ========================
# UI
# ========================

ui <- fluidPage(
  titlePanel("bcftools Shiny 工具"),

  sidebarLayout(
    sidebarPanel(
      fileInput("large_vcf", "上传VCF.gz文件"),
      numericInput(
        "large_vcf_threads",
        "线程数",
        value = 4,
        max = 20,
        min = 1,
        step = 1
      ),
      selectInput("large_vcf_chr", "染色体", choices = NULL),
      numericInput("large_vcf_start", "Start", 1),
      numericInput("large_vcf_end", "End", 10000),
      checkboxInput("large_vcf_header", "保留header", TRUE),
      fluidRow(
        column(width = 6, actionButton("large_vcf_run", "运行")),
        column(width = 6, actionButton("large_vcf_clean", "清理临时文件"))
      )
    ),

    mainPanel(
      verbatimTextOutput("large_vcf_log"),
      downloadButton("large_vcf_download", "下载结果")
    )
  )
)

# ========================
# SERVER
# ========================

server <- function(input, output, session) {
  tools <- get_tools()

  log_val <- reactiveVal("")
  result_file <- reactiveVal(NULL)

  append_log <- function(txt) {
    log_val(paste0(log_val(), txt, "\n"))
  }

  # ===== 读取 contig =====
  observeEvent(input$large_vcf, {
    req(input$large_vcf)

    append_log("读取 contig...")

    contigs <- get_contigs(tools, input$large_vcf$datapath)

    if (!is.null(contigs)) {
      updateSelectInput(
        session,
        "large_vcf_chr",
        choices = contigs,
        selected = contigs[1]
      )
      append_log(paste("检测到染色体:", paste(head(contigs), collapse = ", ")))
    } else {
      append_log("未检测到 contig")
    }
  })

  # ===== 运行 =====
  observeEvent(input$large_vcf_run, {
    req(input$large_vcf, input$large_vcf_chr)
    # 锁定按钮
    shinyjs::disable("large_vcf_run")
    on.exit(shinyjs::enable("large_vcf_run")) # 无论成功失败都最终解锁

    infile <- input$large_vcf$datapath
    outfile <- tempfile(fileext = ".vcf")

    infile_use <- infile
    outfile_use <- outfile

    if (tools$is_windows) {
      infile_use <- to_cygwin_path(infile)
      outfile_use <- to_cygwin_path(outfile)
    }

    region <- paste0(
      input$large_vcf_chr,
      ":",
      input$large_vcf_start,
      "-",
      input$large_vcf_end
    )

    append_log(paste("Region:", region))
    append_log("检查 index...")

    # 自动 index
    # tryCatch({
    #   ensure_index(tools, infile)
    #   append_log("index OK")
    # }, error = function(e) {
    #   append_log(paste("index 失败:", e$message))
    #   return(NULL)
    # })

    tryCatch(
      {
        ensure_index(tools, infile, threads = input$large_vcf_threads)
        append_log("index OK")
      },
      error = function(e) {
        append_log(paste("index 失败:", e$message))
        return(NULL)
      }
    )

    # 构建参数
    if (input$large_vcf_header) {
      args <- c("view", "-r", region, "-o", outfile_use, infile_use)
    } else {
      args <- c("view", "-H", "-r", region, "-o", outfile_use, infile_use)
    }

    append_log("运行 bcftools...")

    res <- tryCatch(
      {
        system2(tools$bcftools, args = args, stdout = TRUE, stderr = TRUE)
      },
      error = function(e) {
        paste("ERROR:", e$message)
      }
    )

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
}

shinyApp(ui, server)
