library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(DT)

library(readxl) # to load compound Excel sheet
library(dplyr) # handle dataframes
library(tibble)
library(stringr) # handle strings
library(plotly)

#****************************************************************************************************#
# read files and store them in a list of dataframes #
read.files <- function(files) {
  df.ls <- lapply(files, function (f) {
    df <- read.csv(file = f, sep = '\t', skip = 1, stringsAsFactors = FALSE, fileEncoding = 'latin1') %>%
      select(c(Sample, Hit.1.Name, R.T...min., Quant.Masses, Tailing.Factor, Area)) %>%
      mutate('Date' = as.Date(unlist(strsplit(Sample, '_'))[1], format = '%y%m%d'),
             'QCLevel' = unlist(strsplit(Sample, '_'))[2]) %>%
      select(-Sample)
    colnames(df) <- c('HitName', 'RT', 'Mass', 'TF', 'Area', 'Date', 'QCLevel')
    return(df)
  })
  names(df.ls) <- sapply(df.ls, function(df) {
    paste(df$Date[1], df$QCLevel[1], sep = '_')
  })
  return(df.ls)
}

# read saved data #
read.data <- function() {
  data <- read.csv2('saved_data.csv',
                    colClasses = c('character', 'character', 'numeric', 'integer',
                                   'character', 'numeric', 'character', 'numeric',
                                   'numeric', 'character', 'character', 'numeric',
                                   'numeric', 'Date', 'character')
                    )
  return(data)
}

# not in #
'%!in%' <- function(x,y) {
  !('%in%'(x,y))
}

# remove data that has already been uploaded #
check.ls <- function(data.ls, data) {
  samples <- unique(paste(data$Date, data$QCLevel, sep = '_'))
  return(data.ls[which(names(data.ls) %!in% samples)])
}

# include exceptions for regex #
to.regex.expr <- function(input.str) {
  unwanted.char <- c('(', ')', '[', ']')
  output.str <- input.str
  for (u.c in unwanted.char) {
    output.str <- str_replace_all(output.str,
                                  paste0('\\', u.c),
                                  paste0('\\\\', u.c))
  }
  return(output.str)
}

# print difference with positive/negative sign #
diff.str <- function(diff, digits) {
  sapply(diff, function(d) {
    if (is.na(d)) {
      return(NA)
    } else {
      d <- round(d, digits)
      if (d > 0) {
        return(paste0('+', d))
      } else {
        return(paste(d))
      }
    }
  }) %>% unlist()
}

# match copounds from compounds excel sheet with raw data #
match.compounds <- function(data.raw, compounds.df) {
  data.matched <- lapply(data.raw, function(df) {
    matched.df <- lapply(rownames(compounds.df), function(compound) {
      # find compounds in raw data using Hitnames #
      HitNames <- str_trim(unlist(strsplit(compounds.df[compound, 'HitNames'], '\n')))
      pattern <- paste0('(',
                       HitNames %>% to.regex.expr() %>% paste(collapse = '|'),
                       ')')
      match <- grep(pattern, df$HitName, ignore.case = TRUE)
      HitNames.df <- df[match, ] %>% mutate(Compound = compound, .before = HitName)
      # assign match to corresponding compound #
      compound.grp <- str_replace(compound, '[\\s]*(1st|2nd|3rd|[0-9]+th)', '')
      if (compound.grp %in% compounds.df$grp) {
        compounds.grp.ls <- compounds.df %>%
          filter(grp == compound.grp) %>%
          select(RT) %>%
          rownames_to_column('Compound') %>%
          deframe()
        keep.row <- sapply(rownames(HitNames.df), function(row) {
          RT.diffs <- abs(compounds.grp.ls - HitNames.df[row, 'RT'])
          return(compound == names(RT.diffs[which(RT.diffs == min(RT.diffs))]))
        })
        HitNames.df <- HitNames.df %>%
          filter(isTRUE(keep.row))
      } 
      Hits <- nrow(HitNames.df)
      if (Hits == 0) {
        HitNames.df <- data.frame(Compound = compound,
                                  HitName = NA,
                                  RT = NA,
                                  Date = df$Date[1], 
                                  QCLevel = df$QCLevel[1])
      }
      return(HitNames.df %>% mutate(Hits = Hits))
    }) %>% bind_rows() %>%
      select(c(Compound, HitName, RT, Mass, TF, Area, Date, QCLevel, Hits))
    
    
    # add some extra columns #
    matched.df <- matched.df %>%
      mutate(Calc.RT = sapply(rownames(matched.df), function(row) {
        compound <- matched.df[row, 'Compound']
        this.RT <- matched.df[row, 'RT']
        return(if_else(is.na(this.RT), compounds.df[compound, 'RT'], this.RT))
        })) %>%
      arrange(Calc.RT) %>% 
      mutate(Ref.Rank = compounds.df[Compound, 'Rank'],
             Rank = seq(nrow(matched.df)),
             Rank.Diff = diff.str(Rank - Ref.Rank, 0),
             Ref.Hits = compounds.df[Compound, 'NumberOfMatches'],
             Hits.Diff = diff.str(Hits - Ref.Hits, 0),
             Ref.RT = round(compounds.df[Compound, 'RT'], 3),
             RT.Diff = diff.str(RT - Ref.RT, 3),
             RT = round(RT, 3),
             .after = HitName) %>%
      relocate(RT, .after = Ref.RT) %>%
      select(-c(Hits, Calc.RT))
    return(matched.df)
  }) %>% bind_rows()
  # Column Names:
  # Compound | HitName | Ref.Rank | Rank | Rak.Diff | Ref.Hits | Hits.Diff | 
  # Ref.RT | RT | RT.Diff | Mass | TF | Area | Date | QCLevel
  return(data.matched)
}

# update the data dataframe #
update.data <- function(data.matched, data) {
  if (is.null(data)) {
    return(data.matched)
  } else {
    return(data %>% bind_rows(data.matched))
  }
}

# update the overview dataframe #
overview.calc <- function(data) {
  dates <- unique(data$Date) %>% sort() %>% as.character()
  overview.df <- lapply(dates, function(date) {
    row <- data.frame(matrix(nrow = 1, ncol = 6))
    colnames(row) <- paste0('eQC', seq(6))
    row[1,] <- sapply(colnames(row), function(level) {
      data %>%
        filter(Date == date,
               QCLevel == level,
               !is.na(Area)) %>%
        nrow()
    })
    row[which(row == 0)] <- NA
    return(row)
  }) %>% bind_rows()
  rownames(overview.df) <- dates
  return(overview.df)
}

# get QC levels for which data has been uploaded #
get.levels <- function(overview.df, date) {
  uploaded <- which(!is.na(overview.df[as.character(date), ]))
  QCLevels <- colnames(overview.df)[uploaded]
  return(QCLevels)
}

# generate data for PCA plot #
pca.calc <- function(unique.data, overview.df) {
  lapply(unique(unique.data$Date), function(date) {
    lapply(get.levels(overview.df, date), function(level) {
      tmp <- data.frame(matrix(ncol = length(unique(unique.data$Compound))+1))
      tmp[1, 1] <- level
      row <- sapply(unique(unique.data$Compound), function(compound) {
        unique.data %>% filter(Date == date,
                               QCLevel == level,
                               Compound == compound) %>%
          select(Area)
      })
      row[is.na(row)] <- 0
      tmp[1, 2:ncol(tmp)] <- row
      colnames(tmp) <- c('QCLevel', unique(unique.data$Compound))
      return(tmp)
    }) %>% bind_rows() %>%
      mutate(Date = date, .before = QCLevel)
  }) %>% bind_rows()
}

# generate data for regression graph #
reg.calc <- function(data, compounds, date) {
  lapply(paste0('eQC', seq(6)), function(level) {
    tmp <- data %>%
      filter(QCLevel == level,
             Date == date,
             Compound %in% compounds,
             !is.na(HitName)) %>%
      select(Area, Compound)
    out.df <- data.frame(QC.amount = data.frame(QC.amount = c(10, 25, 50, 100, 250, 500),
                                           row.names = paste0('eQC', seq(6)))[level,],
                         Area.sum = sum(tmp$Area))
    for (compound in compounds) {
      out.df[, compound] <- tmp %>%
        filter(Compound == compound) %>%
        nrow()
    }
    return(out.df)
  }) %>% bind_rows() %>%
    mutate(QCLevel = paste0('eQC', seq(6)), .before = QC.amount)
}

# generate data for control cards
CC.calc <- function(data, compounds, level) {
  lapply(sort(unique(data$Date)), function(date) {
    tmp <- data %>%
      filter(QCLevel == level,
             Date == date,
             Compound %in% compounds,
             !is.na(HitName)) %>%
      select(RT, TF, Area, Compound)
    if (nrow(tmp) == 0) {
      rt.mean <- NA
      rt.sd <- NA
      tf.mean <- NA
      tf.sd <- NA
      area.sum <- NA
    } else {
      rt.mean = mean(tmp$RT)
      tf.mean = mean(tmp$TF)
      area.sum = sum(tmp$Area)
      if (nrow(tmp) == 1) {
        rt.sd <- NA
        tf.sd <- NA
      } else {
        rt.sd = sd(tmp$RT)
        tf.sd = sd(tmp$TF)
      }
    }
    out.df <- data.frame(Date = date,
                         RT.mean = rt.mean,
                         RT.sd = rt.sd,
                         TF.mean = tf.mean,
                         TF.sd = tf.sd,
                         Area.sum = area.sum)
    for (compound in compounds) {
      out.df[, compound] <- tmp %>%
        filter(Compound == compound) %>%
        nrow()
    }
    return(out.df)
  }) %>% bind_rows()
}

# generate data for ratio control cards #
ratios.calc <- function(unique.data, level,
                        r1u, r1d,
                        r2u, r2d,
                        r3u, r3d) {
  lapply(sort(unique(unique.data$Date)), function(date) {
    tmp <- unique.data %>%
      filter(QCLevel == level,
             Date == date) %>%
      select(Compound, Area, TF)
    
    if (nrow(tmp) > 0) {
      ratio1.area <- tmp %>% filter(Compound == r1u) %>% select(Area) %>% deframe() /
        tmp %>% filter(Compound == r1d) %>% select(Area) %>% deframe()
      
      tario2.area <- tmp %>% filter(Compound == r2u) %>% select(Area) %>% deframe() /
        tmp %>% filter(Compound == r2d) %>% select(Area) %>% deframe()
      
      ratio3.area <- tmp %>% filter(Compound == r3u) %>% select(Area) %>% deframe() /
        tmp %>% filter(Compound == r3d) %>% select(Area) %>% deframe()
      
      ratio1.tf <- tmp %>% filter(Compound == r1u) %>% select(TF) %>% deframe() /
        tmp %>% filter(Compound == r1d) %>% select(TF) %>% deframe()
      
      tario2.tf <- tmp %>% filter(Compound == r2u) %>% select(TF) %>% deframe() /
        tmp %>% filter(Compound == r2d) %>% select(TF) %>% deframe()
      
      ratio3.tf <- tmp %>% filter(Compound == r3u) %>% select(TF) %>% deframe() /
        tmp %>% filter(Compound == r3d) %>% select(TF) %>% deframe()
    } else {
      ratio1.area <- NA
      tario2.area <- NA
      ratio3.area <- NA
      ratio1.tf <- NA
      tario2.tf <- NA
      ratio3.tf <- NA
    }
    
    out.df <- data.frame(Date = date,
                         Ratio1.Area = ratio1.area,
                         Ratio2.Area = tario2.area,
                         Ratio3.Area = ratio3.area,
                         Ratio1.TF = ratio1.tf,
                         Ratio2.TF = tario2.tf,
                         Ratio3.TF = ratio3.tf)

    return(out.df)
  }) %>% bind_rows()
}

# ratio control cards #
ratio.CC <- function(ratio.CC.data, dates.na, first.date, last.date, title) {
  mean.ratio <- mean(ratio.CC.data$Values, na.rm = TRUE)
  sd.ratio <- sd(ratio.CC.data$Values, na.rm = TRUE)
  print(ratio.CC.data)
  print(dates.na)
  p <- plot_ly(data = ratio.CC.data,
               x = ~Date, y = ~Values,
               type = 'scatter', mode = 'lines+markers', color = I('black'),
               name = 'Measured')
  if (length(dates.na) > 0) {
    p <- p %>% 
      add_markers(x = dates.na, y = rep(mean.ratio, length(dates.na)),
                  marker = list(symbol = 'x-open', size = 8),
                  error_y = list(array = rep(NA, length(dates.na))),
                  name = 'Not detected')
  }
  p <- p %>%
    layout(yaxis = list(range = get.min.max(values = ratio.CC.data$Values),
                        zeroline = FALSE, title = title),
           xaxis = list(range = c(first.date-3, last.date+3),
                        zeroline = FALSE),
           shapes = CC.lines(mean.ratio,
                             sd.ratio,
                             NA))
  return(p)
}

# generate data where each compound at a certain date and QC level only has a single value
unique.data.calc <- function(data) {
  multi.matches <- data %>%
    filter(as.numeric(Ref.Hits) + as.numeric(Hits.Diff) > 1) %>%
    select(Compound, Date, QCLevel) %>%
    unique()
  multi.match.df <- lapply(rownames(multi.matches), function(row) {
    tmp <- data %>% 
      filter(Compound == multi.matches[row, 'Compound'],
             Date == multi.matches[row, 'Date'],
             QCLevel == multi.matches[row, 'QCLevel']) %>%
      select(RT, TF, Area)
    return(data.frame(Compound = multi.matches[row, 'Compound'],
                      Date = multi.matches[row, 'Date'],
                      QCLevel = multi.matches[row, 'QCLevel'],
                      RT = mean(tmp$RT, na.rm = TRUE),
                      TF = mean(tmp$TF, na.rm = TRUE),
                      Area = sum(tmp$Area, na.rm = TRUE)))
  }) %>% bind_rows()
  unique.data <- data %>%
    filter(!as.numeric(Ref.Hits) + as.numeric(Hits.Diff) > 1) %>%
    select(Compound, Date, QCLevel, RT, TF, Area) %>%
    bind_rows(multi.match.df) %>%
    arrange(Date, QCLevel, Compound)
  return(unique.data)
}

# generate a dataframe that contains all mean and sd values for all compounds and QC level
alert.calc <- function(unique.data) {
  # iterate over all QC levels and compounds
  lapply(paste0('eQC', seq(6)), function(level) {
    level.df <- unique.data %>% filter(QCLevel == level)
    if (nrow(level.df) == 0) {
      return(NULL)
    }
    dates <- unique(level.df$Date) %>% sort() %>% as.character()
    out.df <- lapply(unique(unique.data$Compound), function(compound) {
      tmp <- level.df %>%
        filter(Compound == compound,
               !is.na(RT)) %>%
        select(RT, Area)
      if (nrow(tmp) >= 3) {
        RT.mean <- mean(tmp$RT)
        RT.sd <- sd(tmp$RT)
        Area.mean <- mean(tmp$Area)
        Area.sd <- sd(tmp$Area)
      } else {
        RT.sd <- NA
        Area.sd <- NA
        if (nrow(tmp) > 0) {
          RT.mean <- mean(tmp$RT)
          Area.mean <- mean(tmp$Area)
        } else {
          RT.mean <- NA
          Area.mean <- NA
        }
      }
      return(data.frame(Compound = compound,
                        QCLevel = level,
                        RT.mean = RT.mean,
                        RT.sd = RT.sd,
                        Area.mean = Area.mean,
                        Area.sd = Area.sd))
    }) %>% bind_rows() %>%
      arrange(QCLevel, Compound) %>% # should be sorted but be sure
      sapply(rep.int, times = length(dates)) %>%
      as.data.frame() %>%
      mutate(RT.mean.diff = level.df$RT - as.numeric(RT.mean),
             RT.sd = as.numeric(RT.sd),
             Area.mean.diff = level.df$Area - as.numeric(Area.mean),
             Area.sd = as.numeric(Area.sd),
             Date = lapply(dates, function(date) {
               rep(date, length(unique(level.df$Compound)))
             }) %>% unlist()) %>%
      select(c(QCLevel, Compound, Date, RT.mean.diff, RT.sd, Area.mean.diff, Area.sd)) %>%
      arrange(QCLevel, Compound, Date)
    return(out.df)
  }) %>% bind_rows()
}

# generate a dataframe of Westgard rule violations
# 1:2s, 1:3s, 2:2s, R:4s, 4:1s, 10x
alert <- function(alert.data) {
  alert.list <- lapply(unique(alert.data$Compound), function(compound) {
    lapply(unique(alert.data$QCLevel), function(level) {
      filter(alert.data,
             QCLevel == level,
             Compound == compound)
    })
  }) %>%
    unlist(recursive = FALSE)
  
  alert.df <- lapply(c('RT', 'Area'), function(parameter){
    mean.diff <- paste(parameter, 'mean.diff', sep = '.')
    sd <- paste(parameter, 'sd', sep = '.')
    out <- lapply(alert.list, function(df) {
      if (df %>% filter(!is.na(df[,sd])) %>% nrow() > 0) {
        df <- filter(df, !is.na(df[,mean.diff]))
        # 1:3s #
        error_1.3s <- df %>% filter(abs(df[,mean.diff]) > 3*df[,sd]) %>% nrow() > 0 # any measurement 3sd diff from mean
        # 1:2s # 
        error_1.2s <- df %>% filter(abs(df[,mean.diff]) > 2*df[,sd]) %>% nrow() > 0 # any measurement 2sd diff from mean
        
        if (error_1.2s) {
          # 2:2s #
          result <- rle(df[,mean.diff] > 2*df[,sd])
          error_2.2s <- any(result$lengths >= 2 & result$values == TRUE) # 2 consecutive measurements 2sd above mean
          if (!error_2.2s) {
            result <- rle(df[,mean.diff] < -2*df[,sd])
            error_2.2s <- any(result$lengths >= 2 & result$values == TRUE) # 2 consecutive measurements 2sd below mean
          }
          # R:4s #
          error_R.4s <- max(df[,mean.diff], na.rm = TRUE) - min(df[,mean.diff], na.rm = TRUE) > 4*df[1,sd] # any 2 measurements with 4sd diff
        } else {
          error_2.2s <- FALSE
          error_R.4s <- FALSE
        }
        
        # 4:1s #
        result <- rle(df[,mean.diff] > df[,sd])
        error_4.1s <- any(result$lengths >= 4 & result$values == TRUE) # 4 consecutive measurements 1sd above mean
        if (!error_4.1s) {
          result <- rle(df[,mean.diff] < -df[,sd])
          error_4.1s <- any(result$lengths >= 4 & result$values == TRUE) # 4 consecutive measurements 1sd below mean
        }
        
      } else {
        error_1.3s <- FALSE
        error_1.2s <- FALSE
        error_2.2s <- FALSE
        error_R.4s <- FALSE
        error_4.1s <- FALSE
      }
      
      # 10x #
      if (nrow(df) > 10) {
        result <- rle(df[,mean.diff] > 0)
        error_10x <- any(result$lengths >= 10 & result$values == TRUE) # 10 consecutive measurements above mean
        if (!error_10x) {
          result <- rle(df[,mean.diff] < df[,sd])
          error_10x <- any(result$lengths >= 10 & result$values == TRUE) # 10 consecutive measurements below mean
        }
      } else {
        error_10x <- FALSE
      }
      return(data.frame(Parameter = parameter,
                        Compound = df$Compound[1],
                        QCLevel = df$QCLevel[1],
                        Error_1.2s = error_1.2s,
                        Error_1.3s = error_1.3s,
                        Error_2.2s = error_2.2s,
                        Error_R.4s = error_R.4s,
                        Error_4.1s = error_4.1s,
                        Error_10x = error_10x))
    }) %>% bind_rows()
    return(out)
  }) %>% bind_rows()
  
  return(alert.df)
}

# percentage of matches per QC Level #
match.freq.calc <- function(unique.data) {
  lapply(paste0('eQC', seq(6)), function(level) {
    df <- data.frame(QCLevel = sapply(unique(unique.data$Compound), function(compound) {
      tmp <- unique.data %>% filter(Compound == compound,
                             QCLevel == level)
      freq <- ifelse(nrow(tmp) > 0,
                     sum(!is.na(tmp$RT)) / nrow(tmp) * 100,
                     0)
      return(freq)
    }))
    colnames(df) <- level
    return(df)
  }) %>% bind_cols(data.frame(Compounds = unique(unique.data$Compound)))
}

# scale of control cards # 
get.min.max <- function(values, errors = 0, ref = NA) {
  mean <- mean(values, na.rm = TRUE)
  CLs <- 3.5 * sd(values, na.rm = TRUE)
  if (is.na(ref)) {ref <- mean}
  ymin <- min(c(mean-CLs,
                values - errors,
                ref),
              na.rm = TRUE)
  ymax <- max(c(mean+CLs,
                values + errors,
                ref),
              na.rm = TRUE)
  if (ymin == ref) {
    ymin = ymin - 0.05*(ymax - ymin)
  } else if (ymax == ref) {
    ymax = ymax + 0.05*(ymax - ymin)
  }
  return(c(ymin, ymax))
}

# draw horizintal lines #
hline <- function(y, color = 'red', dash = 'dot', width = 1) {
  list(type='line',
       x0 = as.Date('2000-01-01'), x1 = as.Date('2100-01-01'),
       y0= y, y1 = y,
       line=list(color = color,
                 dash = dash,
                 width = width))
}

# Title, Y-Axis, UCL, UWL, Mean, LWL, LCL #
CC.lines <- function(m, sd, ref) {
  hline.ls <- list(
    hline(m, color = 'green', width = 2) # Mean
  )
  if (!is.na(sd)) {
    hline.ls[[2]] <- hline(m+3*sd, width = 2) # UCL
    hline.ls[[3]] <- hline(m+2*sd) # UWL
    hline.ls[[4]] <- hline(m-2*sd) # LWL
    hline.ls[[5]] <- hline(m-3*sd, width = 2) # LCL
  }
  if (!is.na(ref)) {
    hline.ls[[length(hline.ls)+1]] <- hline(ref, color = 'gray', dash = 'line') # Reference
  }
  return(hline.ls)
}

# get the position of a compound in the compounds.df
get.name<- function(pattern, compounds.df, sort.by) {
  if (sort.by == 'Name') {
    choices <- compounds.df %>% rownames() %>% sort()
  } else if (sort.by == 'Retention time') {
    choices <- compounds.df %>% arrange(RT) %>% rownames()
  }
  index <- grep(pattern, choices, ignore.case = TRUE)
  index <- ifelse(length(index) == 1, index, 1)
  return(choices[index])
}

# Compound Picker #
UIcompoundPicker <- function(sort.by, compounds.df, id, title = '',
                             selected = NULL, multiple = FALSE) {
  if (sort.by == 'Name') {
    choices <- compounds.df %>% rownames() %>% sort()
  } else if (sort.by == 'Retention time') {
    choices <- compounds.df %>% arrange(RT) %>% rownames()
  }
  return(
    pickerInput(
      inputId = id,
      label = title,
      choices = choices,
      selected = selected,
      options = list(
        `actions-box` = TRUE,
        size = 10,
        `selected-text-format` = 'count > 2'
      ),
      multiple = multiple
    )
  )
}

#****************************************************************************************************#
server <- function(input, output, session) {
  
  rv <- reactiveValues(read.saved.data = file.exists('saved_data.csv'))

  # Help --------------------------------------------------------------------
  observeEvent(input$help, {
    output$help.text <- renderUI({
      if (input$help) {
        return(
          fluidRow(
            column(12,
                   HTML('<b>Steps:</b><br>
                        <ol>
                          <li>Upload the Excel file containg all compounds including their hitnames, retention time etc.</li>
                          <li>Uplaod new Quality Control files and add them for calculation with the \"Add files\" button.
                              Processing of all added files can be started with the \"Start Processing\" button.
                              Tipp: Upload all QC data at once. Thereby, the calculation process doesn\'t have to be repeated for every upload.
                              If data from previous analysis has been saved it will automatically be read in when the procesing starts.</li>
                          <li>An overview table will appear after the calculation has finished. Each cell contains the number of matches for each QC sample.
                              This table can additionally be used to directly choose specific QC samples but also dates and QC levels depending on the chosen analysis.</li>
                          <li>All visualisation will be shown in the second box. Each has its own tab that can be switched between manually.
                              Some tabs also require to choose from the compound list. Multiple selection is possible, retention time and tailing factor are averaged while the area is added up.</li>
                        </ol>
                        The requirements for each tab are shown in the table below:'),
                   tags$hr(style='border-color: grey;'),
                   renderTable(data.frame('Tab' = c('PCA',
                                                    'Alert',
                                                    'Match Frequency',
                                                    'Matches',
                                                    'Regression',
                                                    'Control Cards',
                                                    'Ratios',
                                                    'Delete'),
                                          'Explanation' = c('Calculate a PCA where each sample is an individual point with PC1 and PC2 as its 2-dim coordinates.',
                                                            'Calculate all violations of the quality control Westgard rules to identify inaccuracy and imprecision.',
                                                            'Visualise the match frequency of individual compounds. Compounds that are matched more than expected are listed additionally.',
                                                            'List all matches of a single sample in a detailed table. Use the button "Copy", "Excel" and "PDF" to export.',
                                                            'Plot a fitted linear regression model through all six QC levels.',
                                                            'Display control cards of retention time, tailing factor and area.
                                                             The green line represents the mean with two and three standard deviations above and below shown as red lines.
                                                             The grey line at the retention time control cards indicates the given reference retention time.',
                                                            'Show control cards of important ratios.',
                                                            'Remove unwanted samples from the uploaded data.'),
                                          'Input' = c('None',
                                                      'None',
                                                      'None',
                                                      'Date & QC Level',
                                                      'Date & Compound',
                                                      'QC Level & Compound',
                                                      'QC Level',
                                                      'Date & QC Level')))
            )
          )
        )
      }
    })
  })
  
  # Upload Box --------------------------------------------------------------

  # Load compounds #
  observeEvent(input$compounds.file, {
    out.df <- read_excel(input$compounds.file$datapath, sheet = 'Compounds', skip = 2)
    out.df <- out.df[complete.cases(out.df),] %>%
      arrange(Rank) %>%
      mutate(grp = str_replace(Compound, '[\\s]*(1st|2nd|3rd|[0-9]+th)', '')) %>%
      group_by(grp)
    groups <- out.df %>%
      tally() %>%
      as.data.frame() %>%
      filter(n > 1) %>%
      select(grp) %>%
      unlist()
    rv$compounds.df <- out.df %>%
      as.data.frame() %>%
      mutate(grp = sapply(out.df$grp, function(g) {if_else(g %in% groups, g, NULL)})) %>%
      column_to_rownames('Compound')
  })
  
  # Selection of compound #
  output$compoundPicker <- renderUI({
    if (is.null(rv$compounds.df)) {
      return(
        HTML('<b> &#9888; Please upload a <br> compound Excel file! </b>')
      )
    } else {
      return(
        UIcompoundPicker(sort.by = input$sort,
                              compounds.df = rv$compounds.df,
                              id = 'compoundPicker',
                              title = 'Select compounds',
                              multiple = TRUE)
      )
    }
  })
  
  # Add files #
  observeEvent(input$add, {
    new.files <- input$files$datapath
    names(new.files) <- input$files$name
    if (isTruthy(rv$files)) {
      tmp.files <- c(rv$files, new.files)
    } else {
      tmp.files <- new.files
    }
    rv$files <- rev(tmp.files)[unique(names(tmp.files))]
    output$uploadedFiles <- renderUI({
      tags$b(paste0('Number of uploaded files: ', length(rv$files)))
    })
  })
  
  # Read and match new files #
  observeEvent(input$calc, {
    # require compound.df #
    if (!isTruthy(rv$compounds.df)) {
      showNotification('Please upload a compound Excel file!', type = 'warning', duration = 3)
    }
    req(rv$compounds.df)
    read.QC.files <- isTruthy(rv$files)
    # require files to calculate #
    if (!rv$read.saved.data & !read.QC.files) {
      showNotification('Please enter some files!', type = 'warning', duration = 3)
    } else {
      steps <- 1 # creating overview
      if (rv$read.saved.data) {
        steps <- steps + 1 # reading saved data
      }
      if (read.QC.files) {
        steps <- steps + 3 # reading QC files; matching compounds; updating data
      }
      new.data <- FALSE
      
      withProgress(message = 'Loading data: ', value = 0, {
        if (rv$read.saved.data) {
          # reading saved data #
          incProgress(1/steps, detail = 'Reading saved data')
          rv$data <- read.data()
          Sys.sleep(1) # wait 1 sec to read progress
        }
        
        if (read.QC.files) {
          # reading QC files #
          incProgress(1/steps, detail = 'Reading QC files')
          data.raw <- read.files(rv$files) %>%
            check.ls(rv$data)
          rv$files <- NULL # reset after files have been read
          Sys.sleep(1)
          
          # require new data #
          if (length(data.raw) > 0) {
            # match compounds #
            incProgress(1/steps, detail = 'Matching compounds')
            data.matched <- match.compounds(data.raw, rv$compounds.df)
            new.data <- TRUE
            Sys.sleep(1)
            # update data #
            incProgress(1/steps, detail = 'Updating data')
            rv$data <- update.data(data.matched, rv$data)
            Sys.sleep(1)
          } else {
            incProgress(1/steps, detail = '')
            incProgress(1/steps, detail = '')
          }
        }
        
        if (rv$read.saved.data | new.data) {
          # create overview.df #
          incProgress(1/steps, detail = 'Creating overview')
          rv$overview.df <- overview.calc(rv$data)
          rv$read.saved.data <- FALSE # saved data has already been read
          Sys.sleep(1)
        } else {
          incProgress(1/steps, detail = '')
          showNotification('All uploaded files have already been added.', type = 'warning', duration = 3)
        }
      })
    }
  })
  
  # Save data #
  observeEvent(input$save, {
    if (!isTruthy(rv$data)) {
      showNotification('No data to save!', type = 'warning', duration = 3)
    }
    req(rv$data)
    write.csv2(rv$data, 'saved_data.csv', row.names = FALSE)
    if (file.exists('saved_data.csv')) {
      showNotification('Data has been saved successfully.', type = 'message', duration = 3)
    } else {
      showNotification('Saving failed! Pls try again.', type = 'error', duration = 3)
    }
  })
  
  # Selction of Date and/or QC level #
  # PCA, Alert, Match.Frequency, Matches, Regression, Control.Cards, Ratios, Delete
  observeEvent(c(input$vis.tab, rv$data), {
    if (is.null(rv$overview.df)) {
      rv$selection <- 'none'
    } else {
      if (input$vis.tab %in% c('PCA', 'Alert', 'Match Frequency')) {
        rv$selection = 'none'
      } else if (input$vis.tab %in% c('Matches', 'Delete')) {
        selectable <- matrix(lapply(seq(ncol(rv$overview.df)), function(cols) {
          lapply(seq(nrow(rv$overview.df)), function(rows) {
            c(rows,cols)
          }) %>% unlist()
        }) %>% unlist(),
        ncol = 2,
        byrow = TRUE)[which(!is.na(rv$overview.df)), , drop = FALSE]
        if (input$vis.tab == 'Delete') {
          mode <- 'multiple'
          selected <- NULL
        } else {
          mode <- 'single'
          selected <- rv$matches.selected
        }
        rv$selection = list(mode = mode,
                            target = 'cell',
                            selectable = selectable,
                            selected = selected)
        
      } else if (input$vis.tab == 'Regression') {
        rv$selection = list(mode = 'single',
                            target = 'row',
                            selectable = seq(nrow(rv$overview.df)),
                            selected = rv$regression.selected)
      } else if (input$vis.tab %in% c('Control Cards', 'Ratios')) {
        rv$selection = list(mode = 'single',
                            target = 'column',
                            selectable = as.integer(which(colSums(!is.na(rv$overview.df))>0)),
                            selected = rv$cc.selected)
      }
    }
  })
  
  # Overview dataframe #
  output$overview <- DT::renderDataTable(rv$overview.df,
                                         server = FALSE,
                                         class = 'cell-border stripe',
                                         rownames = TRUE,
                                         options = list(pageLength = -1, info = FALSE,
                                                        lengthMenu = list(c(5, 10, 15, 20, 25, 30, -1),
                                                                          c('5', '10', '15','20',
                                                                            '25', '30', 'All'))),
                                         selection = rv$selection)

  # Ratios compound input ---------------------------------------------------
  
  # Default slected Ratios #
  observeEvent(req(rv$compounds.df), {
    rv$ratio1 = c(get.name('Serine,.3TMS', rv$compounds.df, input$sort),
                  get.name('Serine,.2TMS', rv$compounds.df, input$sort))
    
    rv$ratio2 = c(get.name('Alanine,.2TMS', rv$compounds.df, input$sort),
                  get.name('Valine,.2TMS', rv$compounds.df, input$sort))
    
    rv$ratio3 = c(get.name('Glutamate,.3TMS', rv$compounds.df, input$sort),
                  get.name('Oxoproline,.2TMS', rv$compounds.df, input$sort))
  })
  
  output$compoundPickerRatios <- renderUI({
    req(rv$compounds.df)
    if (input$vis.tab == 'Ratios') {
      return(
        box(width = 12,
            fluidRow(
              column(4,
                     UIcompoundPicker(sort.by = input$sort,
                                      compounds.df = rv$compounds.df,
                                      id = 'compoundPickerR1u',
                                      title = 'Ratio 1',
                                      selected = rv$ratio1[1])
              ),
              column(4,
                     UIcompoundPicker(sort.by = input$sort,
                                      compounds.df = rv$compounds.df,
                                      id = 'compoundPickerR2u',
                                      title = 'Ratio 2',
                                      selected = rv$ratio2[1])
              ),
              column(4,
                     UIcompoundPicker(sort.by = input$sort,
                                      compounds.df = rv$compounds.df,
                                      id = 'compoundPickerR3u',
                                      title = 'Ratio 3',
                                      selected = rv$ratio3[1])
              ),
              column(4,
                     UIcompoundPicker(sort.by = input$sort,
                                      compounds.df = rv$compounds.df,
                                      id = 'compoundPickerR1d',
                                      selected = rv$ratio1[2])
              ),
              column(4,
                     UIcompoundPicker(sort.by = input$sort,
                                      compounds.df = rv$compounds.df,
                                      id = 'compoundPickerR2d',
                                      selected = rv$ratio2[2])
              ),
              column(4,
                     UIcompoundPicker(sort.by = input$sort,
                                      compounds.df = rv$compounds.df,
                                      id = 'compoundPickerR3d',
                                      selected = rv$ratio3[2])
              )
            )
        )
      )
    } else {
      return(NULL)
    }
  })
  
  # update selected ratios #
  observe({
    rv$ratio1 = c(input$compoundPickerR1u, input$compoundPickerR1d)
    rv$ratio2 = c(input$compoundPickerR2u, input$compoundPickerR2d)
    rv$ratio3 = c(input$compoundPickerR3u, input$compoundPickerR3d)
  })

  # Visualisation -----------------------------------------------------------
  
  observeEvent(rv$data, {
    req(rv$data)
    withProgress(message = 'Processing data: ', value = 0, {
      incProgress(1/4, detail = 'Merging multiple matches')
      rv$unique.data <- unique.data.calc(rv$data)
      Sys.sleep(1)
    
      # PCA #
      incProgress(1/4, detail = 'Calculating PCA')
      if (TRUE) {
        PCA.points <- sum(
          sapply(rownames(rv$overview.df), function(date) {
            length(get.levels(rv$overview.df, date))
          })
        )
        req(PCA.points > 1)
        data.pca <- pca.calc(rv$unique.data, rv$overview.df)
        pca <- prcomp(t(data.pca[,3:ncol(data.pca)]), scale = FALSE)
        pca.df <- data.frame(pca$rotation,
                             'Cluster' = factor(data.pca$QCLevel),
                             'Date' = data.pca$Date) %>% select(c(PC1, PC2, Cluster, Date))
  
        output$PCA.plot <- renderPlotly({
          plot_ly(data = pca.df,
                  x = ~PC1, y = ~PC2, text = ~Date,
                  type = 'scatter', mode='markers', color = ~Cluster, marker=list(size=11)) %>%
            layout(xaxis=list(title='PC1', range = c(-1,1)),
                   yaxis=list(title='PC2', range = c(-1,1)))
        })
        
        output$PCA <- renderUI({
          plotlyOutput(outputId = 'PCA.plot')
        })
      }
      
      # Alert #
      incProgress(1/4, detail = 'Checking violations of Westgard rules')
      if (TRUE) {
        alert.data <- alert.calc(rv$unique.data)
        alert.df <- alert(alert.data)
        alert.df <- alert.df %>% filter(rowSums(alert.df[,4:9]) > 0)
        alert.df[alert.df == FALSE] <- NA
        alert.df[alert.df == TRUE] <- 'True'
        
        output$alert.RT <- DT::renderDataTable(alert.df %>% filter(Parameter == 'RT') %>% select(-Parameter),
                                               server = FALSE,
                                               class = 'cell-border stripe',
                                               rownames = FALSE,
                                               selection = 'none',
                                               options = list(scrollX = TRUE,
                                                              pageLength = 5, info = FALSE,
                                                              lengthMenu = list(c(5, 10, 15, 20, 25, -1),
                                                                                c('5', '10', '15', '20', '25', 'All')))
        )
        
        output$alert.Area <- DT::renderDataTable(alert.df %>% filter(Parameter == 'Area') %>% select(-Parameter),
                                                 server = FALSE,
                                                 class = 'cell-border stripe',
                                                 rownames = FALSE,
                                                 selection = 'none',
                                                 options = list(scrollX = TRUE,
                                                                pageLength = 5, info = FALSE,
                                                                lengthMenu = list(c(5, 10, 15, 20, 25, -1),
                                                                                  c('5', '10', '15', '20', '25', 'All')))
        )
        
        output$legend <- renderTable(data.frame('Error' = c('1.2s',
                                                            '1.3s',
                                                            '2.2s',
                                                            'R.4s',
                                                            '4.1s',
                                                            '10x'),
                                                'Explanation' = c('One measurement exceeds 2 standard deviations either above or below the mean of the reference range.',
                                                                  'One measurement exceeds 3 standard deviations either above or below the mean of the reference range.',
                                                                  '2 consecutive measurements exceed 2 standard deviations of the reference range, and on the same side of the mean.',
                                                                  'Two measurements in the same run have a 4 standard deviation difference (such as one exceeding 2 standard deviations above the mean, and another exceeding 2 standard deviations below the mean).',
                                                                  '4 consecutive measurements exceed 1 standard deviation on the same side of the mean.',
                                                                  '10 consecutive measurements are on the same side of the mean.'),
                                                'Suspected' = c('Inaccuracy and/or imprecision',
                                                                'Inaccuracy and/or imprecision',
                                                                'Inaccuracy',
                                                                'Imprecision',
                                                                'Inaccuracy',
                                                                'Inaccuracy')))
        
        output$Alert <- renderUI({
          fluidRow(
            column(12,
                   box(width = 12, status = 'primary', title = tags$b('Retention Time'), 
                       DT::dataTableOutput(outputId = 'alert.RT'))
            ),
            column(12,
                   box(width = 12, status = 'primary', title = tags$b('Area'),
                       DT::dataTableOutput(outputId = 'alert.Area'))
            ),
            column(12,
                   box(width = 12, status = 'primary', title = tags$b('Legend'),
                       tableOutput(outputId = 'legend'))
            )
          )
        })
      }
      
      # Match Frequency #
      incProgress(1/4, detail = 'Calculating match frequency')
      if (TRUE) {
        match.freq.df <- match.freq.calc(rv$unique.data)
        output$Match.Frequency.plot <- renderPlotly({
          plot_ly(match.freq.df, x = ~Compounds, y = ~eQC1, type = 'bar', name = 'eQC1') %>%
            add_trace(y = ~eQC2, name = 'eQC2', visible = 'legendonly') %>%
            add_trace(y = ~eQC3, name = 'eQC3', visible = 'legendonly') %>%
            add_trace(y = ~eQC4, name = 'eQC4', visible = 'legendonly') %>%
            add_trace(y = ~eQC5, name = 'eQC5', visible = 'legendonly') %>%
            add_trace(y = ~eQC6, name = 'eQC6', visible = 'legendonly') %>%
            layout(yaxis = list(title = 'Match frequency [%]', range = c(0,100)),
                   xaxis = list(title = '', tickangle = 45, type='category'),
                   barmode = 'group')
        })
        
        unwanted.matches <- rv$data %>%
          filter(as.numeric(Hits.Diff) > 0) %>%
          select(Compound, HitName, RT, QCLevel, Date) %>%
          arrange(Compound, QCLevel, Date)
        output$Unwanted.Matches <- DT::renderDataTable(unwanted.matches,
                                                       server = FALSE,
                                                       class = 'cell-border stripe',
                                                       rownames = FALSE,
                                                       selection = 'none',
                                                       options = list(scrollX = TRUE,
                                                                      pageLength = 20, info = FALSE,
                                                                      lengthMenu = list(c(10, 20, 30, 40, 50, -1),
                                                                                        c('10', '20', '30', '40', '50', 'All'))))
        
        output$Match.Frequency <- renderUI({
          fluidRow(
            column(12,
                   box(width = 12, status = 'primary',
                       plotlyOutput(outputId = 'Match.Frequency.plot')
                   )
            ),
            column(12,
                   box(width = 12, status = 'primary', title = tags$b('Compounds with more matches than expected'),
                       DT::dataTableOutput(outputId = 'Unwanted.Matches')
                   )
            )
          )
        })
      }
      Sys.sleep(1)
    })
  })

  # Matches #
  observeEvent(input$overview_cells_selected, {
    output$Matches <- NULL
    req(input$overview_cells_selected)
    if (input$vis.tab == 'Matches') {
      rv$matches.selected <- input$overview_cells_selected
      match.df <- rv$data %>%
        filter(Date == rownames(rv$overview.df)[input$overview_cells_selected[1]],
               QCLevel == colnames(rv$overview.df)[input$overview_cells_selected[2]]) %>%
        select(-c(Date, QCLevel))

      output$match.DT <- DT::renderDataTable(match.df,
                                             server = FALSE,
                                             class = 'cell-border stripe',
                                             rownames = FALSE,
                                             selection = 'none',
                                             options = list(scrollX = TRUE,
                                                            pageLength = -1, info = FALSE,
                                                            lengthMenu = list(c(25, 50, 75, 100, -1),
                                                                              c('25', '50', '75', 
                                                                                '100', 'All')),
                                                            dom = 'Bfrtip',
                                                            buttons = c('copy', 'excel', 'pdf')),
                                             extensions = 'Buttons')
      output$Matches <- renderUI({
        DT::dataTableOutput(outputId = 'match.DT')
      })
    }
  })

  # Regression #
  observeEvent(c(input$overview_rows_selected, input$compoundPicker), {
    output$Regression <- NULL
    req(input$overview_rows_selected)
    if (input$vis.tab == 'Regression') {
      rv$regression.selected <- input$overview_rows_selected
      req(input$compoundPicker)
      reg.data <- reg.calc(data = rv$data,
                           compounds = input$compoundPicker,
                           date = rownames(rv$overview.df)[input$overview_rows_selected])
      if (reg.data %>% filter(Area.sum > 0) %>% nrow() < 3) {
        showNotification('Not enough matches to plot linear regression (min. 3 data points)!',
                         type = 'warning', duration = 3)
      } else {
        lm <- reg.data %>% filter(Area.sum > 0) %>% lm(Area.sum ~ QC.amount,.)
        m <- lm$coefficients['QC.amount']
        b <- lm$coefficients['(Intercept)']
        output$Regression.plot <- renderPlotly({
            plot_ly(data = reg.data %>% filter(Area.sum > 0),
                    x = ~QC.amount, y = ~Area.sum,
                    type = 'scatter', mode = 'markers', color = I('black')) %>%
            add_markers(y = ~Area.sum) %>%
            layout(shapes = list(list(type='line',
                                      x0 = 0, x1 = 1000,
                                      y0= b,
                                      y1 = m * 1000 + b,
                                      line=list(color = 'black',
                                                dash = 'line',
                                                width = 1))),
                   showlegend = F,
                   yaxis = list(title = 'Area [#counts]'),
                   xaxis = list(title = 'Amount [ng]', range = c(0,600)))
        })
        
        output$reg.matches <- renderTable(reg.data %>%
                                           select(-c(QC.amount, Area.sum)))
        
        output$Regression <- renderUI({
          fluidRow(
            column(6,
                   box(width = 12, status = 'primary',
                       plotlyOutput(outputId = 'Regression.plot'))
            ),
            column(6,
                   box(width = 12, status = 'primary',
                       title = 'Number of matches',
                       tableOutput(outputId = 'reg.matches'),
                       tags$hr(),
                       HTML(paste0('<h5><b>R<sup>2</sup> value: ', round(summary(lm)$r.squared, 4),
                                  ' with ',
                                  reg.data %>% filter(Area.sum > 0) %>% nrow(), 
                                  '/6 data points</b></h5>')))
            )
          )
        })
      }
    }
  })

  # Control Cards #
  observeEvent(c(input$overview_columns_selected, input$compoundPicker), {
    output$Control.Cards <- NULL
    req(input$overview_columns_selected)
    rv$cc.selected <- input$overview_columns_selected
    if (input$vis.tab %in% c('Control Cards', 'Ratios')) {
      req(input$compoundPicker)
      CC.data <- CC.calc(rv$data,
                         input$compoundPicker,
                         colnames(rv$overview.df)[input$overview_columns_selected])
      
      ref <- ifelse(length(input$compoundPicker) == 1,
                    rv$compounds.df[input$compoundPicker, 'RT'],
                    NA)
      
      dates.na <- CC.data[is.na(CC.data$RT.mean), 'Date']
      first.date <- sort(CC.data$Date)[1]
      last.date <- tail(sort(CC.data$Date),1)
      CC.data <- filter(CC.data, !is.na(RT.mean))
      
      if (nrow(CC.data) == 0) {
        showNotification('Not enough matches to plot Control Cards!',
                         type = 'warning', duration = 3)
      } else {
        # RT #
        mean.rt <- mean(CC.data$RT.mean, na.rm = TRUE)
        sd.rt <- sd(CC.data$RT.mean, na.rm = TRUE)
        p.RT <- plot_ly(data = CC.data,
                     x = ~Date, y = ~RT.mean,
                     type = 'scatter', mode = 'lines+markers', color = I('black'),
                     error_y = ~list(array = RT.sd,
                                     color = '#000000'),
                     name = 'Measured')
        if (length(dates.na) > 0) {
          p.RT <- p.RT %>% 
            add_markers(x = dates.na, y = rep(mean.rt, length(dates.na)),
                        marker = list(symbol = 'x-open', size = 8),
                        error_y = list(array = rep(NA, length(dates.na))),
                        name = 'Not detected')
        }
        p.RT <- p.RT %>%
          layout(yaxis = list(range = get.min.max(values = CC.data$RT.mean,
                                                  errors = CC.data$RT.sd,
                                                  ref = ref),
                              zeroline = FALSE, title = 'Minutes'),
                 xaxis = list(range = c(first.date-3, last.date+3),
                              zeroline = FALSE),
                 shapes = CC.lines(mean.rt,
                                   sd.rt,
                                   ref))
        output$CC.RT <- renderPlotly({p.RT})
        # TF #
        mean.tf <- mean(CC.data$TF.mean, na.rm = TRUE)
        sd.tf <- sd(CC.data$TF.mean, na.rm = TRUE)
        p.TF <- plot_ly(data = CC.data,
                        x = ~Date, y = ~TF.mean,
                        type = 'scatter', mode = 'lines+markers', color = I('black'),
                        error_y = ~list(array = TF.sd,
                                        color = '#000000'),
                        name = 'Measured')
        if (length(dates.na) > 0) {
          p.TF <- p.TF %>% 
            add_markers(x = dates.na, y = rep(mean.tf, length(dates.na)),
                        marker = list(symbol = 'x-open', size = 8),
                        error_y = list(array = rep(NA, length(dates.na))),
                        name = 'Not detected')
        }
        p.TF <- p.TF %>%
          layout(yaxis = list(range = get.min.max(values = CC.data$TF.mean,
                                                  errors = CC.data$TF.sd),
                              zeroline = FALSE, title = 'TF'),
                 xaxis = list(range = c(first.date-3, last.date+3),
                              zeroline = FALSE),
                 shapes = CC.lines(mean.tf,
                                   sd.tf,
                                   NA))
        output$CC.TF <- renderPlotly({p.TF})
        # Area #
        mean.area <- mean(CC.data$Area.sum, na.rm = TRUE)
        sd.area <- sd(CC.data$Area.sum, na.rm = TRUE)
        p.Area <- plot_ly(data = CC.data,
                          x = ~Date, y = ~Area.sum,
                          type = 'scatter', mode = 'lines+markers',
                          color = I('black'), name = 'Measured')
        if (length(dates.na) > 0) {
          p.Area <- p.Area %>% 
            add_markers(x = dates.na, y = rep(0, length(dates.na)),
                        marker = list(symbol = 'x-open', size = 8),
                        name = 'Not detected')
        }
        p.Area <- p.Area %>%
          layout(yaxis = list(range = get.min.max(values = CC.data$Area.sum,
                                                  ref = 0),
                              zeroline = FALSE, title = 'Counts'),
                 xaxis = list(range = c(first.date-3, last.date+3),
                              zeroline = FALSE),
                 shapes = CC.lines(mean.area,
                                   sd.area,
                                   NA))
        output$CC.Area <- renderPlotly({p.Area})
        
        # Matched compounds #
        output$CC.matches <- renderTable(CC.data %>%
                                           select(-c(RT.mean, RT.sd, TF.mean, TF.sd, Area.sum)) %>%
                                           mutate(Date = as.character(Date)))
        output$Control.Cards <- renderUI({
          fluidRow(
            column(6,
                   box(width = 12, status = 'primary',
                       title = tags$b('Retention time'),
                       plotlyOutput(outputId = 'CC.RT') %>% withSpinner(color='#428bca'))
            ),
            column(6,
                   box(width = 12, status = 'primary',
                       title = tags$b('Tailing factor'),
                       plotlyOutput(outputId = 'CC.TF') %>% withSpinner(color='#428bca'))
            ),
            column(6,
                   box(width = 12, status = 'primary',
                       title = tags$b('Area'),
                       plotlyOutput(outputId = 'CC.Area') %>% withSpinner(color='#428bca'))
            ),
            column(6,
                   box(width = 12, status = 'primary',
                       title = 'Number of matches',
                       tableOutput(outputId = 'CC.matches') %>% withSpinner(color='#428bca'))
            )
          )
        })
      }
    }
  })
  
  # Ratios #
  observeEvent(c(input$overview_columns_selected,
                 input$compoundPickerR1u, input$compoundPickerR1d,
                 input$compoundPickerR2u, input$compoundPickerR2d,
                 input$compoundPickerR3u, input$compoundPickerR3d), {
    output$Ratios <- NULL
    req(input$overview_columns_selected)
    rv$cc.selected <- input$overview_columns_selected
    req(input$compoundPickerR1u) # if one exists all do
    if (input$vis.tab %in% c('Control Cards', 'Ratios')) {
      ratios.data <- ratios.calc(rv$unique.data,
                                 colnames(rv$overview.df)[input$overview_columns_selected],
                                 input$compoundPickerR1u, input$compoundPickerR1d,
                                 input$compoundPickerR2u, input$compoundPickerR2d,
                                 input$compoundPickerR3u, input$compoundPickerR3d)
      
      print(ratios.data)
      
      first.date <- sort(ratios.data$Date)[1]
      last.date <- tail(sort(ratios.data$Date),1)

      unable.ratios <- FALSE
      
      # Ratio 1 #
      dates.na.ratio1 <- ratios.data[is.na(ratios.data$Ratio1.Area), 'Date']
      tmp.ratios.data1 <- filter(ratios.data, !is.na(Ratio1.Area))
      if (nrow(tmp.ratios.data1) == 0) {
        unable.ratios <- TRUE
        output$Ratio1.Area <- NULL
        output$Ratio1.TF <- NULL
      } else {
        output$Ratio1.Area <- renderPlotly({
          ratio.CC(data.frame(Date = tmp.ratios.data1$Date,
                              Values = tmp.ratios.data1$Ratio1.Area),
                   dates.na.ratio1, first.date, last.date, 'Area ratio')
        })
        output$Ratio1.TF <- renderPlotly({
          ratio.CC(data.frame(Date = tmp.ratios.data1$Date,
                              Values = tmp.ratios.data1$Ratio1.TF),
                   dates.na.ratio1, first.date, last.date, 'TF ratio')
        })
      }
      
      # Ratio 2 #
      dates.na.ratio2 <- ratios.data[is.na(ratios.data$Ratio2.Area), 'Date']
      tmp.ratios.data2 <- filter(ratios.data, !is.na(Ratio2.Area))
      if (nrow(tmp.ratios.data2) == 0) {
        unable.ratios <- TRUE
        output$Ratio2.Area <- NULL
        output$Ratio2.TF <- NULL
      } else {
        output$Ratio2.Area <- renderPlotly({
          ratio.CC(data.frame(Date = tmp.ratios.data2$Date,
                              Values = tmp.ratios.data2$Ratio2.Area),
                   dates.na.ratio2, first.date, last.date, 'Area ratio')
        })
        output$Ratio2.TF <- renderPlotly({
          ratio.CC(data.frame(Date = tmp.ratios.data2$Date,
                              Values = tmp.ratios.data2$Ratio2.TF),
                   dates.na.ratio2, first.date, last.date, 'TF ratio')
        })
      }
      
      # Ratio 3 #
      dates.na.ratio3 <- ratios.data[is.na(ratios.data$Ratio3.Area), 'Date']
      tmp.ratios.data3 <- filter(ratios.data, !is.na(Ratio3.Area))
      if (nrow(tmp.ratios.data3) == 0) {
        unable.ratios <- TRUE
        output$Ratio3.Area <- NULL
        output$Ratio3.TF <- NULL
      } else {
        output$Ratio3.Area <- renderPlotly({
          ratio.CC(data.frame(Date = tmp.ratios.data3$Date,
                              Values = tmp.ratios.data3$Ratio3.Area),
                   dates.na.ratio3, first.date, last.date, 'Area ratio')
        })
        output$Ratio3.TF <- renderPlotly({
          ratio.CC(data.frame(Date = tmp.ratios.data3$Date,
                              Values = tmp.ratios.data3$Ratio3.TF),
                   dates.na.ratio3, first.date, last.date, 'TF ratio')
        })
      }
      
      if (unable.ratios == TRUE) {
        showNotification('Not enough matches to plot all ratios!',
                         type = 'warning', duration = 3)
      }
      output$Ratios <- renderUI({
        fluidRow(
          column(8, offset = 4,
                 tags$h4(tags$b(paste(input$compoundPickerR1u, input$compoundPickerR1d, sep = ' / ')))
          ),
          column(6,
                 box(width = 12, status = 'primary',
                     title = tags$b('Area'),
                     plotlyOutput(outputId = 'Ratio1.Area') %>% withSpinner(color='#428bca'))
          ),
          column(6,
                 box(width = 12, status = 'primary',
                     title = tags$b('Tailing factor'),
                     plotlyOutput(outputId = 'Ratio1.TF') %>% withSpinner(color='#428bca'))
          ),
          column(8, offset = 4,
                 tags$h4(tags$b(paste(input$compoundPickerR2u, input$compoundPickerR2d, sep = '/')))
          ),
          column(6,
                 box(width = 12, status = 'primary',
                     title = tags$b('Area'),
                     plotlyOutput(outputId = 'Ratio2.Area') %>% withSpinner(color='#428bca'))
          ),
          column(6,
                 box(width = 12, status = 'primary',
                     title = tags$b('Tailing factor'),
                     plotlyOutput(outputId = 'Ratio2.TF') %>% withSpinner(color='#428bca'))
          ),
          column(8, offset = 4,
                 tags$h4(tags$b(paste(input$compoundPickerR3u, input$compoundPickerR3d, sep = '/')))
          ),
          column(6,
                 box(width = 12, status = 'primary',
                     title = tags$b('Area'),
                     plotlyOutput(outputId = 'Ratio3.Area') %>% withSpinner(color='#428bca'))
          ),
          column(6,
                 box(width = 12, status = 'primary',
                     title = tags$b('Tailing factor'),
                     plotlyOutput(outputId = 'Ratio3.TF') %>% withSpinner(color='#428bca'))
          )
        )
      })
    }
  })
  
  observeEvent(c(rv$overview.df, input$vis.tab, input$overview_cells_selected), {
    req(rv$overview.df)
    # Delete #
    if (input$vis.tab == 'Delete') {
      delete.df <- rv$overview.df
      if (!is.null(input$overview_cells_selected)) {
        if(nrow(input$overview_cells_selected) > 0) {
          delete.df[input$overview_cells_selected] <- NA
        }
      }
      
      output$delete.df <- renderTable(delete.df, rownames = TRUE)
      
      output$Delete <- renderUI({
        fluidRow(
          column(2, offset = 1,
                 tags$h5(tags$strong('Delete the marked files?')),
                 actionButton(inputId = 'delete', label = 'Yes')
          ),
          column(9,
                 box(width = 12, status = 'primary',
                     title = 'Preview of the new data',
                     tableOutput(outputId = 'delete.df')
                 )
          )
        )
      })
    }
  })
  
  observeEvent(input$delete, {
    req(input$overview_cells_selected)
    rows.remove <- lapply(seq(nrow(input$overview_cells_selected)), function(row) {
      return(which(rv$data$Date == rownames(rv$overview.df)[input$overview_cells_selected[row, 1]] &
                   rv$data$QCLevel == colnames(rv$overview.df)[input$overview_cells_selected[row, 2]]))
    }) %>% unlist()
    rv$data <- rv$data[-rows.remove, ]
    if (nrow(rv$data) == 0) {
      rv$data <- NULL
      rv$overview.df <- NULL
    } else {
      rv$overview.df <- overview.calc(rv$data)
    }
  })
}
