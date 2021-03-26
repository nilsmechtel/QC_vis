shiny::runApp(appDir = paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/app'),
       port = getOption('shiny.port'),
       launch.browser = getOption('shiny.launch.browser', interactive()))