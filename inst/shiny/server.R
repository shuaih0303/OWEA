library(shiny)


sigma_construct <- function(type = c('identity','first order correlation'), row, ...) {
  ty <- match.arg(type)
  if(ty == 'identity'){
    return(diag(1,row))
  }
}



# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
   
  x <- reactiveValues(out = list())

  
  # tab for crossover dropout #
  observeEvent(input$run_drop,{
    
    withProgress(message = 'Calculation in Progress...', value = 0, {
      incProgress(0.5)
      
      
      set.seed(input$seed_drop)
      x$out <- 
        OWEA::design(
          model = input$model, n = input$n_drop, opt = input$opt_drop,
          t = input$t_drop, p = input$p_drop, 
          drop = as.numeric(strsplit(input$drop_drop,split = ',')[[1]]),
          max_iter_drop = input$max_iter)
      
        
      
      if(input$effNeed){
        incProgress(0.7, message = 'calculating lower bound...')
        x$out$efflb_drop <- effLB(x$out)$efficiency.self
      } else {
        x$out$efflb_drop <- 'Not Available'
      }
      
      incProgress(0.9, message = 'gethering output...')
      
      colnames(x$out$exact_design) <- c(paste('period',1:input$p_drop), 'repetition')
      
      colnames(x$out$approx_design) <- c(paste('period',1:input$p_drop), 'weight')
      
      x$out$eff <- c('efficiency to approximate design :',
                     eff(x$out)$efficiency)
      
      x$out$get <- c('verifying general equivalence theorem :',
                     x$out$verify_equivalence[,input$p_drop + 1])
      
      x$out$exact_design[,1:input$p_drop] <- 
        sprintf("%s",x$out$exact_design[,1:input$p_drop])
      
      x$out$approx_design[,1:input$p_drop] <- 
        sprintf("%s",x$out$approx_design[,1:input$p_drop])
      
      output$exact_design_drop <- 
        renderTable({x$out$exact_design}, rownames = T, align = 'c', digits = 0)
      
      output$approximate_design_drop <- 
        renderTable({x$out$approx_design}, rownames = T, align = 'c',digits = 4)
      
      output$get_drop <- renderText({x$out$get})
      
      output$efficiency_drop <- renderText({x$out$eff})
      
      output$efflb <- renderText({c('lower bound of efficiency :', x$out$efflb_drop)})
      
      setProgress(1, message = 'done')
      
      }) 
    
    output$download_drop <- downloadHandler(
      filename = function() {
        paste(input$model,Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.table(x$out$exact_design, file, sep = ',')
        write.table(x$out$approx_design, file, append = T, sep = ',')
      }
    )
    
    }
  )
  

  
  
  # tab for proportional model #
  observeEvent(input$run_prop,{
    
    updateNavbarPage(session, inputId = 'model', selected = input$model)
    
    withProgress(message = 'Calculation in Progress...', value = 0, {
      incProgress(0.5)
      
      sigma <- sigma_construct('identity',input$p_prop)
      
      set.seed(input$seed_prop)
      x$out <- design(model = input$model, n = input$n_prop, opt = input$opt_prop, 
                      t = input$t_prop, p = input$p_prop,
                      sigma = sigma, 
                      tau = as.numeric(strsplit(input$tau_prop,split = ',')[[1]]), 
                      lambda = input$lambda_prop, max_iter = input$max_iter_prop)
      
      incProgress(0.9, message = 'gethering output...')  
      
      colnames(x$out$exact_design) <- c(paste('period',1:input$p_prop), 'repetition')
      colnames(x$out$approx_design) <- c(paste('period',1:input$p_prop), 'weight')
      
      x$out$eff <- c('efficiency to approximate design :',
                     eff(x$out)$efficiency)
      
      x$out$get <- c('verifying general equivalence theorem :',
                     x$out$verify_equivalence[,input$p_prop + 1])
      
      x$out$exact_design[,1:input$p_prop] <- sprintf("%s",x$out$exact_design[,1:input$p_prop])
      x$out$approx_design[,1:input$p_prop] <- sprintf("%s",x$out$approx_design[,1:input$p_prop])
      
      output$exact_design_prop <- 
        renderTable({x$out$exact_design}, rownames = T, align = 'c', digits = 0)
      
      output$approximate_design_prop <- 
        renderTable({x$out$approx_design}, rownames = T, align = 'c',digits = 4)
      
      output$get_prop <- renderText({x$out$get})
      
      output$efficiency_prop <- renderText({x$out$eff})
      
      output$params_prop <- renderPrint({str(x$out)})
      
      setProgress(1, message = 'done')
      
      }) 
    
    output$download_prop <- downloadHandler(
      filename = function() {
        paste(input$model,Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.table(x$out$exact_design, file, sep = ',')
        write.table(x$out$approx_design, file, append = T, sep = ',')
      }
    )
    
    
    }
  )
  
  # tab for interference model #
  observeEvent(input$run_inter,{
    
    updateNavbarPage(session, inputId = 'model', selected = input$model)
    
    withProgress(message = 'Calculation in Progress...', value = 0, {
      incProgress(0.5)
      
      set.seed(input$seed_inter)
      x$out <- design(model = input$model, n = input$n_inter, opt = input$opt_inter, 
                      t = input$t_inter, p = input$p_inter, 
                      sigma = sigma_construct('identity', input$p_inter),
                      max_iter = input$max_iter_inter)
      
      incProgress(0.9, message = 'gethering output...')  
      
      colnames(x$out$exact_design) <- c(paste('period',1:input$p_inter), 'repetition')
      colnames(x$out$approx_design) <- c(paste('period',1:input$p_inter), 'weight')
      
      x$out$eff <- c('efficiency to approximate design :',
                     eff(x$out)$efficiency)
      
      x$out$get <- c('verifying general equivalence theorem :',
                     x$out$verify_equivalence[,input$p_inter + 1])
      
      x$out$exact_design[,1:input$p_inter] <- sprintf("%s",x$out$exact_design[,1:input$p_inter])
      x$out$approx_design[,1:input$p_inter] <- sprintf("%s",x$out$approx_design[,1:input$p_inter])
      
      output$exact_design_inter <-
        renderTable({x$out$exact_design}, rownames = T, align = 'c', digits = 0)
      
      output$approximate_design_inter <- 
        renderTable({x$out$approx_design}, rownames = T, align = 'c',digits = 4)
      
      output$get_inter <- renderText({x$out$get})
      
      output$efficiency_inter <- renderText({x$out$eff})
      
      setProgress(1, message = 'done')

      })
    
    output$download_inter <- downloadHandler(
      filename = function() {
        paste(input$model,Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.table(x$out$exact_design, file, sep = ',')
        write.table(x$out$approx_design, file, append = T, sep = ',')
      }
    )
    
    
    }
  )
  
  

})
