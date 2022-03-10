pred_multinom <- function(x, beta){
  p <- ncol(x)
  K <- length(beta)/(p+1) + 1
  beta_mat <- matrix(beta, K - 1, p + 1 , byrow = TRUE)
  prod_matrices<-as.matrix(cbind(1, x)) %*% t(beta_mat)
  ## metodo 1 para corregir los exp(prod) que se van a Inf por valores prod>709
  exp_prod_matr <- exp(prod_matrices)
  prod_prev_infinite<-prod_matrices[is.infinite(exp_prod_matr)]
  # prod_matrices[match(prod_prev_infinite,prod_matrices)]<-700 ## se anclan los valores mayores a 709 al valor 700
  prod_matrices[prod_matrices%in%prod_prev_infinite]<-700 ## se anclan los valores mayores a 709 al valor 700
  u_beta <- exp(prod_matrices)
  # u_beta <- exp(as.matrix(cbind(1, x)) %*% t(beta_mat))
  Z <- 1 + apply(u_beta, 1, sum)
  p_beta <- cbind(u_beta, 1)/Z
  as.matrix(p_beta)
}

devianza_calc <- function(x, y){
  dev_fun <- function(beta){
    p_beta <- pred_multinom(x, beta)
    p <- sapply(1:nrow(x), function(i) p_beta[i, y[i]])#+1])
    -2*sum(log(p))
  }
  dev_fun
}

grad_calc <- function(x_ent, y_ent){
  p <- ncol(x_ent)
  K <- length(unique(y_ent)) 
  y_fact <- factor(y_ent) 
  # matriz de indicadoras de clase
  y_dummy <-  model.matrix(~-1 + y_fact)
  salida_grad <- function(beta){
    p_beta <-  pred_multinom(x_ent, beta)
    e_mat <-  (y_dummy  - p_beta)[, -K]
    grad_out <- -2*(t(cbind(1,x_ent)) %*% e_mat)
    as.numeric(grad_out)
  }
  salida_grad
}
descenso <- function(n, z_0, eta, h_deriv, dev_fun,NOdevian=T,alldev=T){
  z <- matrix(0,n, length(z_0))
  z[1, ] <- z_0
  for(i in 1:(n-1)){
    z[i+1, ] <- z[i, ] - eta * h_deriv(z[i, ])
    # if(alldev){
      # devianza_siguiente<-dev_fun(z[i+1, ])
      # if(i %% 50 == 0){
        # print(paste0(i, ' Devianza: ', devianza_siguiente))
        # }
      # }else{
        if(i %% 100 == 0){
          print(paste0("iteration:", i, " - betas from ", i+1))
          if(!NOdevian){
            devianza_siguiente<-dev_fun(z[i+1, ])
            print(paste0('Next deviance: ', devianza_siguiente))
          }
          
          
        }
      # }
    }
  z
}

descenso02 <- function(n, z_0, eta, h_deriv, dev_fun,alldev=T){
  z <- matrix(0,n, length(z_0))
  z[1, ] <- z_0
  for(i in 1:(n-1)){
    z[i+1, ] <- z[i, ] - eta * h_deriv(z[i, ])
    # if(alldev){
    # devianza_siguiente<-dev_fun(z[i+1, ])
    # if(i %% 50 == 0){
    # print(paste0(i, ' Devianza: ', devianza_siguiente))
    # }
    # }else{
    if(i %% 50 == 0){
      devianza_siguiente<-dev_fun(z[i+1, ])
      print(paste0(i, ' Deviance: ', devianza_siguiente))
    }
    # }
  }
  z
}


proceso_prev_iter<-function(pob00_train,
                            vars_not_include=c("id_px","Category","dx_category"),
                            var_obj="dx_category",pob00_test){
  
  x_ent <- pob00_train %>% 
    select(-one_of(vars_not_include)) %>% 
    as.matrix
  p<-ncol(x_ent)
  y_ent <- as.factor(as.numeric(
    pob00_train%>%pull(var_obj)
    ))
  K<-length(unique(y_ent))
  x_ent_s <- scale(x_ent)
  medias <- attr(x_ent_s, 'scaled:center')
  sd <- attr(x_ent_s, 'scaled:scale')
  x_pr <- pob00_test %>% 
    select(-one_of(vars_not_include)) %>% 
    as.matrix
  x_pr_s <- scale(x_pr, center = medias, scale = sd)
  y_pr <- as.factor(as.numeric(
    pob00_test%>%pull(var_obj)
    ))
  
  dev_ent <- devianza_calc(x = x_ent_s,y =  y_ent)
  grad <- grad_calc(x_ent = x_ent_s,y_ent =  y_ent)
  
  return(list(
    "p"=p,"K"=K,
    "dev_ent"=dev_ent,"grad"=grad,
    "x_ent_s"=x_ent_s,"y_ent"=y_ent,
    "x_pr_s"=x_pr_s,"y_pr"=y_pr
    ))

}

update_lista_data_acc<-function(df_devianzas_iteraciones_tmp,
                           iteraciones_nva,
                           proceso_prev_iter_tmp,
                           lista_data_acc){
  x_ent_s<-proceso_prev_iter_tmp$x_ent_s
  y_ent<-proceso_prev_iter_tmp$y_ent
  x_pr_s<-proceso_prev_iter_tmp$x_pr_s
  y_pr<-proceso_prev_iter_tmp$y_pr
  
  top5<-head(
    df_devianzas_iteraciones_tmp%>%
      arrange(devianzas)
    )
  data_acc<-data_frame()
  for(id_top in 1: nrow(top5)){
    id_mindev<-top5[id_top,1]
    print("--------")
    print(id_mindev)
    
    print("deviance:")
    print(proceso_prev_iter_tmp$dev_ent(iteraciones_nva[id_mindev,]))
    
    probas <- pred_multinom(x_ent_s, iteraciones_nva[id_mindev,])
    clase <- apply(probas, 1, which.max)
    print("train:")
    print(table(clase, y_ent ))
    acc_train<-1 - mean(clase != y_ent)
    print(acc_train)
    
    # x_pr_s <- scale(x_pr, center = medias, scale = sd)
    probas <- pred_multinom(x_pr_s, iteraciones_nva[id_mindev,])
    clase <- apply(probas, 1, which.max)
    print("test:")
    print(table(clase, y_pr ))
    acc_test<-1 - mean(clase != y_pr)
    print(acc_test)
    
    data_acc<-data_acc%>%
      bind_rows(
        data_frame(
          id=id_mindev,
          dev_train=proceso_prev_iter_tmp$dev_ent(iteraciones_nva[id_mindev,]),
          acc_train=acc_train,
          acc_test=acc_test 
        )
      )
  }
  
  idmin<-data_acc%>%
    filter(acc_test==max(acc_test))%>%
    filter(acc_train==max(acc_train))%>%
    filter(id==min(id))%>%
    pull(id)
  probas <- pred_multinom(x_ent_s, iteraciones_nva[idmin,])
  clase <- apply(probas, 1, which.max)
  # print("train:")
  table_train<-table(clase, y_ent)
  
  # x_pr_s <- scale(x_pr, center = medias, scale = sd)
  probas <- pred_multinom(x_pr_s, iteraciones_nva[idmin,])
  clase <- apply(probas, 1, which.max)
  table_test<-table(clase, y_pr )
  
  
  iters_prev<-length(lista_data_acc)
  lista_data_acc[[paste0("iterations_",(iters_prev+1))]]<-list(
    "data_acc"=data_acc,
    "idmin"=idmin,
    "table_train"=table_train,
    "table_test"=table_test
  )
  
  print("----------------")
  
  print(
    # as.data.frame(
    as.matrix(
      lista_data_acc[[paste0("iterations_",(iters_prev+1))]]$data_acc
      )
    )
  
  print("----------------")
  lista_data_acc_cop<-lista_data_acc
  lista_data_acc<<-lista_data_acc
  # print(length(lista_data_acc))
  # print(length(lista_data_acc_cop))
  for(i in 1:length(lista_data_acc_cop)){
    print(paste0("i :",i))
    for(d in 1:nrow(lista_data_acc_cop[[i]]$table_train)){
      lista_data_acc_cop[[i]]$table_train[d,d]<-0
      lista_data_acc_cop[[i]]$table_test[d,d]<-0
    }
    print(paste0(i,".... train: ",sum(lista_data_acc_cop[[i]]$table_train), " - test: ",sum(lista_data_acc_cop[[i]]$table_test)))
  }
  
  
  # lista_data_acc<<-lista_data_acc
  
}