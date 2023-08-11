#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(rlang)
library(plotly)
library(plyr)
library(minpack.lm)
library(foreach)
library(dplyr)
library(kableExtra)
library(rstan)
library(readxl)
library(xlsx)
library(reactable)
library(MASS)
library(shinybusy)
###functions
reshape.func <- function(df, vec) {
  df1 <- as.data.frame(df)
  df1$index <- rownames(df)
  tmpl <- df1[,-match(vec,colnames(df1))]
  finaldf <- data.frame()
  for(i in 1:length(vec)) {
    tmpd1 <-  df1[,match(vec[i],colnames(df1))]
    tmpd <- cbind(tmpl, tmpd1)
    tmpd$element <- vec[i]
    colnames(tmpd) <- c(colnames(tmpl), "value", "element")
    finaldf <- rbind(finaldf, tmpd)
  }
  
  return(finaldf)
}
stylefunc <- function(data,acc) {
  value <- data[,match(paste0(acc,".sig"),colnames(data))]
  if (value > 0.05) {
    color <- "#e00000"
  } else if (value <= 0.01) {
    color <- "#008000"
  } else {
    color <- "#e8b600"
  }
  list(color = color, fontWeight = "bold")
}
calculator <- function(fomula, variables, params) {
  
  with(c(as.list(params), variables), {
    eval(fomula)
  }, list())
}
generate.log.samples <- function(vec, intvl=0.1) {
  foreach::foreach(a=seq(min(log10(vec[vec>0])),
                         max(log10(vec[vec>0])),
                         intvl),.combine = 'c') %do% 10^a
}
###gate models
formula.not.model <- function(k3,alpha3,K3,n3,R3) {
  k3*(alpha3+(K3^n3)/(K3^n3+R3^n3))
  #(K3^n3)/(K3^n3+R3^n3)
}
formula.not <- function() {
  expr(k3*(alpha3+(K3^n3)/(K3^n3+R3^n3)))
  #(K3^n3)/(K3^n3+R3^n3)
}
formula.activ <- function() {
  expr(k*(alpha+(I^n1)/(K1^n1+I^n1)))
}
formula.activ.model <- function(k,alpha,n1,K1,I) {
  k*(alpha+(I^n1)/(K1^n1+I^n1))
}
formula.and <- function() {
  expr(Gmax*(R/Kr)^nr*(S/Ks)^ns/(1+(R/Kr)^nr)/(1+(S/Ks)^ns))
}
formula.and.model <- function(Gmax,Kr,nr,Ks,ns,R,S) {
  Gmax*(R/Kr)^nr*(S/Ks)^ns/(1+(R/Kr)^nr)/(1+(S/Ks)^ns)
}

formula.not.model.cello <- function(ymin,ymax,K,n,x) {
  ymin + (ymax -ymin)/(1+(x/K)^n)
}
formula.not.cello <- function() {
  expr(ymin + (ymax -ymin)/(1+(x/K)^n))
}

### nsl model fitting functions
nsl.fit.notgate <- function(mytab, k3, alpha3, K3, n3) {
  modeldat <- reshape.func(mytab, colnames(mytab)[2:length(colnames(mytab))])
  #print(paste(k3, alpha3, K3, n3,sep=","))
  nlsmod.1 <- nlsLM(value ~ formula.not.model(k3,alpha3,K3,n3,R3),
                    data = modeldat,
                    algorithm = "port",
                    start=list(k3=k3,
                               alpha3=alpha3,
                               K3=K3,
                               n3=n3),
                    trace=T
  )
  return(nlsmod.1)
}

nsl.fit.notgate.cello <- function(mytab, ymin,ymax,K,n) {
  modeldat <- reshape.func(mytab, colnames(mytab)[2:length(colnames(mytab))])
  #print(paste(k3, alpha3, K3, n3,sep=","))
  nlsmod.1 <- nlsLM(value ~ formula.not.model.cello(ymin,ymax,K,n,x),
                    data = modeldat,
                    algorithm = "port",
                    start=list(ymin=ymin,ymax=ymax,K=K,n=n),
                    trace=T
  )
  return(nlsmod.1)
}

nsl.fit.sensor <- function(mytab, k,alpha,n1,K1) {
  modeldat <- reshape.func(mytab, colnames(mytab)[2:length(colnames(mytab))])
  #print(paste(k, alpha, K1, n1,sep=","))
  #print(modeldat)
  nlsmod.1 <- nlsLM(value ~ formula.activ.model(k,alpha,n1,K1,I),
                    data = modeldat,
                    algorithm = "port",
                    start=list(k=k,
                               alpha=alpha,
                               n1=n1,
                               K1=K1
                    ),trace=T)
  return(nlsmod.1)
  
}

nsl.fit.and <- function(mytab, Gmax,Kr,nr,Ks,ns) {
  modeldat <- mytab
  #print(paste(k, alpha, K1, n1,sep=","))
  #print(modeldat)
  nlsmod.1 <- nlsLM(O ~ formula.and.model(Gmax,Kr,nr,Ks,ns,R,S),
                    data = modeldat,
                    algorithm = "port",
                    start=list(Gmax=Gmax,
                               Kr=Kr,
                               nr=nr,
                               Ks=Ks,
                               ns=ns
                    ),trace=T)
  return(nlsmod.1)
  
}

call.priors <- function(nls.mod) {
  df <- as.data.frame(summary(nls.mod)$coefficients[,c(1,2)])
  mylist <- foreach::foreach(a=1:length(df[,1])) %do% paste(rownames(df)[a],"~ normal(",df$Estimate[a],",",df$`Std. Error`[a],");\n")
  text <- do.call("paste",c(mylist))
  return(text)
}

call.bmcmc <- function(xvalue,yvalue, stancodep1, stancodep2, priors,chains,iter,cores) {
  testdat <- list(N=length(xvalue),
                  L=length(unique(xvalue)),
                  P=length(xvalue)/length(unique(xvalue)),
                  x=matrix(xvalue,nrow = length(unique(xvalue))),
                  Y= matrix(yvalue,nrow = length(unique(xvalue))))
  
  stanmodel <- paste0(stancodep1,priors,stancodep2)
  cat(stanmodel)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
   fitme <- stan(model_code = stanmodel, 
                 chains=chains,
                 iter=iter,
                 cores=cores,
               model_name = "GateRemodelling", 
               data = testdat)
  
  
  return(fitme)
  
}

call.bmcmc.2i <- function(xvalue,zvalue,yvalue, stancodep1, stancodep2, priors,chains,iter,cores) {

  testdat <- list(N=length(xvalue),
                  z=zvalue,
                  x=xvalue,
                  Y= yvalue)
  
  stanmodel <- paste0(stancodep1,priors,stancodep2)
  cat(stanmodel)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  fitme <- stan(model_code = stanmodel, 
                model_name = "GateRemodelling", 
                chains=chains,
                iter=iter,
                cores=cores,
                data = testdat)
  
  return(fitme)
  
}
stancodes <- function() {
list(c(stanmodel.activ.p1 = "
data {
  int N; 
  int L;
  int P;
  real x[L,P]; 
  real Y[L,P]; 
} 
parameters {
  real alpha;
  real<lower=0> k;
  real<lower=0> K1;
  real<lower=0> n1;
  real<lower=0> tau; 
} 
transformed parameters {
  real sigma; 
  real m[L,P];
  for (i in 1:L) {
     for (j in 1:P) {
       m[i,j] = k * (alpha + pow(x[i,j], n1)/(pow(K1, n1) + pow(x[i,j], n1) ) );
     }
  }
 sigma = 1 / sqrt(tau); 
} 
model {
", stanmodel.activ.p3 = "
tau ~ gamma(.1, .1); 
  // likelihood
  for (i in 1:L) {
     for (j in 1:P) {
       Y[i,j] ~ normal(m[i,j], sigma);
     }
  }
}
generated quantities{
  real Y_mean[L,P]; 
  real Y_pred[L,P]; 
  for (i in 1:L) {
     for (j in 1:P) {
       Y_mean[i,j] =  k * (alpha + pow(x[i,j], n1)/(pow(K1, n1) + pow(x[i,j], n1) ) );
       Y_pred[i,j] = normal_rng(Y_mean[i,j], sigma);  
     }
  }

}
"), c(stanmodel.not.p1 = "
data {
  int N; 
  int L;
  int P;
  real x[L,P]; 
  real Y[L,P]; 
} 
parameters {
  real alpha3;
  real<lower=0> k3;
  real<lower=0> K3;
  real<lower=0> n3;
  real<lower=0> tau; 
} 
transformed parameters {
  real sigma; 
  real m[L,P];
  for (i in 1:L) {
     for (j in 1:P) {
       m[i,j] = k3*(alpha3+1/(1+pow(x[i,j]/K3,n3)));
     }
  }
 sigma = 1 / sqrt(tau); 
} 
model {
", stanmodel.not.p3 = "
tau ~ gamma(.1, .1); 
  // likelihood
  for (i in 1:L) {
     for (j in 1:P) {
       Y[i,j] ~ normal(m[i,j], sigma);
     }
  }
}
generated quantities{
  real Y_mean[L,P]; 
  real Y_pred[L,P]; 
  for (i in 1:L) {
     for (j in 1:P) {
       Y_mean[i,j] = k3*(alpha3+1/(1+pow(x[i,j]/K3,n3)));
       Y_pred[i,j] = normal_rng(Y_mean[i,j], sigma);  
     }
  }

}
"), c(stanmodel.and.p1 = "
data {
  int N;
  real x[N]; 
  real z[N]; 
  real Y[N]; 
} 
parameters {
  real<lower=0> Gmax;
  real<lower=0> Kr;
  real<lower=0> nr;
  real<lower=0> Ks;
  real<lower=0> ns; 
  real<lower=0> tau; 
} 
transformed parameters {
  real sigma; 
  real m[N];
  for (i in 1:N) {
       m[i] =  Gmax*pow(x[i]/Kr,nr)*pow(z[i]/Ks,ns)/(1+pow(x[i]/Kr,nr))/(1+pow(z[i]/Ks,ns));
     
  }
 sigma = 1 / sqrt(tau); 
} 
model {
",
      stanmodel.and.p2 = "
tau ~ gamma(.1, .1); 
  // likelihood
  for (i in 1:N) {
       Y[i] ~ normal(m[i], sigma);
  }
}
generated quantities{
  real Y_mean[N]; 
  real Y_pred[N]; 
  for (i in 1:N) {
       Y_mean[i] =  Gmax*pow(x[i]/Kr,nr)*pow(z[i]/Ks,ns)/(1+pow(x[i]/Kr,nr))/(1+pow(z[i]/Ks,ns));
       Y_pred[i] = normal_rng(Y_mean[i], sigma);  
     
  }

}
" ), 
     c(stanmodel.cello.not.p1 = "
data {
  int N; 
  int L;
  int P;
  real x[L,P]; 
  real Y[L,P]; 
} 
parameters {
  real ymin;
  real<lower=0> ymax;
  real<lower=0> K;
  real<lower=0> n;
  real<lower=0> tau; 
} 
transformed parameters {
  real sigma; 
  real m[L,P];
  for (i in 1:L) {
     for (j in 1:P) {
       m[i,j] = ymin + (ymax - ymin)/(1 + pow(x[i,j]/K,n));
     }
  }
 sigma = 1 / sqrt(tau); 
} 
model {
",stanmodel.cello.not.p3 = "
tau ~ gamma(.1, .1); 
  // likelihood
  for (i in 1:L) {
     for (j in 1:P) {
       Y[i,j] ~ normal(m[i,j], sigma);
     }
  }
}
generated quantities{
  real Y_mean[L,P]; 
  real Y_pred[L,P]; 
  for (i in 1:L) {
     for (j in 1:P) {
       Y_mean[i,j] = ymin + (ymax - ymin)/(1 + pow(x[i,j]/K,n));
       Y_pred[i,j] = normal_rng(Y_mean[i,j], sigma);  
     }
  }

}
"))
}

###plot functions
compare.radar <- function(numvec, mytab, compare1, compare2) {
  plotme <- mytab[numvec,]
  
  fig <- plot_ly(
    type = 'scatterpolar',
    fill = 'toself'
  )
  
  #coln1 <- c("k","n1","R2")
  #coln2 <- c("K1","alpha")
  #coln1 <- c("k3","n3","R2")
  #coln2 <- c("K3","alpha3")
  
  coln1 <- compare1
  coln2 <- compare2
  means1 <- foreach::foreach(a=coln1, .combine = "c") %do% mean(as.numeric(mytab[,colnames(mytab)==a]))
  sd1 <- foreach::foreach(a=coln1, .combine = "c") %do% sd(as.numeric(mytab[,colnames(mytab)==a]))
  means2 <- foreach::foreach(a=coln2, .combine = "c") %do% mean(as.numeric(mytab[,colnames(mytab)==a]))
  sd2 <- foreach::foreach(a=coln2, .combine = "c") %do% sd(as.numeric(mytab[,colnames(mytab)==a]))
  
  for(i in 1:length(plotme$element)) {
    num1 <- as.numeric(plotme[i, match(coln1, colnames(plotme))])
    degree1 <- (num1 - means1)/sd1
    num2 <- as.numeric(plotme[i, match(coln2, colnames(plotme))])
    degree2 <- (means2  -num2)/sd2
    degree <- c(degree1, degree2)
    print(degree)
    print(plotme$element[i])
    
    fig <- fig %>%
      add_trace(
        r = c(degree, degree[1]),
        theta = c(coln1,coln2, coln1[1]),
        name = plotme$element[i]
      )
    
  }
  
  fig <- fig %>%
    layout(
      polar = list(
        radialaxis = list(
          visible = T
        )
      )
    )
  
  return(fig)
}

compare.io <- function(numvec,mytab, compare, Invalue, Outvalue, calfunc) {
  plotme <- mytab[numvec,]
  
  xrange <- plotme[,match(Outvalue,colnames(plotme))]
  
  minx <-  min(log10(xrange[xrange>0] / 100))
  
  maxx <- max(log10(xrange[xrange>0] * 100))
  
  x <- foreach::foreach(a=seq(minx,maxx,0.05),.combine = 'c') %do% 10^a
  
  mylist <- list(x)
  names(mylist) <- Invalue
  
  fig <- plot_ly()
  
  for(i in 1:length(plotme$element)) {
    y <- calculator(calfunc, variables = plotme[i,match(compare, colnames(plotme))], params = mylist )
    mytext <- do.call("paste0", foreach::foreach(a=compare, b=as.numeric(plotme[i,match(compare, colnames(plotme))])) %do% paste0(a,"=",b, "\n "))
    fig <- fig %>%
      add_lines(x=x, y = y,
                name = plotme$element[i],
                text = mytext,
                line = list(shape = "spline")) %>%
      layout(xaxis = list(range = c(log10(min(x)), log10(max(x))),
                          type = 'log',
                          zerolinecolor = '#ffff',
                          zerolinewidth = 2,
                          gridcolor = 'ffff',
                          title = 'x'))
  }
  return(fig)
}

plot.pred <- function(obstab, Invalue, calfunc, model.fit, mytext) {
  minx <-  min(log10(obstab$Input[obstab$Input>0]))
  maxx <- max(log10(obstab$Input[obstab$Input>0]))
  
  
  x <- foreach::foreach(a=seq(minx,maxx,0.05),.combine = 'c') %do% 10^a
  
  mylist <- list(x)
  names(mylist) <- Invalue
  
  y <- calculator(calfunc, variables = coef(model.fit), params = mylist )
  
  fig <- plot_ly(x = ~x)
  fig <- fig %>%
    add_lines(y = y,
              name = "prediction",
              line = list(shape = "spline")) %>%
    layout(xaxis = list(range = c(log10(min(x)), log10(max(x))),
                        type = 'log',
                        zerolinecolor = '#ffff',
                        zerolinewidth = 2,
                        gridcolor = 'ffff',
                        title = 'x'))
  fig <- fig %>%
    add_trace(data = obstab, x = ~ Input, y = ~ Output,  type = 'scatter', name = "observation",
              error_y = ~list(array = sd,
                              color = '#8f8f8f'))
  fig <- fig %>%
    add_trace(x = min(x)*10, y = max(y)*1/2, text = mytext, mode = "text", size = I(1000)) %>%
    config(mathjax = "cdn")
  
  return(fig)
}

pred.mcmc <- function(mcmcfit,nlsfit,formula,modeldat) {
  
  pars <- names(coef(nlsfit))
  df <- as.data.frame(summary(mcmcfit)$summary)
  newpars <- t(df[match(pars,rownames(df)),])[1,]
  
  
  Y_mean <- rstan::extract(mcmcfit, "Y_mean")
  Y_mean_cred <- apply(Y_mean$Y_mean, 2, quantile, c(0.05, 0.95))
  Y_mean_mean <- apply(Y_mean$Y_mean, 2, mean)
  
  Y_pred <- rstan::extract(mcmcfit, "Y_pred")
  Y_pred_cred <- apply(Y_pred$Y_pred, 2, quantile, c(0.05, 0.95))
  Y_pred_mean <- apply(Y_pred$Y_pred, 2, mean)
  xlist <- modeldat$x
  
  fig <- plot_ly() %>% 
    add_trace(x =xlist, y =modeldat$value,  type = 'scatter', name = "observation")%>%
    add_lines(x=unique(xlist), y = Y_pred_mean, 
              name = "Mean_mean", 
              line = list(shape = "spline")) %>%
    add_lines(x=unique(xlist), y = calculator(formula, variables = as.list(newpars), params = list(x=unique(xlist)) ), 
              name = "Bayesian_MCMC", 
              line = list(shape = "spline")) %>%
    add_lines(x=unique(xlist), y = Y_mean_cred[1,], 
              name = "Mean_cred_5%", line = list(color = 'rgb(22, 96, 167)', width = 2, dash = 'dot',shape = "spline")
    ) %>%
    add_lines(x=unique(xlist), y = Y_mean_cred[2,], fill ='tonexty',fillcolor='rgba(22, 96, 167,0.1)',
              name = "Mean_cred_95%", line = list(color = 'rgb(22, 96, 167)', width = 2, dash = 'dot',shape = "spline")
    ) %>%
    add_lines(x=unique(xlist), y = Y_pred_cred[1,], 
              name = "Pred_cred_5%", line = list(color = 'rgb(205, 12, 24)', width = 2, dash = 'dot',shape = "spline")
    ) %>%
    add_lines(x=unique(xlist), y = Y_pred_cred[2,], fill ='tonexty', fillcolor='rgba(205, 12, 24,0.1)',
              name = "Pred_cred_95%", line = list(color = 'rgb(205, 12, 24)', width = 2, dash = 'dot',shape = "spline")
    ) %>%
    add_lines(x=unique(xlist), y = calculator(formula, variables = as.list(coef(nlsfit)), params = list(x=unique(xlist)) ), 
              name = "NLS_prediction", 
              line = list(shape = "spline")) %>%
    layout(xaxis = list(range = c(log10(min(xlist[xlist>0])), log10(max(xlist))),
                        type = 'log',
                        zerolinecolor = '#ffff',
                        zerolinewidth = 2,
                        gridcolor = 'ffff',
                        title = 'x')) 
  
  return(fig)
}

reformat <- function(testdata, incol) {
  reshaped <- reshape.func(testdata, colnames(testdata)[colnames(testdata)!=incol])
  colnames(reshaped)[colnames(reshaped)==incol] <- "Input"
  #reshaped$value <- log10(reshaped$value)
  data_mean <- ddply(reshaped, "Input", summarise, Output = mean(value))
  data_sd <- ddply(reshaped, "Input", summarise, Output = sd(value))
  data.frame(data_mean, sd= data_sd$Output)
}

plot.obs <- function(df, selected=NULL) {
  fig <- plot_ly(data = df, x = ~ Input, y = ~ Output,  type = 'scatter', mode = 'lines',
                 error_y = ~list(array = sd,
                                 color = '#000000'))
  
  if(length(selected)>0) {
    newdf <- df[selected,]
    fig <- fig %>% add_trace(x = newdf$Input, y = newdf$Output, name="selected", mode="markers")
    
  }
  
  fig <- fig %>%
    layout(xaxis = list(range = c(log10(min(df$Input[df$Input>0])), log10(max(df$Input))),
                        type = 'log',
                        zerolinecolor = '#ffff',
                        zerolinewidth = 2,
                        gridcolor = 'ffff',
                        title = 'x'))
  
  return(fig)
}

plot.surface <- function(xlist,ylist,zlist,colorscl='Viridis') {
  Input1 <- xlist
  Input2 <- ylist
  Output <- matrix(zlist, nrow =length(xlist))
  
  fig <- plot_ly( x=Input1, y=Input2) %>%
    add_surface(
      z = Output,
      colorscale = colorscl,
      contours = list(
        z = list(
          show=TRUE,
          usecolormap=TRUE,
          highlightcolor="#ff0000",
          project=list(z=TRUE)
        )
      )
    )
  
  fig <- fig %>%
    layout(
      scene = list(
        xaxis = list(type = 'log'),
        yaxis = list(type = 'log'),
        camera = list(eye = list(x = 0, y = -1, z = 0.5)),
        aspectratio = list(x = .9, y = .8, z = 0.6)))
  
  
  return(fig)
}

plot.mcmc.surface <- function(df1,df2,df3){
  fig <- plot_ly( ) %>% 
    add_surface(x=unique(df1$R), y=unique(df1$S),
                z = matrix(df1$O, ncol=length(unique(df1$S))),
                labels="NLM_pred",
                colorbar=list(title='NLM_pred'),
                #  colorscale = "Jet",
                contours = list(
                  z = list(
                    show=TRUE,
                    usecolormap=TRUE,
                    highlightcolor="#ff0000",
                    project=list(z=TRUE)
                  )
                )
    )
  fig <- fig %>% 
    add_surface(x=unique(df2$R), y=unique(df2$S),
                z = matrix(df2$O, ncol=length(unique(df2$S))),
                colorscale = "Jet",
                colorbar=list(title='BMCMC_pred'),
                labels="BMCMC_pred",
                contours = list(
                  z = list(
                    show=TRUE,
                    usecolormap=TRUE,
                    highlightcolor="#ff0000",
                    project=list(z=TRUE)
                  )
                )
    )
  
  
  fig <- fig %>% 
    add_markers(x=df3$R, y=df3$S, z = df3$O,  name = "observation",size=5)
  
  fig <- fig %>% 
    layout(
      scene = list(
        xaxis = list(type = 'log'),
        yaxis = list(type = 'log'),
        camera = list(eye = list(x = 0, y = -1, z = 0.5)),
        aspectratio = list(x = .9, y = .8, z = 0.6)))
  
  
  return(fig)
}

plot.and.curve <- function(anddf, inval) {
  data_mean <- ddply(anddf, inval, summarise, Output = mean(O))
  data_sd <- ddply(anddf, inval, summarise, Output = sd(O))
  data <- data.frame(data_mean, sd=data_sd$Output)
  
  xlist <- data[,match(inval, colnames(data))]
  print(xlist)
  #data$Input <- log10(as.numeric(as.character(data$Input)))
  
  #data <- rename(data, c("data_sd.Output" = "sd"))
  
  fig <- plot_ly(x = xlist)
  fig <- fig %>%
    add_lines(y = data$Output,
              name = "prediction",
              line = list(shape = "spline")) %>%
    layout(xaxis = list(range = c(log10(min(xlist)), log10(max(xlist))),
                        type = 'log',
                        zerolinecolor = '#ffff',
                        zerolinewidth = 2,
                        gridcolor = 'ffff',
                        title = 'x'))
  fig <- fig %>%
    add_trace(x = xlist, y = data$Output,  type = 'scatter', name = "observation",
              error_y = ~list(array = data$sd,
                              color = '#8f8f8f'))
  
  return(fig)
}

plot.violin <- function(testdata,incol) {
  reshapedtab <- reshape.func(testdata, colnames(testdata)[colnames(testdata)!=incol])
  colnames(reshapedtab)[colnames(reshapedtab)==incol] <- "Input"
  fig <- plot_ly(type = 'violin')
  
  for (i in 1:length(unique(reshapedtab$Input))) {
    fig <- add_trace(
      fig,
      y = reshapedtab$value[reshapedtab$Input==unique(reshapedtab$Input)[i]],
      hoveron = "points+kde",
      name = paste0("Input=",round(unique(reshapedtab$Input)[i],digits = 5)),
      side = 'negative',
      box = list(
        visible = T
      ),
      points = 'all',
      jitter = 0,
      scalemode = 'count',
      meanline = list(
        visible = T
      ),
      color = I("#8dd3c7"),
      marker = list(
        line = list(
          width = 2,
          color = "#8dd3c7"
        ),
        symbol = 'line-ns'
      )
    )
  }
  
  fig <- layout(
    fig,
    title = "Sample distribution",
    yaxis = list(
      zeroline = F
    ),
    violingap = 0,
    violingroupgap = 0,
    violinmode = 'overlay',
    legend = list(
      tracegroupgap = 0
    )
  )
  
  return(fig)
}
###resources
notgatedb <- read.table("www/notgate.tsv", header = T, sep = "\t")
celllo.notgatedb <- read.table("www/cello.notgate.tsv", header = T, sep = "\t")
sensordb <- read.table("www/sensor.tsv", header = T, sep = "\t")
andgatedb <- read.table("www/andgate.tsv", header = T, sep = "\t")
###texts
mygates <- c("Sensor","NotGate","AndGate","CelloNotGate")

texter.activgate1 <- TeX("f(I) = k*(\\alpha + \\frac{{I}^{n_\\text{1}}}{{K}^{n_\\text{1}} + {I}^{n_\\text{1}}})")
texter.notgate1 <- TeX("f(R_\\text{3}) = k_\\text{3}*(\\alpha_\\text{3} + \\frac{{K_\\text{3}}^{n_\\text{3}}}{{K_\\text{3}}^{n_\\text{3}} + {R_\\text{3}}^{n_\\text{3}}})")
texter.andgate1 <- TeX("f(R,S) = G_\\text{max}*\\frac{{R/{K_\\text{r}}}^{n_\\text{r}}}{1+{R/{K_\\text{r}}}^{n_\\text{r}}}*\\frac{{S/{K_\\text{s}}}^{n_\\text{s}}}{1+{S/{K_\\text{s}}}^{n_\\text{s}}}")
texter.notgatecello <- TeX("f(x) = y_\\text{min}+ \\frac{y_\\text{max} - y_\\text{min}}{1 + (\\frac{x}{K})^n}")

functionlist <- c(texter.activgate1[[1]],
                  texter.notgate1[[1]],
                  texter.andgate1[[1]],
                  texter.notgatecello[[1]])

###test data sets

load.not.example <- function() {
  pars1 <- c(k3=7.538e4, alpha3=0.0527, K3 = 111e-3, n3 = 7.647)
  pars2 <- c(k3=7.538e4 + 0.234e4, alpha3=0.0527 -0.0144, K3 = 111e-3 + 4.7e-3, n3 = 7.647 +1.379)
  pars3 <- c(k3=7.538e4 + 0.234e4, alpha3=0.0527 +0.0144, K3 = 111e-3 - 4.7e-3, n3 = 7.647 -1.379)
  
  vars <- list(R3= c(0, 3.9e-4, 1.6e-3, 6.3e-3, 2.5e-2,0.05,0.075, 0.1, 0.13,1.6, 6.4,12.8) )
  
  calculator(formula.not(), variables = vars, params = pars3 )
  
  testdata <- data.frame(vars,
                         set1=calculator(formula.not(), variables = vars, params = pars1 ),
                         set2=calculator(formula.not(), variables = vars, params = pars2 ),
                         set3=calculator(formula.not(), variables = vars, params = pars3 ))
  
  return(testdata)
}

load.sensor.example <- function() {
  pars1 <- c(k=9456 , alpha=0.0012 , n1=1.37 , K1=0.228 * 1e-3)
  pars2 <- c(k=9456+487 , alpha=0.0012 +0.0276, n1=1.37-0.0276 , K1=(0.228 + 0.039)* 1e-3)
  pars3 <- c(k=9456-487 , alpha=0.0012 -0.0276, n1=1.37-0.0276 , K1=(0.228 - 0.039)* 1e-3)
  
  vars <- list(I= c(0, 3.9e-4, 1.6e-3, 6.3e-3, 2.5e-2,0.05,0.075, 0.1, 0.13,1.6, 6.4,12.8) * 1e-3 )
  
  testdata <- data.frame(vars,
                         set1=calculator(formula.activ(), variables = vars, params = pars1 ),
                         set2=calculator(formula.activ(), variables = vars, params = pars2 ),
                         set3=calculator(formula.activ(), variables = vars, params = pars3 ))
  return(testdata)
  
}

###shiny server
shinyServer(function(input, output, session) {
  
  ####dynamic UI functions
  output$moreControls <- renderUI({
    req(input$gatemodel)
   # req(input$file1)
    
    if(is.null(input$file1)) {
      mydiv <-   tagList({
        box(title = "Parameterisation", width = 6, solidHeader = T, status = "primary",height = "300px"
        )
      })
    }
    
    else {
    
    flofile <- rv$obsdf()
    
    #modeldat <- reshape.func(flofile, colnames(flofile)[2:length(colnames(flofile))])
    # print(modeldat)
    
    if(as.numeric(input$gatemodel) == 2) {
      modeldat <- reshape.func(flofile, colnames(flofile)[2:length(colnames(flofile))])
        mydiv <-   tagList({
          box(title = "Parameterisation", width = 6, solidHeader = T, status = "primary",height = "300px",
              fluidRow(
                column(6, numericInput('k3', "k3 (a.u.)", 5, min = 0, max = 1e6, value=max(modeldat$value))),
                column(6, numericInput('n3', "n3", 5, min = 0, max = 1e2, value=log10(max(modeldat$value) -min(modeldat$value))))
              ),
              fluidRow(
                column(6, numericInput('K3', "K3 (M)", 5, min = 0, max = 1e5, value=modeldat$R3[abs(modeldat$value - (max(modeldat$value)+ min(modeldat$value))/2) == min(abs(modeldat$value - (max(modeldat$value)+ min(modeldat$value))/2))])),
                column(6, numericInput('alpha3', "α3", 5, min = 0, max = 1, value=0))
              ),
              fluidRow(
                column(4, p(HTML("<b>Step 4: Simulate NLS fitting</b>"),
                            span(shiny::icon("info-circle"), id = "info_uu"),
                            br(),
                       actionBttn("sim", "Fit NLS", color ="primary", style = "jelly"),
                       tippy::tippy_this(elementId = "info_uu",tooltip = "",placement = "right"))),
                column(4, p(HTML("<b>Step 5 (optional): Simulate BayesianMCMC fitting</b>"),
                            span(shiny::icon("info-circle"), id = "info_uu"),
                            br(),
                       actionBttn("simb", "Fit BayesianMCMC", color ="primary", style = "jelly"),
                       tippy::tippy_this(elementId = "info_uu",tooltip = "",placement = "right"))),
                column(4, p(HTML("<b>Step 6: Add to the record board</b>"),
                            span(shiny::icon("info-circle"), id = "info_uu"),
                            br(),
                       actionBttn("add", "Add",color ="warning", style = "jelly"),
                       tippy::tippy_this(elementId = "info_uu",tooltip = "",placement = "right")))
              )
          )
        })
    }
    if(as.numeric(input$gatemodel) == 4) {
      modeldat <- reshape.func(flofile, colnames(flofile)[2:length(colnames(flofile))])
        mydiv <-   tagList({
          box(title = "Parameterisation", width = 6, solidHeader = T, status = "primary",height = "300px",
              fluidRow(
                column(6, numericInput('ymax', "ymax (a.u.)", 5, min = 0, max = 1e6, value=max(modeldat$value))),
                column(6, numericInput('ymin', "ymin", 5, min = 0, max = 1, value=0))
                
              ),
              fluidRow(
                column(6, numericInput('K', "K (M)", 5, min = 0, max = 1e5, value=modeldat$x[abs(modeldat$value - (max(modeldat$value)+ min(modeldat$value))/2) == min(abs(modeldat$value - (max(modeldat$value)+ min(modeldat$value))/2))])),
                column(6, numericInput('n', "n", 5, min = 0, max = 1e2, value=log10(max(modeldat$value) -min(modeldat$value))))
              ),
              fluidRow(
                column(4, p(HTML("<b>Step 4: Simulate NLS fitting</b>"),
                            span(shiny::icon("info-circle"), id = "info_uu"),
                            br(),
                            actionBttn("sim", "Fit NLS", color ="primary", style = "jelly"),
                            tippy::tippy_this(elementId = "info_uu",tooltip = "",placement = "right"))),
                column(4, p(HTML("<b>Step 5 (optional): Simulate BayesianMCMC fitting</b>"),
                            span(shiny::icon("info-circle"), id = "info_uu"),
                            br(),
                            actionBttn("simb", "Fit BayesianMCMC", color ="primary", style = "jelly"),
                            tippy::tippy_this(elementId = "info_uu",tooltip = "",placement = "right"))),
                column(4, p(HTML("<b>Step 6: Add to the record board</b>"),
                            span(shiny::icon("info-circle"), id = "info_uu"),
                            br(),
                            actionBttn("add", "Add",color ="warning", style = "jelly"),
                            tippy::tippy_this(elementId = "info_uu",tooltip = "",placement = "right")))
              )
          )
        })
 
    }
    if(as.numeric(input$gatemodel) == 1) {
      modeldat <- reshape.func(flofile, colnames(flofile)[2:length(colnames(flofile))])
        mydiv <-  tagList(
          box(title = "Parameterisation", width = 6, solidHeader = T, status = "primary",height = "300px",
              fluidRow(
                column(6, numericInput('k', "k (a.u.)", 5, min = 0, max = 1e6, value=max(modeldat$value))),
                column(6, numericInput('n1', "n", 5, min = 0, max = 1e2, value=log10(max(modeldat$value)- min(modeldat$value))))
              ),
              fluidRow(
                column(6, numericInput('K1', "K1 (M)", 5, min = 0, max = 1e5, value=modeldat$I[abs(modeldat$value - (max(modeldat$value)+ min(modeldat$value))/2) == min(abs(modeldat$value - (max(modeldat$value)+ min(modeldat$value))/2))])),
                column(6, numericInput('alpha', "α", 5, min = 0, max = 1, value=0))
              ),
              
              fluidRow(
                column(4, p(HTML("<b>Step 4: Simulate NLS fitting</b>"),
                            span(shiny::icon("info-circle"), id = "info_uu"),
                            br(),
                            actionBttn("sim", "Fit NLS", color ="primary", style = "jelly"),
                            tippy::tippy_this(elementId = "info_uu",tooltip = "",placement = "right"))),
                column(4, p(HTML("<b>Step 5 (optional): Simulate BayesianMCMC fitting</b>"),
                            span(shiny::icon("info-circle"), id = "info_uu"),
                            br(),
                            actionBttn("simb", "Fit BayesianMCMC", color ="primary", style = "jelly"),
                            tippy::tippy_this(elementId = "info_uu",tooltip = "",placement = "right"))),
                column(4, p(HTML("<b>Step 6: Add to the record board</b>"),
                            span(shiny::icon("info-circle"), id = "info_uu"),
                            br(),
                            actionBttn("add", "Add",color ="warning", style = "jelly"),
                            tippy::tippy_this(elementId = "info_uu",tooltip = "",placement = "right")))
              )
          )
        )
      
    }
    if(as.numeric(input$gatemodel) == 3) {
      modeldat <- rv$applysensor()
        mydiv <-   tagList({
          box(title = "Parameterisation", width = 6, solidHeader = T, status = "primary",height = "400px",
              fluidRow(
                column(6, numericInput('Gmax', "Gmax", 5, min = 0, max = 1e6, value=max(modeldat$Output))),
              ),
              fluidRow(
                column(6, numericInput('Kr', "Kr (M)", 5, min = 0, max = 1e6, value=modeldat$R[abs(modeldat$Output - (max(modeldat$Output)+ min(modeldat$Output))/2) == min(abs(modeldat$Output - (max(modeldat$Output)+ min(modeldat$Output))/2))])),
                column(6, numericInput('nr', "nr", 5, min = 0, max = 1e2, value=max(log10(modeldat$R)) - min(log10(modeldat$R))))
              ),
              fluidRow(
                column(6, numericInput('Ks', "Ks (M)", 5, min = 0, max = 1e5, value=modeldat$S[abs(modeldat$Output - (max(modeldat$Output)+ min(modeldat$Output))/2) == min(abs(modeldat$Output - (max(modeldat$Output)+ min(modeldat$Output))/2))])),
                column(6, numericInput('ns', "ns", 5, min = 0, max = 1, value=max(log10(modeldat$S)) - min(log10(modeldat$S))))
              ),
              
              
              fluidRow(
                column(4, p(HTML("<b>Step 4: Simulate NLS fitting</b>"),
                            span(shiny::icon("info-circle"), id = "info_uu"),
                            br(),
                            actionBttn("sim", "Fit NLS", color ="primary", style = "jelly"),
                            tippy::tippy_this(elementId = "info_uu",tooltip = "",placement = "right"))),
                column(4, p(HTML("<b>Step 5 (optional): Simulate BayesianMCMC fitting</b>"),
                            span(shiny::icon("info-circle"), id = "info_uu"),
                            br(),
                            actionBttn("simb", "Fit BayesianMCMC", color ="primary", style = "jelly"),
                            tippy::tippy_this(elementId = "info_uu",tooltip = "",placement = "right"))),
                column(4, p(HTML("<b>Step 6: Add to the record board</b>"),
                            span(shiny::icon("info-circle"), id = "info_uu"),
                            br(),
                            actionBttn("add", "Add",color ="warning", style = "jelly"),
                            tippy::tippy_this(elementId = "info_uu",tooltip = "",placement = "right")))
              )
          )
        })
    }
    }
    mydiv
    
  })
  output$nameControl <- renderUI({
    tagList(
      textInput('acc', "Name your gate here",
                value=paste0(rv$mygate(),".", sample(423423,1))))
  })
  output$ui_dlbtn_csv <- renderUI({
    tagList(downloadButton("dl_data_csv", "Download CSV"))
    
  })
  output$ui_dlbtn_xls <- renderUI({
    tagList(downloadButton("dl_data_xls", "Download EXCEL"))
    
  })
  output$andControl <- renderUI({
    if(as.numeric(input$gatemodel) == 3) {
      mydiv <- tagList({
        fluidRow(
          column(4, selectInput('andin1', 'Select sensor function for input1', choices = sensordb$element, selectize=TRUE)),
          column(4, selectInput('andin2', 'Select sensor function for input2', choices = sensordb$element, selectize=TRUE))
        )
      })
    }else{
      mydiv <- tagList()
    }
    
    mydiv
  })
  output$predControl1 <- renderUI({
    if(as.numeric(input$gatemodel) == 3) {
      mydiv <- tagList({
        tabBox(title = "", width = 6, height = "600px",
               id = "tabset1",
               tabPanel("Observed florenscent-inducer relation",
                        plotlyOutput("obs_plot_surface1",
                                     width = "100%",
                                     height = "600px") %>% withSpinner(type = 5)
               ),
               tabPanel("Observed Input-Output relation",
                        plotlyOutput("obs_plot_surface2",
                                     width = "100%",
                                     height = "600px") %>% withSpinner(type = 5)
               )
        )
        
      })
    }else{
      mydiv <- tagList(
        tabBox(
          title = "", width = 6, height = "600px",
          id = "tabset1",
          tabPanel("Observed Input-Output relation",
                   plotlyOutput("obs_plot") %>% withSpinner(type = 5)
                   
          ),
          tabPanel("Observed Output variation",
                   plotlyOutput("vio_plot") %>% withSpinner(type = 5)
                   
          )
        )
      )
    }
    
    mydiv
  })
  output$predControl2 <- renderUI({
    if(as.numeric(input$gatemodel) == 3) {
      if(is.null(rv$fitbmcmc)) {
      mydiv <- tagList({
        tabBox(title = "", width = 6, height = "600px",
               id = "tabset2",
               tabPanel("Model fitting summary",
                        verbatimTextOutput("fit_summary")
               ),
               tabPanel("Predicted Output-Input relation (x=R,y=S)",
                        plotlyOutput("pred_io",
                                     width = "100%",
                                     height = "600px") %>% withSpinner(type = 5)
               ),
               tabPanel("Predicted Output-Input relation (x=I1,y=I2)",
                        plotlyOutput("pred_io2",
                                     width = "100%",
                                     height = "600px") %>% withSpinner(type = 5)
               ),
               tabPanel("Prediction curve (R)",
                        plotlyOutput("pred_cur1") %>% withSpinner(type = 5)
               ),
               tabPanel("Prediction curve (S)",
                        plotlyOutput("pred_cur2") %>% withSpinner(type = 5)
               )
        )
      })
      }else{
        mydiv <- tagList({
          tabBox(title = "", width = 6, height = "600px",
                 id = "tabset2",
                 tabPanel("Model fitting summary",
                          verbatimTextOutput("fit_summary")
                 ),
                 tabPanel("Model fitting summary BMCMC",
                          verbatimTextOutput("fit_summary_bmcmc")
                 ),
                 tabPanel("Stan Codes",
                          div(style='height:400px;overflow-y: scroll;',
                              uiOutput("showstancodes"))
                 ),
                 tabPanel("Predicted Output-Input relation (x=R,y=S)",
                          plotlyOutput("pred_io",
                                       width = "100%",
                                       height = "600px") %>% withSpinner(type = 5)
                 ),
                 tabPanel("BMCMC Regression Properties",
                          selectInput("properyselect",NULL,
                                      choices = c("Density"=1,"Trace"=2,"Variation"=3)),
                          plotOutput("mcmcplot",
                                       width = "100%",
                                       height = "400px") %>% withSpinner(type = 5)
                 ),
                 tabPanel("Predicted Output-Input relation (x=I1,y=I2)",
                          plotlyOutput("pred_io2",
                                       width = "100%",
                                       height = "600px") %>% withSpinner(type = 5)
                 )
                 
          )
        })
      }
    }else{
      if(is.null(rv$fitbmcmc)) {
        mydiv <- tagList(
          tabBox(title = "", width = 6, height = "600px",
                 id = "tabset2",
                 tabPanel("Model fitting summary nls",
                          verbatimTextOutput("fit_summary")
                 ),
                 tabPanel("Predicted Input-Output relation",
                          plotlyOutput("pred_io") %>% withSpinner(type = 5)
                 )
          )
        )
      }else{
        mydiv <- tagList(
          tabBox(title = "", width = 6, height = "600px",
                 id = "tabset2",
                 tabPanel("Model fitting summary nls",
                          verbatimTextOutput("fit_summary")
                 ),
                 tabPanel("Model fitting summary BMCMC",
                          verbatimTextOutput("fit_summary_bmcmc")
                 ),
                 tabPanel("Stan Codes",
                          div(style='height:400px;overflow-y: scroll;',
                              uiOutput("showstancodes"))
                 ),
                 tabPanel("Predicted Input-Output relation",
                          plotlyOutput("pred_io") %>% withSpinner(type = 5)
                 ),
                 tabPanel("BMCMC Regression Properties",
                          selectInput("properyselect",NULL,
                                      choices = c("Density"=1,"Trace"=2,"Variation"=3)),
                          plotOutput("mcmcplot",
                                       width = "100%",
                                       height = "400px") %>% withSpinner(type = 5)
                 )
          )
        )
      }
      
    }
    
    mydiv
  })
  
  
  #####plot and table output functions
  
  output$gatefunc <- renderUI({
    print(input$gatemodel)
    shiny::withMathJax(
      shiny::helpText(paste0("$",
                             functionlist[as.numeric(input$gatemodel)],
                             "$")))
    
  })
  
  output$obs_plot <- renderPlotly({
    req(input$gatemodel)
    req(rv$obsdf())
    flofile <- rv$obsdf()
    
    s = input$flo_tab_rows_selected
    print(s)
    
    if(as.numeric(input$gatemodel) == 2) {
      obsdata <- reformat(flofile, "R3")
      
    }
    if(as.numeric(input$gatemodel) == 1) {
      obsdata <- reformat(flofile, "I")
      
    }
    if(as.numeric(input$gatemodel) == 4) {
      obsdata <- reformat(flofile, "x")
      
    }
    
    plot.obs(obsdata,s)
    
  })
  
  output$obs_plot_surface1 <- renderPlotly({
    req(rv$applysensor())
    plotme <- rv$applysensor()
    
    if(as.numeric(input$gatemodel) == 3) {
      plot.surface(unique(plotme$I1), unique(plotme$I2), plotme$Output)
    }
  })
  
  output$obs_plot_surface2 <- renderPlotly({
    req(rv$applysensor())
    
    
    plotme <- rv$applysensor()
    if(as.numeric(input$gatemodel) == 3) {
      plot.surface(unique(plotme$R), unique(plotme$S), plotme$Output)
    }
  })
  
  output$fit_summary <- renderPrint({
    req(input$gatemodel)
    req(rv$modelfit)
    
    print(summary(rv$modelfit))
    # plot.obs(obsdata,s)
    
  })
  
  output$priordist <- renderPrint({
    req(rv$modelfit)
    cat(call.priors(rv$modelfit))
  })
  output$fit_summary_bmcmc <- renderPrint({
    req(input$gatemodel)
    req(rv$fitbmcmc)
    
    print(rv$fitbmcmc,c(names(coef(rv$modelfit)),"sigma"),digits_summary = 8)
    
  })
  
  output$showstancodes <- renderUI({
    output$showstan <- renderPrint({
       req(input$gatemodel)
       req(rv$fitbmcmc)
    
       cat(paste0(stancodep1=stancodes()[[as.numeric(input$gatemodel)]][1], 
               call.priors(rv$modelfit),
               stancodep2=stancodes()[[as.numeric(input$gatemodel)]][2]
               ))
    
  })
    verbatimTextOutput("showstan")
  })
  
  output$pred_io <- renderPlotly({
    req(input$gatemodel)
    req(input$file1)
    req(rv$modelfit)
    
    if(is.null(rv$fitbmcmc)){
      if(as.numeric(input$gatemodel) == 2) {
        obsdata <- reformat(rv$obsdf(), "R3")
        myfig <- plot.pred(obstab=obsdata,
                           Invalue="R3",
                           calfunc=formula.not(),
                           model.fit=rv$modelfit,
                           mytext=functionlist[as.numeric(input$gatemodel)])
        
      }
      if(as.numeric(input$gatemodel) == 4) {
        obsdata <- reformat(rv$obsdf(), "x")
        myfig <-  plot.pred(obstab=obsdata,
                            Invalue="x",
                            calfunc=formula.not.cello(),
                            model.fit=rv$modelfit,
                            mytext=functionlist[as.numeric(input$gatemodel)])
        
      }
      if(as.numeric(input$gatemodel) == 1) {
        obsdata <- reformat(rv$obsdf(), "I")
        myfig <-  plot.pred(obstab=obsdata,
                            Invalue="I",
                            calfunc=formula.activ(),
                            model.fit=rv$modelfit,
                            mytext=functionlist[as.numeric(input$gatemodel)])
        
      }
      if(as.numeric(input$gatemodel) == 3) {
      #  print(rv$modelpred)
        myfig <- plot.surface(unique(rv$modelpred$R),unique(rv$modelpred$S), rv$modelpred$O, colorscl = "Jet")
        
      }
    }else{
      modeldat <- reshape.func(rv$obsdf(), colnames(rv$obsdf())[2:length(colnames(rv$obsdf()))])
      colnames(modeldat)[1] <- "x"
      
      if(as.numeric(input$gatemodel) == 2) {
        myfig <- pred.mcmc(rv$fitbmcmc,rv$modelfit,formula.not,modeldat)
      }
      if(as.numeric(input$gatemodel) == 4) {
        myfig <- pred.mcmc(rv$fitbmcmc,rv$modelfit,formula.not.cello,modeldat)
      }
      if(as.numeric(input$gatemodel) == 1) {
        myfig <-  pred.mcmc(rv$fitbmcmc,rv$modelfit,formula.activ,modeldat)
      }
      if(as.numeric(input$gatemodel) == 3) {
      #  print(rv$mcmcpars)
      #  print(rv$mcmcdf)
        myfig <-  plot.mcmc.surface(rv$modelpred,rv$mcmcdf,rv$applymodeland())
      }
    }
    
    myfig
    
  })
  
  output$pred_io2 <- renderPlotly({
    req(input$gatemodel)
    req(input$file1)
    req(rv$modelfit)
    
    if(as.numeric(input$gatemodel) == 3) {
      myfig <- plot.surface(unique(rv$modelpred$I1),unique(rv$modelpred$I2), rv$modelpred$O, colorscl = "Jet")
      
    }
    
    myfig
    
  })
  
  output$recordtab <- renderReactable({
    print(rv$record)
    req(rv$record)
     if(!is.null(rv$record)) {
       data <- unique(rv$record)
       
       reactable(data, details = function(index) {
         
         sub_data <- foreach::foreach(a=grep(data$element[index],names(rv$recordlist)),.combine = "rbind") %do% rv$recordlist[[a]]
         
         coe <- foreach::foreach(a=grep(data$element[index],names(rv$coef)),.combine = "rbind") %do% rv$coef[[a]]
         
         print(coe)
         coe <- unique(coe)
         
         text <- '
  colDef(style=function(value,index) {
      if (sub_data[, "subme"][index] > 0.05) {
        color <- "#e00000"
      } else if (sub_data[, "subme"][index] <= 0.01) {
        color <- "#008000"
      } else {
        color <- "#e8b600"
      }
      list(color = color, fontWeight = "bold")
})'
         
         myl <- foreach::foreach(i=coe) %do% { 
           eval(parse(text=gsub("subme",paste0(i,".sig"),text)))
         }
         
         names(myl) <- coe
         print(myl)
         
         htmltools::div(style = "padding: 1rem",
                        reactable(sub_data[,match(c(coe,paste0(coe,".se")),colnames(sub_data))], 
                                  columns = myl ,
                                  outlined = TRUE)
         )
       },onClick = "select",selection = "multiple")
     }
  })
  
  output$pred_cur1 <- renderPlotly({
    req(rv$modelpred)
    plot.and.curve(rv$modelpred, "R")
  })
  
  output$pred_cur2 <- renderPlotly({
    req(rv$modelpred)
    plot.and.curve(rv$modelpred, "S")
  })
  
  output$mcmcplot <- renderPlot({
    req(rv$modelfit)
    req(rv$fitbmcmc)
    
    if(input$properyselect ==1) {
      myplot <- stan_dens(rv$fitbmcmc, pars=names(coef(rv$modelfit)))
    }
    if(input$properyselect ==2) {
      myplot <- stan_trace(rv$fitbmcmc, pars=names(coef(rv$modelfit)))
    }
    if(input$properyselect ==3) {
      myplot <- stan_plot(rv$fitbmcmc, pars=names(coef(rv$modelfit)))
    }
    
    myplot
  })
  
  output$flo_tab <- DT::renderDataTable({
    req(input$gatemodel)
    req(input$file1)
    
    if(as.numeric(input$gatemodel) == 2) {
      flofile <- rv$obsdf()
      obsdata <- reformat(flofile, "R3")
      
    }
    if(as.numeric(input$gatemodel) == 1) {
      flofile <- rv$obsdf()
      obsdata <- reformat(flofile, "I")
      
    }
    if(as.numeric(input$gatemodel) == 4) {
      flofile <- rv$obsdf()
      obsdata <- reformat(flofile, "x")
      
    }
    if(as.numeric(input$gatemodel) == 3) {
      obsdata <- rv$applysensor()
      
    }
    
    #print(data)
    obsdata
  }, options = list(pageLength = 5))
  
  output$gate_resources <- DT::renderDataTable({
    rv$resource
    
  }, server = FALSE, options = list(pageLength = 10))
  
  output$gate_compare_radar = renderPlotly({
    s = input$gate_resources_rows_selected
    print(s)
    
    if(as.numeric(input$gatetype) == 1) {
      if(length(s) == 0) {
        myfig <- plot_ly(
          type = 'scatterpolar',
          fill = 'toself'
        )
      }else{
        myfig <-  compare.radar(s, rv$resource, c("k","n1","R2"), c("K1","alpha") )
      }
    }
    if(as.numeric(input$gatetype) == 4) {
      if(length(s) == 0) {
        myfig <- plot_ly(
          type = 'scatterpolar',
          fill = 'toself'
        )
      }else{
        myfig <-  compare.radar(s, rv$resource, c("ymax","n"), c("K","ymin") )
      }
    }
    if(as.numeric(input$gatetype) == 2) {
      if(length(s) == 0) {
        myfig <- plot_ly(
          type = 'scatterpolar',
          fill = 'toself'
        )
      }else{
        myfig <-  compare.radar(s, rv$resource, c("k3","n3","R2"), c("K3","alpha3") )
      }
    }
    if(as.numeric(input$gatetype) == 3) {
      if(length(s) == 0) {
        myfig <- plot_ly(
          type = 'scatterpolar',
          fill = 'toself'
        )
      }else{
        myfig <-  compare.radar(s, rv$resource, c("nr","ns","Gmax"), c("Kr","Ks") )
      }
    }
    
    myfig
    
  })
  
  output$gate_compare_io = renderPlotly({
    s = input$gate_resources_rows_selected
    print(s)
    
    if(as.numeric(input$gatetype) == 1) {
      if(length(s) == 0) {
        myfig <- plot_ly()
      }else{
        myfig <-  compare.io(s, rv$resource, c("k","n1", "K1","alpha"),"I","K1", formula.activ())
      }
    }
    if(as.numeric(input$gatetype) == 4) {
      if(length(s) == 0) {
        myfig <- plot_ly()
      }else{
        myfig <-  compare.io(s, rv$resource, c("ymax","ymin","K","n"),"x","K", formula.not.cello())
      }
    }
    if(as.numeric(input$gatetype) == 2) {
      if(length(s) == 0) {
        myfig <- plot_ly()
      }else{
        myfig <-  compare.io(s, rv$resource, c("k3","n3", "K3","alpha3"),"R3","K3", formula.not() )
      }
    }
    
    myfig
    
  })
  
  output$vio_plot <- renderPlotly({
    req(input$gatemodel)
    req(rv$obsdf())
    flofile <- rv$obsdf()
    
    if(as.numeric(input$gatemodel) == 2) {
      In <- "R3"
    }
    if(as.numeric(input$gatemodel) == 1) {
      In <- "I"
    }
    if(as.numeric(input$gatemodel) == 4) {
      In <- "x"
    }
    plot.violin(flofile,In)
  })
  
  output$dl_data_csv <- downloadHandler(
    filename = function() {
      paste0("geneticgates-", format(Sys.time(),"%Y%m%d-%H%M%S"), ".csv")
    },
    content = function(file) {
      state <- getReactableState("recordtab","selected")
      print(state)
      
      data <- unique(rv$record)
      dllist <- list()
      
      dllist <- foreach::foreach(i=state) %do% {foreach::foreach(a=grep(data$element[i],names(rv$recordlist)),.combine = "rbind") %do% rv$recordlist[[a]]}

      write.csv(rv$record, file, quote = F)
        
    })
  
  output$dl_data_xls <- downloadHandler(
    filename = function() {
      paste0("geneticgates-", format(Sys.time(),"%Y%m%d-%H%M%S"), ".xls")
    },
    content = function(file) {
      state <- getReactableState("recordtab","selected")
      print(state)
      
      data <- unique(rv$record)
      dllist <- list()
      
      dllist <- foreach::foreach(i=state) %do% {foreach::foreach(a=grep(data$element[i],names(rv$recordlist)),.combine = "rbind") %do% rv$recordlist[[a]]}
      
      print(dllist)
      lapply(dllist, function(x) {
        print(x)
        xlsx::write.xlsx(x,file,sheetName = x$element[1],append=T)
      })
      
    })
  
  ###reactive shiny data
  
  rv <- reactiveValues(file1=NULL,fitbmcmc=NULL,record=NULL,mcmcpars=NULL,recordlist=list(),coef=list())
  #rv$record <- reactive(NULL)
  
  rv$origin <- reactive({
    if(as.numeric(input$gatetype) == 1) {
      df <- sensordb
    }
    if(as.numeric(input$gatetype) == 2) {
      df <-  notgatedb
    }
    if(as.numeric(input$gatetype) == 4) {
      df <-  celllo.notgatedb
    }
    if(as.numeric(input$gatetype) == 3) {
      df <-  andgatedb
    }
    
    df
  })
  
  #rv$resource <- reactive(rv$origin())
  
  observeEvent(input$addto, {
    
    if(as.numeric(input$gatetype) == as.numeric(input$gatemodel) ) {
      if(nrow(rv$record) > 0){
        df <- rv$record[,colnames(rv$record)!="gate_type"]
        roundlist <- round(df[1,colnames(df)!="element"], digits = 5)
        mydf <- as.data.frame(c(element=df$element[1], roundlist))
        mydf <- mydf[,match(colnames(rv$origin()), colnames(mydf))]
        
        df <- rbind(rv$origin(),mydf)
      }else{
        df <- rv$origin()
      }
      
    }else{
      df <- rv$origin()
    }
    
    rv$resource <- df
  })
  
  observeEvent(input$gatetype, {
    rv$resource <- rv$origin()
  })
  
  rv$obsdf <- eventReactive({
    input$file1
    input$gatemodel
    }, {
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext %in% c("csv","xlsx","xls"), "Please upload a file in either csv or excel format"))
    
    if(ext == "csv") {
      tab <- read.table(file$datapath, header = T, sep=",")
    }
    if(ext %in% c("xlsx","xls")) {
      tab <- read_excel(file$datapath,col_names = T)
    }
    if(input$gatemodel ==1) {
    colnames(tab)[match("I",colnames(tab))] <- "I"
    }
    
    if(input$gatemodel ==2) {
      colnames(tab)[match("I",colnames(tab))] <- "R3"
    }
    if(input$gatemodel ==4) {
      colnames(tab)[match("I",colnames(tab))] <- "x"
    }
    
    tab
    
  })
  
  observeEvent(list(input$file1,input$gatemodel),{
    rv$fitbmcmc <- NULL
    rv$mcmcpars <- NULL
  })
  
  rv$applysensor <- reactive({
    req(input$andin1)
    req(input$andin2)
    req(rv$obsdf())
    
    df <- rv$obsdf()
    
    data_mean <- ddply(df, c("I1","I2"), summarise, Output = round(mean(O), digits = 4))
    data_sd <- ddply(df, c("I1","I2"), summarise, Output = round(sd(O), digits = 4))
    data <- data.frame(data_mean,sd= data_sd$Output)
    
    R <- calculator(formula.activ(),
                    variables = list(I=data$I2),
                    params = sensordb[match(input$andin1, sensordb$element),match(c("k","n1","K1", "alpha"),colnames(sensordb))] )
    data$R <- round(R, digits = 4)
    S <- calculator(formula.activ(),
                    variables = list(I=data$I2),
                    params = sensordb[match(input$andin2, sensordb$element),match(c("k","n1","K1", "alpha"),colnames(sensordb))] )
    data$S <- round(S, digits = 4)
    
    data[, match(c("I1","I2","R", "S", "Output","sd"), colnames(data))]
    
  })
  
  rv$applymodeland <- reactive({
    req(input$andin1)
    req(input$andin2)
    req(rv$obsdf())
    
    df <- rv$obsdf()
    
    R <- calculator(formula.activ(),
                    variables = list(I=df$I1),
                    params = sensordb[match(input$andin1, sensordb$element),match(c("k","n1","K1", "alpha"),colnames(sensordb))] )
    df$R <- round(R, digits = 4)
    S <- calculator(formula.activ(),
                    variables = list(I=df$I2),
                    params = sensordb[match(input$andin2, sensordb$element),match(c("k","n1","K1", "alpha"),colnames(sensordb))] )
    df$S <- round(S, digits = 4)
    
    df
    
    
    
  })
  
  rv$mygate <- reactive({
    mygates[as.numeric(input$gatemodel)]
  })
  
  observeEvent(input$sim, {
    req(input$gatemodel)
    req(input$file1)
    
    
    if(as.numeric(input$gatemodel) == 2) {
      flofile <- rv$obsdf()
      try( nsl.fit <- nsl.fit.notgate(flofile, k3=input$k3, alpha3=input$alpha3, K3=input$K3, n3=input$n3))
      modeldata <- reshape.func(flofile, colnames(flofile)[2:length(colnames(flofile))])
      rv$modelr2 <- cor(fitted(nsl.fit), modeldata$value)
    }
    if(as.numeric(input$gatemodel) == 4) {
      flofile <- rv$obsdf()
      try( nsl.fit <- nsl.fit.notgate.cello(flofile, ymax=input$ymax, K=input$K, n=input$n, ymin=input$ymin))
      modeldata <- reshape.func(flofile, colnames(flofile)[2:length(colnames(flofile))])
      rv$modelr2 <- cor(fitted(nsl.fit), modeldata$value)
    }
    if(as.numeric(input$gatemodel) == 1) {
      flofile <- rv$obsdf()
      try( nsl.fit <- nsl.fit.sensor(flofile, k=input$k, alpha=input$alpha,  n1=input$n1, K1=input$K1))
      modeldata <- reshape.func(flofile, colnames(flofile)[2:length(colnames(flofile))])
      rv$modelr2 <- cor(fitted(nsl.fit), modeldata$value)
    }
    if(as.numeric(input$gatemodel) == 3) {
      flofile <- rv$applymodeland()
      try(  nsl.fit <- nsl.fit.and(flofile, Gmax=input$Gmax,Kr=input$Kr,nr=input$nr,Ks=input$Ks,ns=input$ns))
      
      rv$modelr2 <- cor(fitted(nsl.fit), flofile$O)
    }
    
    rv$modelfit <- nsl.fit
    
    print(rv$modelr2)
    
    fitdf <- as.data.frame(matrix(round(as.numeric(summary(rv$modelfit)$coefficients[,c(1,2,4)]),digits = 12),nrow = 1))
    colnames(fitdf) <- c(names(coef(rv$modelfit)), paste0(names(coef(rv$modelfit)),".se"),paste0(names(coef(rv$modelfit)),".sig"))
    fitdf$R2 <- rv$modelr2
    fitdf$element <- input$acc
    fitdf$gate_type <- mygates[as.numeric(input$gatemodel)]
    fitdf$fit_method <- "Nonlinear_least_Square"
    fitdf$formula <- paste0("$",functionlist[1],"$")
    fitdf$R2 <- rv$modelr2
    
    rownames(fitdf) <- paste0(input$acc,".NLS")
    rv$fitdf <- fitdf
  })
  
  observeEvent(input$gomcmc,{
    req(rv$modelfit)
    
    show_modal_spinner(
      spin = "orbit",
      # color = "firebrick",
      color = "royalblue",
      text = "MCMC SAMPLING FOR MODEL ..."
    )
    
    if(as.numeric(input$gatemodel) %in% c(1,2,4)) {
      
      if(as.numeric(input$gatemodel) == 2) {
        df <- reshape.func(rv$obsdf(), colnames(rv$obsdf())[2:length(colnames(rv$obsdf()))])
        xvalue <- df$R3
        yvalue <- df$value
        
      }
      if(as.numeric(input$gatemodel) == 1) {
        df <- reshape.func(rv$obsdf(), colnames(rv$obsdf())[2:length(colnames(rv$obsdf()))])
        xvalue <- df$I
        yvalue <- df$value
      }
      if(as.numeric(input$gatemodel) == 4) {
        df <- reshape.func(rv$obsdf(), colnames(rv$obsdf())[2:length(colnames(rv$obsdf()))])
        xvalue <- df$x
        yvalue <- df$value
      }
      
      print(xvalue)
      print(yvalue)
      
      rv$fitbmcmc <- call.bmcmc(xvalue=xvalue,
                                yvalue=yvalue, 
                                stancodep1=stancodes()[[as.numeric(input$gatemodel)]][1], 
                                stancodep2=stancodes()[[as.numeric(input$gatemodel)]][2], 
                                call.priors(rv$modelfit),
                                chains=input$chains,
                                iter=input$iter,
                                cores=input$cores)
      
      print(rv$fitbmcmc,names(coef(rv$modelfit)),digits_summary = 8)
      
      df <- as.data.frame(summary(rv$fitbmcmc)$summary)
      
      rv$mcmcpars <- t(df[match(names(coef(rv$modelfit)),rownames(df)),])[1,]
      
    }else{
      if(as.numeric(input$gatemodel) %in% c(3)) {
        flofile <- rv$applymodeland()
        print(flofile)
        xvalue <- flofile$R
        zvalue <- flofile$S
        yvalue <- flofile$O
      }
      rv$fitbmcmc <- call.bmcmc.2i(xvalue=xvalue,zvalue=zvalue,yvalue=yvalue, 
                                   stancodep1=stancodes()[[as.numeric(input$gatemodel)]][1], 
                                   stancodep2=stancodes()[[as.numeric(input$gatemodel)]][2], 
                                   priors= call.priors(rv$modelfit),
                                   chains=input$chains,
                                   iter=input$iter,
                                   cores=input$cores)
      
      df <- as.data.frame(summary(rv$fitbmcmc)$summary)
      
      rv$mcmcpars <- t(df[match(names(coef(rv$modelfit)),rownames(df)),])[1,]
    }
    
    coe <- names(coef(rv$modelfit))
    fitdf <- as.data.frame(matrix(round(c(summary(rv$fitbmcmc,coe)$summary[,1],
                                          summary(rv$fitbmcmc,coe)$summary[,2], 
                                          abs(summary(rv$fitbmcmc,coe)$summary[,10]-1)),
                                        digits = 12),nrow = 1))
    colnames(fitdf) <- c(coe, paste0(coe,".se"),paste0(coe,".sig"))
    fitdf$R2 <- cor(rep(apply(rstan::extract(rv$fitbmcmc, "Y_mean")$Y_mean, 2, mean),
                        length(yvalue)/length(apply(rstan::extract(rv$fitbmcmc, "Y_mean")$Y_mean, 2, mean))),
                    yvalue)
    fitdf$element <- input$acc
    fitdf$gate_type <- mygates[as.numeric(input$gatemodel)]
    fitdf$fit_method <- "Bayesian_MCMC"
    fitdf$formula <- paste0("$",
                            functionlist[as.numeric(input$gatemodel)],
                            "$")
    
    rownames(fitdf) <- paste0(input$acc,".BMCMC")
    rv$fitdf <- fitdf
    
    remove_modal_spinner()
    removeModal()
  })
  
  observeEvent(input$add, {
    req(rv$modelfit)
    
    mydf <- data.frame(element=input$acc,
                       gate_type=mygates[as.numeric(input$gatemodel)],
                       formula=paste0("$",functionlist[as.numeric(input$gatemodel)],"$"))
    
    rv$record <- rbind(rv$record,mydf)
    
    coe <- names(coef(rv$modelfit))
    rec2 <- c("element","fit_method", coe, paste0(coe,".se"), paste0(coe,".sig"), "R2")
    oldnames <- names(rv$recordlist)
    rv$recordlist <- c(rv$recordlist, list(rv$fitdf[,match(rec2,colnames(rv$fitdf))]))
    names(rv$recordlist) <- c(oldnames,rownames(rv$fitdf)[1])
    
    oldnames <- names(rv$coef)
    rv$coef <- c(rv$coef, list(names(coef(rv$modelfit))))
    names(rv$coef) <- c(oldnames,rownames(rv$fitdf)[1])
  })
  
  observeEvent(list(input$sim, input$gomcmc), {
    req(input$gatemodel)
    req(rv$modelfit)
    
    if(as.numeric(input$gatemodel) == 3) {
      flofile <- rv$obsdf()
      
      x <- generate.log.samples(unique(flofile$I1))
      y <- generate.log.samples(unique(flofile$I2))
      
      tmpvar <- list(R=calculator(formula.activ(),
                                  variables = list(I=x),
                                  params = sensordb[match(input$andin1, sensordb$element),match(c("k","n1","K1", "alpha"),colnames(sensordb))] ),
                     S=calculator(formula.activ(),
                                  variables = list(I=y),
                                  params = sensordb[match(input$andin2, sensordb$element),match(c("k","n1","K1", "alpha"),colnames(sensordb))] ))
      
      
      anddf <- data.frame(I1= x[rep(1:length(x), time=length(y))],
                          I2=y[rep(1:length(y), each=length(x))],
                          R=tmpvar[[1]][rep(1:length(tmpvar[[1]]), time=length(tmpvar[[2]]))],
                          S=tmpvar[[2]][rep(1:length(tmpvar[[2]]), each=length(tmpvar[[1]]))],
                          O=calculator(formula.and(),
                                       variables = list(R=tmpvar[[1]][rep(1:length(tmpvar[[1]]), time=length(tmpvar[[2]]))],
                                                        S=tmpvar[[2]][rep(1:length(tmpvar[[2]]), each=length(tmpvar[[1]]))]),
                                       params =  coef(rv$modelfit)))
      
      rv$modelpred <- anddf
      
      if(!is.null(rv$mcmcpars)){
        print(rv$mcmcpars)
        rv$mcmcdf
        anddf <- data.frame(I1= x[rep(1:length(x), time=length(y))],
                                I2=y[rep(1:length(y), each=length(x))],
                                R=tmpvar[[1]][rep(1:length(tmpvar[[1]]), time=length(tmpvar[[2]]))],
                                S=tmpvar[[2]][rep(1:length(tmpvar[[2]]), each=length(tmpvar[[1]]))],
                                O=calculator(formula.and(), 
                                             variables = list(R=tmpvar[[1]][rep(1:length(tmpvar[[1]]), time=length(tmpvar[[2]]))],
                                                              S=tmpvar[[2]][rep(1:length(tmpvar[[2]]), each=length(tmpvar[[1]]))]), 
                                             params =  rv$mcmcpars))
        
        rv$mcmcdf <- anddf
        print(rv$mcmcdf)
      }
    }
  })
  
  observeEvent(input$btn_remove,{
    if(nrow(rv$record) > 0){
      rv$record <- dplyr::slice(rv$record, 1:nrow(rv$record)-1)
    }
  })
  
  dataModal <- function(failed = FALSE) {
    modalDialog(size="m",title="Bayesian MCMC Regression Options",
      helpText("Run Bayesian MCMC regression with priors below:"
      ),
      verbatimTextOutput("priordist"),
      numericInput("chains","Number of chains:",
                   value = 4,
                   min = 1,
                   max = 16),
      numericInput("iter","Number of iterations:",
                   value = 2000,
                   min = 0,
                   max = 10000),
      numericInput("cores","Number of cores for sampling:",
                   value = parallel::detectCores(),
                   min = 1,
                   max = 8),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("gomcmc", "GO")
      )
    )
  }
  
  # Show modal when button is clicked.
  observeEvent(input$simb, {
    showModal(dataModal())
  })
})
