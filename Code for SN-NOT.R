#Main functions

cumsum_plinear_constrast <- function(ts, type, start_t, trim_size, total_length){
  ### this function calculates the normalizer L or R depending on the type 
  n <- length(ts)
  end_t <- start_t+n-1
  time_index <- (start_t:end_t)/total_length
  a_result <- b_result <- rep(0,n)
  if(type=='L'){
    # first sum
    ydm <- ts
    xdm <- time_index
    y_bar <- cumsum(ydm)/(1:n)
    x_bar <- cumsum(xdm)/(1:n)
    xy <- cumsum(ydm*xdm)
    # y2 <- cumsum(ydm^2)
    x2 <- cumsum(xdm^2)
    b_cumsum_f <- (xy-(1:n)*x_bar*y_bar)/(x2-(1:n)*x_bar^2) # slope
    a_cumsum_f <- y_bar-x_bar*b_cumsum_f # intercept
    
    # second sum
    ydm <- rev(ts)
    xdm <- rev(time_index)
    y_bar <- cumsum(ydm)/(1:n)
    x_bar <- cumsum(xdm)/(1:n)
    xy <- cumsum(ydm*xdm)
    # y2 <- cumsum(ydm^2)
    x2 <- cumsum(xdm^2)
    b_cumsum_b <- (xy-(1:n)*x_bar*y_bar)/(x2-(1:n)*x_bar^2)
    a_cumsum_b <- y_bar-x_bar*b_cumsum_b
    b_cumsum_b <- rev(b_cumsum_b)
    a_cumsum_b <- rev(a_cumsum_b)
    
    # contrast
    b_result[2:(n-2)] <- b_cumsum_f[2:(n-2)]-b_cumsum_b[3:(n-1)]
    a_result[2:(n-2)] <- a_cumsum_f[2:(n-2)]-a_cumsum_b[3:(n-1)]
  }else{
    # first sum
    ydm <- ts
    xdm <- time_index
    y_bar <- cumsum(ydm)/(1:n)
    x_bar <- cumsum(xdm)/(1:n)
    xy <- cumsum(ydm*xdm)
    # y2 <- cumsum(ydm^2)
    x2 <- cumsum(xdm^2)
    b_cumsum_f <- (xy-(1:n)*x_bar*y_bar)/(x2-(1:n)*x_bar^2) # slope
    a_cumsum_f <- y_bar-x_bar*b_cumsum_f # intercept
    
    # second sum
    ydm <- rev(ts)
    xdm <- rev(time_index)
    y_bar <- cumsum(ydm)/(1:n)
    x_bar <- cumsum(xdm)/(1:n)
    xy <- cumsum(ydm*xdm)
    # y2 <- cumsum(ydm^2)
    x2 <- cumsum(xdm^2)
    b_cumsum_b <- (xy-(1:n)*x_bar*y_bar)/(x2-(1:n)*x_bar^2)
    a_cumsum_b <- y_bar-x_bar*b_cumsum_b
    b_cumsum_b <- rev(b_cumsum_b)
    a_cumsum_b <- rev(a_cumsum_b)
    
    # contrast
    b_result[3:(n-1)] <- b_cumsum_b[3:(n-1)]-b_cumsum_f[2:(n-2)]
    a_result[3:(n-1)] <- a_cumsum_b[3:(n-1)]-a_cumsum_f[2:(n-2)]
  }
  # trimming contrast
  if(trim_size>0){
    b_result[1:trim_size] <- 0
    b_result[(n-trim_size+1):n] <- 0
    a_result[1:trim_size] <- 0
    a_result[(n-trim_size+1):n] <- 0
  }
  return(cbind(a_result, b_result))
}


####################################

SN_test_plinear <- function(ts, k, pre_grid_position, trim_size, total_length){
  ###this function calculates the test statistic at each k with given subsample data 
  n <- length(ts)
  start_t <- pre_grid_position
  end_t <- start_t+n-1
  time_index <- (start_t:end_t)/total_length
  no_stat <- 2 # two statistics
  
  ### the numerator D
  ### for single change point
  ### eq 2.1 in paper: t1=1, t2=n, k is the cpt pos.
  data_x <- cbind(1, time_index)
  tmp1 <- data_x[1:k,]
  D1 <- as.numeric(solve(t(tmp1)%*%tmp1)%*%t(tmp1)%*%ts[1:k])
  tmp2 <- data_x[(k+1):n,]
  D2 <- as.numeric(solve(t(tmp2)%*%tmp2)%*%t(tmp2)%*%ts[(k+1):n])
  
  D <- k*(n-k)/n^1.5*(D1-D2)
  
  #### the denominator V 
  inter1 <- cumsum_plinear_constrast(ts[1:k], 'L', pre_grid_position, trim_size, total_length)
  inter2 <- cumsum_plinear_constrast(ts[(k+1):n], 'R', pre_grid_position+k, trim_size, total_length)
  multiplier1 <- ((1:k)*((k-1):0))^2/n^2/k^2
  multiplier2 <- ((0:(n-k-1))*((n-k):1))^2/n^2/(n-k)^2
  M1 <- M2 <- matrix(0, no_stat, no_stat)
  for(index1 in 1:(no_stat-1)){
    for(index2 in (index1+1):no_stat){
      M1[index1, index2] <- M1[index2, index1] <- sum(inter1[,index1]*inter1[,index2]*multiplier1)
      M2[index1, index2] <- M2[index2, index1] <- sum(inter2[,index1]*inter2[,index2]*multiplier2)
    }
  }
  for(index1 in 1:no_stat){
    M1[index1, index1] <- sum(inter1[,index1]^2*multiplier1)
    M2[index1, index1] <- sum(inter2[,index1]^2*multiplier2)
  }
  test_SN <- t(D)%*%solve(M1+M2)%*%D
  return(test_SN)
}

###############################

notsn <- function(s,e,y,M,grid_size,trim_size, total_length, critic){
  ### this function is the NOT-SN, s is the start of the interval, e is the end, M is the random sample numbers, total_length 
  ###is the full data length 
  ### set.seed(8) since NOT has some randomness
  cpt <- NULL
  if ( (e-s) < (2 * grid_size)){
    return(cpt)
  }
  flag <- e-s   ## used to find the the narrowest interval 
  seg <- matrix(0,nrow=2,ncol=M)  
  for (m in 1:M)
  {
    
    temp<- c()
    seg[,m] <- sort(floor(runif(2,s,e))) ## random subinterval 
    a <- seg[1,m]
    b <- seg[2,m]
    if ( (b-a) <2*grid_size) {
      next;
    }
    for (k in (a+grid_size) : (b-grid_size)){
      temp <- c(temp,SN_test_plinear(y[a:b], k-a+1, a, trim_size, total_length))
    }
    if( ( max(temp) > critic ) & ( (b-a) <flag ) ){
      flag <- b-a+1 # update the narrowest interval 
      cpt <- which.max(temp)+grid_size+a-1  # update change-point location in subsample
    }
  }
  if (!is.null(cpt)){
    cpt <- c(notsn(s,cpt,y,M,grid_size,trim_size,total_length,critic),cpt,notsn(cpt+1,e,y,M,grid_size,trim_size,total_length,critic))
  }
  return(cpt)
}

########################################################################
#Threshold calculation functions

WBSset=function(M,total_length)
{
  seg=matrix(0,nrow=2,ncol=M)
  for (m in 1:M){
    seg[,m]=sort(floor(runif(2,1,total_length+1)))
  }
  return(seg)
}

bootcritic=function(B,M,total_length,grid_size,trim_size){
  bootval=c()
  for (b in 1:B)
  { 
    y=rnorm(total_length,0,1)
    seg=WBSset(M,total_length)
    cvalb=c()
    for (m in 1:M)
    {
      a=seg[1,m]
      b=seg[2,m]
      if ((b-a)<2*grid_size) {
        next;
      }
      temp=c()
      for (k in (a+grid_size):(b-grid_size)){
        temp=c(temp,SN_test_plinear(y[a:b], k-a+1, a, trim_size, total_length))
      }
      cvalb=c(cvalb,max(temp))
    }
    bootval=c(bootval,max(cvalb))
    print(max(cvalb))
  }
  return(bootval)
}

########################################################################
#Applying change-point model to cumulative cases

library(latex2exp)

set.seed(8)

#Imported dataset named Malta

n <- nrow(Malta)
date <- as.Date(Malta$date)
tot_cases <- Malta$total_cases
plot(tot_cases,xaxt = "n",xlab="Date",ylab="Cumulative Total Cases",cex.lab=1,type="l")
axis(1,c(1,floor(n/5)*(1:4),n),date[c(1,floor(n/5)*(1:4),n)],cex.axis=1)
y <- log(tot_cases)
plot(y,xaxt = "n",xlab="Date",ylab="Log of Cumulative Total Cases",cex.lab=1,type="l")
axis(1,c(1,floor(n/5)*(1:4),n),date[c(1,floor(n/5)*(1:4),n)],cex.axis=1)

grid_size=floor(n*0.045) #h in SN-NOT where 0.05=epsilon
trim_size=floor(n*0.02) #d in SN-NOT where 0.02=delta

#Finding the threshold of the SN-NOT for cumulative cases
B=1000
M=1000
bootvalue=bootcritic(B,M,n,grid_size,trim_size)
critic=quantile(bootvalue,0.95)

#Applying SN-NOT
x1 <- notsn(1,n,y,M, grid_size,trim_size,n,critic)

cpt <- c(0,x1,length(y))
m <- length(cpt)-2 #Number of estimated change-points
begseg <- date[cpt[1:(m+1)]+1] #Beginning of segments
endseg <- date[cpt] #End of segments
num_seg <- rep(0,length(cpt)-1)
slopesn <- rep(0,length(cpt)-1)
interceptsn <- rep(0,length(cpt)-1)
res <- c()
fit_log <- c()
plot(y,xaxt = "n",xlab="Date",ylab="Log of Cumulative Total Cases",cex.lab=1)
axis(1,c(1,floor(n/5)*(1:4),n),date[c(1,floor(n/5)*(1:4),n)],cex.axis=1)
for (i in 1:(length(cpt)-1)){
  x <-((cpt[i]+1):(cpt[i+1]))/length(y)
  num_seg[i] <- length(y[(cpt[i]+1):(cpt[i+1])])
  reg <- lm(y[(cpt[i]+1):(cpt[i+1])]~x)
  slopesn[i] <- reg$coefficients[2]
  interceptsn[i] <- reg$coefficients[1]
  res <- c(res,reg$residuals)
  fit_log <- c(fit_log, reg$fitted.values)
  fix <- ((cpt[i]+1):(cpt[i+1]))
  #Fitted model of piecewise linear of log-scale cumulative cases
  lines(fix,reg$fitted.values,col="cyan",lwd=2)
}
mean(num_seg) #Average number of days of segments
normslope <- slopesn/n

change <-c()
for (i in 2:(length(cpt)-1)){
  points(cpt[i],y[cpt[i]],pch=16,cex=1.5)
  change <- c(change,slopesn[i]/slopesn[i-1]-1)
}
change <- c(NA,change)

#Fitting the exponential of the linear model
plot(tot_cases,xaxt = "n",xlab="Date",ylab="Cumulative Total Cases",cex.lab=1)
axis(1,c(1,floor(n/5)*(1:4),n),date[c(1,floor(n/5)*(1:4),n)],cex.axis=1)
lines(exp(fit_log),col="cyan",lwd=2)
for (i in 2:(length(cpt)-1)){
  points(cpt[i],tot_cases[cpt[i]],pch=16,cex=1.5)
}


#Data frame for summary
df <- data.frame(begseg, endseg, num_seg, interceptsn, normslope, change)
write.csv(df,"summary statistics cases.csv")

#Plotting residuals
plot(res,xaxt = "n",xlab="Date",ylab="Residuals",cex.lab=1,type='l')
axis(1,c(1,floor(n/5)*(1:4),n),date[c(1,floor(n/5)*(1:4),n)],cex.axis=1)
acf(res,lag.max = 40, cex.lab=1)
mean(res)

########################################################################
#Applying change-point model to cumulative deaths

#Subsetting the dataset from when cumul. deaths exceed 10
Maltasub <- subset(Malta,total_deaths>10)
n2 <- nrow(Maltasub)
date2 <- as.Date(Maltasub$date)
tot_deaths <- Maltasub$total_deaths
plot(tot_deaths,xaxt = "n",xlab="Date",ylab="Cumulative Total Deaths",cex.lab=1,type="l")
axis(1,c(1,floor(n2/5)*(1:4),n2),date2[c(1,floor(n2/5)*(1:4),n2)],cex.axis=1)
z <- log(tot_deaths)
plot(z,xaxt = "n",xlab="Date",ylab="Log of Cumulative Total Cases",cex.lab=1,type="l")
axis(1,c(1,floor(n2/5)*(1:4),n2),date2[c(1,floor(n2/5)*(1:4),n2)],cex.axis=1)
grid_size2=floor(n2*0.095) #h in SN-NOT where 0.05=epsilon
trim_size2=floor(n2*0.02) #d in SN-NOT where 0.02=delta

#Finding the threshold of the SN-NOT for cumulative deaths
bootvalue2=bootcritic(B,M,n2,grid_size2,trim_size2)
critic2=quantile(bootvalue2,0.95)

#Applying SN-NOT
x2 <- notsn(1,n2,z,M, grid_size2,trim_size2,n2,critic2)

cpt2 <- c(0,x2,length(z))
m2 <- length(cpt2)-2 #Number of estimated change-points
begseg2 <- date2[cpt2[1:(m2+1)]+1] #Beginning of segments
endseg2 <- date2[cpt2] #End of segments
num_seg2 <- rep(0,length(cpt2)-1)
slopesn2 <- rep(0,length(cpt2)-1)
interceptsn2 <- rep(0,length(cpt2)-1)
res2 <- c()
fit_log2 <- c()
plot(z,xaxt = "n",xlab="Date",ylab="Log of Cumulative Total Deaths",cex.lab=1)
axis(1,c(1,floor(n2/5)*(1:4),n2),date2[c(1,floor(n2/5)*(1:4),n2)],cex.axis=1)
for (i in 1:(length(cpt2)-1)){
  x <-((cpt2[i]+1):(cpt2[i+1]))/length(z)
  num_seg2[i] <- length(z[(cpt2[i]+1):(cpt2[i+1])])
  reg <- lm(z[(cpt2[i]+1):(cpt2[i+1])]~x)
  slopesn2[i] <- reg$coefficients[2]
  interceptsn2[i] <- reg$coefficients[1]
  res2 <- c(res2,reg$residuals)
  fit_log2 <- c(fit_log2, reg$fitted.values)
  fix <- ((cpt2[i]+1):(cpt2[i+1]))
  #Fitted model of piecewise linear of log-scale cumulative Deaths
  lines(fix,reg$fitted.values,col="orange",lwd=2)
}
mean(num_seg2) #Average number of days of segments
normslope2 <- slopesn2/n2

change2 <-c()
for (i in 2:(length(cpt2)-1)){
  points(cpt2[i],z[cpt2[i]],pch=16,cex=1.5)
  change2 <- c(change2,slopesn2[i]/slopesn2[i-1]-1)
}
change2 <- c(NA,change2)

#Fitting the exponential of the linear model
plot(tot_deaths,xaxt = "n",xlab="Date",ylab="Cumulative Total Deaths",cex.lab=1)
axis(1,c(1,floor(n2/5)*(1:4),n2),date2[c(1,floor(n2/5)*(1:4),n2)],cex.axis=1)
lines(exp(fit_log2),col="orange",lwd=2)
for (i in 2:(length(cpt2)-1)){
  points(cpt2[i],tot_deaths[cpt2[i]],pch=16,cex=1.5)
}


#Data frame for summary
df2 <- data.frame(begseg2, endseg2, num_seg2, interceptsn2, normslope2, change2)
write.csv(df2,"summary statistics deaths.csv")

#Plotting residuals
plot(res2,xaxt = "n",xlab="Date",ylab="Residuals",cex.lab=1,type='l')
axis(1,c(1,floor(n2/5)*(1:4),n2),date2[c(1,floor(n2/5)*(1:4),n2)],cex.axis=1)
acf(res2,lag.max = 40, cex.lab=1)
mean(res2)
