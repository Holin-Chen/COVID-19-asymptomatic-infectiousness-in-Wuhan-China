library(EpiModel)
library(ggplot2)

covid19_SEIR <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    ## Derived values
    num <- s.num + e.num + i.num + r.num
    R0 <- 2.2
    if (t >= 55) {
      R0 <- 1.05
    }
    if (t >= 70) {
      R0 <- 0.41
    }
    
    omega <- 200000/14000000
    kappa <- 0.6*omega
    if (t >= 55) {
      omega <- 0
    }
    beta <- R0/(num*i.dur) 
    lambda <- beta * (i.num + thita*e.num)
    gamma <- 1/e.dur
    rho <- 1/i.dur
    severe_prop <- 0.35
    if (t >= 77) {
      severe_prop <- 0.2
    }
    
    ## Differential equations
    dS <- -lambda*s.num - omega*s.num + kappa*s.num
    dE <- lambda*s.num - gamma*e.num - omega*e.num
    dI <- gamma*e.num - rho*i.num - severe_prop*muI*i.num
    dR <- rho*i.num
    dD <- severe_prop*muI*i.num
    
    ## Output
    list(c(dS, dE, dI, dR, dD,
           se.flow = lambda * s.num,
           ei.flow = gamma * e.num,
           ir.flow = rho * i.num))
  })
}

param <- param.dcm(e.dur = 6.4, i.dur = 12, thita = c(0.25,0.35,0.45,0.55,0.65,0.75,0.85), muI = 1/16.1)
init <- init.dcm(s.num = 14000000, e.num = 150, i.num =50, r.num = 0, d.num=0,
                 se.flow = 0, ei.flow = 0, ir.flow = 0)
control <- control.dcm(nsteps = 150, dt=1, new.mod = covid19_SEIR)

mod <- dcm(param, init, control)
mod
df_covid19 <- as.data.frame(mod)
#Adding the variable of date into the dataframe
df_covid19$date <- format(seq(as.Date("2019-12-01"), as.Date("2020-04-28"), by = "1 day"))
df_covid19$date <- as.Date(df_covid19$date)
#Adding the variable of thita into the dataframe
df_covid19$thita[df_covid19$run == 1] <- 0.25
df_covid19$thita[df_covid19$run == 2] <- 0.35
df_covid19$thita[df_covid19$run == 3] <- 0.45
df_covid19$thita[df_covid19$run == 4] <- 0.55
df_covid19$thita[df_covid19$run == 5] <- 0.65
df_covid19$thita[df_covid19$run == 6] <- 0.75
df_covid19$thita[df_covid19$run == 7] <- 0.85
df_covid19$thita <- as.factor(df_covid19$thita)

##Adding the variable of cumulated cases into dataframe
df_covid19$cum_cases[1]<- df_covid19$i.num[1] + df_covid19$ei.flow[1]
df_covid19$cum_cases[151]<- df_covid19$i.num[1] + df_covid19$ei.flow[151]
df_covid19$cum_cases[301]<- df_covid19$i.num[1] + df_covid19$ei.flow[301]
df_covid19$cum_cases[451]<- df_covid19$i.num[1] + df_covid19$ei.flow[451]
df_covid19$cum_cases[601]<- df_covid19$i.num[1] + df_covid19$ei.flow[601]
df_covid19$cum_cases[751]<- df_covid19$i.num[1] + df_covid19$ei.flow[751]
df_covid19$cum_cases[901]<- df_covid19$i.num[1] + df_covid19$ei.flow[901]

for (i in (1:1049)){
  if (df_covid19$thita[i] == df_covid19$thita[i+1]){
    df_covid19$cum_cases[i+1] <- df_covid19$cum_cases[i] + df_covid19$ei.flow[i+1]}
  else if (df_covid19$thita[i] != df_covid19$thita[i+1]){
    df_covid19$cum_cases[i+1] <- df_covid19$i.num[1] + df_covid19$ei.flow[i+1]
  }}

##Draw the pictures of the COVID-19 dynamics from the 7 simulations by ggplot2
p1 <- ggplot(df_covid19, aes(date, i.num, group=thita, color = thita)) +
  geom_line(size = 0.8) +
  geom_point(pch = 21, size = 1)  +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = c(0.9,0.85),
    plot.margin = margin(3,15,3,3,"mm")
  ) +
  labs(x = "Date", y= "Infected Number",
       subtitle = "Real-Time Infected Population with Different Estimated Thita")
p1

p2 <- ggplot(df_covid19, aes(date, ei.flow, group=thita, color = thita)) +
  geom_line(size = 0.8) +
  geom_point(pch = 21, size = 1)  +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = c(0.9,0.85),
    plot.margin = margin(3,15,3,3,"mm")
  ) +
  labs(x = "Date", y= "Incidence",
       subtitle = "Incident Cases with Different Estimated Thita")
p2

p3 <- ggplot(df_covid19, aes(date, cum_cases, group = thita, color = thita)) +
  geom_line(size = 0.8) +
  geom_point(pch = 21, size = 1)  +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = c(0.1,0.5),
    plot.margin = margin(3,15,3,3,"mm")
  ) +
  labs(x = "Date", y= "Cumulated Cases",
       subtitle = "Cumulated Cases with Different Estimated Thita")
p3

p4 <- ggplot(df_covid19, aes(date, d.num, group = thita, color = thita)) +
  geom_line(size = 0.8) +
  geom_point(pch = 21, size = 1)  +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = c(0.1,0.5),
    plot.margin = margin(3,15,3,3,"mm")
  ) +
  labs(x = "Date", y= "Death Number",
       subtitle = "Mortality with Different Estimated Thita")
p4

##Define the dates when infected number and incidence reach the peaks
df_covid19$date[which.max(df_covid19$i.num[df_covid19$thita==0.25])] #2020-02-06
df_covid19$date[which.max(df_covid19$i.num[df_covid19$thita==0.35])] #2020-02-08
df_covid19$date[which.max(df_covid19$i.num[df_covid19$thita==0.45])] #2020-02-09
df_covid19$date[which.max(df_covid19$i.num[df_covid19$thita==0.55])] #2020-02-09
df_covid19$date[which.max(df_covid19$i.num[df_covid19$thita==0.65])] #2020-02-09
df_covid19$date[which.max(df_covid19$i.num[df_covid19$thita==0.75])] #2020-02-10
df_covid19$date[which.max(df_covid19$i.num[df_covid19$thita==0.85])] #2020-02-10


df_covid19$date[which.max(df_covid19$ei.flow[df_covid19$thita==0.25])] #2020-01-23
df_covid19$date[which.max(df_covid19$ei.flow[df_covid19$thita==0.35])] #2020-01-23
df_covid19$date[which.max(df_covid19$ei.flow[df_covid19$thita==0.45])] #2020-01-23
df_covid19$date[which.max(df_covid19$ei.flow[df_covid19$thita==0.55])] #2020-01-24
df_covid19$date[which.max(df_covid19$ei.flow[df_covid19$thita==0.65])] #2020-02-07
df_covid19$date[which.max(df_covid19$ei.flow[df_covid19$thita==0.75])] #2020-02-07
df_covid19$date[which.max(df_covid19$ei.flow[df_covid19$thita==0.25])] #2020-02-07

#Add the vertical lines of dates when reaching the peaks
p1 + geom_vline(xintercept = as.Date("2020-02-09"), col="hotpink")
p2 + geom_vline(xintercept = as.Date("2020-01-23"), col="hotpink") + 
  geom_vline(xintercept= as.Date("2020-02-07"), col="indianred")

library(nCov2019)#a new package that can crawl the updated COVID-19 reported data
x <- load_nCov2019()
wh <- subset(x['Hubei'], city == 'Wuhan') #creat the dataframe of reported data in Wuhan

#Draw the curves of cumulated cases calibrated with reported data
max(df_covid19$cum_cases[df_covid19$thita==0.55]) #help define the ylim 
par(mfrow = c(1,1), mar = c(3,3,1,1), mgp = c(2,1,0))
plot(df_covid19$date[df_covid19$thita==0.45],df_covid19$cum_cases[df_covid19$thita==0.45], 
     ylim = c(0,61800), type="l", col="red",lwd=2, xlab="Date", ylab="Cumulated Cases")
title(main="Cumulated Cases under Thita of 0.45, 0.55 and Reported Cumulated Cases", cex.main=0.78)
lines(df_covid19$date[df_covid19$thita==0.55],df_covid19$cum_cases[df_covid19$thita==0.55], type="l", col="cadetblue",lwd=2)
lines(wh$time[1:150],wh$cum_confirm[1:150], type="l", col="green", lwd=2)
legend("bottomright", legend = c("thita=0.55","reported cases","thita=0.45"), 
       col=c("red","green","cadetblue"),lty = 1, lwd = 2)

#Draw the curves of death calibrated with reported data
max(df_covid19$d.num[df_covid19$thita==0.55]) #help define the ylim 
par(mfrow = c(1,1), mar = c(3,3,1,1), mgp = c(2,1,0))
plot(df_covid19$date[df_covid19$thita==0.45],df_covid19$d.num[df_covid19$thita==0.45],ylim=c(0,10700), 
     type="l", col="red",lwd=2,xlab="Date", ylab="Cumulated Death")
title(main="Cumulated Death under Thita of 0.45, 0.55 and Reported Cumulated Death", cex.main=0.78)
lines(df_covid19$date[df_covid19$thita==0.55],df_covid19$d.num[df_covid19$thita==0.55], type="l", col="cadetblue",lwd=2)
lines(wh$time[1:150],wh$cum_dead[1:150], type="l", col="green", lwd=2)
legend("topleft", legend = c("thita=0.55","reported cases","thita=0.45"), col=c("red","green","cadetblue"),lty = 1, lwd = 2)

##Change the simulating days to 350 
control <- control.dcm(nsteps = 350, dt=1, new.mod = covid19_SEIR)

mod <- dcm(param, init, control)
mod
covid19_SEIR_350 <- as.data.frame(mod)
covid19_SEIR_350$date <- format(seq(as.Date("2019-12-01"), as.Date("2020-11-14"), by = "1 day"))
covid19_SEIR_350$date <- as.Date(covid19_SEIR_350$date)
covid19_SEIR_350$thita[covid19_SEIR_350$run == 1] <- 0.25
covid19_SEIR_350$thita[covid19_SEIR_350$run == 2] <- 0.35
covid19_SEIR_350$thita[covid19_SEIR_350$run == 3] <- 0.45
covid19_SEIR_350$thita[covid19_SEIR_350$run == 4] <- 0.55
covid19_SEIR_350$thita[covid19_SEIR_350$run == 5] <- 0.65
covid19_SEIR_350$thita[covid19_SEIR_350$run == 6] <- 0.75
covid19_SEIR_350$thita[covid19_SEIR_350$run == 7] <- 0.85

#See when the COVID-19 would be eliminated
covid19_SEIR_350$date[which(covid19_SEIR_350$i.num[covid19_SEIR_350$thita==0.45] < 1)[1]]
covid19_SEIR_350$date[which(covid19_SEIR_350$i.num[covid19_SEIR_350$thita==0.55] < 1)[1]]

max(covid19_SEIR_350$i.num[covid19_SEIR_350$thita==0.55])
plot(covid19_SEIR_350$date[covid19_SEIR_350$thita==0.45], covid19_SEIR_350$i.num[covid19_SEIR_350$thita==0.45], 
     ylim=c(0,12000), type="l", col="cadetblue",lwd=2,xlab="Date", ylab="Infected Number", 
     main="Infected Number with Different thita Estimation")
lines(covid19_SEIR_350$date[covid19_SEIR_350$thita==0.55], covid19_SEIR_350$i.num[covid19_SEIR_350$thita==0.55], col="red", lwd=2)
abline(v= as.Date("2020-09-12"), col="hotpink")
abline(v= as.Date("2020-10-03"), col="indianred")
legend("topright", legend = c("thita=0.45", "thita=0.55","elimination date 2020-09-12 for thita=0.45", 
                              "elimination date 2020-10-03 for thita=0.55"  ), 
       col = c("cadetblue", "red", "hotpink", "indianred"),  lty = 1, lwd = 2,cex = 0.75)

#Build the counterfactual model without any intervention
covid19_SEIR <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    ## Derived values
    num <- s.num + e.num + i.num + r.num
    R0 <- 2.2
    
    
    omega <- 200000/14000000
    kappa <- 0.6*omega
    
    beta <- R0/(num*i.dur) 
    lambda <- beta * (i.num + thita*e.num)
    gamma <- 1/e.dur
    rho <- 1/i.dur
    severe_prop <- 0.35
    
    ## Differential equations
    dS <- -lambda*s.num - omega*s.num + kappa*s.num
    dE <- lambda*s.num - gamma*e.num - omega*e.num
    dI <- gamma*e.num - rho*i.num - severe_prop*muI*i.num
    dR <- rho*i.num
    dD <- severe_prop*muI*i.num
    
    ## Output
    list(c(dS, dE, dI, dR, dD,
           se.flow = lambda * s.num,
           ei.flow = gamma * e.num,
           ir.flow = rho * i.num))
  })
}

param <- param.dcm(e.dur = 6.4, i.dur = 12, thita = c(0.25,0.35,0.45,0.55,0.65,0.75,0.85), muI = 1/16.1)
init <- init.dcm(s.num = 14000000, e.num = 150, i.num =50, r.num = 0, d.num=0,
                 se.flow = 0, ei.flow = 0, ir.flow = 0)
control <- control.dcm(nsteps = 350, dt=1, new.mod = covid19_SEIR)


mod <- dcm(param, init, control)
mod
df_nointervention <- as.data.frame(mod)
#Adding the variable of date into the dataframe
df_nointervention$date <- format(seq(as.Date("2019-12-01"), as.Date("2020-11-14"), by = "1 day"))
df_nointervention$date <- as.Date(df_nointervention$date)
#Adding the variable of thita into the dataframe
df_nointervention$thita[df_nointervention$run == 1] <- 0.25
df_nointervention$thita[df_nointervention$run == 2] <- 0.35
df_nointervention$thita[df_nointervention$run == 3] <- 0.45
df_nointervention$thita[df_nointervention$run == 4] <- 0.55
df_nointervention$thita[df_nointervention$run == 5] <- 0.65
df_nointervention$thita[df_nointervention$run == 6] <- 0.75
df_nointervention$thita[df_nointervention$run == 7] <- 0.85

##Adding the variable of cumulated cases into dataframe
df_nointervention$cum_cases[1]<- df_nointervention$i.num[1] + df_nointervention$ei.flow[1]
df_nointervention$cum_cases[351]<- df_nointervention$i.num[1] + df_nointervention$ei.flow[351]
df_nointervention$cum_cases[701]<- df_nointervention$i.num[1] + df_nointervention$ei.flow[701]
df_nointervention$cum_cases[1051]<- df_nointervention$i.num[1] + df_nointervention$ei.flow[1051]
df_nointervention$cum_cases[1401]<- df_nointervention$i.num[1] + df_nointervention$ei.flow[1401]
df_nointervention$cum_cases[1751]<- df_nointervention$i.num[1] + df_nointervention$ei.flow[1751]
df_nointervention$cum_cases[2101]<- df_nointervention$i.num[1] + df_nointervention$ei.flow[2101]

for (i in (1:2449)){
  if (df_nointervention$thita[i] == df_nointervention$thita[i+1]){
    df_nointervention$cum_cases[i+1] <- df_nointervention$cum_cases[i] + df_nointervention$ei.flow[i+1]}
  else if (df_nointervention$thita[i] != df_nointervention$thita[i+1]){
    df_nointervention$cum_cases[i+1] <- df_nointervention$i.num[1] + df_nointervention$ei.flow[i+1]
  }}

##Draw the pictures of the COVID-19 dynamics from the 7 simulations without any intervention
max(df_nointervention$i.num)
plot(df_nointervention$date[df_nointervention$run==1], df_nointervention$i.num[df_nointervention$run==1], ylim=c(0, 1500000), 
     type="l", col="red",lwd=2,xlab="Date", ylab="Infected number", 
     main="Infected Number with Different thita Estimation without any Intervention")
lines(df_nointervention$date[df_nointervention$run==2], df_nointervention$i.num[df_nointervention$run==2], col="orange", lwd=2)
lines(df_nointervention$date[df_nointervention$run==3], df_nointervention$i.num[df_nointervention$run==3], col="brown", lwd=2)
lines(df_nointervention$date[df_nointervention$run==4], df_nointervention$i.num[df_nointervention$run==4], col="green", lwd=2)
lines(df_nointervention$date[df_nointervention$run==5], df_nointervention$i.num[df_nointervention$run==5], col="cyan", lwd=2)
lines(df_nointervention$date[df_nointervention$run==6], df_nointervention$i.num[df_nointervention$run==6], col="blue", lwd=2)
lines(df_nointervention$date[df_nointervention$run==7], df_nointervention$i.num[df_nointervention$run==7], col="purple", lwd=2)
legend("topright", legend = c("thita=0.25","thita=0.35","thita=0.45", "thita=0.55","thita=0.65","thita=0.75","thita=0.85"), 
       col = c("red","orange","brown","green","cyan","blue","purple"),  lty = 1, lwd = 2, cex=0.75)

max(df_nointervention$ei.flow)
plot(df_nointervention$date[df_nointervention$run==1], df_nointervention$ei.flow[df_nointervention$run==1], ylim=c(0, 190000), 
     type="l", col="red",lwd=2,xlab="Date", ylab="Incidence", 
     main="Incident Cases with Different thita Estimation without any Intervention")
lines(df_nointervention$date[df_nointervention$run==2], df_nointervention$ei.flow[df_nointervention$run==2], col="orange", lwd=2)
lines(df_nointervention$date[df_nointervention$run==3], df_nointervention$ei.flow[df_nointervention$run==3], col="brown", lwd=2)
lines(df_nointervention$date[df_nointervention$run==4], df_nointervention$ei.flow[df_nointervention$run==4], col="green", lwd=2)
lines(df_nointervention$date[df_nointervention$run==5], df_nointervention$ei.flow[df_nointervention$run==5], col="cyan", lwd=2)
lines(df_nointervention$date[df_nointervention$run==6], df_nointervention$ei.flow[df_nointervention$run==6], col="blue", lwd=2)
lines(df_nointervention$date[df_nointervention$run==7], df_nointervention$ei.flow[df_nointervention$run==7], col="purple", lwd=2)
legend("topright", legend = c("thita=0.25","thita=0.35","thita=0.45", "thita=0.55","thita=0.65","thita=0.75","thita=0.85"), 
       col = c("red","orange","brown","green","cyan","blue","purple"),  lty = 1, lwd = 2, cex=0.75)

max(df_nointervention$cum_cases)
plot(df_nointervention$date[df_nointervention$run==1], df_nointervention$cum_cases[df_nointervention$run==1], ylim=c(0, 7151400),
     type="l", col="red",lwd=2,xlab="Date", ylab="Cumulated Cases", 
     main="Cumulated Cases with Different thita Estimation without any Intervention")
lines(df_nointervention$date[df_nointervention$run==2], df_nointervention$cum_cases[df_nointervention$run==2], col="orange", lwd=2)
lines(df_nointervention$date[df_nointervention$run==3], df_nointervention$cum_cases[df_nointervention$run==3], col="brown", lwd=2)
lines(df_nointervention$date[df_nointervention$run==4], df_nointervention$cum_cases[df_nointervention$run==4], col="green", lwd=2)
lines(df_nointervention$date[df_nointervention$run==5], df_nointervention$cum_cases[df_nointervention$run==5], col="cyan", lwd=2)
lines(df_nointervention$date[df_nointervention$run==6], df_nointervention$cum_cases[df_nointervention$run==6], col="blue", lwd=2)
lines(df_nointervention$date[df_nointervention$run==7], df_nointervention$cum_cases[df_nointervention$run==7], col="purple", lwd=2)
legend("topleft", legend = c("thita=0.25","thita=0.35","thita=0.45", "thita=0.55","thita=0.65","thita=0.75","thita=0.85"), 
       col = c("red","orange","brown","green","cyan","blue","purple"),  lty = 1, lwd = 2, cex=0.75)

max(df_nointervention$d.num)
plot(df_nointervention$date[df_nointervention$run==1], df_nointervention$d.num[df_nointervention$run==1], ylim=c(0, 1480000), 
     type="l", col="red",lwd=2,xlab="Date", ylab="Death", 
     main="Death Number with Different thita Estimation without any Intervention")
lines(df_nointervention$date[df_nointervention$run==2], df_nointervention$d.num[df_nointervention$run==2], col="orange", lwd=2)
lines(df_nointervention$date[df_nointervention$run==3], df_nointervention$d.num[df_nointervention$run==3], col="brown", lwd=2)
lines(df_nointervention$date[df_nointervention$run==4], df_nointervention$d.num[df_nointervention$run==4], col="green", lwd=2)
lines(df_nointervention$date[df_nointervention$run==5], df_nointervention$d.num[df_nointervention$run==5], col="cyan", lwd=2)
lines(df_nointervention$date[df_nointervention$run==6], df_nointervention$d.num[df_nointervention$run==6], col="blue", lwd=2)
lines(df_nointervention$date[df_nointervention$run==7], df_nointervention$d.num[df_nointervention$run==7], col="purple", lwd=2)
legend("topleft", legend = c("thita=0.25","thita=0.35","thita=0.45", "thita=0.55","thita=0.65","thita=0.75","thita=0.85"), 
       col = c("red","orange","brown","green","cyan","blue","purple"),  lty = 1, lwd = 2, cex=0.75)

##Calculate the NIA and PIA of COVID-19 comparing with having and no intervention
cum_inc_0.45 <- sum(covid19_SEIR_350$ei.flow[covid19_SEIR_350$thita == 0.45]) + covid19_SEIR_350$i.num[1]
cum_inc_0.45_nointv <- sum(df_nointervention$ei.flow[df_nointervention$thita == 0.45]) + df_nointervention$i.num[1]
NIA_0.45 <- cum_inc_0.45_nointv - cum_inc_0.45
PIA_0.45 <- NIA_0.45 / cum_inc_0.45_nointv
cum_inc_0.45 #35215
cum_inc_0.45_nointv #5192639
NIA_0.45 #5157424
PIA_0.45 #0.9932184

cum_inc_0.55 <- sum(covid19_SEIR_350$ei.flow[covid19_SEIR_350$thita==0.55]) + covid19_SEIR_350$i.num[1]
cum_inc_0.55_nointv <- sum(df_nointervention$ei.flow[df_nointervention$thita==0.55]) + df_nointervention$i.num[1]
NIA_0.55 <- cum_inc_0.55_nointv - cum_inc_0.55 
PIA_0.55 <- NIA_0.55 / cum_inc_0.55_nointv
cum_inc_0.55 #62548
cum_inc_0.55_nointv #5746607
NIA_0.55 #5684059
PIA_0.55 #0.9891157
