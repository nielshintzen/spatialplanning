# NEED TO READ IN PREDGAM.CSV FILE !



k<-0.85
Pam<-390
pm<-19.4
Em<-2500
Eg<-5600
TA<-7000
Tref<-283
TAlow<-50000
TL<-277
TAhigh<-100000
TH<-297
m<-0.219
foodh <- 0.000069


Growth_seasonal_lim <- function(temp, food, L1) {

                f_2 <- 0.000075 # calibrated
                V1 <- (m*L1)^3
                f <- ((46*f_2*food)/(foodh+(46*f_2*food))) * f_1
                G <- Growth <- ((k*(f)*(Pam*
                (exp((TA/Tref)-(TA/(273+temp)))
                *((1+exp((TAlow/Tref)-(TAlow/TL))+exp((TAhigh/(TH-(0.1*L1)))-(TAhigh/Tref)))
                /(1+exp((TAlow/(273+temp))-(TAlow/TL))+exp((TAhigh/(TH-(0.1*L1)))-(TAhigh/(273+temp))))))))
                *((m*L1)^3)^(2/3)
                -((pm*(exp((TA/Tref)-(TA/(273+(temp+0))))))*((m*L1)^3)))/
                ((k*(f)*Em)+Eg)

                V2<- ((m*L1)^3)+G         # Volume achieved after day of growth
                L2<- (V2^(1/3))/m         # Backcalculating volume to Length
                G2<- L2-L1                # difference in length after weekly growth
                return(G2)
}

#Looping Growth through length matrix by time
# f1 defined by stomach fullness percentages - gam output

for(i in 1:z) {

  f_1 <- predgam$f1[i]

  Length_summer[,,i+1]<-Length_summer[,,i]+(Growth_seasonal_lim(temp_summer[,,i],food_summer[,,i],Length_summer[,,i]))
  }
