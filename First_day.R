rnorm(100)
my_rn<-rnorm(100)
mean(my_rn)
sd(my_rn)
my_rn+5
mean(my_rn+5)
speakers<-sample(c(1:25))
speakers
speakers<-sample(c(1:25),3)
speakers
speakers<-sample(c(1:25),3)
speakers
speakers<-sample(c(1:25),3)
speakers<-sample(c(1:25))
day1<-speakers[1:13]
day2<-speakers[14:25]
day1
rep("A",10)
my_rn[sd(my_rn)>1]
my_rn[c(sd(my_rn)>1)]
sd(my_rn)>1
my_rn>(mean(my_rn)-sd(my_rn))
my_rn[my_rn>(mean(my_rn)-sd(my_rn))]
my_rn.onesd<-my_rn[my_rn>(mean(my_rn)-sd(my_rn))]
hist(my_rn)
abline(1,0)
abline(0,1)
hist(my_rn)
abline(-1)
abline(h=-1)
abline(h=-1.1)
hist(my_rn)
abline(h=-1)
abline(h=-2)
abline(v = -1)
abline(v = c(-1,1))
abline(v= c((mean(my_rn)-sd(my_rn)),(mean(my_rn)+sd(my_rn)))
)
