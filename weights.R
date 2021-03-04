weights<-function(prey_sizes,niche_center,niche_min,niche_max){
  W<-dnorm(prey_sizes,mean=niche_center,
      sd=sd(seq(from=niche_min,to=niche_max,
                by=((niche_max-niche_min)/100))))/
  max(dnorm(prey_sizes,mean=niche_center,
            sd=sd(seq(from=niche_min,to=niche_max,
                      by=((niche_max-niche_min)/100)))))
W/sum(W)}
