get_niche_attributes<-function(consumer_size, consumer_category) {
  
  if(consumer_category=="vertebrate"){
     # Unwrap the input parameters
  qrsup = Param_regvert[[2]]
  qrinf = Param_regvert[[3]]
  
  # Estimate the parameters for the niche model
  n = consumer_size						# The niche n
  low = qrinf[1] + qrinf[2]*consumer_size	# The lower limit of the range
  high = qrsup[1] + qrsup[2]*consumer_size	# The higher limit of the range
  c = low+(high-low)/2			# The centroid c
  }
  
  if(consumer_category=="invertebrate"){
    # Unwrap the input parameters
    qrsup = Param_reginvert[[2]]
    qrinf = Param_reginvert[[3]]
    
    # Estimate the parameters for the niche model
    n = consumer_size						# The niche n
    low = qrinf[1] + qrinf[2]*consumer_size	# The lower limit of the range
    high = qrsup[1] + qrsup[2]*consumer_size	# The higher limit of the range
    c = low+(high-low)/2			# The centroid c
  }
  
  return(cbind(n,c,low,high))	
}
