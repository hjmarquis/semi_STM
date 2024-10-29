wgt.sd = function(x,w,na.rm=T)
{
  sqrt(mean((x-mean(x*w, na.rm = na.rm))^2*w,na.rm = na.rm))
}

feature.descr = function(x,...,nsmall = 0)
{
  out = rep("", ...length())
  
  if(is.logical(x))
  {
    for(i in 1:...length())
    {
      out[i] = paste(format(sum(x * ...elt(i), na.rm = T),digits = 0,scientific =F, nsmall = nsmall),' (',
                     format(mean(x * ...elt(i), na.rm = T)*100, digits = 0, nsmall = nsmall),
                     "\\%)", sep='')
    }
  }else
  {
    for(i in 1:...length())
    {
      out[i] = paste(format(mean(x * ...elt(i), na.rm = T),digits = 0, nsmall = nsmall),' (',
                     format(wgt.sd(x,...elt(i), na.rm = T), digits = 0, nsmall = nsmall),
                     ")", sep='')
    }
  }
  return(out)
}

feature.descr.simple = function(x,...,nsmall = 0)
{
  out = rep("", ...length())
  
  if(is.logical(x))
  {
    for(i in 1:...length())
    {
      out[i] = paste(format(mean(x * ...elt(i), na.rm = T)*100, digits = 0, nsmall = nsmall),
                     "\\%", sep='')
    }
  }else
  {
    for(i in 1:...length())
    {
      out[i] = format(mean(x * ...elt(i), na.rm = T),digits = 0, nsmall = nsmall)
    }
  }
  return(out)
}