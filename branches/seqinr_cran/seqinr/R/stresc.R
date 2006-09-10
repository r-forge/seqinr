stresc <- function(string = "->\\<->{<->}<->$<->^<->_<->%<->#<->&<->~<-")
{
  fromchar <- s2c("\\{}$^_%#&~")
  tochar <- c("$\\backslash$", "\\{", "\\}", "\\$", "\\^{}", "\\_", "\\%", "\\#", "\\&", "\\~{}")
  c2s( sapply(s2c(string), function(x) ifelse(x %in% fromchar, tochar[which(x == fromchar)], x)))
}
