BEGIN {
  OFS="\t"
  print "loop", "time_index", "lambda"
}
{
  split($3, a, ",")
  for(i=1; i <= length(a); i++) {
    print NR, i, a[i]
  }
}