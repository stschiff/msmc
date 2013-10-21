BEGIN {
  OFS="\t"
  print "loop", "time_index", "lambda"
}
{
  split($4, a, ",")
  T = length(a) / 3
  for(i=1; i <= T; i++) {
    print NR, i, 2 * a[3 * i - 1] / (a[3 * i - 2] + a[3 * i])
  }
}