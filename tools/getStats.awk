function isSegSite(str) {
  split(str, c, "")
  for(i = 2; i <= length(c); i += 1) {
    if(c[i] != c[1])
      return 1
  }
  return 0
}

{
  called_sites += $3
  split($4, allele_str, ",")
  if(NR==1)
    nr_haplotypes = length(allele_str[1])
  for(a in allele_str) {
    if(isSegSite(allele_str[a])) {
      segsites += 1
      break
    }
  }
}

END {
  watterson_fac = 0
  for(i = 1; i < nr_haplotypes; i++)
    watterson_fac += 1.0 / i
  printf("Segregating Sites:\t%s\n", segsites)
  printf("Called Sites:\t%s\n", called_sites)
  printf("Watterson's Theta:\t%s\n", (segsites / called_sites) / watterson_fac)
}