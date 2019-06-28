yeastFasta = function (orf) {
  print(orf)
  url=paste0(urlpieces[1],orf,urlpieces[2])
  print(paste0("url=",url))
  readDNAStringSet(url)
}