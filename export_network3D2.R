ExportToNetwork3D2 <- function(community, spmet,dir, fname_root=CP(community,'title'))
{
links <- cbind(PredatorID=NodeNameIndices(community, TLPS(community)[,'consumer']),
PreyID=NodeNameIndices(community, TLPS(community)[,'resource']),Weights=c(round(((community$trophic.links$weight)^2),2)))

write.table(links, file.path(dir, paste(fname_root,'_links.web', sep='')),
row.names=FALSE, sep=' ')
species <- cbind(ID=1:NumberOfNodes(community),
CommonName=NP(community,'node'))
write.table(links, file.path(dir, paste(fname_root,'_species.txt', sep='')),row.names=FALSE, sep=' ')

  Info<-data.frame(Index=c(seq(from=1,to=nrow(spmet),by=1)),Sp=spmet[,1],Name=spmet[,2],Size=c(round(spmet[,4],2)))#
  write.table(Info, file.path(dir, paste(fname_root,'_speciesInfo.txt', sep='')),row.names=FALSE, sep=' ')
  
}