
iref14=read.table("~/Work/Data/iRefindex_14.0/9606.mitab.04072015.txt",sep="\t",header = T,as.is = T,fill = T)
human_human_list = data.frame(iref14$taxa,iref14$taxb)
tmp = do.call(`paste`, c(unname(human_human_list), list(sep=".")))
iref14_human = iref14[tmp == "taxid:9606(Homo sapiens).taxid:9606(Homo sapiens)" | tmp == "-.taxid:9606(Homo sapiens)",] #Filters human-human PPIs
filter_UP=intersect(grep("uniprotkb:",iref14_human$uidA),grep("uniprotkb:",iref14_human$uidB)) # To filter interactions with known Uniprot IDs
selfLoops=which(iref14_human_UP$uidA==iref14_human_UP$uidB) # To remove 12163 self loops and edges
iref14_human_UP=iref14_human[filter_UP,c(1:4,15,39:40)][-selfLoops,] # Applying all filters
iref14_ppi=graph.data.frame(iref14_human_UP,directed = F)
iref14_np_confidence=as.numeric(gsub("np:","",lapply(strsplit(grep("np:[[:digit:]]{1,}",iref14_human_UP$confidence,value = T),split="\\|"),`[[`,3)))
# Preserve nodes with known PDB ids
hs_UP=UniProt.ws(taxId = 9606)
columns=c(columns(hs_UP)[c(116,86,37,56)])
kt=(columns(hs_UP)[116])
V(iref14_ppi)$name=gsub(pattern = "uniprotkb:",replacement = "",V(iref14_ppi)$name)
keys=V(iref14_ppi)$name
res=select(hs_UP,keys,columns,kt)
iref14_ppi_PDB=delete_vertices(iref14_ppi,V(iref14_ppi)[-which(V(iref14_ppi)$name%in%res[(is.na(res$PDB)==T),1])])

