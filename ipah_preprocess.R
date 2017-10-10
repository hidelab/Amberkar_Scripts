library(data.table)
library(tximport)
library(org.Hs.eg.db)

system("cp --parents `find ./final-full/ -name 'quant.sf'` /shared/hidelab2/user/md4zsa/Work/Data/IPAH/")
quant_countsList=list()
quant_files=system("find . -name 'quant.sf'",intern=T)
for(i in 1:length(quant_files)){
   quant_countsList[[i]]=fread(quant_files[i],header=T,sep="\t",showProgress=T,data.table=F,stringsAsFactors=F)
}
quant_countsList2=lapply(quant_countsList,function(x){x$Name<-unlist(lapply(strsplit(x=x$Name,split="\\."),`[[`,1));x})
for(i in 1:length(quant_files)){
  fwrite(x=quant_countsList2[[i]],file=quant_files[i],sep="\t",col.names=T,row.names=F,buffMB=100,nThread=4,showProgress=T)
}
