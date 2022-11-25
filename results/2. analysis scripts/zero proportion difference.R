remove(list=ls())

dt=readRDS("pbmc_10x-v2.rds")

dt=dt[,-(which(colSums(dt[,-1])==0)+1)]

celltype=names(table(dt[,1]))


dt_1=dt[dt[,1]==celltype[1],-1]
dt_2=dt[dt[,1]==celltype[2],-1]
dt_3=dt[dt[,1]==celltype[3],-1]
dt_4=dt[dt[,1]==celltype[4],-1]
dt_5=dt[dt[,1]==celltype[5],-1]
dt_6=dt[dt[,1]==celltype[6],-1]
dt_7=dt[dt[,1]==celltype[7],-1]
dt_8=dt[dt[,1]==celltype[8],-1]
dt_9=dt[dt[,1]==celltype[9],-1]


dt_p1=as.vector(colSums(dt_1==0))
dt_p2=as.vector(colSums(dt_2==0))
dt_p3=as.vector(colSums(dt_3==0))
dt_p4=as.vector(colSums(dt_4==0))
dt_p5=as.vector(colSums(dt_5==0))
dt_p6=as.vector(colSums(dt_6==0))
dt_p7=as.vector(colSums(dt_7==0))
dt_p8=as.vector(colSums(dt_8==0))
dt_p9=as.vector(colSums(dt_9==0))



p.matrix=matrix(NA,nrow =((8+1)*8/2),ncol = ncol(dt_1))
for(i in 1:ncol(dt_1)){
  p.matrix[1,i]=prop.test(x=c(dt_p1[i],dt_p2[i]),n=c(nrow(dt_1),nrow(dt_2)))$p.value
  p.matrix[2,i]=prop.test(x=c(dt_p1[i],dt_p3[i]),n=c(nrow(dt_1),nrow(dt_3)))$p.value
  p.matrix[3,i]=prop.test(x=c(dt_p1[i],dt_p4[i]),n=c(nrow(dt_1),nrow(dt_4)))$p.value
  p.matrix[4,i]=prop.test(x=c(dt_p1[i],dt_p5[i]),n=c(nrow(dt_1),nrow(dt_5)))$p.value
  p.matrix[5,i]=prop.test(x=c(dt_p1[i],dt_p6[i]),n=c(nrow(dt_1),nrow(dt_6)))$p.value
  p.matrix[6,i]=prop.test(x=c(dt_p1[i],dt_p7[i]),n=c(nrow(dt_1),nrow(dt_7)))$p.value
  p.matrix[7,i]=prop.test(x=c(dt_p1[i],dt_p8[i]),n=c(nrow(dt_1),nrow(dt_8)))$p.value
  p.matrix[8,i]=prop.test(x=c(dt_p1[i],dt_p9[i]),n=c(nrow(dt_1),nrow(dt_9)))$p.value

  p.matrix[9,i]=prop.test(x=c(dt_p2[i],dt_p3[i]),n=c(nrow(dt_2),nrow(dt_3)))$p.value
  p.matrix[10,i]=prop.test(x=c(dt_p2[i],dt_p4[i]),n=c(nrow(dt_2),nrow(dt_4)))$p.value
  p.matrix[11,i]=prop.test(x=c(dt_p2[i],dt_p5[i]),n=c(nrow(dt_2),nrow(dt_5)))$p.value
  p.matrix[12,i]=prop.test(x=c(dt_p2[i],dt_p6[i]),n=c(nrow(dt_2),nrow(dt_6)))$p.value
  p.matrix[13,i]=prop.test(x=c(dt_p2[i],dt_p7[i]),n=c(nrow(dt_2),nrow(dt_7)))$p.value
  p.matrix[14,i]=prop.test(x=c(dt_p2[i],dt_p8[i]),n=c(nrow(dt_2),nrow(dt_8)))$p.value
  p.matrix[15,i]=prop.test(x=c(dt_p2[i],dt_p9[i]),n=c(nrow(dt_2),nrow(dt_9)))$p.value

  p.matrix[16,i]=prop.test(x=c(dt_p3[i],dt_p4[i]),n=c(nrow(dt_3),nrow(dt_4)))$p.value
  p.matrix[17,i]=prop.test(x=c(dt_p3[i],dt_p5[i]),n=c(nrow(dt_3),nrow(dt_5)))$p.value
  p.matrix[18,i]=prop.test(x=c(dt_p3[i],dt_p6[i]),n=c(nrow(dt_3),nrow(dt_6)))$p.value
  p.matrix[19,i]=prop.test(x=c(dt_p3[i],dt_p7[i]),n=c(nrow(dt_3),nrow(dt_7)))$p.value
  p.matrix[20,i]=prop.test(x=c(dt_p3[i],dt_p8[i]),n=c(nrow(dt_3),nrow(dt_8)))$p.value
  p.matrix[21,i]=prop.test(x=c(dt_p3[i],dt_p9[i]),n=c(nrow(dt_3),nrow(dt_9)))$p.value
 
  p.matrix[22,i]=prop.test(x=c(dt_p4[i],dt_p5[i]),n=c(nrow(dt_4),nrow(dt_5)))$p.value
  p.matrix[23,i]=prop.test(x=c(dt_p4[i],dt_p6[i]),n=c(nrow(dt_4),nrow(dt_6)))$p.value
  p.matrix[24,i]=prop.test(x=c(dt_p4[i],dt_p7[i]),n=c(nrow(dt_4),nrow(dt_7)))$p.value
  p.matrix[25,i]=prop.test(x=c(dt_p4[i],dt_p8[i]),n=c(nrow(dt_4),nrow(dt_8)))$p.value
  p.matrix[26,i]=prop.test(x=c(dt_p4[i],dt_p9[i]),n=c(nrow(dt_4),nrow(dt_9)))$p.value

  p.matrix[27,i]=prop.test(x=c(dt_p5[i],dt_p6[i]),n=c(nrow(dt_5),nrow(dt_6)))$p.value
  p.matrix[28,i]=prop.test(x=c(dt_p5[i],dt_p7[i]),n=c(nrow(dt_5),nrow(dt_7)))$p.value
  p.matrix[29,i]=prop.test(x=c(dt_p5[i],dt_p8[i]),n=c(nrow(dt_5),nrow(dt_8)))$p.value
  p.matrix[30,i]=prop.test(x=c(dt_p5[i],dt_p9[i]),n=c(nrow(dt_5),nrow(dt_9)))$p.value

  p.matrix[31,i]=prop.test(x=c(dt_p6[i],dt_p7[i]),n=c(nrow(dt_6),nrow(dt_7)))$p.value
  p.matrix[32,i]=prop.test(x=c(dt_p6[i],dt_p8[i]),n=c(nrow(dt_6),nrow(dt_8)))$p.value
  p.matrix[33,i]=prop.test(x=c(dt_p6[i],dt_p9[i]),n=c(nrow(dt_6),nrow(dt_9)))$p.value

  p.matrix[34,i]=prop.test(x=c(dt_p7[i],dt_p8[i]),n=c(nrow(dt_7),nrow(dt_8)))$p.value
  p.matrix[35,i]=prop.test(x=c(dt_p7[i],dt_p9[i]),n=c(nrow(dt_7),nrow(dt_9)))$p.value
 
  p.matrix[36,i]=prop.test(x=c(dt_p8[i],dt_p9[i]),n=c(nrow(dt_8),nrow(dt_9)))$p.value
 

}

p.matrix[is.na(p.matrix)]=1
saveRDS(p.matrix,"p.matrix.rds")



right_v=celltype

for(ii in 1:(length(celltype)-1)){
  left_v=rep(celltype[ii],(length(celltype)-ii))
  right_v=right_v[-1]

  temp=paste0(left_v,",",right_v)
  
  if(ii==1){
    all_v=temp
  }else{
    all_v=c(all_v,temp)
  }
  
}
View(all_v)
rownames(p.matrix)=all_v

celltype
for(ii in 1:nrow(p.matrix)){
  print(table(p.matrix[ii,]<0.05))
}


