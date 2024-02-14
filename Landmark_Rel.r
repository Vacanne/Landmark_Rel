### In my defense, this was one of my Internship projects and i didn't name those variables
sampleDf1 <- data.frame(x=c(160.604990592718, 183.899441767906, 208.419916689156, 372.707098661533, 533.316209395722, 550.480541840598, 566.41885053941, 201.063774212781, 236.618462848594, 329.796267549345, 413.165882281596, 514.925853204785, 543.124399364222, 376.38516989972, 376.38516989972, 324.892172565095, 370.255051169408, 424.200095996158, 231.714367864344, 278.303270214719, 324.892172565095, 280.755317706844, 413.165882281596, 472.015022092597, 513.699829458722, 472.015022092597),
y=c(430.240182431859, 293.80792938545, 208.537771231444, 87.8477012288511, 219.032559927321, 304.302718081327, 432.863879605828, 427.616485257889, 453.853456997583, 428.928333844874, 422.369090909951, 456.477154171553, 436.799425366782, 402.69136210518, 275.442049167664, 272.818351993694, 251.828774601939, 266.259109058771, 388.261027648348, 413.186150801058, 386.949179061364, 377.766238952471, 381.701784713425, 411.874302214073, 396.132119170256, 375.142541778501))
sampleDf2 <- data.frame(x=c(150.892982933491, 181.822442306398, 209.658955742015, 375.647054376618, 526.170423324768, 548.8520268649, 569.471666446839, 198.318153971949, 237.495469177632, 330.283847296354, 408.63847770772, 513.798639575605, 541.635153011222, 375.647054376618, 375.647054376618, 319.974027505385, 377.709018334812, 419.979279477786, 234.402523240341, 280.796712299702, 319.974027505385, 279.765730320605, 426.165171352367, 467.404450516244, 508.643729680121, 475.652306349019),
y=c(432.177345760942, 299.073725809619, 217.247729937904, 101.600322439214, 216.156716659614, 300.164739087908, 436.5413988741, 426.722279369494, 452.906598048443, 428.904305926073, 422.358226256336, 451.815584770153, 435.45038559581, 409.266066916862, 276.162446965539, 272.88940713067, 260.888261069485, 269.616367295802, 391.809854464229, 413.63012003002, 390.71884118594, 374.353642011597, 388.53681462936, 411.448093473441, 391.809854464229, 378.717695124755))

reliability<-function(db,n,k,A,B){

#db <- as.integer(dlgInput("NUMBER OF DIMENSION ", Sys.info()[""])$res)
#db <- as.integer(readline(prompt="NUMBER OF DIMENSION "))
#n <- as.integer(dlgInput("NUMBER OF SUBJECTS ", Sys.info()[""])$res)
#n <- as.integer(readline(prompt="NUMBER OF SUBJECTS"))
#k <- as.integer(dlgInput("NUMBER OF LANDMARKS ", Sys.info()[""])$res)
#k <- readline(prompt="NUMBER OF LANDMARKS")
#baslangic={'','',''};
a=0; rr=0; ka=0; kb=0; sn1=1; tlsr=0; tlr2=0; tr2=0; tr1a=0; tr1b=0;x=0
tl2=0; ts2=0; tlsb2=0; td=0; tlre=0; tsre=0; sn2=0; d=0; trs2=0; srkt=0;i=0;
C = choose(k,2); sn2=(n*(C)); s=C*n;
R<-matrix(nrow=2*n*C,ncol=4)

P<-matrix(nrow=n*C,ncol=4)
V<-matrix(nrow=n*C,ncol=4)
print(C)
for (a in 1:2){
  for (b in 1:n){
    for (c in 1:C){
          d=d+1; R[d,1]=a; R[d,2]=b; R[d,3]=c;
    }
  }
}



# Calculation of the Euclidean distances between landmark pairs for raters
if (db < 3){
# calculations for two dimensions
  for(p in seq(from = 1, to = n*k, by = k)){
    t=k+p; m=t-k;
    for (j in m:(t-1)){
      for (i in m:(t-1)){
        if ((j<i) & (i!=j)){
          acia1= sqrt(((A[j,1]-A[i,1])^2)+((A[j,2]-A[i,2])^2));
          R[sn1,4]= acia1; sn2=sn2+1;
          acib1=sqrt(((B[j,1]-B[i,1])^2)+((B[j,2]-B[i,2])^2));
          R[sn2,4]= acib1;
#calculation of sum of squares for lxrxs
          tlsr= tlsr + (acia1^2) + (acib1^2);
#calculation of sum of squares for r (cont. Outside the loop)
          tr1a= tr1a + acia1; tr1b = tr1b + acib1;
#calculation of adjustment factor (cont. Outside the loop)
          td= td + acia1 + acib1;
          sn1= sn1+1;
        }
      }
    }
  }
}
if (db>2){
#calculations for three dimensions
  for(p in seq(from = 1, to = n*k, by = k)){
    t=k+p; m=t-k;
    for (j in m:(t-1)){
      for (i in m:(t-1)){
        if ((j<i) & (i!=j)){
          acia1=sqrt(((A[j,1] - A[i,1])^2) + ((A[j,2] - A[i,2])^2) + ((A[j,3] - A[i,3])^2));
          R[sn1,4]=as.integer(acia1); sn2=sn2+1;
          acib1 = sqrt(((B[j,1]- B[i,1])^2) + ((B[j,2]-B[i,2])^2)+((B[j,3]-B[i,3])^2));
          R[sn2,4] =as.integer(acib1);
#calculation of sum of squares for lxrxs
          tlsr= tlsr + (acia1^2) + (acib1^2);
#calculation of sum of squares for r (cont. Outside the loop)
          tr1a= tr1a+ acia1; tr1b=tr1b + acib1;
#calculation of adjustment factor (cont. Outside the loop)
          td= td + acia1 + acib1;
          sn1= sn1+1;
        }
      }
    }
  }
}

#calculation of adjustment factor
df=(td^2)/(2*C*n)
#calculation of sum of square of r and mean of square of r
tr2=(tr1a^2)+(tr1b^2); rkt=(tr2/(C*n))-df; rko=((tr2/(C*n))-df)/1; ms_r=rko;
#Division of R matrix into 2 parts for analysis (P for rater 1, V for rater 2)
for (i in 1:s){
  for (j in 1:4){
    P[i,j]= R[i,j];
  }
}
v=0
xyz=s+1
for (i in xyz:(2*s)){
  v=v+1;
  for (j in 1:4){
    V[v,j]=R[i,j];
  }
}
#calculation of sum of squares for l
for (l in 1:C){
  tl1=0;
  for (i in 1:s){
    if (l>P[i,3]-1 & (l <P[i,3]+1)){
      tl1=tl1+P[i,4];}
    if (l>(V[i,3]-1) & (l<V[i,3]+1)){
      tl1=tl1+V[i,4];}
}
tl2=tl2+(tl1^2);
}
lkt=(tl2/(n*2))-df
#calculation of sum of squares for s
for (sb in 1:n){
  ts1=0;
  for (i in 1:s){
    if ((sb>(P[i,2]-1) & sb<(P[i,2]+1))){
      ts1=ts1+P[i,4];}
    if ((sb>(V[i,2]-1) & sb<(V[i,2]+1))){
      ts1=ts1+V[i,4];}
  }
  ts2=ts2+(ts1^2);
}

skt=(ts2/(C*2))-df;

#calculation of sum of squares for lxs
for (l in 1:C){
  for (sb in 1:n){
    tls1=0;
    for (i in 1:s){
      if ((l>(P[i,3]-1)) & (l<(P[i,3]+1))){
        if ((sb>(P[i,2]-1)) & (sb<(P[i,2]+1))){
          tls1=tls1+P[i,4];}
      }

      if ((l>(V[i,3]-1) & l<(V[i,3]+1))){
        if ((sb>(V[i,2]-1) & sb<(V[i,2]+1))){
          tls1=tls1+V[i,4];}
      }
    }
    tlsb2=tlsb2+(tls1^2);
  }
}
lskt=(tlsb2/2)-df-lkt-skt;

# calculation of sum of squares for lxr
for (l in 1:C){
  tlrp=0;
  tlrv=0;
  for (i in 1:s){
    if ((l>(P[i,3]-1) & l<(P[i,3]+1))){
      tlrp=tlrp+P[i,4]; tlrv=tlrv+V[i,4];}
    }
  tlre=tlre+(tlrp^2+tlrv^2);
}
lrkt=(tlre/n)-df-lkt-rkt;
# calculation of sum of squares for rxs
for (r in 1:2){
  for (sb in 1:n){
    trsp=0; trsv=0;
    for (i in 1:s){
      if ((r>(P[i,1]-1) & r<(P[i,1]+1))){
        if ((sb>(P[i,2]-1) & sb<(P[i,2]+1))){
          trsp=trsp+P[i,4];}
      }
      if ((r>(V[i,1]-1) & r<(V[i,1]+1))){
        if ((sb>(V[i,2]-1) & sb<(V[i,2]+1))){
          trsv=trsv+V[i,4];}
      }
    }
    trs2=trs2+(trsp^2)+(trsv^2);
  }
}
srkt=(trs2/C)-df-rkt-skt;

#calculation of sum of squares for lxrxs
tlsrkt=tlsr-df-lkt-rkt-skt-lrkt-lskt-srkt
#calculation of Mean Squares
ms_lrs=tlsrkt/((n-1)*(C-1));
ms_lr=lrkt/(C-1);
ms_ls=lskt/((n-1)*(C-1));
ms_rs=srkt/(n-1);
ms_l=lkt/(C-1);
ms_s=skt/(n-1);
#estimates of the variance components
var_lrs=ms_lrs;
var_lr=(ms_lr-ms_lrs)/n;
var_ls=(ms_ls-ms_lrs)/2;
var_rs=(ms_rs-ms_lrs)/C;
var_l=(ms_l-ms_lr-ms_ls+ms_lrs)/(2*n);
var_r=(ms_r-ms_lr-ms_rs+ms_lrs)/(C*n);
var_s=(ms_s-ms_ls-ms_rs+ms_lrs)/(2*C);
if (var_lrs<0){
  var_lrs=0;}

if (var_lr<0){
  var_lr=0;}

if (var_ls<0){
  var_ls=0;}

if (var_rs<0){
  var_rs=0;}

if (var_l<0){
  var_l=0;}

if (var_r<0){
  var_r=0;}

if (var_s<0){
var_s=0;}

var_rel=(var_lr/2)+(var_ls/n)+(var_lrs/(2*n));

#calculation of G coefficient

G=var_l/(var_l+var_rel);
cat('Rater-->r, Subject-->s, Landmark eslestirmesi-->l' )
cat(' \n r s l Euclidean distances\n ');print(R)
cat('\n r sum of square',rkt)
cat('\n l sum of square', lkt)
cat('\n s sum of square', skt)
cat('\n lxs sum of square', lskt)
cat('\n lxr sum of square', lrkt)
cat('\n rxs sum of square', srkt)
cat('\n lxsxr sum of square', tlsrkt)
cat('\n r mean square', ms_r)
cat('\n l mean square', ms_l)
cat('\n s mean square', ms_s)
cat('\n lxs mean square', ms_ls)
cat('\n lxr mean square', ms_lr)
cat('\n rxs mean square', ms_rs)
cat('\n lxrxs mean square', ms_lrs)
cat('\n variance of r', var_r)
cat('\n variance of l', var_l)
cat('\n variance of s', var_s)
cat('\n variance of lxs',var_ls)
cat('\n variance of lxr', var_lr)
cat('\n variance of rxs', var_rs)
cat('\n variance of lxrxs', var_lrs)
cat('\n variance of rel', var_rel)
cat(G);cat('\n G COEFFICIENT')
}
reliability(2,2,26,sampleDf1,sampleDf2)