`RiboSort` <-
function(data, dataformat="standard", dye="B", output="proportions", zerorows=FALSE, repeats=1, mergerepeats="none")
{
# ------------SECTION 1----------------

n=seq(0,1500,by=1)
dataribo=array(n)
dataabund=array(n)


# ------------SECTION 2----------------

# ------------STEP 1----------------

nfiles=length(data)

if (dataformat != "abimultiple")
{
namelengths=nchar(data)-4
trimmednames=strtrim(data,namelengths)
}

origdetections=c()

step3merges=c()
step3percentmerges=c()

step4merges=c()
step4percentmerges=c()
step4OutofInt=c()
step4percentOutofInt=c()

types=identical(dataformat,"standard") || identical(dataformat,"abisingle") || identical(dataformat,"abimultiple") || identical(dataformat,"beckman")
if (types=="FALSE")
{
print("The dataformat argument supplied is invalid. It can take one")
print("of four values: standard, abisingle, abimultiple or beckman.")
print("Any value supplied must be within inverted commas. Remember ")
print("that R is a case-sensitive environment.")
}else
{

types=identical(mergerepeats,"none") || identical(mergerepeats,"presentinall") || identical(mergerepeats,"presentintwo") || identical(mergerepeats,"presentinone")
if (types=="FALSE")
{
print("The mergerepeats argument supplied is invalid. It can take")
print("one of four values: none, presentinall, presentintwo or ")
print("presentinone. Any value supplied must be within inverted commas.")
print("Remember that R is a case-sensitive environment.")
}else
{

types=identical(output,"abundances") || identical(output,"proportions")
if (types=="FALSE")
{
print("The output argument supplied is invalid. It can take")
print("one of two values: abundances or proportions.")
print("Any value supplied must be within inverted commas.")
print("Remember that R is a case-sensitive environment.")
}else
{

for (k in 1:nfiles)
{

#----------------------------------------------------------------------------------------------------------------------------
#---------------STEP 2-------------- Reading original file and extracting the Bacterial profile only.


if (dataformat=="standard")
{
profile=read.table(data[k])
p=as.matrix(profile)
}

if (dataformat=="abisingle")
{
orig=read.table(data[k],skip=1,fill=TRUE)
onlybs=orig[substr(orig[,1],1,1)==dye,]
onlybs2=onlybs[,3:4]
p=as.matrix(onlybs2)
#write.table(p,"profile.xls",col.names=FALSE,row.names=FALSE)
}

if (dataformat=="beckman")
{
orig=read.table(data[k],header=TRUE,skip=1,fill=TRUE,comment.char="",)#Skip first line, fill=True means fill blank spaces
origm=as.matrix(orig)

r=nrow(origm)#r is the number of rows read in from text file
c=ncol(origm)
#c is the number of cols read in from text file
for (i in 1:r)#In each row containing a *, the standard indicator
{#basically delete the first column point and move 
if (identical(origm[i,1],"*"))#other points in this row back a column.
{
for (j in 1:(c-1))
{
origm[i,j]=origm[i,j+1]
}
}
}

colswanted=c(2,5,6)#Datavector indicating the 3 columns required from
origm1=origm[,colswanted]
origm2=origm1[origm1[,1]==dye,]
profile=origm2[,-1]
mode(profile)="numeric"#Profile created containing only ribotype numbers and
p=as.matrix(profile)#their respective abundances. Converted to matrix p.

#write.table(p,"profile.xls",col.names=FALSE,row.names=FALSE)
}

if(dataformat=="abimultiple")
{
orig=read.csv(data[k],skip=1,fill=TRUE)
colorig=ncol(orig)

#----Removing empty columns from the original file 
colsums=c()
for (u in 2:colorig)
{
	una=sum(orig[,u], na.rm=TRUE)
	colsums=c(colsums,una)
}

del=c()
for (v in 1:length(colsums))
{
	if (identical(all.equal(colsums[v],0,tolerance=0), TRUE))
	{
	del=c(del,-(v+1))
	}
}

ld=length(del)
if (identical(all.equal(ld,0,tolerance=0), TRUE))
{
	orig2=orig
}else
{	
	orig2=orig[,del]
}
#----empty columns removed


roworig2=nrow(orig2)
colorig2=ncol(orig2)
obs=roworig2/2
names=colnames(orig2)[2:colorig2]

origa=orig2[1:obs,]
origb=orig2[(obs+1):roworig2,]

for (i in 2:colorig2)
{
cola=as.matrix(origa[,(i)])
colb=as.matrix(origb[,(i)])
profile=matrix(data=c(cola,colb),nrow=obs)

rid=c()
for (w in 1:obs)
{ 
if (is.na(profile[w,1]))
{
rid=c(rid,-w)
}
}
p=profile[rid,]
p=matrix(p,ncol=2)


#----------------------------------------------------------------------------------------------------------------------------
#---------------STEP 3-------------- Merging ribotypes that are the same; those within 0.5 of each other.

r=nrow(p)#r is the number of ribotypes in sample profile#null vector created for use when omitting rows later
origdetections=c(origdetections,r)

diff=c()
mergecount=0

if (identical(all.equal(r,1,tolerance=0), TRUE))
{
diff=c(diff,1)#A vector with element greater than 0.5 to avoid the min(diff) statement below producing a warning.
}else
{
for (i in 1:(r-1))#For all rows excluding the last one, the difference 
{#between that ribotype and the next one is measured.
diff=c(diff,(p[i+1,1]-p[i,1]))#diff1 is the vector of these differences.
}
}

while (min(diff)<0.5)#Merge the associated pair
{
i=which.min(diff)

p[i,1]=((p[i+1,2]*p[i+1,1])+(p[i,2]*p[i,1]))/(p[i+1,2]+p[i,2])
p[i,2]=p[i+1,2]+p[i,2]#Now we proceed to merge the closest pair (i & i+1)
p=p[-(i+1),]
mergecount=mergecount+1

if (identical(all.equal(i,1,tolerance=0), TRUE))#The smallest diff is between the first two detections.
{
diff[i]=p[i+1,1]-p[i,1]
diff=diff[-(i+1)]
}
else
{
if (identical(all.equal(i,length(diff),tolerance=0),TRUE))#The smallest diff is between the last two detections.
{
diff[i-1]=p[i,1]-p[i-1,1]
diff=diff[-i]
}
else#The smallest diff is neither between the first or the 
{#last two detections, but lies somewhere in the middle.
diff[i-1]=p[i,1]-p[i-1,1]
diff[i]=p[i+1,1]-p[i,1]
diff=diff[-(i+1)]
}
}
}
percent=(mergecount*2*100)/r
step3merges=c(step3merges, mergecount)
step3percentmerges=c(step3percentmerges, percent)

newprofile=round(p,digits=2)


#ASIDE: Below is code to create an excel file containing the newprofile.If such is desired,remove # from code and run.
#write.table(newprofile,"revisedprofile.xls",col.names=FALSE,row.names=FALSE)


#----------------------------------------------------------------------------------------------------------------------------
#---------------STEP 4-------------- Ordering and assigning ribotypes to whole numbers.

s=newprofile#Renaming our new profile as "s" for sample.

n=seq(0,1500,by=1)#An array n is created containing 0 to 1500.
a=matrix(0,nrow=length(n),ncol=2)#a is a zero matrix of length n with 2 columns to 
#which we will assign our sample profile in order.
mergecount=0
OutofIntcount=0

for (j in 1:nrow(s))#For each recorded ribotype number in the sample...
{
b=c(0)#b vector used to break 'i' For Loop later.
for (i in 1:length(n)) #For each possible whole number the ribotype could be 
{#assigned to...
if (((n[i]-0.5) <= s[j,1]) & (s[j,1] < (n[i]+0.5)))#Check where the sample ribo lies.
{
b=b+1
if (identical(all.equal(a[i,1],0,tolerance=0), TRUE))#If this position in 'a'matrix is still empty:assign.
{
a[i,]=s[j,]
}
else#The a position is already filled:must compare value 
{#in it with the value wanting 2go there &see who wins.
x=a[i,1]-n[i]#How close original is to n[i]
y=s[j,1]-n[i]#How close new one is to n[i]
if (abs(y)>abs(x))#Original a[i] is closer and so we don't change it. 
{
if (identical(all.equal(i,length(n),tolerance=0), TRUE)) #If i=1501 the length of n,then row i+1 doesn't exist, 
{#so we'll just merge the original a[i] and new s[j].
a[i,1]=((a[i,2]*a[i,1])+(s[j,2]*s[j,1]))/(a[i,2]+s[j,2])
a[i,2]=a[i,2]+s[j,2]
mergecount=mergecount+1
}
else
{
a[i+1,]=s[j,]#assigning above:a[i+1]always empty @this stage;>a[i].
OutofIntcount=OutofIntcount+1
}
}
else#New s[j] is closer: must compare orig a[i] & a[i-1].
{
if (identical(all.equal(a[i-1,1],0,tolerance=0), TRUE))#If this position in 'a'matrix is empty:
{
a[i-1,]=a[i,]#Assign below.
a[i,]=s[j,]#Assign s[j] to a[i] where it wanted to go.
OutofIntcount=OutofIntcount+1
}
else#a[i-1] is filled,we must see if a[i] lies within n[i] 
{#or was it pushed there?
if (((n[i]-0.5) <= a[i,1]) & (a[i,1] < (n[i]+0.5)))#a[i] lies in n[i];wasn't pushed
{
if ((s[j,1]) < (n[i]+0.25))#If s[j] is more than 0.25 away from the next interval,
{#its not close enough 2b pushed up,merge with a[i].
a[i,1]=((a[i,2]*a[i,1])+(s[j,2]*s[j,1]))/(a[i,2]+s[j,2])
a[i,2]=a[i,2]+s[j,2]
mergecount=mergecount+1
}
else#It is within 0.25 of upper interval,thus push up.
{
a[i+1,]=s[j,]
OutofIntcount=OutofIntcount+1
}
}
else#a[i] doesn't lie in n[i]: must have been pushed out
#of n[i-1]. Check if a[i-2] is empty,if so push down  
{#a[i-1] if within .25 of n[i-2].
if ((a[i-2,1]==0)&(a[i-1,1]<=(n[i-1]-0.25)))
{
a[i-2,]=a[i-1,]
a[i-1,]=a[i,]
a[i,]=s[j,]
}
else#either a[i-2]is filled or a[i-1] is not within 0.25 
{#of n[i-2] interval.Merge a[i]+a[i-1],store in a[i-1].
a[i-1,1]=((a[i-1,2]*a[i-1,1])+(a[i,2]*a[i,1]))/(a[i-1,2]+a[i,2])
a[i-1,2]=a[i-1,2]+a[i,2]
a[i,]=s[j,]
mergecount=mergecount+1
OutofIntcount=OutofIntcount-1
}
}
}
}
}
}
if(b!=0) break#This breaks the 'i' For Loop when s[j] has found the 
}#number in n which it should be assigned to. 
}

percentmerges=(mergecount*2*100)/r
step4merges=c(step4merges, mergecount)
step4percentmerges=c(step4percentmerges, percentmerges)

percentOutofInt=(OutofIntcount*100)/r
step4OutofInt=c(step4OutofInt, OutofIntcount)
step4percentOutofInt=c(step4percentOutofInt, percentOutofInt)

a=round(a,digits=2)#Rounds merged values in 'a' to 2 decimel places.
finalsample=matrix(data=c(n,a),nrow=length(n),ncol=3)#Creates matrix with 1st col=n and 'a' in 2nd & 3rd.
#write.table(finalsample,"sorted.xls",row.names=FALSE,col.names=FALSE)#Code available to create a text file of this matrix.


#----------------------------------------------------------------------------------------------------------------------------
#---------------STEP 5-------------- Combining all datasets together into the final output.

dataribo=c(dataribo,a[,1])#Adding first column (ribotypes) to the data array.
dataabund=c(dataabund,a[,2])#Adding first column (abundances) to the data array.

}#end of for abimultiple (i in 2:colorig2) loop, going thru individual samples
}else#END OF if(dataformat=="abimultiple") TRUE statement
{
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#---------------STEP 3-------------- Merging ribotypes that are the same; those within 0.5 of each other.

r=nrow(p)#r is the number of ribotypes in sample profile#null vector created for use when omitting rows later
origdetections=c(origdetections,r)

diff=c()
mergecount=0

if (identical(all.equal(r,1,tolerance=0), TRUE))
{
diff=c(diff,1)#A vector with element greater than 0.5 to avoid the min(diff) statement below producing a warning.
}else
{
for (i in 1:(r-1))#For all rows excluding the last one, the difference 
{#between that ribotype and the next one is measured.
diff=c(diff,(p[i+1,1]-p[i,1]))#diff1 is the vector of these differences.
}
}

while (min(diff)<0.5)#Merge the associated pair
{
i=which.min(diff)

p[i,1]=((p[i+1,2]*p[i+1,1])+(p[i,2]*p[i,1]))/(p[i+1,2]+p[i,2])
p[i,2]=p[i+1,2]+p[i,2]#Now we proceed to merge the closest pair (i & i+1)
p=p[-(i+1),]
mergecount=mergecount+1

if (identical(all.equal(i,1,tolerance=0), TRUE))#The smallest diff is between the first two detections.
{
diff[i]=p[i+1,1]-p[i,1]
diff=diff[-(i+1)]
}
else
{
if (identical(all.equal(i,length(diff),tolerance=0),TRUE))#The smallest diff is between the last two detections.
{
diff[i-1]=p[i,1]-p[i-1,1]
diff=diff[-i]
}
else#The smallest diff is neither between the first or the 
{#last two detections, but lies somewhere in the middle.
diff[i-1]=p[i,1]-p[i-1,1]
diff[i]=p[i+1,1]-p[i,1]
diff=diff[-(i+1)]
}
}
}
percent=(mergecount*2*100)/r
step3merges=c(step3merges, mergecount)
step3percentmerges=c(step3percentmerges, percent)

newprofile=round(p,digits=2)


#ASIDE: Below is code to create an excel file containing the newprofile.If such is desired,remove # from code and run.
#write.table(newprofile,"revisedprofile.xls",col.names=FALSE,row.names=FALSE)


#----------------------------------------------------------------------------------------------------------------------------
#---------------STEP 4-------------- Ordering and assigning ribotypes to whole numbers.

s=newprofile#Renaming our new profile as "s" for sample.

n=seq(0,1500,by=1)#An array n is created containing 0 to 1500.
a=matrix(0,nrow=length(n),ncol=2)#a is a zero matrix of length n with 2 columns to 
#which we will assign our sample profile in order.
mergecount=0
OutofIntcount=0

for (j in 1:nrow(s))#For each recorded ribotype number in the sample...
{
b=c(0)#b vector used to break 'i' For Loop later.
for (i in 1:length(n)) #For each possible whole number the ribotype could be 
{#assigned to...
if (((n[i]-0.5) <= s[j,1]) & (s[j,1] < (n[i]+0.5)))#Check where the sample ribo lies.
{
b=b+1
if (identical(all.equal(a[i,1],0,tolerance=0), TRUE))#If this position in 'a'matrix is still empty:assign.
{
a[i,]=s[j,]
}
else#The a position is already filled:must compare value 
{#in it with the value wanting 2go there &see who wins.
x=a[i,1]-n[i]#How close original is to n[i]
y=s[j,1]-n[i]#How close new one is to n[i]
if (abs(y)>abs(x))#Original a[i] is closer and so we don't change it. 
{
if (identical(all.equal(i,length(n),tolerance=0), TRUE)) #If i=1501 the length of n,then row i+1 doesn't exist, 
{#so we'll just merge the original a[i] and new s[j].
a[i,1]=((a[i,2]*a[i,1])+(s[j,2]*s[j,1]))/(a[i,2]+s[j,2])
a[i,2]=a[i,2]+s[j,2]
mergecount=mergecount+1
}
else
{
a[i+1,]=s[j,]#assigning above:a[i+1]always empty @this stage;>a[i].
OutofIntcount=OutofIntcount+1
}
}
else#New s[j] is closer: must compare orig a[i] & a[i-1].
{
if (identical(all.equal(a[i-1,1],0,tolerance=0), TRUE))#If this position in 'a'matrix is empty:
{
a[i-1,]=a[i,]#Assign below.
a[i,]=s[j,]#Assign s[j] to a[i] where it wanted to go.
OutofIntcount=OutofIntcount+1
}
else#a[i-1] is filled,we must see if a[i] lies within n[i] 
{#or was it pushed there?
if (((n[i]-0.5) <= a[i,1]) & (a[i,1] < (n[i]+0.5)))#a[i] lies in n[i];wasn't pushed
{
if ((s[j,1]) < (n[i]+0.25))#If s[j] is more than 0.25 away from the next interval,
{#its not close enough 2b pushed up,merge with a[i].
a[i,1]=((a[i,2]*a[i,1])+(s[j,2]*s[j,1]))/(a[i,2]+s[j,2])
a[i,2]=a[i,2]+s[j,2]
mergecount=mergecount+1
}
else#It is within 0.25 of upper interval,thus push up.
{
a[i+1,]=s[j,]
OutofIntcount=OutofIntcount+1
}
}
else#a[i] doesn't lie in n[i]: must have been pushed out
#of n[i-1]. Check if a[i-2] is empty,if so push down  
{#a[i-1] if within .25 of n[i-2].
if ((a[i-2,1]==0)&(a[i-1,1]<=(n[i-1]-0.25)))
{
a[i-2,]=a[i-1,]
a[i-1,]=a[i,]
a[i,]=s[j,]
}
else#either a[i-2]is filled or a[i-1] is not within 0.25 
{#of n[i-2] interval.Merge a[i]+a[i-1],store in a[i-1].
a[i-1,1]=((a[i-1,2]*a[i-1,1])+(a[i,2]*a[i,1]))/(a[i-1,2]+a[i,2])
a[i-1,2]=a[i-1,2]+a[i,2]
a[i,]=s[j,]
mergecount=mergecount+1
OutofIntcount=OutofIntcount-1
}
}
}
}
}
}
if(b!=0) break#This breaks the 'i' For Loop when s[j] has found the 
}#number in n which it should be assigned to. 
}

percentmerges=(mergecount*2*100)/r
step4merges=c(step4merges, mergecount)
step4percentmerges=c(step4percentmerges, percentmerges)

percentOutofInt=(OutofIntcount*100)/r
step4OutofInt=c(step4OutofInt, OutofIntcount)
step4percentOutofInt=c(step4percentOutofInt, percentOutofInt)

a=round(a,digits=2)#Rounds merged values in 'a' to 2 decimel places.
finalsample=matrix(data=c(n,a),nrow=length(n),ncol=3)#Creates matrix with 1st col=n and 'a' in 2nd & 3rd.
#write.table(finalsample,"sorted.xls",row.names=FALSE,col.names=FALSE)#Code available to create a text file of this matrix.

#----------------------------------------------------------------------------------------------------------------------------
#---------------STEP 5-------------- Combining all datasets together into the final output.

dataribo=c(dataribo,a[,1])#Adding first column (ribotypes) to the data array.
dataabund=c(dataabund,a[,2])#Adding first column (abundances) to the data array.

}#END OF if(dataformat=="abimultiple") ELSE statement
}#end of for (k in 1:nfiles) loop

if (dataformat=="abimultiple")
{
names=c()
for (k in 1:nfiles)
{
orig=read.csv(data[k],skip=1,fill=TRUE)
colorig=ncol(orig)

#----Removing empty columns from the original file 
colsums=c()
for (u in 2:colorig)
{
	una=sum(orig[,u], na.rm=TRUE)
	colsums=c(colsums,una)
}

del=c()
for (v in 1:length(colsums))
{
	if (identical(all.equal(colsums[v],0,tolerance=0), TRUE))
	{
	del=c(del,-(v+1))
	}
}

ld=length(del)
if (identical(all.equal(ld,0,tolerance=0), TRUE))
{
	orig2=orig
}else
{	
	orig2=orig[,del]
}
#----empty columns removed

roworig2=nrow(orig2)
colorig2=ncol(orig2)
names=c(names,colnames(orig2)[2:colorig2])
}
nfiles=length(names)
}

#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------



#-------------SECTION 3----------------

outputribo=matrix(dataribo,nrow=1501)
outputabund=matrix(dataabund,nrow=1501)

#---------Labelling the Sample Columns and Rows---------

if(dataformat=="abimultiple")
{
colnames(outputribo)=c("Samples:",names)
rownames(outputribo)=outputribo[,1]

colnames(outputabund)=c("Samples:",names)
rownames(outputabund)=outputabund[,1]
}else
{
colnames(outputribo)=c("Samples:",trimmednames)
rownames(outputribo)=outputribo[,1]

colnames(outputabund)=c("Samples:",trimmednames)
rownames(outputabund)=outputabund[,1]
}

#----------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------

#----Making adjustments required by ARGUMENTS 4 & 5: repeats, mergerepeats


#-----------------------------NONE

if (mergerepeats=="none")
{
finalabund=outputabund
finalribo=outputribo
}

#-------------ASIDE-just dealing with the issue when repeats!=1 and it doesn't divide samples evenly.

q=nrow(outputabund)
g=nfiles/repeats#number of groups of replicates
roundg=round(g)
gf=c((1:roundg)*repeats)#list of group-finishing column indices

if ((g-roundg)!=0 && mergerepeats!="none")
{
print("The number of repeats specified does not divide evenly")
print("into the total number of profiles submitted. This suggests")
print("that the number of repeats taken was not constant for")
print("all samples, a condition required by the RiboSort")
print("function.")
print("Suggestion: If your number of repeats is not constant")
print("across samples, then submit repeats=1 to obtain a sorted")
print("output file containing all profiles, and proceed to")
print("manually merge the repeats.")
}else
{


#-----------------------------PRESENT-IN-ALL

if (mergerepeats=="presentinall")
{
if (repeats!=1)
{

namedoutputribo=outputribo[,-1]#Removing first column containing ribotype numbers
namedoutputabund=outputabund[,-1]#Removing first column containing ribotype numbers

outabundmerged=c()
for (i in gf)#creating a loop that identifies each group of repeats
{#and produces a composite sample for each group
group=namedoutputabund[,(i-repeats+1):i]
abundtot=c()
for(c in 1:repeats)
{
abundtot=c(abundtot,sum(group[,c]))
}
avtotal=sum(abundtot)/repeats
for (j in 1:q)#looking in each row
{
notempty=0
for (k in 1:repeats)#looking in each column of each row
{
if (group[j,k]!=0)
{
notempty=notempty+1
}
}
if (notempty==repeats)
{
abundprop=c()
for(a in 1:repeats)
{
abundprop=c(abundprop, group[j,a]/abundtot[a])
}
avabundprop=sum(abundprop)/repeats
group[j,1]=avtotal*avabundprop
u=group[j,1]
u=round(u)
group[j,1]=u
}
else
{
group[j,1]=0
}
}
outabundmerged=matrix(data=c(outabundmerged,group[,1]),nrow=q)
}

outribomerged=c()
for (i in gf)#creating a loop that identifies each group of repeats
{#and produces a composite sample for each group 
group=namedoutputribo[,(i-repeats+1):i]
for (j in 1:q)#looking in each row
{
notempty=0
for (k in 1:repeats)#looking in each column of each row
{
if (group[j,k]!=0)
{
notempty=notempty+1
}
}
if (notempty==repeats)
{
group[j,1]=sum(group[j,])/repeats
u=group[j,1]
u=round(u,digits=2)
group[j,1]=u
}
else
{
group[j,1]=0
}
}
outribomerged=matrix(data=c(outribomerged,group[,1]),nrow=q)
}

#-------Naming the rows & columns in the two new merged outputs

outabundmerged=matrix(data=c(outputabund[,1],outabundmerged),nrow=q)#Inserting rownames as the first column of the matrix
rownames(outabundmerged)=rownames(namedoutputabund)#Formally labelling rows-won't print them in output:they are 1st col.

outribomerged=matrix(data=c(outputribo[,1],outribomerged),nrow=q)#Same for ribo matrices
rownames(outribomerged)=rownames(namedoutputribo)

d=colnames(namedoutputabund)
cmerge=c("Samples:")
for(i in gf)
{
cmerge=c(cmerge,d[(i-repeats+1)])
}

colnames(outabundmerged)=cmerge
colnames(outribomerged)=cmerge

finalabund=outabundmerged
finalribo=outribomerged

}else
{#end of if(reps!=1) statement, thus following lines executed
finalabund=outputabund#when reps=1 and mergerepeats="presentinall"
finalribo=outputribo
}
}#end of if(mergerepeats="presentinall") statement


#---------------------------PRESENT-IN-TWO

if (mergerepeats=="presentintwo")
{
if(repeats!=1)
{

namedoutputribo=outputribo[,-1]#Removing first column containing ribotype numbers
namedoutputabund=outputabund[,-1]#Removing first column containing ribotype numbers

outabundmerged=c()
for (i in gf)#creating a loop that identifies each group of repeats
{#and produces a composite sample for each group
group=namedoutputabund[,(i-repeats+1):i]
abundtot=c()
for(c in 1:repeats)
{
abundtot=c(abundtot,sum(group[,c]))
}
avtotal=sum(abundtot)/repeats
for (j in 1:q)#looking in each row
{
notempty=0
for (k in 1:repeats)#looking in each column of each row
{
if (group[j,k]!=0)
{
notempty=notempty+1
}
}
if (notempty>=2)
{
abundprop=c()
for(a in 1:repeats)
{
abundprop=c(abundprop, group[j,a]/abundtot[a])
}
avabundprop=sum(abundprop)/repeats
group[j,1]=avtotal*avabundprop
u=group[j,1]
u=round(u)
group[j,1]=u
}
else
{
group[j,1]=0
}
}
outabundmerged=matrix(data=c(outabundmerged,group[,1]),nrow=q)
}

outribomerged=c()
for (i in gf)#creating a loop that identifies each group of repeats
{#and produces a composite sample for each group
group=namedoutputribo[,(i-repeats+1):i]
for (j in 1:q)#looking in each row
{
notempty=0
for (k in 1:repeats)#looking in each column of each row
{
if (group[j,k]!=0)
{
notempty=notempty+1
}
}
if (notempty>=2)
{
group[j,1]=sum(group[j,])/repeats
u=group[j,1]
u=round(u,digits=2)
group[j,1]=u
}
else
{
group[j,1]=0
}
}
outribomerged=matrix(data=c(outribomerged,group[,1]),nrow=q)
}

#-------Naming the rows & columns in the two new merged outputs

outabundmerged=matrix(data=c(outputabund[,1],outabundmerged),nrow=q)#Inserting rownames as the first column of the matrix
rownames(outabundmerged)=rownames(namedoutputabund)#Formally labelling rows-won't print them in output:they are 1st col.

outribomerged=matrix(data=c(outputribo[,1],outribomerged),nrow=q)
rownames(outribomerged)=rownames(namedoutputribo)#Same for ribo matrices

d=colnames(namedoutputabund)
cmerge=c("Samples:")
for(i in gf)
{
cmerge=c(cmerge,d[(i-repeats+1)])
}

colnames(outabundmerged)=cmerge
colnames(outribomerged)=cmerge

finalabund=outabundmerged
finalribo=outribomerged

}else
{#end of if(reps!=1) statement, thus following lines executed
finalabund=outputabund#when reps=1 and mergerepeats="presentintwo"
finalribo=outputribo
}
}#end of if(mergerepeats="presentintwo") statement



#---------------------------PRESENT-IN-ONE


if (mergerepeats=="presentinone")
{
if(repeats!=1)
{

namedoutputribo=outputribo[,-1]#Removing first column containing ribotype numbers
namedoutputabund=outputabund[,-1]#Removing first column containing ribotype numbers

outabundmerged=c()
for (i in gf)#creating a loop that identifies each group of repeats
{#and produces a composite sample for each group
group=namedoutputabund[,(i-repeats+1):i]
abundtot=c()
for(c in 1:repeats)
{
abundtot=c(abundtot,sum(group[,c]))
}
avtotal=sum(abundtot)/repeats
for (j in 1:q)#looking in each row
{
notempty=0
for (k in 1:repeats)#looking in each column of each row
{
if (group[j,k]!=0)
{
notempty=notempty+1
}
}
if (notempty>=1)
{
abundprop=c()
for(a in 1:repeats)
{
abundprop=c(abundprop, group[j,a]/abundtot[a])
}
avabundprop=sum(abundprop)/repeats
group[j,1]=avtotal*avabundprop
u=group[j,1]
u=round(u)
group[j,1]=u
}
else
{
group[j,1]=0
}
}
outabundmerged=matrix(data=c(outabundmerged,group[,1]),nrow=q)
}

outribomerged=c()
for (i in gf)#creating a loop that identifies each group of repeats
{#and produces a composite sample for each group
group=namedoutputribo[,(i-repeats+1):i]
for (j in 1:q)#looking in each row
{
notempty=0
for (k in 1:repeats)#looking in each column of each row
{
if (group[j,k]!=0)
{
notempty=notempty+1
}
}
if (notempty>=1)
{
group[j,1]=sum(group[j,])/repeats
u=group[j,1]
u=round(u,digits=2)
group[j,1]=u
}
else
{
group[j,1]=0
}
}
outribomerged=matrix(data=c(outribomerged,group[,1]),nrow=q)
}

#-------Naming the rows & columns in the two new merged outputs

outabundmerged=matrix(data=c(outputabund[,1],outabundmerged),nrow=q)#Inserting rownames as the first column of the matrix
rownames(outabundmerged)=rownames(namedoutputabund)#Formally labelling rows-won't print them in output:they are 1st col.

outribomerged=matrix(data=c(outputribo[,1],outribomerged),nrow=q)
rownames(outribomerged)=rownames(namedoutputribo)#Same for ribo matrices

d=colnames(namedoutputabund)
cmerge=c("Samples:")
for(i in gf)
{
cmerge=c(cmerge,d[(i-repeats+1)])
}

colnames(outabundmerged)=cmerge
colnames(outribomerged)=cmerge

finalabund=outabundmerged
finalribo=outribomerged

}else
{#end of if(reps!=1) statement, thus following lines executed
finalabund=outputabund#when reps=1 and mergerepeats="presentinone"
finalribo=outputribo
}
}#end of if(mergerepeats="presentinone") statement


#----------------------------------------------------------------------------------------------------------------------------


#-----REMOVING ZEROROWS----from top and bottom from all files, and all zerorows if requested to do so.


#Removing zero rows from top of list.
delete=c()
for (i in 1:(nrow(finalribo)))
{

if (identical(all.equal((sum(finalribo[i,])-finalribo[i,1]),0,tolerance=0), TRUE)) 
{
delete=c(delete, -i)
}
else
{
if (identical(all.equal(length(delete),0,tolerance=0), TRUE)) 
{break}
finalribo=finalribo[delete,]
finalabund=finalabund[delete,]
break 
}
} 

#Removing zero rows from bottom of list.
delete2=c()
for (i in (nrow(finalribo)):1)
{
if (identical(all.equal((sum(finalribo[i,])-finalribo[i,1]),0,tolerance=0), TRUE)) 
{
delete2=c(delete2, -i)
}
else
{
if (identical(all.equal(length(delete2),0,tolerance=0), TRUE)) 
{break}
finalribo=finalribo[delete2,]
finalabund=finalabund[delete2,]
break 
}
} 

#----------------------------------------------------------------------------------------------------------------------------
#------Making adjustments requested by ARGUMENT 3:zerorows.

if (zerorows==FALSE)
{
#Remove all remaining zerorows.
delete3=c()
for (i in 1:(nrow(finalribo)))
{
if (identical(all.equal((sum(finalribo[i,])-finalribo[i,1]),0,tolerance=0), TRUE)) 
{
delete3=c(delete3, -i)
}
}
finalribo=finalribo[delete3,]
finalabund=finalabund[delete3,]
} 

#----------------------------------------------------------------------------------------------------------------------------
#---Creating the proportions files

col1=as.matrix(finalabund[,1])
onlysamps=as.matrix(finalabund[,-1])
coltots=apply(onlysamps,2,sum)
props=sweep(onlysamps,2,coltots,"/")

finalprop=matrix(data=c(col1,props),nrow=nrow(finalabund))
rownames(finalprop)=rownames(finalabund)
colnames(finalprop)=colnames(finalabund)


#----------------------------------------------------------------------------------------------------------------------------
#---------Creating the output datafiles in Excel

write.table(finalribo,"Ribotypes Output.xls",col.names=TRUE,row.names=FALSE)
write.table(finalabund,"Abundances Output.xls",col.names=TRUE,row.names=FALSE)
write.table(finalprop,"Proportions Output.xls",col.names=TRUE,row.names=FALSE)

#----------------------------------------------------------------------------------------------------------------------------
#---------Creating Files to inform users of the extent of changes made to samples

step3percentmerges=round(step3percentmerges, digits=0)
step4percentmerges=round(step4percentmerges, digits=0)
step4percentOutofInt=round(step4percentOutofInt, digits=0)
totalpercentchanged=step3percentmerges+step4percentmerges+step4percentOutofInt

blank=rep(" ",times=(nfiles+1))

if(dataformat=="abimultiple")
{
inform=c("Samples:",names,blank,"Number of original sequencer detections:",origdetections,"Number of mergings during Pre-Assignment Stage:",step3merges,"Number of mergings during Assignment Stage:",step4merges,"Number of detections assigned to intervals (n+1) or (n-1) when they lie in interval n:",step4OutofInt,blank,"Percentage detections involved in mergings during the Pre-Assignment Stage:",step3percentmerges,"Percentage detections involved in mergings during the Assignment Stage:",step4percentmerges,"Percentage detections assigned to interval (n-1) or (n+1) when they lie in interval n:",step4percentOutofInt,"Total percentage of data adjusted by RiboSort:",totalpercentchanged)
}else
{
inform=c("Samples:",trimmednames,blank,"Number of original sequencer detections:",origdetections,"Number of mergings during Pre-Assignment Stage:",step3merges,"Number of mergings during Assignment Stage:",step4merges,"Number of detections assigned to intervals (n+1) or (n-1) when they lie in interval n:",step4OutofInt,blank,"Percentage detections involved in mergings during the Pre-Assignment Stage:",step3percentmerges,"Percentage detections involved in mergings during the Assignment Stage:",step4percentmerges,"Percentage detections assigned to interval (n-1) or (n+1) when they lie in interval n:",step4percentOutofInt,"Total percentage of data adjusted by RiboSort:",totalpercentchanged)
}
informmatrix=matrix(inform,ncol=(nfiles+1),byrow=TRUE)

write.table(informmatrix,"Information File.xls",col.names=FALSE,row.names=FALSE)

#----------------------------------------------------------------------------------------------------------------------------
finalabund2=finalabund[,-1]#Removing 1st col containing ribo numbers for output matrix

if(class(finalabund2)!="matrix")
{
f=finalabund
rownames(f)=c()
f2=f[,-1]
f2m=matrix(f2)
rownames(f2m)=rownames(finalabund)
colnames(f2m)=colnames(finalabund)[2]
finalabund2=f2m
}

finalprop2=finalprop[,-1]#Removing 1st col containing ribo numbers for output matrix

if(class(finalprop2)!="matrix")
{
f=finalprop
rownames(f)=c()
f2=f[,-1]
f2m=matrix(f2)
rownames(f2m)=rownames(finalprop)
colnames(f2m)=colnames(finalprop)[2]
finalprop2=f2m
}

if (output=="abundances")
{
return(finalabund2)
}else
{
return(finalprop2)
}


}#end of ifelse stating whether or not nfiles/repeats divides evenly
}#end of ifelse stating whether or not output argument has been supplied correctly
}#end of ifelse stating whether or not mergerepeats argument has been supplied correctly
}#end of ifelse stating whether or not dataformat argument has been supplied correctly

}#end of the RiboSort Function body

