rm(list=ls(all=T))

# Time to run ~6-7 min

setwd("/Users/maxrubino/Library/CloudStorage/OneDrive-AuburnUniversity/Agent based models/Simple model code") #set working directory
directory = getwd()
#dir.create("output")
outdir    = paste(directory,"/output/",sep="")                    #directory to save model output  
source(paste(directory, "/source/FunctionSource.R", sep = ''))   #source functions and set source directory

#herbivores
K = 500                   #carrying capacity
maxsize    = 10           #maximum inital herbivore size, lets say this is lbs
searchnum  = 20          #Number of search's per period
G = 0.15 # metabolic growth parameter
CTF = 0.2 #Metabolic cost of foraging
C=0.2 #Metabolic cost per period of foraging, related to size

#food resources
lowsuccess.V  = c(0.8, 0.9 )         #how likely an individual is in finding low quality food 
highsuccess.V = c(0.1, 0.2)        #how likely an individual is in finding high quality food

qual.lo.V=c(1.2,1.5) #mean caloric value (lbs of growth) of lo quality food
qual.hi.V=c(2.5,3.5)#mean caloric value (lbs of growth) of hi quality food
sd.lo=0.15 #sd caloric value (lbs of growth) of lo quality food
sd.hi=0.2 #sd caloric value (lbs of growth) of hi quality food

 

runvars = Replicates(lowsuccess.V, highsuccess.V,qual.lo.V,qual.hi.V)

#summarizing and output
nreps   = 100
towrite = c("r",'reps','K','maxsize','lowsuccess','highsuccess','meandenied','meanhi','num.dead',
            'meansize','meangrowth','denial_size_slope','denial_size_R2', 'growth_hi_slope','growth_hi_R2') 
#list of headers for everything that you will include
write.table(t(towrite), paste(outdir,"outputsummary.csv", sep=""), append=F, sep=",", row.names=F, col.names=F)

#run model iterating over parameters in Replicates
Start=Sys.time()
for(r in 1:nrow(runvars)){
  #set parameter values
  lowsuccess = runvars$lowsuccess[r] 
  highsuccess = runvars$highsuccess[r]
  
  qual.lo=runvars$qual.lo #mean caloric value (lbs of growth) of lo quality food
  qual.hi=runvars$qual.hi #mean caloric value (lbs of growth) of high quality food
  
  for(reps in 1:nreps){
    #initialize population of herbivores
    herbs = data.frame(id = seq(1,K,1), init.size = runif(K,0,10)) #data.frame
    #set up column to store food outcome
    herbs$denied.lo=rep(NA, nrow(herbs)) #number of lo quality patches skipped
    herbs$final.size=rep(NA, nrow(herbs)) #final.size
    
    #iterate over individuals in herbs
    for(i in 1:nrow(herbs)){
      Debt=rep(0,searchnum) #Calori Debt lbs
      Size=c(herbs$init.size[i],rep(NA,searchnum-1)) #animals can grow
      lo.denied=NULL #number of denied food patches
      lo.A=NULL #number of accepted food patches
      hi.found=NULL #number of hi patches found
      for(rr in 1:searchnum){
        lo = rbinom(1, 1, prob = lowsuccess) #did an animal find low quality food?
        hi = rbinom(1, 1, prob = highsuccess) #did an animal find high quality food?
        lo.qual=rnorm(1,lo*qual.lo,lo*sd.lo) #Quality of lo food if failed qual is 0 beacause lo would = 0
        hi.qual=rnorm(1,hi*qual.hi,hi*sd.hi) #Quality of hi food
        Cost=Size[rr]*C+CTF #What is the cost of foraging this patch
        Benefit.lo=lo.qual-Cost #benefit of the foraging
        Benefit.hi=qual.hi-Cost
        A_lo = rbinom(1,1,(1-hi)*plogis(Benefit.lo-Benefit.hi*highsuccess)) #here is where the time between foraging events comes into play, larger animals need more food and therfore are more likely to skip a foraging event
        A_hi = hi #high foods always accepted
        foodobtained= A_lo*lo.qual + A_hi*hi.qual #food that the animal got if based on search results and adaptive choice
        Cost2=ifelse(foodobtained==0,Size[rr]*C,Size[rr]*C+CTF) #take out CTF if they dont forage
        Growth=foodobtained-Cost2-Debt[rr] #Growth of an animal
        Debt[rr+1]=ifelse(Growth<0,Growth+Debt[rr]-CTF,0) #take out foraging cost because they didnt forage
        Growth.A=ifelse(Growth>0,Growth,0) #if the animals dont grow, there growth is 0
        Size[rr+1]=Size[rr]+Growth.A #increase in size
        temp.denied=lo-A_lo #denials
        temp.A=A_lo*lo #accepted
        lo.denied=c(lo.denied,temp.denied) 
        lo.A=c(lo.A,temp.A)
        hi.found=c(hi.found,A_hi) #found high quality patches
        
      }
      
      #record food collected 
      herbs$denied.lo[i]=sum(lo.denied) #number denied
      herbs$final.size[i]=Size[searchnum] #final size of animal
      herbs$hi.found[i]=sum(hi.found) #number of hi quality patches found
      
    }
    herbs$growth=herbs$final.size-herbs$init.size #growth of animals

    #determine denial and size relationship
    lm.denial = lm(denied.lo~init.size, data=herbs) 
    denial_size_slope = summary(lm.denial)$coeff[2,1]
    denial_size_R2 = summary(lm.denial)$r.squared
    
    #determine growth and high food founds relationship
    lm.growth = lm(growth~hi.found,data=herbs)
    growth_hi_slope = summary(lm.growth)$coeff[2,1]
    growth_hi_R2 = summary(lm.growth)$r.squared
    #calculated some summaries
    meansize   = mean(herbs$final.size, na.rm=T) #mean final size of animals
    meangrowth   = mean(herbs$growth, na.rm=T) #mean growth of animals
    meandenied   = mean(herbs$denied.lo, na.rm=T) #mean lo denied
    meanhi   = mean(herbs$hi.found, na.rm=T) #mean number of hi patches found
    num.dead=length(herbs$growth[herbs$growth==0]) #if you dont grow you die
    #record data to file
    towrite = c(r,reps,K,maxsize,lowsuccess,highsuccess,meandenied,meanhi,num.dead,
                meansize,meangrowth,denial_size_slope,denial_size_R2, growth_hi_slope,growth_hi_R2)
    
    #write data, append = F
    write.table(t(towrite), paste(outdir,"outputsummary.csv", sep=""), append=T, sep=",", row.names=F, col.names=F)
    
  }
}
time2run=Sys.time()-Start #just a way to measure how long the model takes to run
#create a folder with date and time so it is unique to new folder names
folder = gsub(" ", "", paste("output_", Sys.time(),time2run, ""), fixed = TRUE)
#dir.create(paste(folder))

#copy output into folder to use later
data = read.table("output/outputsummary.csv", header=T, sep=",")
write.table(data, paste(folder, "/outputsummary.csv", sep=""), row.names=F, col.names=T, sep=",")  

data

#end