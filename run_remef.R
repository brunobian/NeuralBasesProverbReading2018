tStart     = 1
tEnd       = 1
nIter      = 1
outPath    = '/media/brunobian/ExtraDrive1/Co-registro_2018/Data/LMM/results/'
modType    = 'lmm'
fixEf      = 'pos:sntType + pred:sntType + freq'
fixEf      = 'pos:sntType + pred'
ranEf      = '(1|pal) + (1|suj_id)'
perPath    = '/media/brunobian/ExtraDrive1/Co-registro_2018/Data/LMM/permutations/'
perVar     = 'across'


cstPath    = '/home/brunobian/Documents/Repos/CuBaPeTo/example//'
rPath      = '/home/brunobian/Documents/Repos/CuBaPeTo/R_functions/'
inPath     = '/media/brunobian/DATABRUNO/csvLMM/'

# Load libraries and functions from CuBaPeTo
require(lme4)
source(paste0(rPath, '/generateIterMatrix.R'))
source(paste0(rPath, '/lmmElecTime.R'))
source(paste0(rPath, '/lmElecTime.R'))
source(paste0(rPath, '/lmmTimeWindow.R'))
source(paste0(cstPath, '/processData.R'))
library(remef)
# source('/home/brunobian/Dropbox/Labo Juan/gerardo_juan_diego/oraciones120/mfiles/Codigos_LMER/functions/remef.v0.6.10.R')

fileName = paste(inPath, 't1.csv', sep = '')
todo <- read.csv(fileName, comment.char = "", sep=";")
todo$lngth <-todo$length  
todo <- processData(todo)
todo <- todo[,-which(names(todo) == "E1")] 
todo <- todo[,-which(names(todo) == "time")]
todo <- todo[,-which(names(todo) == "pal")]
todo <- todo[,-which(names(todo) == "MaxJump")]

for (iTime in c(1:103)){
  fileName = paste(inPath, 't', iTime, '.csv', sep = '')
  print(paste0('Loading file T', iTime))
  tmp <- read.csv(fileName, comment.char = "", sep=";") 
  tmp$lngth <-tmp$length  
  dataSet <- processData(tmp)
  
  summary(tmp)
  electrode = 'E1'
  variables = 'freq + tipo:palnum + tipo:pred + (1|pal) + (1|suj_id)'
  model = paste0(electrode,' ~ ' , variables)
  thisLmm <- lmer(model , data = dataSet, REML = FALSE) 
  
  # If it don't converge, continue iterating
  # https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
  tt <- getME(thisLmm,"theta")
  ll <- getME(thisLmm,"lower")
  singularityCheck = min(tt[ll==0])
  if (singularityCheck < 10e-8){
    ss <- getME(thisLmm,c("theta","fixef"))
    system.time(thisLmm <- update(thisLmm, start = ss, 
                                  control = lmerControl(optCtrl=list(maxfun=2e4))))
  }

  var <- paste0('t', iTime)
  todo[var] <- remef(thisLmm, fix=c('tipo0:palnum', 'tipo1:palnum'))  
}

write.csv(todo, '/home/brunobian/Dropbox/Labo Juan/gerardo_juan_diego/oraciones120/mfiles/newEEG.csv')
