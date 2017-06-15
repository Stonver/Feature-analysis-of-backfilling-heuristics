#!/usr/bin/env Rscript

library(docopt)
library(ggplot2)



'usage: tool.R <input> ... [--out=<output>] 
#[-b <o2>]
        tool.R -h | --help

options:
 <input>        The input data
 --out=<output>  Output file [Default: default.features] 
 #-b <o2>                Output file in case of pdf/tikz/png output.
#			[Default: tmp.pdf]

#--output=FILE    Output file [default: test.txt]
 ' -> doc

args<- docopt(doc)


#reading swf file
swf_read <- function(f)
{
    	df <- read.table(f,comment.char=';')
    	names(df) <- c('job_id','submit_time','wait_time','run_time','proc_alloc','cpu_time_used','mem_used','proc_req','time_req','mem_req','status','user_id','group_id','exec_id','queue_id','partition_id','previous_job_id','think_time')
	return(df)
}


#reading multiple files
#temp = list.files(path = args$input, pattern="*week.swf")
#temp = paste ("o/zymakefile", temp, sep = "/")
#myfiles = lapply(temp, swf_read)

#reading one file
data <- swf_read(args$input)

	
RunFeatures <- c(mean(data$run_time),max(data$run_time),min(data$run_time), sum(data$run_time), quantile(data$run_time, c(0.1,0.5,0.9), names = F))
#max(data$run_time)-min(data$run_time)

TimeReqFeatures <- c(mean(data$time_req),max(data$time_req),min(data$time_req), sum(data$time_req), quantile(data$time_req, c(0.1,0.5,0.9), names = F))	

JobSizeFeatures <- c(mean(data$proc_alloc),max(data$proc_alloc),min(data$proc_alloc), sum(data$proc_alloc), quantile(data$proc_alloc, c(0.1,0.5,0.9), names = F))

JobArea <- data$run_time * data$proc_alloc
JobAreaFeatures <- c(mean(JobArea),max(JobArea),min(JobArea), sum(JobArea), quantile(JobArea, c(0.1,0.5,0.9), names = F))

RatioRunReq <- data$run_time / data$time_req
RatioRunReqFeatures <- c(mean(RatioRunReq),max(RatioRunReq),min(RatioRunReq), sum(RatioRunReq), quantile(RatioRunReq, c(0.1,0.5,0.9), names = F))

RatioRunJob <- data$run_time / data$proc_alloc
RatioRunJobFeatures <- c(mean(RatioRunJob),max(RatioRunJob),min(RatioRunJob), sum(RatioRunJob), quantile(RatioRunJob, c(0.1,0.5,0.9), names = F))

secperday <- 60*60*24
CosinusTime <- cos(2*pi*(data$submit_time %% secperday)/secperday)
SinusTime <- sin(2*pi*(data$submit_time %% secperday)/secperday)
CosinusTimeFeatures <- c(mean(CosinusTime),max(CosinusTime),min(CosinusTime), sum(CosinusTime), quantile(CosinusTime, c(0.1,0.5,0.9), names = F))
SinusTimeFeatures <- c(mean(SinusTime),max(SinusTime),min(SinusTime), sum(SinusTime), quantile(SinusTime, c(0.1,0.5,0.9), names = F))


week <- cbind(t(RunFeatures),t(TimeReqFeatures), 
t(JobSizeFeatures), t(JobAreaFeatures), 
t(RatioRunReqFeatures), t(RatioRunJobFeatures),
t(CosinusTimeFeatures),t(SinusTimeFeatures))

#rownames(week)<-1:length(week[,1])
#week <- as.data.frame(week)
#colnames(week) <- c('AvgRun','MaxRun','MinRun','SumRun', '10Run', 'MedRun', '90Run', 'AvgTReq','MaxTReq','MinTReq','SumTReq', '10TReq', 'MedTReq', '90TReq', 'AvgSize','MaxSize','MinSize','SumSize', '10Size', 'MedSize', '90Size', 'AvgArea','MaxArea','MinArea','SumArea', '10Area', 'MedArea', '90Area', 'AvgRR','MaxRR','MinRR','SumRR', '10RR', 'MedRR', '90RR', 'AvgRJ','MaxRJ','MinRJ','SumRJ', '10RJ', 'MedRJ', '90RJ', 'AvgCos','MaxCos','MinCos','SumCos', '10Cos', 'MedCos', '90Cos', 'AvgSin','MaxSin','MinSin','SumSin', '10Sin', 'MedSin', '90Sin')

cat(week)

