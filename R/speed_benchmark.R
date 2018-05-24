library(calibrate)


netEMDSpeedTest <- function()
{
    ##load the data
    source_dir <- system.file(file.path("extdata", "random"), package = "netdist")
    print(source_dir)
    edge_format = "ncol"
    file_pattern = ""
#    source_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")
#    edge_format = "ncol"
#    file_pattern = ".txt"
    graphs <- read_simple_graphs(source_dir = source_dir, format = edge_format, pattern = file_pattern)
    n1=names(graphs)
    lab1=c()
    gddBuildTime=c()
    netEMDtime=c()
    for (i in 1:length(graphs))
    {
        for (j in 1:(i))
        {
            g1=graphs[[i]]
            g2=graphs[[j]]
            lab1=append(lab1,paste(n1[i],n1[j],sep=','))
            print(paste(n1[i],n1[j],sep=','))
            fulltimeStart=Sys.time()
            gdd1=gdd(g1)
            gdd2=gdd(g2)
            netEMDStart=Sys.time()
            net_emd(gdd1,gdd2)
            endTime=Sys.time()
            gddBuildTime=append(gddBuildTime,as.double(netEMDStart-fulltimeStart))
            netEMDtime=append(netEMDtime,as.double(endTime-netEMDStart))
        }
    }

    result1=list(lab1,gddBuildTime,netEMDtime)
#    textxy(gddBuildTime,netEMDtime,lab1)
    combinedTimes=c()
    gddTimes=c()
    netEMDTimes=c()
    indices=c()
    netEMDFullTime=c()
    for (i in 2:length(graphs))
    {
        graphs1=graphs[1:i]
        fulltimeStart=Sys.time()
        gdds=lapply(graphs1,gdd)
        netEMDStart=Sys.time()
        net_emds_for_all_graphs(gdds)
        endTime=Sys.time()
        combinedTimes=append(combinedTimes,as.double(endTime-fulltimeStart))
        gddTimes=append(gddTimes,as.double(netEMDStart-fulltimeStart))
        netEMDTimes=append(netEMDTimes,as.double(endTime-netEMDStart))
        netEMDFullTime=append(netEMDFullTime,as.double(endTime-fulltimeStart))
        indices=append(indices,i)
    }
    result2=list(indices,netEMDFullTime)
    result3=list(result1,result2)
    result3
}


# This function runs the speed tests and makes a table
#' @export
netEMDSpeedTestPlot <- function()
{
    ##load the data
    source_dir <- system.file(file.path("extdata", "random"), package = "netdist")
    print(source_dir)
    edge_format = "ncol"
    file_pattern = ""
#    source_dir <- system.file(file.path("extdata", "VRPINS"), package = "netdist")
#    edge_format = "ncol"
#    file_pattern = ".txt"
    graphs <- read_simple_graphs(source_dir = source_dir, format = edge_format, pattern = file_pattern)
    n1=names(graphs)
    lab1=c()
    gddBuildTime=c()
    netEMDtime=c()
    for (i in 1:length(graphs))
    {
        for (j in 1:(i))
        {
            g1=graphs[[i]]
            g2=graphs[[j]]
            lab1=append(lab1,paste(n1[i],n1[j],sep=','))
            print(paste(n1[i],n1[j],sep=','))
            fulltimeStart=Sys.time()
            gdd1=gdd(g1)
            gdd2=gdd(g2)
            netEMDStart=Sys.time()
            net_emd(gdd1,gdd2)
            endTime=Sys.time()
            gddBuildTime=append(gddBuildTime,as.double(netEMDStart-fulltimeStart))
            netEMDtime=append(netEMDtime,as.double(endTime-netEMDStart))
        }
    }
quartz()
    plot(gddBuildTime,netEMDtime)
#    textxy(gddBuildTime,netEMDtime,lab1)
    combinedTimes=c()
    gddTimes=c()
    netEMDTimes=c()
    for (i in 2:length(graphs))
    {
        graphs1=graphs[1:i]
        fulltimeStart=Sys.time()
        gdds=lapply(graphs1,gdd)
        netEMDStart=Sys.time()
        net_emds_for_all_graphs(gdds)
        endTime=Sys.time()
        combinedTimes=append(combinedTimes,as.double(endTime-fulltimeStart))
        gddTimes=append(gddTimes,as.double(netEMDStart-fulltimeStart))
        netEMDTimes=append(netEMDTimes,as.double(endTime-netEMDStart))
    }
quartz()

    plot(combinedTimes)
quartz()

    plot(gddTimes,netEMDTimes)
}
