####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                  LOAD PACKAGES AND FUNCTIONS                           ####
library(xml2)       # to read simulatorConfig xml files

####  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ####
####                           SET PATHS                                    ####
# WARNING: Docker desktop must be running in Windows and user signed in
path_simulation     <- 'G:/GRUPO_BIGDATA/Proyecto_ESSNet Big Data II/Simulations/template'
path_configFiles    <- file.path(path_simulation, 'data/simulatorConfig')
path_resources      <- file.path(path_simulation, 'param/resources')
path_groundTruth    <- file.path(path_simulation, 'param/groundTruth')
path_events         <- file.path(path_simulation, 'data/networkEvents')
path_docker         <- "C:/simulator-master/docker"
#path_docker         <- "C:/simulator/simulator-master-docker/docker"
path_docker_WIN     <- normalizePath(path_docker)
path_docker_output  <- file.path(path_docker, "output")
path_docker_output_WIN <- normalizePath(path_docker_output)

###                            :::::::::::::::                              ####
#####                      SIMULATION PARAMETERS                           #####
simulatorConfigFileNames <- c('antennas.xml', 'map.wkt', 'persons.xml', 
                              'probabilities.xml', 'simulation.xml')

antennas.xml <- as_list(read_xml(file.path(path_configFiles, 'antennas.xml')))$antennas
MNOname <- antennas.xml$antenna$mno_name[[1]]

simulatorOutputFileNames <- list(
  resources = c(paste0('AntennaCells_', MNOname, '.csv'),
                'antennas.csv', 
                paste0('SignalMeasure_', MNOname, '.csv'),
                'grid.csv'),
  groundTruth = c('persons.csv'),
  events = c(paste0('AntennaInfo_MNO_', MNOname, '.csv'))
  
  )

###                            :::::::::::::::                              ####
#####                COPY CONFIG FILES IN DOCKER FOLDER STR                #####
for (fn in simulatorConfigFileNames){
  
  if (file.exists(file.path(path_docker_output, fn))) {
    
    file.remove(file.path(path_docker_output, fn))
    
  }
  file.copy(file.path(path_configFiles, fn),
            file.path(path_docker_output, fn))
  
}

###                            :::::::::::::::                              ####
#####                        EXECUTE SIMULATION                            #####
fileConn <- file(file.path(path_docker, "exe_simulator_docker.cmd"))
commands_console <- paste0(
  "cd ", path_docker_WIN, " \n",
  "docker load -i ", path_docker_WIN, "\\simulator.tar \n", 
  "docker run --rm -v ", path_docker_output_WIN, ":/repo/output ",
  "-t -i simulator -m output/map.wkt -s output/simulation.xml ",
  "-p output/persons.xml -a output/antennas.xml ",
  "-v -o -pb output/probabilities.xml")
writeLines(commands_console, fileConn)
close(fileConn)

shell.exec(file.path(path_docker, "exe_simulator_docker.cmd"))



###                            :::::::::::::::                              ####
#####             MOVE OUTPUT FILES TO SIMULATOR FOLDER STR                #####
for (fn in simulatorOutputFileNames$resources){
  
  file.copy(file.path(path_docker_output, fn),
            file.path(path_resources, fn), overwrite = TRUE)
  
}

for (fn in simulatorOutputFileNames$groundTruth){
  
  file.copy(file.path(path_docker_output, fn),
            file.path(path_groundTruth, fn), overwrite = TRUE)
  
}

for (fn in simulatorOutputFileNames$events){
  
  file.copy(file.path(path_docker_output, fn),
            file.path(path_events, fn), overwrite = TRUE)
  
}


###                            :::::::::::::::                              ####
#####               REMOVE ALL FILES FROM DOCKER FOLDER STR                #####
for (fn in list.files(path_docker_output)){
  
  file.remove(file.path(path_docker_output, fn))
  
}
