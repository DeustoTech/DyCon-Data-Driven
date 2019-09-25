%%
% This lines add all dependences of Data Drive 

% add DyCon Toolbox
addpath(genpath('/home/djoroya/Documentos/GitHub/DyCon-toolbox'))

% Add folders 
File = 'StartDriveData.m';

path = which(File);
path = replace(path,File,'');

addpath(genpath(path))