%%
% This lines add all dependences of Data Drive 

% add DyCon Toolbox
addpath(genpath('C:\Users\JorgeJuan\Documents\GitHub\DyCon-toolbox'))
% Add folders 
File = 'StartDriveData.m';

path = which(File);
path = replace(path,File,'');

addpath(genpath(path))

%%
unzip('https://github.com/rezaahmadzadeh/Expectation-Maximization/archive/master.zip')