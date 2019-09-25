function matrixplots(mats)
%MATRIXPLOTS Summary of this function goes here
%   Detailed explanation goes here

NumPlots = length(mats);
for i=1:NumPlots
    
    subplot(NumPlots,4,4*(i-1) + 1)
    surf(mats(i).A)
    
    subplot(NumPlots,4,4*(i-1) + 1)
    surf(mats(i).B)
    
    subplot(NumPlots,4,4*(i-1) + 1)
    surf(mats(i).C)

    
    subplot(NumPlots,4,4*(i-1) + 1)
    surf(mats(i).D)

    
end

