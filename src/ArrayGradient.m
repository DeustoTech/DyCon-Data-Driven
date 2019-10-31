function dY = ArrayGradient(xline,Y)
    % Y is a matrix [ Nx x Nms ]
    % Nx: Number of Discretization
    % Nms: Number of Measurements 

    [Nx   ,Nms  ] = size(Y);
    [xrows,xcols] = size(xline);
    
    
    dY = zeros(Nx,Nms);
    for j = 1:Nms
        dY(:,j) = gradient(Y(:,j)',xline)';
    end
end

