function PV = ProbabilityVector(Data,States)

    N = length(Data);
    %
    iter = 0;
    PV = zeros(length(States),1);
    for istate = States'
        iter = iter + 1;
        PV(iter) = sum(Data == istate)/N;
    end
end