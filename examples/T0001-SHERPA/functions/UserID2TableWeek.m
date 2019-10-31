function TableOfProbabilities = UserID2TableWeek(data,uid)

    %%
    states = unique(data.O_tag);
    days   = ["lunes","martes","miércoles","jueves","viernes","sábado","domingo"]';
    hours  = 0:23;
    %%
    DateByUser = data(data.userid ==uid,1:3);
    %%
    iter = 0;

    for iday = days'

        iter = iter + 1;
        DataByUserByDate = DateByUser(DateByUser.Day == iday,[1 3]);

        TableOfProbabilities(iter).Day = iday; 

        StateData = {};
        Probabilities = {};
        gp = [];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ----------------------- Take Probabilities by Day and Hour ----------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ihour = hours                                                                        %
            StateData{ihour+1} = DataByUserByDate(DataByUserByDate.Hour == ihour,2);             %
            if ~isempty(StateData{ihour+1})                                                      %
                Probabilities{ihour+1} = ProbabilityVector(StateData{ihour+1}.O_tag,states);     %
                [~ ,gp(ihour+1)]= max(Probabilities{ihour+1});
            else                                                                                 %
                Probabilities{ihour+1} = []; 
                gp(ihour+1) = nan;%
            end                                                                                  %
        end                                                                                      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gptime = [gp; 1:24];
        gptime = gptime(:,~isnan(gp));
        try
        gpi = interp1(gptime(2,:),gptime(1,:),1:24,'nearest','extrap');
        catch
            'hola'
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TableOfProbabilities(iter).Hours         = StateData;
        TableOfProbabilities(iter).Probabilities = Probabilities;
        TableOfProbabilities(iter).MP = gp;
        %% Interpolate
        TableOfProbabilities(iter).MPi = gpi;
    end


end



