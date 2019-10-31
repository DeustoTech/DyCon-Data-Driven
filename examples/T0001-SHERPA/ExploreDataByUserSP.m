clear all
ImportGolden_users_SP
data = goldenusersSP;
data = data(:,[3 4 6 8]);
%%
usersids = unique(data.userid);
states = unique(data.O_tag);
days   = ["lunes","martes","miércoles","jueves","viernes","sábado","domingo"]';
hours  = 0:23;
%%
uid = usersids(3);
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
    gpi = interp1(gptime(2,:),gptime(1,:),1:24,'nearest','extrap');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    TableOfProbabilities(iter).Hours         = StateData;
    TableOfProbabilities(iter).Probabilities = Probabilities;
    TableOfProbabilities(iter).MP = gp;
    %% Interpolate
    TableOfProbabilities(iter).MPi = gpi;
end
%%
surf(reshape([TableOfProbabilities(:).MPi],24,7))
view(0,-90)


%%
DH = [];
for i = 1:7
stringdays  = repmat(days(i),24,1);
stringhours = num2str((0:23)','%.2d');
DH = [DH;strcat(stringdays,stringhours)];
end

