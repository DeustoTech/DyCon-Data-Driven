clear all
ImportGolden_users_SP
data = goldenusersSP;
data = data(:,[3 4 6 8]);
%%
usersids = unique(data.userid);

iter = 0;
for iuser = usersids'
    iter = iter + 1;
    TableOfProbabilities{iter} = UserID2TableWeek(data,usersids(iter));
end