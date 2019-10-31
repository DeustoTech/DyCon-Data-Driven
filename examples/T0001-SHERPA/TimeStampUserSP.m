clear all
ImportGolden_users_SP
dataSP = goldenusersSP;
data   = 
userid = dataSP.userid(1);

datauser = dataSP(dataSP.userid==userid,:);