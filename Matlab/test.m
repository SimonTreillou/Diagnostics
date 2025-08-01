T = readtable('/Users/simon/Desktop/param_global.csv');

T(ismember(T.iribarren,999.0),:)=[];
T(ismember(T.neg_frac,999.0),:)=[];
T(ismember(T.NomRun,*"KUMAR"*),:)=[];

irB = T.iribarren;
NF  = T.neg_frac;
LD  = T.ld_std;
OWK = T.IsOWOk_;

scatter(irB,NF);

scatter(irB,LD);

scatter(OWK,NF);

scatter(OWK,LD);