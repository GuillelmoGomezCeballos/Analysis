#/bin/sh

export SEL=$1;
export NJETS=$2;

if [ $SEL == 3 ] && [ $NJETS == 0 ]; then
root -l -q -b optimalCutszh_53x.C+'(1105,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh105inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,0,0,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1115,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh115inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,0,0,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1125,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh125inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,0,0,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1135,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh135inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,0,0,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1145,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh145inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,0,0,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1175,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh175inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,0,0,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1200,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh200inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,0,0,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1300,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh300inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,0,0,120,2.7,0.25,180)';

root -l -q -b optimalCutszh_53x.C+'(1105,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh105inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,1,0,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1115,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh115inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,1,0,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1125,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh125inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,1,0,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1135,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh135inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,1,0,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1145,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh145inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,1,0,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1175,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh175inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,1,0,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1200,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh200inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,1,0,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1300,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh300inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,1,0,120,2.7,0.25,180)';

elif [ $SEL == 3 ] && [ $NJETS == 1 ]; then
root -l -q -b optimalCutszh_53x.C+'(1105,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh105inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,0,1,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1115,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh115inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,0,1,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1125,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh125inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,0,1,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1135,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh135inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,0,1,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1145,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh145inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,0,1,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1175,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh175inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,0,1,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1200,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh200inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,0,1,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1300,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh300inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,0,1,120,2.7,0.25,180)';

root -l -q -b optimalCutszh_53x.C+'(1105,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh105inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,1,1,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1115,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh115inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,1,1,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1125,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh125inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,1,1,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1135,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh135inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,1,1,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1145,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh145inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,1,1,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1175,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh175inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,1,1,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1200,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh200inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,1,1,120,2.7,0.25,180)';
root -l -q -b optimalCutszh_53x.C+'(1300,17,"ntuples_zh_53x/backgroundA_skim10.root","ntuples_zh_53x/zh300inv.root","ntuples_zh_53x/data_skim10.root","ntuples_zh_53x/hww_syst_skim10_gj.root",3,1,1,120,2.7,0.25,180)';

elif [ $SEL == 4 ] && [ $NJETS == 0 ]; then
root -l -q -b optimalCutszh_53x.C+'(1105,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh105inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,0,0,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1115,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh115inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,0,0,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1125,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh125inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,0,0,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1135,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh135inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,0,0,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1145,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh145inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,0,0,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1175,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh175inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,0,0,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1200,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh200inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,0,0,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1300,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh300inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,0,0,120,2.7,0.25,180)'

root -l -q -b optimalCutszh_53x.C+'(1105,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh105inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,1,0,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1115,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh115inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,1,0,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1125,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh125inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,1,0,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1135,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh135inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,1,0,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1145,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh145inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,1,0,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1175,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh175inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,1,0,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1200,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh200inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,1,0,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1300,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh300inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,1,0,120,2.7,0.25,180)'

elif [ $SEL == 4 ] && [ $NJETS == 1 ]; then
root -l -q -b optimalCutszh_53x.C+'(1105,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh105inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,0,1,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1115,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh115inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,0,1,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1125,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh125inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,0,1,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1135,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh135inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,0,1,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1145,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh145inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,0,1,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1175,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh175inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,0,1,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1200,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh200inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,0,1,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1300,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh300inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,0,1,120,2.7,0.25,180)'

root -l -q -b optimalCutszh_53x.C+'(1105,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh105inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,1,1,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1115,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh115inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,1,1,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1125,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh125inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,1,1,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1135,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh135inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,1,1,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1145,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh145inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,1,1,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1175,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh175inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,1,1,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1200,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh200inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,1,1,120,2.7,0.25,180)'
root -l -q -b optimalCutszh_53x.C+'(1300,17,"ntuples_zh_42x/backgroundZH_skim10.root","ntuples_zh_42x/zh300inv.root","ntuples_zh_42x/data_skim10.root","ntuples_zh_42x/zh_syst_skim10_gj.root",4,1,1,120,2.7,0.25,180)'

fi
