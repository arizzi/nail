import collections

nbins=30
#ordered with increasing priority
binningrules=[
("n.*" , "%s , 0 , 30"%nbins),
("N.*" , "%s , 0 , 30"%nbins),
(".*_pt(_|$).*" , "%s , 0 , 300"%nbins),
(".*HT.*" , "%s , 0 , 300"%nbins),
(".*_eta(_|$).*" , "%s , -5 , 5"%nbins),
(".*_phi(_|$).*" , "%s , -3.1415 , +3.1415"%nbins),
(".*_m(_|$).*" , "%s , 0,500"%nbins),
(".*_M(_|$).*" , "%s , 0,500"%nbins),
(".*_Mass(_|$).*" , "%s , 0,500"%nbins),
(".*Higgs_m.*" , "%s , 70,150"%(nbins*3)),
(".*_.*tag.*" , "%s , 0,1"%nbins),
(".*Class.*" , "%s , 0,1"%nbins),
(".*Atan.*" , "%s , 0,2"%(nbins*2)),
(".*Soft.*" , "10 , -0.5,9.5"),
(".*Mqq_log.*" , "%s , 0,10"%nbins),
(".*zstar.*" , "%s , -2,2"%nbins),
(".*mmjj_pz.*" , "%s , 0,5000"%nbins),
(".*mmjj_pz_logabs.*" , "%s , 0,10"%nbins),
(".*_pt_log.*" , "%s , 0,10"%nbins),
(".*DeltaEta.*" , "%s , 0,10"%nbins),
(".*theta.*" , "%s , -1,1"%nbins),
(".*AbsEta.*" , "%s , 0,5"%nbins),
("PV_npvs.*" , "%s , 0,60"%nbins),
("LeadingSAJet_pt.*" , "%s , -10 , 150"%16),

]

