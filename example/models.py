background={
"DY0J":["DY0J"],
"DY1J":["DY1J"],
"DY2J":["DY2J"],
"EWKZ":["EWKZ","EWKZint"],
"Top":["STs","STwtbar","STwt","STtbar","STt","TTlep","TTsemi"],
"Other":["W2J","W0J","WWdps","WW2l2n","WZ1l1n2q","WZ1l3n","WZ2l2q","WZ3l1n",
#"WZ",
##"ZZ2l2n", #manca
"ZZ2l2q",
#"ZZ2q2n", missing xsec
#"ZZ4l"
],
#"ZZ"
}

#sorting
backgroundSorted=["Other","Top","DY0J","DY1J","DY2J","EWKZ"]
backgroundSorted+=[x for x in background if x not in backgroundSorted]


signal={
"VBF H":["vbfHmm"],
"gg H":["ggHmm"]
}

data={
"2017":["data"]
}

import ROOT
fillcolor={
"DY0J": ROOT.kOrange+2,
"DY1J": ROOT.kOrange+1,
"DY2J": ROOT.kOrange,
"EWKZ": ROOT.kViolet,
"Top": ROOT.kGreen,
"Other" : ROOT.kTeal,
"VBF H":ROOT.kRed,
"gg H":ROOT.kRed+4,
}

systematicsToPlot=["MuScaleUp"]

linecolor=fillcolor
markercolor=fillcolor

toblind={}
