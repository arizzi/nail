from analysisflow import flow
from yields import yields
from samples import signal #,background,data

a=Analysis(flow,yields,[signal])
a.Compute(yields=".*",samples=".*",systematics=".*")


