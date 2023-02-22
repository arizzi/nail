from .nail import *
import ROOT

flow = SampleProcessing(
    "Simple Test", "/scratch/arizzi/0088F3A1-0457-AB4D-836B-AC3022A0E34F.root")


flow.SubCollection("GenMuon", "GenPart", "abs(GenPart_pdgId)==13")
flow.MatchDeltaR("GenMuon", "GenJet")
flow.Define("GenJet_goodGenMuonIdx",
            "ROOT::VecOps::Where(GenJet_GenMuonDr<0.4,GenJet_GenMuonIdx,-1)")

nthreads = 10
histos = {}
targets = ["nGenJet", "GenJet_pt", "GenJet_eta", "GenJet_phi", "GenJet_GenMuonDr",
           "GenJet_GenMuonIdx", "nGenMuon", "GenMuon_pt", "GenMuon_eta", "GenMuon_phi", "GenJet_goodGenMuonIdx"]
processor = flow.CreateProcessor(
    "eventProcessor", targets, histos, [], "", nthreads)


rdf = ROOT.RDataFrame(
    "Events", "/scratch/arizzi/0088F3A1-0457-AB4D-836B-AC3022A0E34F.root")
result = processor(rdf)

processed_rdf = result.rdf.find("").second
processed_rdf.Snapshot("Events", "out.root", targets)
