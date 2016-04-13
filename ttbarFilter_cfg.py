import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("opendata")

goodJSON = 'Cert_136033-149442_7TeV_Apr21ReReco_Collisions10_JSON_v2.txt'
myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')

import FWCore.Utilities.FileUtils as FileUtils

singleMuFiles = FileUtils.loadListFromFile('CMS_Run2011A_SingleMu_AOD_12Oct2013-v1_10000_file_index.txt')
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(*singleMuFiles)
                            )

process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
process.source.lumisToProcess.extend(myLumis)

process.load("TTBarFilter.TTBarFilter.TTBarFilter_cfi")

process.TTBarFilter.csvFileName = cms.string("SingleMu_Run2011A.csv")
process.TTBarFilter.minPt = cms.double(20.0)
process.TTBarFilter.maxEta = cms.double(2.1)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))

process.mypath = cms.Path(process.TTBarFilter)
process.schedule = cms.Schedule(process.mypath)
