Steps
* cmsrel CMSSW_10_1_7
* cd CMSSW_10_1_7/src
* cmsenv
* voms-proxy-init --voms cms
* export XRD_NETWORKSTACK=IPv4
* git clone https://github.com/pkalbhor/SusyAnalyzer.git
* scram b
* cmsRun python/Conf_cfg.py
