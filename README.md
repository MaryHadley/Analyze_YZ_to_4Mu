#Analyze_YZ_to_4Mu


mkdir < workArea >  
cd < workArea >  
cmsrel CMSSW_10_6_27  
cd CMSSW_10_6_27/src  
cmsenv  
git clone git@github.com:MaryHadley/Analyze_YZ_to_4Mu.git  
cd Analyze_YZ_to_4Mu  
mv * ..  
cd ..  
rm -rf Analyze_YZ_to_4Mu  

#To run code  
root -l -b Fit_goodToUse_basedOn2016MC_fromSL.C++



