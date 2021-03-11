%% Convert DICOM to plan C object (which should be performed at 2020a)
% Generated plan C matfiles can be used in 'Deep Learning Model.m'.


%Specify source DICOM folder and destination folder
srcDir = 'Dicom object files generated abve code';    
dstDir = 'Destination Folder';

%Set flags for compression and merging
zipFlag = 'No'; %Set to 'Yes' for compression to bz2 zip 
mergeFlag = 'No'; %Set to 'Yes' to merge all scans into a single series

%Batch import
init_ML_DICOM;
batchConvert (srcDir,dstDir,zipFlag,mergeFlag);

%% Load configuration json file and get radiogenomics profile
%Detail JSON structure are found in https://github.com/cerr/CERR/wiki/Radiomics

paramFileName = 'configuration.json';
dirName = 'XXXX';

paramS = getRadiomicsParamTemplate(paramFileName);
outXlsFile = 'results_1120_threshould.xlsx';

featureS = batchExtractRadiomics(dirName,paramFileName,outXlsFile);

%% Write output to Excel file if outXlsFile exists
fileNamC = {featureS.fileName};
combinedFieldNamC = {};
combinedFeatureM = [];
featureForStructS = [featureS.('struct')];
imgC = fieldnames(featureForStructS);
        
% Exclude 'Shape' property due to the almost round shape of 
for iImg = 2:length(imgC)
    [featureM,allFieldC] = featureStructToMat([featureForStructS.(imgC{iImg})]);
    combinedFieldNamC = [combinedFieldNamC; strcat(allFieldC,'_',imgC{iImg})];
    combinedFeatureM = [combinedFeatureM, featureM];
end

xlswrite(outXlsFile, ['File Name';combinedFieldNamC]', 'All', 'A1');
xlswrite(outXlsFile, fileNamC(:), 'All', 'A2');
xlswrite(outXlsFile, combinedFeatureM, 'All', 'B2');
       
 


