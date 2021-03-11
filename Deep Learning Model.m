%% Convert DICOM Images and RS structures to Matfiles

% The place where image mat file to be saved, please change this with your folder
imageDir = 'XXXXX/Images';
% The place where label mat file (CTV structure) to be saved, please change this with your folder
labelDir = 'XXXXX/Labels';


% Replace XXXXX with your DICOM folder
dicomDir = 'XXXXXX';
collection = dicomCollection(dicomDir,'IncludeSubfolders',true);

% Convert DICOM images to MatFiles
CTcollection = collection(collection.Modality=='CT',:);
% CT images to Matfiles (volumes)
for idx = 1:numel(CTcollection.Row)
    dicomFileName = CTcollection.Filenames{idx};
    if length(dicomFileName) > 1
        matFileName = fileparts(dicomFileName(1));
        matFileName = split(matFileName,filesep);
        matFileName = replace(strtrim(matFileName(end))," ","_");
    else
        [~,matFileName] = fileparts(dicomFileName);
    end
    matFileName = fullfile(imageDir+'/images',matFileName);
    
    try
        V = dicomreadVolume(CTcollection,CTcollection.Row{idx});
    catch ME
        V = dicomread(dicomFileName);
        if ndims(V)<4
            % Skip files that are not volumes
            continue;
        end
    end
    
    % For multi-file DICOM, dicomFileName is a string array.    
    save(matFileName,'V','dicomFileName');
    
end

% Convert DIcom RS structure to MatFiles
paramS = getRadiomicsParamTemplate('JBS_radiomics_settings2.json');
RScollection = collection(collection.Modality=='RTSTRUCT',:);
for idx = 1:numel(RScollection.Row)
    [currentFolder,name,ext] = fileparts(RScollection.Filenames{idx});
    RSfile = RScollection.Filenames{idx};
    
    % Load PlanC matfiles to find struct number
    pathparts = strsplit(RSfile,'\');
    ptno = pathparts(5);
    % Please chagge planC_Matfile_Folder with your folder containing
    % plan C matfile generated from 'Radiogenomics Score Model.m'
    matfile = fullfile('./planC_Matfile_Folder',strcat(ptno,'.mat'));
    planC = loadPlanC(matfile, tempdir);    
    strC = {planC{4}.structureName};
    % Please include the RF CTV names such as 'GTV2' or 'RF' etc..
    structures = {"GTV2", "GTV-HD", "GTVp", "RF", "GTV_54.6"};
    for iStr = 1:length(structures)
        structNum = getMatchingIndex(structures{iStr},strC,'exact');
        if ~isempty(structNum)
            break;
        end
    end
        
    indexS = planC{end};
    maskStruct3M = getUniformStr(structNum,planC); 
    mask3Dgtv = maskStruct3M(:,:,:) == 1;
    currentFile = fullfile(labelDir,strcat(ptno,'.mat'));
    save(currentFile,'mask3Dgtv');
    
end

%% To gernerate CTV volume matfile (3D). 

% Preprocessing
    volds = imageDatastore(imageDir,'FileExtensions','.mat','ReadFcn',@matRead);
    classNames = ["background","tumor"];
    pixelLabelID = [0 1];
    pxds = pixelLabelDatastore(labelDir, classNames, pixelLabelID,'FileExtensions','.mat','ReadFcn',@matLabelRead);
    reset(volds);
    reset(pxds);
    
% Crop relevant region
    NumFiles = length(pxds.Files);
    id = 1;
    reset(volds);
    reset(pxds);
    while hasdata(pxds)
        outL = readNumeric(pxds);
        outV = read(volds);
        temp = outL>0;
        sz = size(outL);
        reg = regionprops3(temp,'BoundingBox');
        tol = 132;
        ROI = ceil(reg.BoundingBox(1,:));
        ROIst = ROI(1:3) - tol;
        ROIend = ROI(1:3) + ROI(4:6) + tol;
        
        ROIst(ROIst<1)=1;
        ROIend(ROIend>sz)=sz(ROIend>sz);          
        
        tumorRows = ROIst(2):ROIend(2);
        tumorCols = ROIst(1):ROIend(1);
        tumorPlanes = ROIst(3):ROIend(3);
        
        tcropVol = outV(tumorRows,tumorCols, tumorPlanes,:);
        tcropLabel = outL(tumorRows,tumorCols, tumorPlanes);

        % Data set with a valid size for 3-D U-Net (multiple of 8).
        ind = floor(size(tcropVol)/8)*8;
        incropVol = tcropVol(1:ind(1),1:ind(2),1:ind(3),:);
        mask = incropVol == 0;
        cropVol = channelWisePreProcess(incropVol);
        
        % Set the nonbrain region to 0.
        cropVol(mask) = 0;
        cropLabel = tcropLabel(1:ind(1),1:ind(2),1:ind(3));    
        [filepath,name,ext] = fileparts(volds.Files{id});
        
        cropVol(cropLabel==0)=0;
        save([preDir '\' name '.mat'],'cropVol');
       
        id=id+1;       
    end
    
% Preview
volume = preview(volds);
label = preview(pxds);

viewPnl = uipanel(figure,'Title','Labeled Training Volume');
hPred = labelvolshow(label,volume(:,:,:,1),'Parent',viewPnl, ...
    'LabelColor',[0 0 0;1 0 0],'VolumeOpacity',0.01, 'VolumeThreshold', 0.1);
hPred.LabelVisibility(1) = 0;

% Mask segmentation
volume(label=='background')=0;
cropVol(cropLabel==0)=0;
mat
temp=preview(dsTest);
stats = regionprops3(temp,'all')

% After this, files should be divided according to KRAS or MSI labels.
% Subfolder name represents label.
% Example of folder structure:
% /KRAS -|-- 0_Wild    (K-ras wild type CTV matfiles generated from above)
%        |-- 1_Mutant  (K-ras mutant type CTV matfiles generated from above)
 

%% Deep learning Model structure 

Inputsize = [160 160 80];

layers = [image3dInputLayer(Inputsize,'Name','inputLayer','Normalization','none'),...
    convolution3dLayer(8,32,'Stride',4,'Name','Conv3'),...
    leakyReluLayer(0.1,'Name','leakyRelu3'),...
    convolution3dLayer(4,32,'Stride',2,'Name','Conv1'),...
    leakyReluLayer(0.1,'Name','leakyRelu1'),...
    convolution3dLayer(2,32,'Stride',1,'Name','Conv2'),...
    leakyReluLayer(0.1,'Name','leakyRulu2'),...
    maxPooling3dLayer(2,'Stride',2,'Name','maxPool'),...
    fullyConnectedLayer(400,'Name','fc1'),...
    reluLayer('Name','relu'),...
    dropoutLayer(0.5,'Name','dropout1'),...
    fullyConnectedLayer(2,'Name','fc2'),...
    softmaxLayer('Name','softmax'),...
    classificationLayer('Name','crossEntropyLoss')];
voxnet = layerGraph(layers);
voxnet



%% 5-Fold testing

inputDir = 'XXXX/Kras';
% MSI prediction 
inputDir = 'XXXX/MSI';

ds = imageDatastore(inputDir,'FileExtensions','.mat','IncludeSubfolders',true, 'LabelSource', 'foldernames','ReadFcn',@matRead);

% For 10-fold
[imd1 imd2 imd3 imd4 imd5 imd6 imd7 imd8 imd9 imd10] = splitEachLabel(ds,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,'randomize');
% For 2-fold
[imd1 imd2 ] = splitEachLabel(ds,0.5,'randomize');

partStores{1} = imd1.Files ;
partStores{2} = imd2.Files ;
partStores{3} = imd3.Files ;
partStores{4} = imd4.Files ;
partStores{5} = imd5.Files ;
partStores{6} = imd6.Files ;
partStores{7} = imd7.Files ;
partStores{8} = imd8.Files ;
partStores{9} = imd9.Files ;
partStores{10} = imd10.Files; 

% Please chage the K-value according to fold number
k = 10;
idx = crossvalind('Kfold', k, k);
intervals= linspace(0, 1, 100);
cla;
pause(2);

for i = 1:k
    i
    test_idx = (idx == i);
    train_idx = ~test_idx;
    
    imdsTest = imageDatastore(partStores{test_idx}, 'IncludeSubfolders', true, 'FileExtensions','.mat','LabelSource', 'foldernames','ReadFcn',@matRead);
    imdsTrain = imageDatastore(cat(1, partStores{train_idx}), 'IncludeSubfolders', true,'FileExtensions','.mat', 'LabelSource', 'foldernames','ReadFcn',@matRead);

    
    miniBatchSize = 1;
    numObservations = numel(imdsTrain.Labels);
    numIterationsPerEpoch = floor(numObservations / miniBatchSize);

    options = trainingOptions('rmsprop', ...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',5, ...
    'InitialLearnRate',1e-4, ...
    'Shuffle','never', ...
    'Verbose', true, ...
    'Plots','none', ...
    'ResetInputNormalization',false,...
    'ValidationData',imdsTest,...
    'ValidationFrequency', numIterationsPerEpoch)

    trained_net = trainNetwork(imdsTrain,voxnet,options);
        
    [predicted_labels,posterior]=classify(trained_net, imdsTest, 'MiniBatchSize',miniBatchSize);
    actual_labels=imdsTest.Labels;
    test_labels=double(nominal(imdsTest.Labels));
    accuracy = mean(predicted_labels == actual_labels)
    [fp_rate,tp_rate,T,AUC_fold(i), OPTROCPT]=perfcurve(test_labels,posterior(:,1),1);
    plot(fp_rate, tp_rate, 'LineWidth', 1.5); legends{i}= sprintf('Fold %d (AUC = %.2f)', i, AUC_fold(i)); hold on
    x_adj= adjust_unique_points(fp_rate); %interp1 requires unique points
    if i==1 %if is the first fold 
        mean_curve= (interp1(x_adj, tp_rate, intervals))/k; 
    else
        mean_curve= mean_curve+ (interp1(x_adj, tp_rate, intervals))/k; 
    end
end

    grid on;
    figure(1); plot(intervals, mean_curve, 'Color', 'Black', 'LineWidth', 3.0); 
    legends{i+1}= sprintf('Average (AUC = %.2f)', mean(AUC_fold));
    legend(legends)

