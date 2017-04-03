% Demo for Structured Edge Detector (please see readme.txt first).
addpath(genpath('toolbox-master/'))
% addpath('/media/New_Volume/Research/Print_PDF/')
% addpath(genpath('/media/New_Volume/Research/Project_contour/third_order/TO-edge-detectorToolbox'))
addpath ../TO-edge-detectorToolbox/
%% set opts for training (see edgesTrain.m)
opts=edgesTrain();                % default options (good settings)
opts.modelDir='models/';          % model will be in models/forest
opts.modelFnm='modelBsds';        % model name
opts.nPos=5e5; opts.nNeg=5e5;     % decrease to speedup training
opts.useParfor=0;                 % parallelize if sufficient memory

%% train edge detector (~20m/8Gb per tree, proportional to nPos/nNeg)
tic, model=edgesTrain(opts); toc; % will load model if already trained

%% set detection parameters (can set after training)
model.opts.multiscale=0;          % for top accuracy set multiscale=1
model.opts.sharpen=2;             % for top speed set sharpen=0
model.opts.nTreesEval=4;          % for top speed set nTreesEval=1
model.opts.nThreads=4;            % max number threads for evaluation
model.opts.nms=true;                 % set to true to enable nms

%% evaluate edge detector on BSDS500 (see edgesEval.m)
if(0), edgesEval( model, 'show',1, 'name','' ); end

%% detect edge and visualize results
th = 0.07;
sigma = 0.7;
src_path = '../../Data/CFGD_img/';
% src_path = '/media/New_Volume/Research/Datasets/capital_dataset_resized/jpeg/';
% src_path = '/media/New_Volume/Research/Datasets/Illumination_selected/SET003/';
% src_path = '/media/New_Volume/Research/Datasets/random_noise/sigma5/';
% src_path = 'Data/CUB14/all_imgs/';
img_files = dir([src_path '*.jpg']);
% out_path = 'results/CUB14/';
% out_path = 'results/capital_dataset_resized_test/';
% out_path = 'results/Illumination_selected/SET003/';
% out_path = 'results/random_noise/sigma5/';
out_path = '../../Data/SE_SEL_CFGD/edges/';
mkdir(out_path);
for i = 1:length(img_files)
    i
    imgFile = [src_path img_files(i).name(1:end-4) '.jpg'];
    outFile = [out_path img_files(i).name(1:end-4) '.edg'];
    I = imread(imgFile);
%     tic, [E,O,~,~]=edgesDetect(I,model); toc
    tic,[E,O,edginfo,inds,segs] = edgesDetect_TO(I,model, th, sigma); toc
%     [xx,yy,mag] = find(E~=0);
%     Angle = zeros(length(xx),1); 
%     for j=1:length(xx)
%         Angle(j) = O(xx(j),yy(j));
%     end
%     edge_map  = [xx yy Angle mag];
%     figure(1); imshow(I, 'border', 'tight');
%     figure(2); imshow(1-E, 'border', 'tight');
%     pause(.2)
    
    outImg = [out_path img_files(i).name(1:end-4) '.png'];
    imwrite(E, outImg, 'PNG');
    
    % change edges to cxx coodinates
    edginfo(:,1) = edginfo(:,1)-1;
    edginfo(:,2) = edginfo(:,2)-1;
    save_edg(outFile, edginfo, [size(I,2) size(I,1)]);
   
end