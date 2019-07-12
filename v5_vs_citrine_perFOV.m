%   =======================================================================================
%   Copyright (C) 2013  Erlend Hodneland
%   Email: erlend.hodneland@biomed.uib.no 
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   =======================================================================================

% -- FIRST, combine all 10 FOVs horizonally in ImageJ/Fiji -- %
% -- This script assumes 512x512 size images -- %

clear all
close all

% load the image tif files data
filepath='C:/Users/Tina Kim/Dropbox/PAPERS/Nanoluc lov/Data/';
filefolder='20180909_spark2_hatag/';
filename_base='15minlight_iso'; 
filename_v5=strcat(filename_base,'_v5.tif'); 
imsegm=imread(strcat(filepath,filefolder,filename_v5));
imsegm=double(imsegm);
filename_cit=strcat(filename_base,'_cit.tif');
imsegm_cit=imread(strcat(filepath,filefolder,filename_cit));
imsegm_cit=double(imsegm_cit);

% Segmentation by adaptive thresholding
prm.method = 'adth';
prm.adth.th = 0.1;

[cellbw1,wat,imsegmout,prmout] = cellsegm.segmct(imsegm,0,500,'prm',prm);
cellsegm.show(imsegmout,1);title('Raw image');axis off;
cellsegm.show(cellbw1,1);title('Cell segmentation by ADTH');axis off;

% % improving the results by splitting of cells
% splitth = 1;
% plane = 1;
% % cells above this threshold are split (all cells here)
% n = 100;
% h = [0.5 0.5 1.5];
% cellbw2 = cellsegm.splitcells(cellbw1,splitth,n,h);
% cellsegm.show(cellbw2,2);title('Cell segmentation by ADTH with splitting');axis off;
cellbw2=cellbw1;

% create cellmasks
[L,num_cells]=bwlabel(cellbw2);
segcentroids=zeros(num_cells,2);
meanfluo=zeros(num_cells,1);
for a=1:num_cells
    [tempx,tempy]=find(L==a);
    fluo=[]; % calculate each cell's v5 fluorescence
    fluo_cit=[]; % calculate each cell's citrine fluroescence
    cellmasktemp=zeros(size(cellbw2,1),size(cellbw2,2));
    for b=1:length(tempx)
        cellmasktemp(tempx(b),tempy(b))=1;
        fluo(b)=imsegm(tempx(b),tempy(b));
        fluo_cit(b)=imsegm_cit(tempx(b),tempy(b));
    end
    temp=regionprops(cellmasktemp);
    segcentroids(a,:)=temp.Centroid;
    meanfluo(a)=mean(fluo);
    meanfluo_cit(a)=mean(fluo_cit);
    meanratio(a)=mean(fluo_cit)./mean(fluo); % calculate average citrine/v5 ratio
end

% save data 
savename=strcat(filename_base,'.mat');
save(strcat(filepath,filefolder,'results/',savename),'num_cells','meanfluo','meanfluo_cit','meanratio');

%% Calculate the average ratios per FOV

% This assumes FOVs are combined in 1 row x 10 columns.

imsize=512; % pixel width/height of each individual square image.
for z=1 % number of rows in combined image
    temp=find(segcentroids(:,2)>=(z-1)*imsize+1 & segcentroids(:,2)<(z)*imsize);
    for a=1:10 % number of columns in combined image
        tempx=segcentroids(temp,1);
        tempa=find(tempx>=(a-1)*imsize+1 & tempx<(a)*imsize);
        
        v5=mean(meanfluo(tempa));
        cit=mean(meanfluo_cit(tempa));
        FOV_ratios(z,a)=cit/v5;
    end
end

savename=strcat(filename_base,'_FOV.mat');
save(strcat(filepath,filefolder,'results/',savename),'FOV_ratios');
