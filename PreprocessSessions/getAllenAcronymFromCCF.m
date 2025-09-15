function [ann, name, acr, layer] = getAllenAcronymFromCCF(coordsArray)

% coordsArray is an N * 3 array where each row is the CCF coordinate of a
% single unit (from allen neuropixels dataset).
% Each row is: [AP, DV, ML]

% example data can be generated with:
% coordsArray = [vertcat(units(1:10).anterior_posterior_ccf_coordinate), vertcat(units(1:10).dorsal_ventral_ccf_coordinate), vertcat(units(1:10).left_right_ccf_coordinate)]

%% location of allen annotation volume and sructure tree 

annotation_volume_location = 'C:\Users\edward.horrocks\Documents\Code\allenCCF\annotation_volume_10um_by_index.npy';
structure_tree_location = 'C:\Users\edward.horrocks\Documents\Code\allenCCF\structure_tree_safe_2017.csv';

%% load the reference brain annotations
tic
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
else
    error('Please specify file locations of annotation_volume_location and structure_tree_location')
end
toc
%% divide the coordsArray by 10 to match annotation volume
coordsArray = round(coordsArray./10);

%% get annotation

nUnits = size(coordsArray,1);

for iunit = 1:nUnits
    
% divide coordinates by 10 to match annotation volume
ap = coordsArray(iunit,1);
dv = coordsArray(iunit,2);
ml = coordsArray(iunit,3);

if ap > 0 && dv > 0 && ml > 0 % check a valid coordinate (e.g. not 'grey')
ann_temp = av(ap,dv,ml); 
ann{iunit} = ann_temp;
name{iunit} = st.safe_name{ann_temp};
acr{iunit} = st.acronym{ann_temp};
layer{iunit} = extractAfter(name{iunit}, 'layer ');

else
ann{iunit} = nan;
name{iunit} = nan;
acr{iunit} = nan;
layer{iunit} = nan;

end

end
