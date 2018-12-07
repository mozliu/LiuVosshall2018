% This MATLAB script analyzes images from the heat-seeking assay. It analyzes each image to record the position of each mosquito and count the number of mosquitoes in the Peltier area. A very similar script is used to count mosquitoes in the heat-seeking choice assay.

%Requires the Local Adaptive Thresholding script available in the MATLAB
%file exchange:
%http://www.mathworks.com/matlabcentral/fileexchange/8647-local-adaptive-thresholding

function [count_data, filelist] = heatseeking_assay_count(cur_dir)

%prompts selection of directory
if nargin < 1
    cur_dir = uigetdir('D:\Dropbox\labwork\behavior\heatseeking-pilot\');
end
filelist = dir(cur_dir);
dir_name = regexp(cur_dir,'\','split');
dir_name = dir_name{end};

%create waitbar for processing
%countWait = waitbar(0,['processing ' num2str(length(filelist)) ' images']);

% finds background image in folder within current directory
BR = imread([cur_dir '\' dir_name '_background\R_' dir_name '_background.tif']);

count_data = struct('right',[],'right_pelt',[],'on_dot',0,'pos',[]);

for i=1:length(filelist)
    
    %countWait = waitbar(i/length(filelist),countWait);
    
    if strncmp(filelist(i).name,'R',1) %checks if Right
        
        %extracts picture number
        split_name = regexp(filelist(i).name, '_','split');
        split_number = regexp(split_name(end), '\.','split');
        split_number_str = split_number{1}(1);
        number_str = split_number_str;
        number = str2num(number_str{:});
        
        %sets the right entry to the count for the right picture of this number
        I = imread([cur_dir '/' filelist(i).name]);
        [count_data(number).right, count_data(number).right_pelt, count_data(number).pos] = count_pelt(I,BR);
        count_data(number).on_dot = 0;
    end
end

%close(countWait);
%clear countWait

% save for MATLAB
save([cur_dir '\' dir_name '_count']);

% save for Python
load([cur_dir '\' dir_name]);
Rtot = ([count_data.on_dot] + [count_data.right_pelt])/mosq;
Rtemp = temparray(:,2)/100;
save([cur_dir '\' dir_name '_topy'],'Rtot','Rtemp');

% plot the results
figure;
subplot(2,1,1);
plot(Rtemp);
ylabel('Pelt T (C)');
subplot(2,1,2);
plot(Rtot*100);
xlabel('time (s)');
ylabel('% on Pelt');

end