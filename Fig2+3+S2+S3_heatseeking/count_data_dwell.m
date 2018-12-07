function count_data2 = count_data_dwell()

path = 'D:\Dropbox\labwork\behavior\heatseeking-pilot\';
load([path 'dwell36.mat'],'dwell');
dates = unique({dwell.dir_name});
temp_i = unique([dwell.fn]);

% blocking out the dot
% load([path 'ellipse.mat'],'polypos');
% [x,y] = meshgrid(1:640,1:480);
% block = inpolygon(x,y,polypos(:,1),polypos(:,2));
% clear x y polypos

% if you don't have a dot to block out (untested):
block = meshgrid(1:5,1:5);

% load model
load('D:\Dropbox\labwork\scripts\heatseeking\hsfit.mat','model');

count_data2 = struct('dir_name',{},'fn',[],'right',[],'pos',[],'landing',[],'takeoff',[]);
count = 1;

% selecting the data
for i=1:length(dates)
    cur_dir = [path dates{i}];
    dir_name = dates{i};
    disp(['Loading ' dir_name '...']);

    % get files
    load(fullfile(cur_dir,[dir_name '_count']),'filelist','count_data');

    % sort filelist according to the order of the images
    fn = {filelist.name};
    fn = fn(contains(fn,'R_'));
    file_order = regexp(fn,'\d+.tif$','match');
    file_ind = ones(1,length(file_order));
    for k=1:length(file_ind)
        file_ind(k) = sscanf(char(file_order{k}),'%i');
    end
    [~,ix] = sort(file_ind);

    BR = imread([cur_dir '\' dir_name '_background\R_' dir_name '_background.tif']);

    for k=temp_i
        next = ix(k);
        count_data2(count).dir_name = dir_name;
        count_data2(count).fn = k;
        count_data2(count).landing = dwell(count).landing;
        count_data2(count).takeoff = dwell(count).takeoff;
        
        I = imread([cur_dir '/' fn{next}]);
        [count_data2(count).right,~,count_data2(count).pos] = count_pelt_model(I,BR,block,model);
        count_data2(count).right = count_data2(count).right + count_data(k).on_dot;

        count = count+1;
    end
end

end
