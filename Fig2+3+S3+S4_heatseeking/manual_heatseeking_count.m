function [count_data,filelist] = manual_heatseeking_count(cur_dir)

set(0,'DefaultFigureWindowStyle','docked');

% get list of files
if nargin < 1
    cur_dir = uigetdir('D:\Dropbox\labwork\behavior\heatseeking-pilot\');
end
dir_name = regexp(cur_dir,'\','split');
dir_name = dir_name{end};

% load counted data
load([cur_dir '\' dir_name '_count'],'count_data','filelist');

% sort filelist according to the order of the images
fn = {filelist.name};
fn = fn(contains(fn,'R_'));
file_order = regexp(fn,'\d+.tif$','match');
file_ind = ones(1,length(file_order));
for i=1:length(file_ind)
    file_ind(i) = sscanf(char(file_order{i}),'%i');
end
[~,ix] = sort(file_ind);

% display image and count mosquitoes on the dot
figure;
for i=ix
    n = file_ind(i); % "true" number
    imshow(imread([cur_dir '\' fn{i}]));
    title(num2str(n));
    hold on;
    if ~isempty(count_data(n).pos)
        scatter(count_data(n).pos(:,1),count_data(n).pos(:,2),'ro');
    end
    hold off;
    count_data(n).on_dot = input('> ');
    
    if isempty(count_data(n).on_dot)
        count_data(n).on_dot = count_data(n-1).on_dot; % assumes empty vector is same as last
    elseif count_data(n).on_dot >= 10
        count_data(n).on_dot = mod(count_data(n).on_dot,10); % over 10 is unrealistically high, so take the last digit
    end
end

% save for MATLAB
save([cur_dir '\' dir_name '_count']);

% save for Python
load([cur_dir '\' dir_name],'mosq','temparray');
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
