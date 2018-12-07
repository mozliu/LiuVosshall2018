function dwell = manual_dwell(name)

path = 'D:\Dropbox\labwork\behavior\heatseeking-pilot\';
%path = '/home/mliu/Dropbox/labwork/behavior/heatseeking-pilot/';
dates = {'0627','0628','0629','0711','0712','0713'};
% Peltier = [1.369011799410029e+02,1.512964601769910e+02;...
%     4.682286135693216e+02,1.352492625368730e+02;...
%     4.710604719764012e+02,4.118274336283185e+02;...
%     1.425648967551622e+02,3.957802359882005e+02];
set(0,'DefaultFigureWindowStyle','docked');

% find which temperature to sample
templabels = [400,26,28.5,31,33.5,36,38.5,40,45,50,55,60];
temp = 36;
temp_n = find(templabels == temp);
temp_i = (540+720*(temp_n-1)+1):(540+720*(temp_n-1)+180);

try
    load([path name '.mat'],'dwell');
    count = length(dwell)+1;
catch
    dwell = struct('dir_name',{},'fn',[],'landing',[],'takeoff',[],'dmoz',[]);
    count = 1;
end

% selecting the data
for i=1:length(dates)
    order = randperm(3);
    for j=order
        cur_dir = [path '2018' dates{i} '_' num2str(j)];
        dir_name = ['2018' dates{i} '_' num2str(j)];

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

        gcf;
        KEY_IS_PRESSED = 0;
        for k=temp_i
            prev = ix(k-1);
            next = ix(k);
            dwell(count).dir_name = dir_name;
            dwell(count).fn = k;
            dwell(count).dmoz = (count_data(k).right+count_data(k).on_dot) - (count_data(k-1).right+count_data(k-1).on_dot);

            set(gcf, 'KeyPressFcn', @myKeyPressFcn)
            imshowpair(imread(fullfile(cur_dir,fn{prev})),imread(fullfile(cur_dir,fn{next})));
            hold on;
            scatter(count_data(k-1).pos(:,1),count_data(k-1).pos(:,2),'mo');
            scatter(count_data(k).pos(:,1),count_data(k).pos(:,2),'go');
            title([num2str(k) ' landing, dmoz = ' num2str(dwell(count).dmoz)]);
            %impoly(gca,Peltier);

            while ~KEY_IS_PRESSED
                %disp(KEY_IS_PRESSED);
                h = impoint();
                setColor(h,'g');
                if KEY_IS_PRESSED
                    break;
                end
                wait(h);
                dwell(count).landing = [dwell(count).landing ; getPosition(h)];
            end
            %disp(KEY_IS_PRESSED);

            KEY_IS_PRESSED = 0;
            title([num2str(k) ' takeoff, dmoz = ' num2str(dwell(count).dmoz)]);
            while ~KEY_IS_PRESSED
                %disp(KEY_IS_PRESSED);
                h = impoint();
                setColor(h,'m');
                if KEY_IS_PRESSED
                    break;
                end
                wait(h);
                dwell(count).takeoff = [dwell(count).takeoff ; getPosition(h)];
            end
            %disp(KEY_IS_PRESSED);

            KEY_IS_PRESSED = 0;
            count = count+1;
            clf;
            
            save([path name],'dwell');
        end
        
    end
end

function myKeyPressFcn(hObject, event)
    KEY_IS_PRESSED  = 1;
    %disp('next frame...')
end

end

