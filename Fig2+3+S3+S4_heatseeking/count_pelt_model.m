function [count, pelt_count, pos] = count_pelt_model(I,BR,block,model)
% updated version using a classifier instead of manual parameters

% declaring threshold parameters
ws = 10;
C = .05;

Peltier = [1.369011799410029e+02,1.512964601769910e+02;...
    4.682286135693216e+02,1.352492625368730e+02;...
    4.710604719764012e+02,4.118274336283185e+02;...
    1.425648967551622e+02,3.957802359882005e+02];

% processing image
IB = (BR - I);
IBinv = imcomplement(IB);
IBinv(block)=0;
IBcbw = ~(adaptivethreshold(IBinv,ws,C,0));

% properties analysis
props = regionprops(IBcbw, 'Eccentricity','Perimeter','Centroid','Area','Orientation');
propstable = struct2table(props);

% filtering by Perimeter and Eccentricity to get mosquitoes
fit = model.predictFcn(propstable);
count = sum(strcmp(fit,'landed'))+sum(strcmp(fit,'collision'))*2;
f = strcmp(fit,'landed') | strcmp(fit,'collision');
c = strcmp(fit,'collision');

if count == 0
    pos = [];
    pelt_count = 0;
else
    % obtaining positions of all mosquitoes
    moz = props(f);
    pos = zeros(sum(f),2);
    centroid = [moz.Centroid];
    pos(:,1) = centroid(1:2:sum(f)*2);
    pos(:,2) = centroid(2:2:sum(f)*2);

    % counting points that are inside the Peltier
    p = inpolygon(pos(:,1),pos(:,2),Peltier(:,1),Peltier(:,2));
    pelt_count = sum(p);
    
    if sum(c) > 0
        centroid = [props(c).Centroid];
        colpos(:,1) = centroid(1:2:sum(c)*2);
        colpos(:,2) = centroid(2:2:sum(c)*2);
        p = inpolygon(colpos(:,1),colpos(:,2),Peltier(:,1),Peltier(:,2));
        pelt_count = pelt_count+p;
    end
end

%debugging: seeing image
% figure, hold on;
% imshowpair(I,IBcbw);
% %h = impoly(gca,Peltier);
% %setColor(h,'r');
% %n = 0;
% for k=1:length(props)
%     c = 'c'; % not a mosquito
%     if f(k)
%         c = 'w'; % mosquito
%         %n = n+1;
%         %if p(n)
%         %    c = 'r';
%         %end
%     end
%     text(props(k).Centroid(1),props(k).Centroid(2),num2str(k),'Color',c);
% end
% hold off;

end