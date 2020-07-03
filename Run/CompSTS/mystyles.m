  function [linestyles,MarkerEdgeColors,Markers,mycolor]=mystyles(n)
	linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',':','-',':',...
    '-.','--','-',':','-.','--','-',':','-.'));
     
    %MarkerEdgeColors=jet(n);  % n is the number of different items you have
    MarkerEdgeColors=hsv(n);  % n is the number of different items you have
    Markers=['o','x','+','*','s','d','v','^','<','>','p','h','.',...
    '+','*','o','x','^','<','h','.','>','p','s','d','v',...
    'o','x','+','*','s','d','v','^','<','>','p','h','.'];
    mycolor=jet(n); 
    % [...]
     
%    figure
%    hold on
%    for i=1:n
%      plot(X(i,:), Y(i,:),[linestyles{i} Markers(i)],'Color',MarkerEdgeColors(i,:));
%    end
