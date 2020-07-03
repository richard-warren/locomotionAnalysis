function [segment,tc,wc] = pcaseg(data,num_segments,q,flag)
% The basic algorithm works by creating
% a fine segmented representation then merging 
% the lowest cost segments until only 'num_segments' remain.


minres=floor(length(data)/200); %minimal resolution
degree=1;
left_x       = [1 : minres : size(data,1)-1];    % Find the left x values vector for the "fine segmented representation".
right_x      = left_x + minres;                  % Find the right x values vector for the "fine segmented representation".
right_x(end) = size(data,1);                     % Special case, the rightmost endpoint.
number_of_segments = length(left_x );            % Calculate the number of segments in the initial "fine segmented representation".
% Initialize the segments in the "fine segmented representation". 
for i = 1 : number_of_segments 
   segment(i).lx = left_x(i);
   segment(i).rx = right_x(i);
   segment(i).mc = inf;
   segment(i).c = inf;
end;
tc=[];
wc=[];

% Initialize the merge cost of the segments in the "fine segmented representation". 
for i = 1 : number_of_segments -1
   sx=data(segment(i).lx :segment(i+1).rx,:);
   segment(i).mc = pcaresid(sx,q,flag);
   sx=data(segment(i).lx :segment(i).rx,:);
   segment(i).c = pcaresid(sx,q,flag);
end;
   sx=data(segment(i+1).lx :segment(i+1).rx,:);
   segment(i+1).c = pcaresid(sx,q,flag);

% Keep merging the lowest cost segments until only 'num_segments' remain. 
while  length(segment) > num_segments  
    
   [value, i ] = min([segment(:).mc]);                              % Find the location "i", of the cheapest merge.
   if i > 1 & i < length(segment) -1								% The typical case, neither of the two segments to be merged are end segments
           segment(i).c=segment(i).mc;
           sx=data(segment(i).lx :segment(i+2).rx,:);
           segment(i).mc = pcaresid(sx,q,flag);
    	   segment(i).rx = segment(i+1).rx;
           segment(i+1) = [];
            
           i = i - 1;
           sx=data(segment(i).lx :segment(i+1).rx,:);
           segment(i).mc = pcaresid(sx,q,flag);
           
           
   elseif i == 1                                                    % Special case: The leftmost segment must be merged.
           segment(i).c=segment(i).mc;
           sx=data(segment(i).lx :segment(i+2).rx,:);
           segment(i).mc = pcaresid(sx,q,flag);
	 	   segment(i).rx = segment(i+1).rx;
           segment(i+1) = [];
           sx=data(segment(i).lx :segment(i).rx,:);
           segment(i).c = pcaresid(sx,q,flag);

              
   else                                                             % Special case: The rightmost segment must be merged.
          segment(i).rx = segment(i+1).rx;
          segment(i).c=segment(i).mc;
          segment(i).mc = inf;
          segment(i+1) = [];
 
          i = i - 1;
       
          sx=data(segment(i).lx :segment(i+1).rx,:);
          segment(i).mc = pcaresid(sx,q,flag);


   end; % end of if i > 1 & i < length(segment) -1

% total cost 
   tc=[tc; sum([segment.c])];

% weighted cost
   for j=1:num_segments
        seglength(j) = segment(j).rx - segment(j).lx;
        norm_seglength(j) = seglength(j)/right_x(end);
        segcost(j) = segment(j).c;
   end
   wc = [wc; sum([norm_seglength .* segcost])];


end; % end of while  length(segment) > num_segments



for i = 1 : size(segment,2)
    sx=data(segment(i).lx :segment(i).rx,:);
    [segment(i).mc,segment(i).pc,segment(i).avg,segment(i).latent] = pcaresid(sx,q,flag);
end;

% Figure
%     plot(data)
%     hold on
%     b=axis;
%     b=b(3:4);
%     for i=1:size(segment,2)
%         line([segment(i).lx  segment(i).lx], b);
%     end

