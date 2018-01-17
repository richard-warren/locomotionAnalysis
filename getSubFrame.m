function [subFrame, isPadded] = getSubFrame(frame, yx, subFrameSize)

% !!! need to document
% given a frame, pulls out a subframe of size subFrameSize (h, w), centered at coordinates yx
% if the subframe is not fully contained within frame, frame is padded with zeros and flag isPadded is set to true


% get inds for top left corner of image
y1 = yx(1) - floor(subFrameSize(1)/2);
x1 = yx(2) - floor(subFrameSize(2)/2);


% get ind vectors
yInds = y1:y1+subFrameSize(1)-1;
xInds = x1:x1+subFrameSize(1)-1;


try
% if image falls outside of frame, pad with zeros and then get subframe
if any(yInds<=0) || any(yInds>size(frame,1)) || any(xInds<=0) || any(xInds>size(frame,2))
    
    frame = [zeros(subFrameSize(1), size(frame,2)); frame; zeros(subFrameSize(1), size(frame,2))]; % pad top and bottom
    frame = [zeros(size(frame,1), subFrameSize(2)), frame, zeros(size(frame,1), subFrameSize(2))]; % pad sides
    subFrame = frame(yInds+subFrameSize(1), xInds+subFrameSize(2));
    isPadded = true;
    
% otherwise just take the friggin subframe
else
    subFrame = frame(yInds, xInds);
    isPadded = false;
end
catch; keyboard; end