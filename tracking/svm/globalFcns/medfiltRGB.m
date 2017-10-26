function medFiltFrame = medfiltRGB(frame, kernel)
    % simply computes medfilt2 on each channel of an rbg (or other color
    % scheme) frame
    
    medFiltFrame = uint8(nan(size(frame)));
    
    for i=1:size(frame,3)
        medFiltFrame(:,:,i) = medfilt2((frame(:,:,i)), kernel);
    end
end

