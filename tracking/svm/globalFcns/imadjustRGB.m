function imadjustFrame = imadjustRGB(frame)
    imadjustFrame = uint8(nan(size(frame)));
    
    for i=1:size(frame,3)
        imadjustFrame(:,:,i) = imadjust(frame(:,:,i));
    end
end

