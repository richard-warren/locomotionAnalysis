function wheelMask = getWheelMask(circRoiPts, frameSize)
        
[wheelRadius, wheelCenter] = fitCircle(circRoiPts);
wheelMask = uint8(ones(frameSize));

for y=1:frameSize(1)
    for x=1:frameSize(2)
        if norm(wheelCenter-[x;y]) < wheelRadius
            wheelMask(y,x) = 0;
        end
    end
end
