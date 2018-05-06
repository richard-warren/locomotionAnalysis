session = '180122_000';
frame = 189567;

vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
frame = read(vid, frame);
imwrite(frame, 'C:\Users\rick\Desktop\github\DeepLabCut\matlab.png')
