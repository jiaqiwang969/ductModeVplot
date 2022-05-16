  %Example 1, Basic Surf Point Detection
  %close all
  % Load image
    I=imread('C:\Users\wjq\Desktop\OpenSURF_version1c\TestImages\新建文件夹\Test1.png');
  % Set this option to true if you want to see more information
    Options.verbose=false; 
  % Get the Key Points
    Ipts=OpenSurf(I,Options);
  % Draw points on the image
    PaintSURF(I, Ipts);
 