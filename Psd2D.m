function [impsd,corearea,samplearea] =Psd2D(t1b,graphicscale)
%An improved algorithm is used to measure the PSD distribution of porous geomaterials
%Input：imURL:The absolute path URL or variable name of the binary two-dimensional pore-matrix image to be measured 
%             (white corresponds to pore pixels, black corresponds to skeleton pixels);
%      imagescale：Image pixel size, unit: micron-meter;；
% Output：impsd:PSD measurement results, where the first column radius, 
%              the second column area, the third column proportion.
%      coreAREA：Total pore area（Square micron-meter）
%      samplearea：Total area of section（Square micron-meter）
% Before running, you need to set the folder path where the code is stored as the search path of the software.
load 2DTEst.mat;%
t1b=~logical(t1b);%The binary image to be measured is introduced so that 1 is the matrix and 0 black is the pore
labelmap=zeros(size(t1b,1),size(t1b,2));

corearea=numel(find(t1b==0));%Pore pixel count
samplearea=numel(t1b)*graphicscale^2;%Calculate the cross-sectional area of the sample (according to the number of image pixels)

for i=1:8000%The number of fill attempts can be varied
    Dtb1=bwdist(t1b);
    if box(end).area>1 
        break
    end
    [x,y]=find(Dtb1==max(Dtb1,[],'all'),1);
    for j=ceil(max(Dtb1,[],'all')):-1:1
        xlu=max(1,x-j);
        xld=min(x+j,size(t1b,1));
        yll=max(1,y-j);
        ylr=min(y+j,size(t1b,2));
        block=t1b(xlu:xld,yll:ylr)+box(j).cell(1+j-x+xlu:1+j+xld-x,1+j-y+yll:1+j-y+ylr);%Fill area + original image
        if isempty(find(block==2,1))%If 2 appears, it means that the filled circle interferes with the skeleton pixel, otherwise no interference occurs
            t1b(xlu:xld,yll:ylr)=block;% If this fill does not interfere, update the image
            box(j).area=box(j).area+sum(box(j).cell(1+j-x+xlu:1+j+xld-x,1+j-y+yll:1+j-y+ylr),'all');%The dimension of the filling circle
            labelmap(xlu:xld,yll:ylr)=labelmap(xlu:xld,yll:ylr)+box(j).cell(1+j-x+xlu:1+j+xld-x,1+j-y+yll:1+j-y+ylr)*j;
            break % Break out of the current cycle of trying to reduce the size
        end
        if j==1
            t1b(x,y)=1;%If the radius is 1 and still cannot be filled, then fill with the radius 0.5
            labelmap(x,y)=0.5;
            box(end).area=box(end).area+1;
        end
    end
end
BW2=bwpropfilt(~t1b,'Area',[1 2],4);
BW3 = bwpropfilt(~t1b,'MinorAxisLength',[0 2],4);
labelmap(BW2==1|BW3==1)=0.5;%The pores that meet these conditions are considered to have a radius of 0.5
t1b(BW2==1|BW3==1)=1;
box(end).area=box(end).area+sum(logical(BW2+BW3),'all');
for i=1:numel(find(t1b==0))  %The rest start by pressing radius 1 to try filling
    [x,y]=find(t1b==0,1);
    xlu=max(1,x-1);
    xld=min(x+1,size(t1b,1));
    yll=max(1,y-1);
    ylr=min(y+1,size(t1b,2));
    block=t1b(xlu:xld,yll:ylr)+box(1).cell(1+1-x+xlu:1+1+xld-x,1+1-y+yll:1+1-y+ylr);%The binary image is updated according to the filling result
    if isempty(find(block==2,1))% If the filling is successful, the radius is 1
        t1b(xlu:xld,yll:ylr)=block;% If not, update the image
        box(1).area=box(1).area+sum(box(1).cell(1+1-x+xlu:1+1+xld-x,1+1-y+yll:1+1-y+ylr),'all');
        labelmap(xlu:xld,yll:ylr)=labelmap(xlu:xld,yll:ylr)+box(1).cell(1+1-x+xlu:1+1+xld-x,1+1-y+yll:1+1-y+ylr);
    else
        t1b(x,y)=1;% If the radius is 1 and still cannot be filled, then fill with the radius 0.5
        labelmap(x,y)=0.5;
        box(end).area=box(end).area+1;
    end
end


%% Results summary and unit conversion
impsd=zeros(length(box),3);
impsd(1,:)=[0.5*graphicscale,numel(find(labelmap==0.5))*graphicscale^2,numel(find(labelmap==0.5))*100/corearea];
for k=2:length(impsd)
    impsd(k,:)=[(k-1)*graphicscale,box(k-1).area*graphicscale^2,box(k-1).area*100/corearea];
end
impsd(impsd(:,3)==0,:)=[];
end