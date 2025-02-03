function [impsd,corevol,samplevol] = Psd3D(t1b,graphicscale)
%An improved algorithm is used to measure the PSD distribution of porous geomaterials
%Input：imURL:The absolute path URL or variable name of the binary 3D pore-matrix image to be measured 
%             (white corresponds to pore pixels, black corresponds to skeleton pixels);
%      imagescale：Image pixel size, unit: micron-meter;；
% Output：impsd:PSD measurement results, where the first column radius, 
%              the second column area, the third column proportion.
%      coreAREA：Total pore area（Square micron-meter）
%      samplearea：Total area of section（Square micron-meter）
% Before running, you need to set the folder path where the code is stored as the search path of the software.
%%
load 3DTEst.mat;
t1b=~logical(t1b);%The binary image to be measured is introduced so that 1 is the matrix and 0 black is the pore
sz=size(t1b);
labelmap=zeros(sz);
corevol=numel(find(t1b==0));
samplevol=numel(t1b)*graphicscale^3;% Sample volume calculation (number of image pixels and converted)

for i=1:8000%The number of fill attempts can be varied

    Dtb1=bwdist(t1b);%

    if box(end).count>5 
        break
    end
    [x,y,z]=ind2sub(sz,find(Dtb1==max(Dtb1,[],'all'),1));
    % This loop is used to keep trying to reduce the size in the original position until it is filled in, updating the image if successful
    for j=floor(max(Dtb1,[],'all')):-1:1
        xlu=max(1,x-j);
        xld=min(x+j,size(t1b,1));
        yll=max(1,y-j);
        ylr=min(y+j,size(t1b,2));
        zls=max(1,z-j);
        zlx=min(z+j,size(t1b,3));
        block=t1b(xlu:xld,yll:ylr,zls:zlx)+box(j).cell(1+j-x+xlu:1+j+xld-x,1+j-y+yll:1+j-y+ylr,1+j-z+zls:1+j-z+zlx);
        if isempty(find(block==2,1))%If 2 occurs, interference between the filling circle and the skeleton occurs, otherwise no interference occurs
            t1b(xlu:xld,yll:ylr,zls:zlx)=block;% If this fill does not interfere, update the image
            box(j).count=box(j).count+sum(box(j).cell(1+j-x+xlu:1+j+xld-x,1+j-y+yll:1+j-y+ylr,1+j-z+zls:1+j-z+zlx),'all');
          
            break %Get out of the current cycle of trying to downsize
        end

        if j==1 && numel(find(t1b_sub==0))>=2
            t1b_sub=t1b(xlu:xld,yll:ylr,zls:zlx);
            t1b_sub(2,2,2)=1;
            t1b_sub(find(t1b_sub==0,2))=1;
            t1b(xlu:xld,yll:ylr,zls:zlx)=t1b_sub;%If left alone, update the image
            box(1).count=box(1).count+3;
        elseif j==1
            t1b(x,y,z)=1;%If the radius is 1 and still cannot be filled, then fill with the radius 0.5
            labelmap(x,y,z)=0.5;
            box(end).count=box(end).count+1;
        end
    end
end

BW2=false(size(t1b));
BW3=false(size(t1b));
CC =bwconncomp(~t1b,26);% 
stats = regionprops3(CC,'PrincipalAxisLength',"Volume");%
small_components = stats .Volume <= 3;% 
short_components = stats .PrincipalAxisLength <= 3;% 
small_indices = ismember(labelmatrix(CC), find(small_components));% 
short_indices = ismember(labelmatrix(CC), find(short_components));% 
BW2(small_indices) = true;% 
BW3(short_indices) = true;% 

labelmap(BW2==1 | BW3==1)=0.5;%
t1b(BW2==1|BW3==1)=1;
box(end).count=box(end).count+sum(logical(BW2+BW3),'all');

while (sum(t1b,'all')/numel(t1b)*100)<=97 && box(end).count<=3000 %Start with the rest of the radius 1 try filling (until 97% filled)
    [x,y,z]=ind2sub(sz,find(t1b==0,1));
    xlu=max(1,x-1);
    xld=min(x+1,size(t1b,1));
    yll=max(1,y-1);
    ylr=min(y+1,size(t1b,2));
    zls=max(1,z-j);
    zlx=min(z+j,size(t1b,3));
    t1b_sub=t1b(xlu:xld,yll:ylr,zls:zlx);
    t1b_sub(2,2,2)=1;
    if numel(find(t1b_sub==0))>=2
        t1b_sub(find(t1b_sub==0,2))=1;
        t1b(xlu:xld,yll:ylr,zls:zlx)=t1b_sub;%If this fill does not interfere, update the image
        box(1).count=box(1).count+3;
    else
        t1b(xlu:xld,yll:ylr,zls:zlx)=t1b_sub;%If the radius is 1 and still cannot be filled, then fill with the radius 0.5
        labelmap(x,y,z)=0.5;
        box(end).count=box(end).count+1;
    end
end
% Fill the rest with radius 0.5
box(end).count=box(end).count+numel(find(t1b==0));
labelmap(t1b==0)=0.5;
%% Results summary and unit conversion
impsd=zeros(length(box),3);
impsd(1,:)=[0.5*graphicscale,numel(find(labelmap==0.5))*graphicscale^3,numel(find(labelmap==0.5))*100/corevol];%
for k=2:length(impsd)
    impsd(k,:)=[(k-1)*graphicscale,box(k-1).count*graphicscale^3,box(k-1).count*100/corevol];
end
impsd(impsd(:,3)==0,:)=[];
end