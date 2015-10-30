function [] = KymographAnalysis(tseries, z_proj,n)

% Script written by Diego Ramirez in collaboration with Daniela Garcia
% MPI Martisnried, 2015

%proj_stack, z-projection of stack, 10 with average intensity is more than
%enough
%stack, adquisition treated with filter mean = 2 pixels

stack = tseries;
proj_stack = z_proj;

time=3;
off=1;

%DistancePix=37.940;
%TIRF-corefacility, 1/14.779
%Confocal, 1/37.940 Zoom10X
%Confocal, 1/22.764 Zoom6X

%Check number!!!

%pixel_size=1/DistancePix;

pixel_size=0.0423;

tsStack = tiffread(stack);


Im=imread(proj_stack);

%Open figure to select ring
figure
imshow(Im);
[y0,x0,p0] = impixel;
close;



for j=1:1:length(x0)/2
    

    
I=imread(proj_stack);

r=sqrt((x0(2*j-1)-x0(2*j))^2+(y0(2*j-1)-y0(2*j))^2)+2;

B=I(x0(2*j-1)-r:x0(2*j-1)+r, y0(2*j-1)-r:y0(2*j-1)+r);

%B=Im(x0(2*j-1)-(r+offset):x0(2*j-1)+(r+offset), y0(2*j-1)-(r+offset):y0(2*j-1)+(r+offset));

thresh = multithresh(B);

%imshow(B);

tam=size(B);

A1=B(:,1:round(tam(2)/2));
A2=B(:,round(tam(2)/2)+1:round(tam(2)));
tam1=size(A1);
% tam2=size(A2);

xy3=zeros(tam(1),2);
xy4=zeros(tam(1),2);
 
 
for i=1:1:tam(1)
     
    max_value1=max(A1(i,:));
    if  max_value1(1)>thresh
    ind1=find(A1(i,:)==max_value1);
    xy3(i,:)=[mean(ind1),i];
    end
    
    max_value2=max(A2(i,:));
    if  max_value1(1)>thresh
    ind2=find(A2(i,:)==max_value2);
    xy4(i,:)=[mean(ind2)+tam1(2),i];
    end 
             
     
end    
       
xy_fit = cat(1, xy3,xy4);
xy_fit(all(xy_fit==0,2),:)=[];
%scatter(xy_fit(:,1), xy_fit(:,2));
initial=[round((tam(1)/2)) round((tam(2)/2)) 0.8*round((tam(2)/2))];
out=LM(xy_fit,initial,1);

%viscircles([out(1) out(2)], out(3),'EdgeColor','b');
th = 0:pi/50:2*pi;
xunit_1 = round((out(3)) * cos(th) + out(1));
yunit_1 = round((out(3)) * sin(th) + out(2));
xunit_2 = round((out(3)-off) * cos(th) + out(1));
yunit_2 = round((out(3)-off) * sin(th) + out(2));
xunit_3 = round((out(3)+off) * cos(th) + out(1));
yunit_3 = round((out(3)+off) * sin(th) + out(2));
xunit_4 = round((r) * cos(th)) + y0(2*j-1);
yunit_4 = round((r) * sin(th)) + x0(2*j-1);
%Im=imread(proj_stack);

%Get unique pair of values between the vectors xunit_n, yunit_n
unit_1 = unique(table(xunit_1', yunit_1'), 'stable');
xunit_1 = table2array(unit_1(:,1));
yunit_1 = table2array(unit_1(:,2));

unit_2 = unique(table(xunit_2', yunit_2'), 'stable');
xunit_2 = table2array(unit_2(:,1));
yunit_2 = table2array(unit_2(:,2));

unit_3 = unique(table(xunit_3', yunit_3'), 'stable');
xunit_3 = table2array(unit_3(:,1));
yunit_3 = table2array(unit_3(:,2));

unit_4 = unique(table(xunit_4', yunit_4'), 'stable');
xunit_4 = table2array(unit_4(:,1));
yunit_4 = table2array(unit_4(:,2));


%Area selected plot
figure
subplot(1,2,1),
imshow(B);

hold on 
plot(xunit_1,yunit_1,xunit_2, yunit_2, xunit_3, yunit_3,'-');
hold off

subplot(1,2,2)
imshow(I);

hold on 
plot(xunit_4,yunit_4,'-');
hold off

%Save plot
s=num2str(j);
print(s,'-dpng');
close;

%Same size as unit_1 since xunit_1 and yunit_1 are being printed
diameter=zeros(1,length(xunit_1));
vel=zeros(1,length(xunit_1));
kymograp = zeros(size(tsStack,3),length(xunit_1));
kymograp2 = zeros(size(tsStack,3),length(xunit_2));
kymograp3 = zeros(size(tsStack,3),length(xunit_3));
kymograp_mean = zeros(size(tsStack,3),length(xunit_1));


for i = 1:size(tsStack,2)
   
A=tsStack(i).data;

%background = imopen(A,strel('disk',20));

%A = A - background;

%Im = mat2gray(A);
 try
 Im=A(x0(1)-r:x0(1)+r, y0(1)-r:y0(1)+r);
 catch exception  
 continue 
 end
     try  
    lin_ind = sub2ind(size(Im),yunit_1,xunit_1);
    kymograp(i,:)= Im(lin_ind);
    %kymograp(i,:)=(kymograp(i,:)-min(min(kymograp(i,:))))/(max(max(kymograp(i,:)))-min(min(kymograp(i,:))));

     catch exception  
     continue       
     end
     try
    lin_ind = sub2ind(size(Im),yunit_2,xunit_2);
    kymograp2(i,:)= Im(lin_ind);
    %kymograp2(i,:)=(kymograp2(i,:)-min(min(kymograp2(i,:))))/(max(max(kymograp2(i,:)))-min(min(kymograp2(i,:))));
     catch exception  
     continue       
     end
     try
    lin_ind = sub2ind(size(Im),yunit_3,xunit_3);
    kymograp3(i,:)= Im(lin_ind);
    %kymograp3(i,:)=(kymograp3(i,:)-min(min(kymograp3(i,:))))/(max(max(kymograp3(i,:)))-min(min(kymograp3(i,:))));
    catch exception
    end    
     
end
    try 
%spacer=mat2gray(ones(size(tsStack,2),length(th)));
kymograp=mat2gray(kymograp);
kymograp2=mat2gray(kymograp2);
kymograp2=imresize(kymograp2,'OutputSize',size(kymograp));
kymograp3=mat2gray(kymograp3);
kymograp3=imresize(kymograp3,'OutputSize',size(kymograp));
diameter(1)=2*out(3)*pixel_size;

kymograp_mean=(kymograp+kymograp2+kymograp3)/3;

kymograp_mean = imadjust(kymograp_mean);
kymograp = imadjust(kymograp);
kymograp2 = imadjust(kymograp2);
kymograp3 = imadjust(kymograp3);
l=size(tsStack,2)+7;
%Save kymograph
figure
subplot(1,3,1)
imshow(kymograp);

subplot(1,3,2)
imshow(kymograp2);
 
subplot(1,3,3)
imshow(kymograp3);

print(strcat(n,s,'kymo_mean'),'-dpng');
close;

figure
%imshow(kymograp);
imshow(kymograp_mean);
print(strcat(n,s,'kymo_means'),'-dpng');
%imshow([kymograp,spacer,kymograp2,spacer,kymograp3]);
p = impoly(gca);
pos = wait(p);

if (isempty(pos) == 0)
m= abs((pos(2,1) - pos(1,1)))/abs((pos(2,2) - pos(1,2)));

%Turn pix/seg to um/min
v = 60*pixel_size*m/(time);
vel(1) = v;

tx = text(50,l,sprintf('V = %d um/min', v));
set(tx, 'FontSize', 10);
print(strcat(n,s,'kymo_meanV'),'-dpng');
close;

T=table(xunit_1, yunit_1, diameter', vel');
writetable(T,strcat(n,s,'kymo_mean.txt'),'Delimiter','\t');

else
    
close;
T=table([xunit_1, yunit_1, diameter']);
writetable(T,strcat(n,s,'kymo_mean.txt'),'Delimiter','\t');
end

    catch exception  
    close all;
    continue       
    end

end

end
