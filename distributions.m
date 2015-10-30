% Script written by Diego Ramirez in collaboration with Daniela Garcia
% MPI Martisnried, 2015

clear all;
close all;
names=dir('*.txt');

vel=zeros(length(names),1);
diam=zeros(length(names),1);
%p='dist_velocities_FtsZ-MTS.txt';
 
for j=1:1:length(names)

     a=importdata(names(j).name);
     diam(j)=a.data(1,3);
     vel(j)=a.data(1,4);


      
end

m_d=mean(diam);
s_d=std(diam);
m_v=mean(vel);
s_v=std(vel);

bins1=5;
bins2=5;
figure
hist(diam,bins1);
xlabel('\mum','FontSize',24);
ylabel('Frequency','FontSize',24);
axis([0.5,1.5,0.1,50]);
set(gca,'fontsize',20)
legend({['N=' num2str(length(names)) '  Mean=' num2str(m_d) '  Std=' num2str(s_d)]},'FontSize',14,'Location','north');

figure
hist(vel,bins2);
xlabel('\mum/min','FontSize',24);
ylabel('Frequency','FontSize',24);
axis([1.5,4.5,0.1,60]);
set(gca,'fontsize',20)
legend({['N=' num2str(length(names)) '  Mean=' num2str(m_v) '  Std=' num2str(s_v)]},'FontSize',14,'Location','north');

%p='dist_vel_FtsZ-mts.txt';

% fileID = fopen(p,'w');
% fprintf(fileID,'%6f\n',data);
% fclose(fileID);

