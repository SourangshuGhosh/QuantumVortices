clear;clc;clf;

N = 128;

for t = 0:1:10;

fid=fopen(['data/wf/wf',num2str(t),'p',num2str,'.dat'],'r');
dum = fread(fid,1,'float32');
wf = fread(fid,[N/4 N],'float64');
dum = fread(fid,1,'float32');
fclose(fid);
