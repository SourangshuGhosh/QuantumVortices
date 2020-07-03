clear;clc;

N=128;

a=load('Pot1p1.dat');
b=load('Pot1p2.dat');
c=load('Pot1p3.dat');
d=load('Pot1p4.dat');
a1=reshape(a,N,N/4);
b1=reshape(b,N,N/4);
c1=reshape(c,N,N/4);
d1=reshape(d,N,N/4);
V = a1;
V(1:N,N/4+1:N/2)=b1;
V(1:N,N/2+1:3*N/4)=c1;
V(1:N,3*N/4+1:N)=d1;

V = transpose(V);
