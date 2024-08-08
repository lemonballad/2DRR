%UPTOWN_FUNK   uptown_funk
%  uptown_funk is the wicked fast (x,y,z) PDB coordinate grabber
%  
%  [X,Y,Z]=uptown_funk(file,natom,ntime)    takes the path and filename to
%  be read as file, the number atoms per model natom, and the number of
%  models/times to be read and returns the matrices X, Y, and Z
%  size(NAtoms,NTime) where NAtoms and NTime are actual dimension read in
%  from the file.
%  
% FUCK YEAH MOFOS!

function [X,Y,Z]=uptown_funk(file,natom,ntime)
formatspec='%*s %*d %*s %*s %*d %f %f %f %*f %*f %*s';
fid=fopen(file,'r');
it=0;
X=zeros(natom,ntime,'double');
Y=zeros(natom,ntime,'double');
Z=zeros(natom,ntime,'double');
while ~feof(fid)
    data=[];
    while isempty(data)&&~feof(fid)
        data=fscanf(fid,formatspec,[3 inf]);
    end
    nd2=size(data,2);
    if ~feof(fid)&&nd2>1
        it=it+1;
        X(1:nd2,it)=data(1,:);
        Y(1:nd2,it)=data(2,:);
        Z(1:nd2,it)=data(3,:);
    end
end
fclose(fid);
X(:,it+1:end)=[];
Y(:,it+1:end)=[];
Z(:,it+1:end)=[];
end