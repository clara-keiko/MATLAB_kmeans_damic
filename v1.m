%input 
data_raw=load('out.csv');

%charge
charge = log(data_raw(:,1));
charge2 = data_raw(:,1);

%track lenght
len=data_raw(:,7);

%"de/dx" (not phically, just charge divided by lenght)
dedx=log(charge2./len);

%number of pixels 
npix = log(data_raw(:,2));
npix2 = data_raw(:,2);

%initial position
xi=data_raw(:,3);
xf=data_raw(:,5);
yi=data_raw(:,4);
yf=data_raw(:,6);

%distance between first and last pixel
s=sqrt(((xf-xi).^2) + ((yf-yi).^2));

%ocupation ratio
ocup=log(data_raw(:,8));

%center parameter (distance between most energetic pixel and the center of the image)
cp=log(data_raw(:,10));
cp2=data_raw(:,11);

%varition between track lenght and distance between first and last pixel
delta = log(abs(len-s));

%data
data = [charge,dedx,npix,ocup,cp,delta,s];

%number of clusters
[idx,C] = kmeans(data,6);
gplotmatrix(data,[],idx,[],'.',15);

fileID = fopen('grupo.txt','w');
fprintf(fileID,'%d\n',idx);
fclose(fileID);
