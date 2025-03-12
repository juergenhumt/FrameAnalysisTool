%
% Test various randomly generated biarc
%
close all force

jX = -1
if jX > 0
    nn   = 2;
    npts = nn*nn;

    ra     = 0.8*pi;
    x0     = rand(1,npts);
    y0     = rand(1,npts);
    theta0 = rand(1,npts)*2*ra-ra;
    x1     = rand(1,npts)+3;
    y1     = rand(1,npts);
    theta1 = rand(1,npts)*2*ra-ra;

    fmt1 = {'Color','blue','Linewidth',3};
    fmt2 = {'Color','red','Linewidth',3};
    for k=1:npts
      xx0 = x0(k);
      yy0 = y0(k);
      th0 = theta0(k);
      xx1 = x1(k);
      yy1 = y1(k);
      th1 = theta1(k);
      tx  = mod(k-1,nn)*3;
      ty  = floor((k-1)/nn)*3;
      xx0 = xx0+tx;
      yy0 = yy0+ty;
      xx1 = xx1+tx;
      yy1 = yy1+ty;
      [l0,th0,k0,l1,th1,k1,xs,ys,thetas] = biarc(xx0,yy0,th0,xx1,yy1,th1);
      plot(xx0,yy0,'o');
      hold on;
      plot(xx1,yy1,'o');
      biarc_plot(xx0,yy0,l0,th0,k0,...
                 xx1,yy1,l1,th1,k1,...
                 fmt1,fmt2);
    end

else 
    
   
 xKar = [0.0, -425, -970, -1355, -1730, -1810];
 zKar = [1240, 685, 420, 338, 140, -400];
 thAr = [pi, deg2rad(211), deg2rad(211-16.5), deg2rad(211-14.5-2.5),  deg2rad(226), 1.03*0.5*pi*3]; 
 
 xKar = [0.0, -180, -445, -970, -1355, -1700, -1810];
 zKar = [1240, 1030,  605,  400,   295,   90, -400];
 
 thAr = [pi, 1.5*pi, deg2rad(211), deg2rad(211-16.5), deg2rad(211-14.5-0.5),  deg2rad(226), 1.03*0.5*pi*3];
 thN = thAr(1);
 
 j1=1; j2=2
 
 for k=1:6

   xN = xKar(j1);
   zN = zKar(j1);     
     
   x1 = xKar(j2);
   z1 = xKar(j2);
     
   x1 = xKar(j2);
   z1 = zKar(j2);
   th1= thAr(k+1);  ;

   k=0;

   xx0=xN; yy0=zN; th0= thN;
   xx1=x1; yy1=z1; th1= th1;

   [l0,th0,k0,l1,th1,k1,xs,ys,thetas] = biarc(xx0,yy0,th0,xx1,yy1,th1);
   plot(xx0,yy0,'o');
   axis equal   
   hold on;
   
   plot(xx1,yy1,'o');
   biarc_plot(xx0,yy0,l0,th0,k0,...
       xx1,yy1,l1,th1,k1,...
       fmt1,fmt2);
   
   thN = th1; 
   j1=j1 + 1; j2=j2+1;
     
 end            
             
             
end    


grid