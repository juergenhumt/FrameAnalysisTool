%
% Plot a biarc:
%
%  x0, y0 = initial point for the first circle arc 
%  theta0 = initial angle for the first circle arc
%  l0     = length of the first circle arc
%  kappa0 = curvature of the first circle arc
%  x1, y1 = final point for the second circle arc 
%  theta1 = final angle for the second circle arc
%  l1     = length of the second circle arc
%  kappa1 = curvature of the second circle arc
%
%  fmt1   = cell array with the format commmand 
%           for the first cicle arc
%  fmt2   = cell array with the format commmand 
%           for the second cicle arc
%  for example
%    fmt1 = {'Color','blue','Linewidth',3};
%    fmt2 = {'Color','red','Linewidth',3};
%
function biarc_plot(x0,y0,l0,theta0,kappa0,...
                    x1,y1,l1,theta1,kappa1,...
                    fmt1,fmt2)
  %
  ell = 0:l0/100:l0;
  tmp = (kappa0/2)*ell;
  S   = Sinc(tmp);
  x   = x0 + ell.*S.*cos(theta0+tmp);
  y   = y0 + ell.*S.*sin(theta0+tmp);  
  plot(x,y,fmt1{:});
  hold on;
  ell = 0:l1/100:l1;
  tmp = (kappa1/2)*ell;
  S   = Sinc(tmp);
  x   = x1 - ell.*S.*cos(theta1-tmp);
  y   = y1 - ell.*S.*sin(theta1-tmp);  
  plot(x,y,fmt2{:});
end

function r = Sinc( x )
  r = sinc(x/pi);
end
