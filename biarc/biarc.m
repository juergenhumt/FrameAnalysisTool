%
%  Compute the biarc passig from point (x0,y0) angle th0
%  to point (x1,y1) angle th1
%
%  On Input
%  x0, y0 = initial point for the first circle arc 
%  th0    = initial angle for the first circle arc
%  x1, y1 = final point for the second circle arc 
%  th1    = final angle for the second circle arc
%
%  On Output
%  l0     = length of the first circle arc
%  theta0 = initial angle for the first circle arc
%            (same as th0 but +2*k*pi if necessary)
%  kappa0 = curvature of the first circle arc
%
%  l1     = length of the second circle arc
%  theta1 = final angle for the second circle arc
%            (same as th1 but +2*k*pi if necessary)
%  kappa1 = curvature of the second circle arc
%
%  The solution algorithm is described in
%  A Note on Robust Biarc Computation
%  by Enrico Bertolazzi and Marco Frego
%  Computer-Aided Design & Applications, 16(5), 2019, 822-835
%  https://doi.org/10.14733/cadaps.2019.822-835
%
function [ l0, theta0, kappa0, l1, theta1, kappa1, xs, ys, thetas ] = biarc( x0, y0, th0, x1, y1, th1 )
  dx     = x1-x0;
  dy     = y1-y0;
  d      = hypot(dx,dy);
  omega  = atan2(dy,dx);
  theta0 = omega+Range(th0-omega);
  theta1 = omega+Range(th1-omega);
  dt     = (theta1-theta0)/2;
  t      = d * Sinc( dt/2 ) / Sinc( dt );
  thetas = 2*omega-(theta0+theta1)/2;
  dt0    = (thetas-theta0)/2;
  dt1    = (thetas-theta1)/2;
  l0     = t/( 2*Sinc( dt0 ));
  l1     = t/( 2*Sinc( dt1 ));
  %kappa0 = 2*dt0/l0;
  %kappa1 = -2*dt1/l1;
  kappa0 = 4*sin(dt0)/t;
  kappa1 = -4*sin(dt1)/t;
  xs     = x0+(t/2)*cos( (thetas+theta0)/2 );
  ys     = y0+(t/2)*sin( (thetas+theta0)/2 );
end

function r = Sinc( x )
  r = sinc(x/pi);
end

function theta = Range( theta )
  while theta > pi
    theta = theta - 2*pi;
  end
  while theta < -pi
    theta = theta + 2*pi;
  end
end
