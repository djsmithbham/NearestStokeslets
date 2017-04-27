function [X1r,X2r,X3r]=RotateMesh(X1,X2,X3,X0,theta,axx)

 % X1 X2 X3 is a mesh for visualiation
 % X0 is a 3-vector (centre of rotation)
 % theta is an angle
 % axx is a 3-vector
 
  R=RotationMatrix(theta,axx);
  
  X1r=R(1,1)*X1+R(1,2)*X2+R(1,3)*X3;
  X2r=R(2,1)*X1+R(2,2)*X2+R(2,3)*X3;
  X3r=R(3,1)*X1+R(3,2)*X2+R(3,3)*X3;
  
