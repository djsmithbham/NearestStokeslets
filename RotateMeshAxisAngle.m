function [X1r,X2r,X3r]=RotateMeshAxisAngle(X1,X2,X3,rotvec)

% takes general axis vector and angle as 4-vector rotvec, generates rotation matrix and performs rotation on grid vector X
% requires simulink for vrrotvec2mat

% works with a mesh X1 X2 X3 for visualisation purposes

  R=vrrotvec2mat(rotvec);
  
  X1r=R(1,1)*X1+R(1,2)*X2+R(1,3)*X3;
  X2r=R(2,1)*X1+R(2,2)*X2+R(2,3)*X3;
  X3r=R(3,1)*X1+R(3,2)*X2+R(3,3)*X3;
 
