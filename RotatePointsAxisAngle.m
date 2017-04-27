function X=RotatePointsAxisAngle(X,rotvec)

% takes general axis vector and angle as 4-vector rotvec, generates rotation matrix and performs rotation on grid vector X
% requires simulink for vrrotvec2mat

  [X1 X2 X3]=ExtractComponents(X);
  R=vrrotvec2mat(rotvec);
  X=R*[(X1)';(X2)';(X3)'];
  X=reshape(X',[],1);
