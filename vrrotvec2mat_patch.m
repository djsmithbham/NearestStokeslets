function R=vrrotvec2mat_patch(rotvec)

% replaces simulink 3D animation function vrrotvec2mat
% takes rotvec (3-vector axis and scalar angle)
% and produces equivalent rotation matrix

b1=rotvec(1:3)/norm(rotvec(1:3));

v=[1 0 0];
if norm(cross(b1,v))<1e-10
    v=[1 1 0];
end
b2=cross(b1,v);
b2=b2/norm(b2);
b3=cross(b1,b2);b3=b3/norm(b3);

A=[b1(:) b2(:) b3(:)];
theta=rotvec(4);
M=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
R=A*M*A';

