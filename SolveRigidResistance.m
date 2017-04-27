function [FF,MM,f,condno]=SolveRigidResistance(x,X,X0,U,Om,ep,domain,blockSize)

  % Solves the rigid body resistance problem with nearest neighbour regularised stokeslets, force nodes = coll nodes
  % Input:  x,         collocation/traction points with all x1 coords first, then all x2...
  %         X,         stokeslet points
  %         X0,        origin of rotation
  %         U, Om,     vectors with rigid body velocity and angular velocity
  %         ep,       regularisation parameter
  %         domain,    'i' infinite fluid (regularised stokeslet) or 'h' halfspace x3 > 0 (regularised blakelet)
  %         blockSize, maximum memory storage for stokeslet matrix in GB
  % Output: FF,        vector of force on body
  %         MM,        vector of moment on body
  %         f,         stokeslet strengths 

  M=length(x)/3;
  Q=length(X)/3;

  % lhs assembly
  [AS, NN]=AssembleStokesletMatrix(x,X,x,ep,domain,blockSize);  
  % rhs assembly
  
  [x1, x2, x3]=ExtractComponents(x);ze=0*x1;
  u1 = U(1)*ones(M,1);
  u2 = U(2)*ones(M,1);
  u3 = U(3)*ones(M,1);

  x1=x1-X0(1);
  x2=x2-X0(2);
  x3=x3-X0(3);
  rhs=[u1;u2;u3]+[Om(2)*x3-Om(3)*x2;Om(3)*x1-Om(1)*x3;Om(1)*x2-Om(2)*x1];
  
  condno=cond(AS);
  
  f=AS\rhs;

  % calculate force
  [F1 F2 F3]=ExtractComponents(NN*f);
  FF=[sum(F1);sum(F2);sum(F3)];
  
  % calculate moment
  [X1 X2 X3]=ExtractComponents(X'*NN);Ze=0*X1;
  X1=X1-X0(1);
  X2=X2-X0(2);
  X3=X3-X0(3);
  AM=[Ze -X3 X2; X3 Ze -X1; -X2 X1 Ze];
  MM=AM*f;

