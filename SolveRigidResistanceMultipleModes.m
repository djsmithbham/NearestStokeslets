function [FF,MM,f,condno]=SolveRigidResistanceMultipleModes(x,X,X0,U,Om,ep,domain,blockSize,condOpt)

  % Solves the rigid body resistance problem with nearest neighbour regularised stokeslets, force nodes = coll nodes
  %         optional calculation of condition number
  % Input:  x,         collocation/traction points with all x1 coords first, then all x2...
  %         X,         stokeslet points
  %         X0,        origin of rotation
  %         U, Om,     matrices with column vectors for each mode's rigid body velocity and angular velocity
  %         ep,       regularisation parameter
  %         domain,    'i' infinite fluid (regularised stokeslet) or 'h' halfspace x3 > 0 (regularised blakelet)
  %         blockSize, maximum memory storage for stokeslet matrix in GB
  %         condOpt    'y' calculate condition number, otherwise don't
  % Output: FF,        matrix; each column is vector of force on body for a given mode
  %         MM,        matrix; each column is vector of moment on body for a given mode
  %         f,         stokeslet strengths for each mode

  M=length(x)/3;
  Q=length(X)/3;

  % lhs assembly
  [AS, NN]=AssembleStokesletMatrix(x,X,x,ep,domain,blockSize);  

  switch condOpt
    case 'y'
  	condno=cond(AS);
  end

  % rhs assembly
  rhs = zeros(3*M,size(U,2));
  for j = 1:size(U,2) % loop over modes to set up right hand sides
    [x1, x2, x3]=ExtractComponents(x);ze=0*x1;
    x1c=x1-X0(1);
    x2c=x2-X0(2);
    x3c=x3-X0(3);
    u1 = U(1,j)*ones(M,1);
    u2 = U(2,j)*ones(M,1);
    u3 = U(3,j)*ones(M,1);
    rhs(:,j)=[u1;u2;u3]+[Om(2,j)*x3c-Om(3,j)*x2c;Om(3,j)*x1c-Om(1,j)*x3c;Om(1,j)*x2c-Om(2,j)*x1c];
  end

  f=AS\rhs;

  for j = 1:size(U,2)
    % calculate force
    [F1 F2 F3]=ExtractComponents(NN*f(:,j));
    FF(:,j)=[sum(F1);sum(F2);sum(F3)];
  
    % calculate moment
    [X1 X2 X3]=ExtractComponents(X'*NN);Ze=0*X1;
    X1C=X1-X0(1);
    X2C=X2-X0(2);
    X3C=X3-X0(3);
    AM=[Ze -X3C X2C; X3C Ze -X1C; -X2C X1C Ze];
    MM(:,j)=AM*f(:,j);
  end
