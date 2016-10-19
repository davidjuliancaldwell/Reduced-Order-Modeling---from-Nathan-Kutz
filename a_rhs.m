function rhs=a_rhs(t,a,dummy,phi,Lr)

rhs= Lr*a + i*phi'*( (abs(phi*a).^2).*(phi*a) ); % here you have to do the complex conjugate transpose - dont put dot!! if you do, messes it all up 