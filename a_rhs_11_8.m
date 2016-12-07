function rhs=a_rhs_11_8(t,a,dummy,phi,Lr,PNL)

rhs= Lr*a + i*PNL'*( (abs(PNL*a).^2).*(PNL*a) ); % here you have to do the complex conjugate transpose - dont put dot!! if you do, messes it all up 