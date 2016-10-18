function rhs=a_rhs(t,a,dummy,phi,Lr)

rhs= Lr*a + i*phi'*( (abs(phi*a).^2).*(phi*a) );