function ra_e = ForwardKin(l1,l2,th1,th2)
%FORWARDKIN
%    RA_E = FORWARDKIN(L1,L2,TH1,TH2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    08-Apr-2021 14:50:34

t2 = th1+th2;
ra_e = [-l2.*sin(t2)-l1.*sin(th1);l2.*cos(t2)+l1.*cos(th1);0.0];
