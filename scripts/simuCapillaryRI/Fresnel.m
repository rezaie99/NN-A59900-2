function [Thetat, Transmission] = Fresnel(ni,nt,Thetai)
if abs(ni*sin(Thetai)/nt)>=1
    Thetat = pi-Thetai;
    Transmission = 0;
else
    Thetat=asin(ni*sin(Thetai)/nt);
    R_s = ((ni*cos(Thetai)-nt*cos(Thetat))/(ni*cos(Thetai)+nt*cos(Thetat)))^2;
    R_p = ((ni*cos(Thetat)-nt*cos(Thetai))/(ni*cos(Thetat)+nt*cos(Thetai)))^2;
    R=0.5*(R_s+R_p);
    Transmission = 1-R;
end
if Transmission<0
    error('Ouch')
end
end