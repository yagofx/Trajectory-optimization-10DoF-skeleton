    function CF=ComputationContactForces(pen,dpen,dtan)
%pen is penetration
%dpen is the velocity of penetration (positive pointing downwards)
%dtan is the tangential velocity, positive to the right
%CF(1) is the tangential contact forces, CF(2) is the normal contact force
import casadi.*
if isnumeric(pen)
    CF=zeros(1,2);
else
    CF=MX.zeros(1,2);
end

stiffness=1e7;
k = (1./2)*(stiffness.^(2/3));
contactSphereRadius=0.01;
cf=1e-8;
bd=1000;
bv=1000;
fh_pos = (4./3)*k*sqrt(contactSphereRadius*k)*(sqrt(pen*pen+cf).^(3/2));
fh_smooth=fh_pos*(1./2.+(1./2.)*tanh(bd*pen));
dissipation=0.8;
c = dissipation;
vnormal=dpen;
fhc_pos = fh_smooth*(1.+(3./2)*c*vnormal);
fhc_smooth = fhc_pos*(1/2+(1/2)*tanh(bv*(vnormal+(2/(3*c)))));
CF(2)=fhc_smooth;

% Calculate the friction force.
vtangent=dtan;
aux = vtangent+cf;
vslip = sqrt(aux.^2);
vt=0.1;
vrel = vslip / vt;
us=0.82;
ud=0.75;
uv=0.8;
ff = fhc_smooth*(min(vrel,1)*(ud+2*(us-ud)/(1+vrel*vrel))+uv*vslip);
force_tan = -ff*(vtangent) / vslip;

CF(1)=force_tan;


end