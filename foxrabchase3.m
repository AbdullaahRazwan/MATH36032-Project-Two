function foxrabchase3
sr0 = 13;
sf0 = 16;
rz0 = [0 0 0];
fz0 = [250 -550 0];
ts = [0 norm(fz0)/(sf0-sr0)];
mindist = 0.1;
burrow = [600 600];
options = odeset('Events', @(t,z)foxrab2(t,z,sr0, mindist, burrow));
[t, r] = ode45(@(t,r)rabode(t, r, sr0), ts, rz0);
[t, z, te, ze, zi] = ode45(@(t,z)foxode2(t, z, sf0), ts, fz0, options);
te, ze, zi
plot(r(:,1), r(:,2), z(:,1), z(:,2))
legend('Rabbit', 'Fox', 'Location', 'Best')

function drdt = rabode(t,r, s_r)
%definition of ODE for rabbit - calculate distance travlled by rabbit
    drdt = zeros(3,1);
    drdt(1) = s_r*cos(pi/4)*exp(-0.0008*r(3));
    drdt(2) = s_r*sin(pi/4)*exp(-0.0008*r(3));
    drdt(3) = sqrt((drdt(1))^2 + (drdt(2))^2);
end

function dzdt = foxode2(t, z, s_f)
%definition of ODE for the Fox
        dist = sqrt((r(1)-z(1))^2 + (r(2)-z(2))^2);%distance between fox and rabbit
        dzdt = zeros(3,1); %2x1 column vector
        dzdt(1) = s_f*exp(-0.0002*z(3))*(r(1)-z(1))/dist; %horizantal velocity
        dzdt(2) = s_f*exp(-0.0002*z(3))*(r(2)-z(2))/dist; %vertical velocity
        dzdt(3) = sqrt(dzdt(1)^2 + dzdt(2)^2);
end

function [value, isterminal, direction] = foxrab2(t, z, s_r, mindist, burrow)
        burrow = [600, 600]; 
        value(1) = sqrt((r(1)-z(1))^2 + (r(2)-z(2))^2) - mindist; %event to 
                              %check if fox has caught the rabbit
        isterminal(1) = 1; %stop solution
        direction(1) = -1;
        value(2) = (burrow(1) - r(1) > 0 && burrow(2) - r(2) > 0);
        isterminal(2) = 1;
        direction(2) = -1;
end
end
