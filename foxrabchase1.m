function foxrabchase1 
%A fox rabbit pursuit simulation based on the values and criteria of 
%question 1

sr = 13; %speed of rabbit
sf = 16; % "     " fox 
z0 = [250 -550]; %ic for fox
ts = [0 norm(z0)/(sf-sr)]; % timespan
mindist = 0.1;
burrow = [600 600];
options = odeset('Events', @(t,z)foxrab1(t,z,sr, mindist, burrow));
[t, z, te, ze, zi] = ode45(@(t,z)foxode1(t, z, sr, sf), ts, z0, options);
te, ze, zi
plot(z(:,1), z(:,2), sr*cos(pi/4)*t, sr*sin(pi/4)*t)

function dzdt = foxode1(t, z, s_r, s_f)
%definition of ODE for the Fox-Rabbit chase

r = [s_r*cos(pi/4)*t s_r*sin(pi/4)*t]; %the position of the rabbit

dist = sqrt((r(1)-z(1))^2 + (r(2)-z(2))^2);%distance between fox and rabbit
dzdt = zeros(2,1); %2x1 column vector
dzdt(1) = s_f*(r(1)-z(1))/dist; %horizantal velocity
dzdt(2) = s_f*(r(2)-z(2))/dist; %vertical velocity
end

function [value, isterminal, direction] = foxrab1(t, z, s_r,mindist,burrow)
r = [s_r*cos(pi/4)*t s_r*sin(pi/4)*t];
burrow = [600, 600]; 
value(1) = sqrt((r(1)-z(1))^2 + (r(2)-z(2))^2) - mindist; %fox catches 
                                                          %rabbit
isterminal(1) = 1;
direction(1) = -1;
value(2) = (burrow(1) - r(1) > 0 && burrow(2) - r(2) > 0); %if rabbit 
                                                           %reaches burrow
isterminal(2) = 1;
direction(2) = -1;
end
end
