function [ b_v ] = rot_srt( a_v, theta, lambda )
%SRT a = original vector, lambda = rotation axis unit vector
% theta = rotation angle (degrees)

   C11 = cosd(theta)+lambda(1,:).^2*(1-cosd(theta));
   C12 = -lambda(3,:)*sind(theta) + lambda(1,:)*lambda(2,:)*(1-cosd(theta));
   C13 = lambda(2,:)*sind(theta) + lambda(3,:)*lambda(1,:)*(1-cosd(theta));
   C21 = lambda(3,:)*sind(theta) + lambda(1,:)*lambda(2,:)*(1-cosd(theta));
   C22 = cosd(theta) + lambda(2,:).^2*(1-cosd(theta));
   C23 = -lambda(1,:)*sind(theta) + lambda(2,:)*lambda(3,:)*(1-cosd(theta));
   C31 = -lambda(2,:)*sind(theta) + lambda(3,:)*lambda(1,:)*(1-cosd(theta));
   C32 = lambda(1,:)*sind(theta) + lambda(2,:)*lambda(3,:)*(1-cosd(theta));
   C33 = cosd(theta) + lambda(3,:).^2*(1-cosd(theta));
   
   R = [C11,C12,C13;C21,C22,C23;C31,C32,C33];
   b_v = R*a_v;
   
end

