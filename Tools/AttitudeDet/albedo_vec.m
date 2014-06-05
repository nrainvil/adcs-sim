function [ albedo_v] = albedo_vec( a, sat_v )
%ALBEDO_VEC 
CONST.EMR = 6371.01e3; % m 

k = 1;
[sy sx] = size(a);
for i=1:sy
	j=1:sx;
	[grid_theta grid_phi] = idx2rad(i,j,sy,sx);
	[grid(1,:) grid(2,:) grid(3,:)] = sph2cart(grid_theta,pi/2-grid_phi,CONST.EMR);
	satpoint_v = sat_v*ones(1,length(grid(1,:)))-grid;
	satdist_v = sqrt(sum(satpoint_v.^2));
	satpoint_norm_v = satpoint_v./(ones(3,1)*satdist_v);
	if ~(exist('albedo_m','var'))
		albedo_v = (ones(3,1)*a(i,:)).*satpoint_norm_v;
	else
		albedo_v = [albedo_m ,(ones(3,1)*a(i,:)).*satpoint_norm_v];
	end

%    	for j=1:sx
%		[grid_theta grid_phi] = idx2rad(i,j,sy,sx);
%		[grid(1) grid(2) grid(3)] = sph2cart(grid_theta,pi/2-grid_phi,CONST.EMR);
%		satdist = norm(sat_v-grid');
%		albedo_v(:,k) = a(i,j)*(sat_v-grid')/satdist;
%		k = k+1;
%	end
end
end
