%Plot atomic structure
function Aplot(title,data)

% Parameters
Atom_res = 20; % resolution for your spheres
Atom_rad = .35; % radius of spheres
Atom_bond = 3.2; % length of bond

ps = data;
% Load the data
%load('ps.mat');

% Make base sphere data
[xb, yb, zb] = sphere(Atom_res);

figure('color','w','name',title);
hold on;
grid on;

% Plot spheres
for i=1:1:size(ps,1)
    surf(Atom_rad*xb+ps(i,1), Atom_rad*yb+ps(i,2), Atom_rad*zb+ps(i,3), 'facecolor', 'b', 'edgealpha', 0);
end

for i=1:3
    xleng(1,i) = max(ps(:,i));
    xleng(2,i) = min(ps(:,i));
end

ax = gca;
ax.DataAspectRatio = [1 1 1];
ax.XLabel.String = 'x';
ax.YLabel.String = 'y';
ax.ZLabel.String = 'z';
axis([xleng(2,1)-5 xleng(1,1)+5 xleng(2,2)-5 xleng(1,2)+5 xleng(2,3)-5 xleng(1,3)+5]);

% Make sure they're smooth and shaded
light;
lighting gouraud;

% Plot bonds
for i=1:size(ps,1)
    for j=1:size(ps,1)
        bl = norm(ps(i,1:3) - ps(j,1:3));
        if i~= j & bl < Atom_bond;
           [xc,yc,zc]=cylinder2P([0.08,0.08],10,ps(i,1:3),ps(j,1:3));
           surf(xc,yc,zc,'facecolor', 'b', 'edgealpha', 0);
        end
    end
end

hold off