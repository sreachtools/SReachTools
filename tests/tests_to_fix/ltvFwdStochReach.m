%% Script for sanity checking getConcatMat for LtvSystem --- TO BE DELETED
T = 0.25;
decay_rate = 2;
cov_dist = 0.1;
B_mat = @(t) exp(-decay_rate *t)*[T^2/2;
                                  T];
k_gain = 3;

 %sys = LtvSystem('StateMatrix', @(t) [1,t;
                                      %0,1],...
                 %'DisturbanceMatrix',B_mat, ...
                 %'Disturbance', RandomVector('Gaussian', 0, cov_dist));
sys = LtvSystem('StateMatrix', [1,T;
                                0,1],...
                'InputMatrix',[T^2/2;
                               T], ...
                 'InputSpace', Polyhedron('lb',-[1],'ub',[1]),...
                'DisturbanceMatrix',[T^2/2;
                                     T], ...
                'Disturbance', RandomVector('Gaussian', 0, cov_dist));
                           
lsys = LtiSystem('StateMatrix', [1,T;
                                 0,1],...
                 'InputMatrix',[T^2/2;
                                T], ...
                 'InputSpace', Polyhedron('lb',-[1],'ub',[1]),...
                 'DisturbanceMatrix',[T^2/2;
                                      T], ...
                 'Disturbance', RandomVector('Gaussian', 0, cov_dist));
[Z,H,G]=getConcatMats(sys,3);                 
[Zlti,Hlti,Glti]=getConcatMats(lsys,3);                 

if all(all(abs(Z-Zlti)<1e-8)) && all(all(abs(H-Hlti)<1e-8)) && all(all(abs(G-Glti)<1e-8))
    disp('Matching!')
else
    disp('Not matching!')
end

[Z,H,G]=sys.getConcatMats(10);                 
[Zlti,Hlti,Glti]=lsys.getConcatMats(10);                 
if all(all(abs(Z-Zlti)<1e-8)) && all(all(abs(H-Hlti)<1e-8)) && all(all(abs(G-Glti)<1e-8))
    disp('Matching!')
else
    disp('Not matching!')
end
sys = LtvSystem('StateMatrix', @(t) [1,t;
                                0,1],...
                'InputMatrix', @(t) [1+t^2/2;
                                     t], ...
                'InputSpace', Polyhedron('lb',-[1],'ub',[1]),...
                'DisturbanceMatrix', @(t) [1+t^2/2;
                                           t], ...
                'Disturbance', RandomVector('Gaussian', 0, cov_dist));
[Z,H,G]=getConcatMats(sys,3);     
