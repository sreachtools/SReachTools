cvx_clear

%% 4D case
clear
A = Polyhedron('lb',-ones(4,1),'ub',ones(4,1));
B = A;
C = A + B;

disp('>>> C_inner = minkSumInner(A, B)');
tic
C_inner = minkSumInner(A, B);
toc;
if C.contains(C_inner) && size(C_inner.V,1)>1
    disp('Containment rules satisfied');
else
    disp('### Containment rules not satisfied');keyboard
end    

disp('>>> D_inner = minkSumInner(A, B, D)');
D = Polyhedron('V',allcomb([-1,1],[-1,1],[-1,1],[-1,1]));
tic
D_inner = minkSumInner(A, B, D);
toc;
if C.contains(D_inner) && size(D_inner.V,1)>1
    disp('Containment rules satisfied');
else
    disp('### Containment rules not satisfied');keyboard
end    

%% 2D case
clear
disp(' ');
disp('2D case');
A = Polyhedron('lb',-ones(2,1),'ub',ones(2,1));
B = [cos(pi/3), -sin(pi/3);
     sin(pi/3), cos(pi/3)] * A + 10 *ones(2,1);
C = A + B;

disp('>>> C_inner = minkSumInner(A, B)');
tic;
C_inner = minkSumInner(A,B);
toc;
if C.contains(C_inner) && size(C_inner.V,1)>1
    disp('Containment rules satisfied');
else
    disp('### Containment rules not satisfied');keyboard        
end

disp('>>> D_inner = minkSumInner(A, B, D)');
D = Polyhedron('V',[0,1;1,0;-1,0;0,-1]);
tic;
D_inner = minkSumInner(A, B, D);
toc;
if C.contains(D_inner) && size(D_inner.V,1)>1
    disp('Containment rules satisfied');
else
    disp('### Containment rules not satisfied');keyboard        
end

disp('>>> E_inner = minkSumInner(A, B, B)');
tic;
E_inner = minkSumInner(A, B, B - 10*ones(2,1));
toc;
if C.contains(E_inner) && size(E_inner.V,1)>1
    disp('Containment rules satisfied');
else
    disp('### Containment rules not satisfied');keyboard        
end

figure();
plot(C)
hold on
plot(C_inner,'alpha',0.6,'color','g')
plot(D_inner,'alpha',0.6,'color','b')
plot(E_inner,'alpha',0.6,'color','y')
plot(A,'color','c');
plot(B,'color','m');
leg=legend('A+B','minkSumInner(.,.,A)','minkSumInner(.,.,D)','minkSumInner(.,.,B)','A','B');
set(leg,'Location','SouthEast');
axis square
box on;
% disp('Testing cvx based containment');
% for vertex_indx = 1:size(C_inner.V,1)
%     test_vertex = C_inner.V(vertex_indx,:);
%     cvx_begin quiet
%         variable theta(size(A.V,1),1)
%         variable delta(size(B.V,1),1)
%         minimize 0
%         subject to
%             theta'*A.V + delta'*B.V == test_vertex;
%             0<=theta<=1
%             0<=delta<=1
%             sum(delta) == 1;
%             sum(theta) == 1;
%     cvx_end
%     if ~strcmpi(cvx_status,'Solved')
%         fprintf('Failed in vertex_indx: %d at C_inner\n',vertex_indx);
%         disp('### Containment rules not satisfied');keyboard        
%         break;
%     end
% end
