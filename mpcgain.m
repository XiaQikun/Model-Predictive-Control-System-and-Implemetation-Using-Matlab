function [Phi_Phi ,Phi_F,Phi_R,A_e,B_e,C_e] = mpcgain(Ap,Bp,Cp,Nc,Np);
[m1,n1] = size(Cp);
[n1,n_in] = size(Bp);
A_e = eye(n1+m1,n1+m1);
A_e(1:n1,1:n1) = Ap;
A_e(n1+1:n1+m1,1:n1) = Cp*Ap;
B_e = zeros(n1+m1,n_in);
B_e(1:n1,:) = Bp;
B_e(n1+1:n1+m1,:) = Cp*Bp;
C_e = zeros(m1,n1+m1);
C_e(:,n1+1:n1+m1) = eye(m1,m1);
% 由于F和Phi特殊的矩阵格式，利用该特点构造矩阵
%先构造矩阵F
n = n1+m1;
h(1,:) = Ce;
F(1,:) = C_e*A_e;
for kk = 2:Np
    h(kk,:) = h(kk-1,:)*A_e;
    F(kk,:) = F(kk-1,:)*A_e;
end
% 构造矩阵Phi
v = h*B_e;
Phi = zeros(Np,Nc); % Phi矩阵的维度
Phi(:,1) = v; %矩阵第一列
for i = 2:Nc
    Phi(:,i) = [zeros(i-1,1);v(1:Np-i+1,1)]; %Toeplitz 矩阵
end
BarRs = ones(Np,1);
Phi_Phi = Phi'*Phi;
Phi_F = Phi'*F;
Phi_R = Phi'*BarRs;
end
