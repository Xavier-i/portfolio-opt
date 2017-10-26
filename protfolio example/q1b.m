load('hist_data_ps2.mat')
r_p = [1, 1.023, 1.045, 1.13];
x1 = zeros(1,12);
x1(1:4) = 1;
x2 = zeros(1,12);
x2(10:12) = 1;

r_bar = mean(R,1);

b = [1;0.25;0.2;r_p(1)];%r_P(i) i from 1 to 4


A = [ones(1,12);x1; x2;r_bar];% is the A2 matrix in part a
[m,n]=size(A);

%partition

A1 = A(:,1:m);
A2 = A(:,m+1:end);
%Lu method find x_0 and Z
%Can't directly use it since A1 is singular
[L,U,P]=lu(transpose(A));

AP = A*P; 
A1p = AP(:,1:m);
x1 = A1p\b;
x0 = transpose(P)*[x1;zeros(n-m,1)];
%

L1 = L(1:m,:);
L2 = L(m+1:end,:);
%becasue P^T * P = I
Z = transpose(P)*[-transpose(L1)\transpose(L2); eye(n-m)];

% solve w using lu
Aone = r_bar - R;%is tge A1 matrix in part a
A_bar = Aone*Z;
b_bar = Aone*x0;
[Y_bar,R_bar]=qr(A_bar,0);
w_star = - R_bar\(transpose(Y_bar))*b_bar;
x_star = x0+ Z * w_star

risk = (norm(Aone * x_star)^2)/100;
disp('risk is');
disp(risk);

bias = (norm(A * x_star-b));
disp('bias is');
disp(bias);