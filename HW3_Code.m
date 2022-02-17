%Connor Hughes
%CH E 152B HW3

%% Exercise 7: RGA Analysis
K = [-0.04, 0.0005; 0.22, -0.02];
K_inv = inv(K)
RGA = K.*(inv(K)')

%% Exercise 8: RGA Analysis
K = [3.38, 2.50, 0.953; 3.20, 2.38, 0.986; 3.13, 2.33, 1.054];
K_inv = inv(K)
RGA = K.*(inv(K)')  
 

