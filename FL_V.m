function [A,B,K] = FL_V()
    clc;clear;close all
    A1 = [0 1 0 0;
          0 0 1 0;
          0 0 0 1;
          0 0 0 0];
    A2 = [0 1;
          0 0];
    A = [A1 zeros([4 4]) zeros([4 4]) zeros([4 2]);
         zeros([4 4]) A1 zeros([4 4]) zeros([4 2]);
         zeros([4 4]) zeros([4 4]) A1 zeros([4 2]);
         zeros([2 4]) zeros([2 4]) zeros([2 4]) A2];

    B1 = [0;
          0;
          0;
          1];
    B = [B1 zeros([4 1]) zeros([4 1]) zeros([4 1]);
         zeros([4 1]) B1 zeros([4 1]) zeros([4 1]);
         zeros([4 1]) zeros([4 1]) B1 zeros([4 1]);
         zeros([2 4])];
    
    Q = diag([1 1 1 1 1 1 1 1 1 1 1 1 1 1]);
    R = diag([1 1 1 1]);

    [K,S,P] = lqr(A,B,Q,R);
end