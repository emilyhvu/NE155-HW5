function [SOR, SOR_iter] = mySOR(n,error,tol,w)

A=zeros(n,n+2); %A matrix made of zeros with n rows and n+2 columns
for i=1:n %for every row in A
    A(i,i:i+2)=[-1 3 -1]; %replace corresponding indices with recurring coefficients
end
A=A(:,2:n+1); %get rid of first and last columns
b=100*ones(n,1); %b vector with n 100s
x_0=zeros(n,1); %initial guess of all zeros

D=zeros(n,n);
for i=1:n %constructs diagonal matrix D
    D(i,i)=A(i,i);
end

L=zeros(n,n+2); %L matrix made of zeros with n rows and n+2 columns
for i=1:n %for every row in A
    L(i,i:i+2)=[-1 0 0]; %replace corresponding indices with recurring coefficients (similar to A)
end
L=L(:,2:n+1); %get rid of first and last columns
U=A-L-D; %construct upper matrix

P=(D+w*L)\((1-w)*D-w*U); %multiply -U and L+D

SOR=P*x_0+(D+w*L)\(w*b); %first iteration with backslash to take inverse of L+D
SOR_iter=1;
x=A\b; %actual solution
e_vector=SOR-x; 
e=norm(e_vector,2); %normalizing error vector

while e>error %tolerance
    SOR_prev=SOR; %store previous vector
    SOR=P*SOR_prev+(D+w*L)\(w*b); %uses last guess
    SOR_iter=SOR_iter+1; %counts iterations
    
    e_vector=SOR-SOR_prev; %error vector
    if strcmpi(tol,'absolute')
        e=norm(e_vector,2); %normalize to get absolute error 
    elseif strcmpi(tol,'relative')
        e=norm(e_vector,2)/norm(SOR,2); %finds relative error
    end
end
end


