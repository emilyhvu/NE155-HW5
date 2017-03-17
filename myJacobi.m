function [jacobi, J_iter] = myJacobi(n,error,tol)

A=zeros(n,n+2); %A matrix made of zeros with n rows and n+2 columns
for i=1:n %for every row in A
    A(i,i:i+2)=[-1 3 -1]; %replace corresponding indices with recurring coefficients
end
A=A(:,2:n+1); %get rid of first and last columns
b=100*ones(n,1); %b vector with n 100s
x_0=zeros(n,1); %initial guess of all zeros

%jacobi method
D=zeros(n,n);
for i=1:n %extracts elements from diagonal of A and creates D matrix
    D(i,i)=A(i,i);
end

D_inv=zeros(n,n);
for i=1:n %inverts D by taking reciprocal of diagonal elements
    D_inv(i,i)=1/D(i,i);
end
P=D_inv*(D-A); %combines D inverse and (D-A)

jacobi=P*x_0+D_inv*b; %first iteration
J_iter=1;
x=A\b; %actual solution
e_vector=jacobi-x; 
e=norm(e_vector,2); %normalizing error vector

while e>error %tolerance
    jacobi_prev=jacobi; %store previous vector
    jacobi=P*jacobi+D_inv*b; %uses last guess
    J_iter=J_iter+1; %counts iterations
    
    e_vector=jacobi-jacobi_prev; %error vector
    
    if strcmpi(tol,'absolute')
        e=norm(e_vector,2);  %absolute error
    elseif strcmpi(tol,'relative')
        e=norm(e_vector,2)/norm(jacobi,2); %relative error
    end
end
end

