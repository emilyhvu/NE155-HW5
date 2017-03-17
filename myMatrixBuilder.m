function [A,b]=myMatrixBuilder(n)

b=[0,1:n-1];
b=b'; %transpose b to be a column vector

A=zeros(n,n+2); %zero matrix with n rows and n+2 columns

for i=1:n %for every row in A
    A(i,i:i+2)=[-1 2 -1]; %replace corresponding indices with recurring coefficients
end

A=A(:,2:n+1); %get rid of first and last columns

end
