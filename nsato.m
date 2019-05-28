function nsato(A,B,n)
clf
assert(isequal(size(B),size(A)), 'The matrices do not have the same dimensions or are not square')
D=eye(length(A));%Create a matrix to check the conditions of A and C
idx=zeros(1,length(A)); %Build a permutation array to switch rows of D
idx(1)=length(A);
for j=2:length(A)
    idx(j)=j-1;
end
D=D(idx,:);%permute the columns of matrix
e=zeros(1,length(A));%Make a list to grab values for the magnitude
F=zeros(length(A)); %Make matrix to fill with values for which eigenvalues will be pulled
for a=1:length(A)
    for b=1:length(A)
        if D(a,b)==0
            assert(A(a,b)==0, 'The first input matrix is not a weighted shift matrix of desired form')
        end
        if D(a,b)==1 %for correct spots in weighted shift, we build F from products and collect abs terms 
            e(a)=abs(A(a,b)*B(b,a));
            F(a,b)=A(a,b)*B(b,a);
        end
    end
end
Y=diagrange(A,B,n);
G=eig(F); %get our eigenvalues
H=zeros(1,length(G));
for k=1:length(G)
    H(k)=G(k)/abs(G(k)); %build H as normalized eigenvalues
end
H=H*sum(e); %scale H by sum of norms

hold on
plot(real(Y),imag(Y))
plot(real(H),imag(H),'r*')
plot(real(H),imag(H),'ro')
axis equal
end

