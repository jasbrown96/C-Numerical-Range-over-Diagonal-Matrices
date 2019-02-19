% The function diagrange below takes 3 inputs:C, A, and n. C and A are the
% two matrices for which we wish to approximate W_diag(C,A). The input n is
% a parameter to increase the precision of the approximation. For kxk
% matrices, this function will create n^k diagonal unitary matrices, which
% have entries evenly spaced via the roots of unity. The output of the
% function, f, will be an array containing the n^k values of tr(CU*AU). This
% can be plotted via the command "plot(real(f),imag(f))"


function f = diagrange(C,A,n) %C fixed, finding W(A), taking n^dim points
assert(isequal(size(C),size(A)), 'The matrices do not have the same dimensions or are not square')
L=zeros(1,n); %Initialize list for roots of unity
for x=1:n
    L(x)=exp(1i*2*pi*x/n); %Fills list with roots of unity
end
B=permn(L,length(C)); %Generates all permutations allowing repetition of roots of unity as a mxn matrix where m is the number of permutation and n is the length of each permutation
[D,E]=size(B); 
F=eye(E); %Initialize Identity matrix of dimension equal to the length of each permutation
G=zeros(1,D); %Initialize a list of length equal to the number of permutations
for x=1:D %For each slot in G, we compute and fill it with the trace of the corresponding diagonal matrix of permutations of the roots of unity
    for y=1:E %Builds the diagonal matrix with permutations of roots of unity grabbed from the matrix B. F(y,y) changes the intialized matrix's diagonal entry and we replace it with B(x,y), which is the corresponding permutation character
        F(y,y)=B(x,y);
    end
    G(x)=trace(C*F'*A*F);
end
f=G; %output the list. plotting this is ideal
    
